#!/usr/bin/env python
"""
ncp_identification.py

Identifies histone octamers (NCPs) and key DNA features from a PDB file.
This version uses the 'biotite' library for robust structure analysis and base pair identification.
"""

import argparse
import os
import subprocess
import numpy as np
import biotite.structure as struc
import biotite.structure.io as strucio
from biotite.structure.io.pdb import PDBFile
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from collections import defaultdict
import warnings
import logging

warnings.filterwarnings("ignore", message=".*has less than 3 base atoms.*")

# --- Constants ---
DNA_RESIDUES = ['DA', 'DC', 'DG', 'DT']
INTERACTION_DISTANCE = 10.0
NCP_DNA_LENGTH = 147

HISTONE_CORE_RANGES = {
    'H2A': (21, 116),
    'H2B': (34, 121),
    'H3': (44, 135),
    'H4': (24, 102),
}

# --- Helper Functions ---
def get_center_of_mass(atom_array):
    """Calculates the center of mass for a biotite AtomArray."""
    if atom_array.array_length() == 0:
        return None
    return np.mean(atom_array.coord, axis=0)

def get_core_center_of_mass(chain, histone_type):
    """Get center of mass using only the structured core region of a histone chain."""
    if histone_type not in HISTONE_CORE_RANGES:
        return get_center_of_mass(chain)
    
    core_start, core_end = HISTONE_CORE_RANGES[histone_type]
    # Filter atoms based on residue ID range
    core_filter = (chain.res_id >= core_start) & (chain.res_id <= core_end)
    core_atoms = chain[core_filter]
    
    return get_center_of_mass(core_atoms) if core_atoms.array_length() > 0 else get_center_of_mass(chain)

# --- Core Logic ---
def get_chain_sequences(structure):
    """Extracts protein sequences from a biotite AtomArray."""
    sequences = {}
    for chain in struc.chain_iter(structure):
        protein_chain = chain[struc.filter_amino_acids(chain)]
        if protein_chain.array_length() > 0:
            sequence_obj = struc.to_sequence(protein_chain)
            sequence_str = str(sequence_obj)
            chain_id = protein_chain.chain_id[0]
            sequences[chain_id] = SeqRecord(Seq(sequence_str), id=chain_id, description="")
    return sequences

def identify_histones_by_blast(sequences, db_path):
    histone_map = {}
    temp_fasta_query = f"temp_query_{os.getpid()}.fasta"
    SeqIO.write(list(sequences.values()), temp_fasta_query, "fasta")
    cmd = ["blastp", "-query", temp_fasta_query, "-db", db_path, "-outfmt", "6 qseqid sseqid pident", "-max_target_seqs", "1"]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        for line in result.stdout.strip().split('\n'):
            if not line: continue
            parts = line.split()
            chain_id, subject_id, pident = parts[0], parts[1], float(parts[2])
            histone_type = subject_id.split('_')[0]
            if pident > 50: histone_map[chain_id] = histone_type
    except Exception as e:
        logging.error(f"BLAST error: {e}")
        return None
    finally:
        if os.path.exists(temp_fasta_query): os.remove(temp_fasta_query)
    return histone_map

def cluster_chains_into_ncps(structure, histone_map):
    histone_chains = {chain_id: structure[structure.chain_id == chain_id] 
                      for chain_id in histone_map.keys()}
    if not histone_chains: return []

    # Use CORE center of mass for more stable clustering
    chain_coms = {cid: get_core_center_of_mass(chain, histone_map.get(cid)) 
                  for cid, chain in histone_chains.items()}
    
    chains_by_type = defaultdict(list)
    for cid, h_type in histone_map.items(): chains_by_type[h_type].append(histone_chains[cid])
    
    ncp_octamers, assigned_ids = [], set()
    if not chains_by_type.get('H3'): return []

    for h3 in chains_by_type['H3']:
        h3_id = h3.chain_id[0]
        if h3_id in assigned_ids: continue
        
        unassigned = [(cid, com) for cid, com in chain_coms.items() if cid not in assigned_ids and cid != h3_id]
        unassigned.sort(key=lambda x: np.linalg.norm(chain_coms[h3_id] - x[1]))
        
        potential_ids = [h3_id] + [cid for cid, com in unassigned[:7]]
        stoichiometry = defaultdict(int)
        for cid in potential_ids: stoichiometry[histone_map[cid]] += 1
        
        if all(stoichiometry.get(ht, 0) == 2 for ht in ['H2A', 'H2B', 'H3', 'H4']):
            mask = np.isin(structure.chain_id, potential_ids)
            octamer_atoms = structure[mask]
            ncp_octamers.append(octamer_atoms)
            assigned_ids.update(potential_ids)
            
    # Sorting is now handled in generate_config_file
    return ncp_octamers

def find_ncp_features(octamer_atoms, histone_map, full_structure, args):
    logging.info(f"--- Processing NCP ---")
    octamer_chain_ids = np.unique(octamer_atoms.chain_id)
    logging.info(f"Histone chains: {octamer_chain_ids}")

    dna_structure = full_structure[struc.filter_nucleotides(full_structure)]

    # Find interacting DNA residues using NumPy for guaranteed correctness.
    octamer_coords = octamer_atoms.coord
    dna_coords = dna_structure.coord
    diffs = dna_coords - octamer_coords[:, np.newaxis]
    dist_sq = np.sum(diffs**2, axis=-1)
    min_dist_sq = np.min(dist_sq, axis=0)
    interacting_dna_indices = np.where(min_dist_sq < INTERACTION_DISTANCE**2)[0]
    interacting_dna = dna_structure[interacting_dna_indices]
    
    if interacting_dna.array_length() == 0:
        return {'error': 'No interacting DNA residues found.'}
    logging.info(f"Found {interacting_dna.array_length()} interacting DNA atoms.")

    base_pairs = struc.base_pairs(interacting_dna)
    if base_pairs.shape[0] == 0:
        return {'error': 'biotite found no canonical base pairs in the interacting DNA.'}
    logging.info(f"biotite found {base_pairs.shape[0]} base pairs.")

    dyad_bp_info = None
    h3_chains = [octamer_atoms[octamer_atoms.chain_id == cid] for cid in octamer_chain_ids if histone_map[cid] == 'H3']
    logging.info(f"H3 chains for dyad calculation: {[c.chain_id[0] for c in h3_chains]}")

    h3_116_atoms = [c[c.res_id == 116] for c in h3_chains]
    h3_116_atoms = [a for a in h3_116_atoms if a.array_length() > 0]

    if len(h3_116_atoms) == 2:
        r116_midpoint = np.mean([get_center_of_mass(a) for a in h3_116_atoms], axis=0)
        logging.info(f"Using H3 residue 116 restraint. Midpoint at: {r116_midpoint}")
        
        min_dist = float('inf')
        for bp_indices in base_pairs:
            res1_id = interacting_dna.res_id[bp_indices[0]]
            res1_chain_id = interacting_dna.chain_id[bp_indices[0]]
            res2_id = interacting_dna.res_id[bp_indices[1]]
            res2_chain_id = interacting_dna.chain_id[bp_indices[1]]
            mask = (
                (interacting_dna.chain_id == res1_chain_id) & (interacting_dna.res_id == res1_id) |
                (interacting_dna.chain_id == res2_chain_id) & (interacting_dna.res_id == res2_id)
            )
            bp_atoms = interacting_dna[mask]
            com = get_center_of_mass(bp_atoms)
            if com is None: continue
            dist = np.linalg.norm(com - r116_midpoint)
            if dist < min_dist:
                min_dist = dist
                dyad_bp_info = {
                    'res1': {'chain_id': res1_chain_id, 'res_id': res1_id},
                    'res2': {'chain_id': res2_chain_id, 'res_id': res2_id}
                }
    else:
        missing_chains = [c.chain_id[0] for c in h3_chains if c[c.res_id == 116].array_length() == 0]
        logging.warning(f"Could not find H3 residue 116 in H3 chain(s): {missing_chains}. Falling back to geometric method.")
        middle_bp_indices = base_pairs[len(base_pairs) // 2]
        res1_id = interacting_dna.res_id[middle_bp_indices[0]]
        res1_chain_id = interacting_dna.chain_id[middle_bp_indices[0]]
        res2_id = interacting_dna.res_id[middle_bp_indices[1]]
        res2_chain_id = interacting_dna.chain_id[middle_bp_indices[1]]
        dyad_bp_info = {
            'res1': {'chain_id': res1_chain_id, 'res_id': res1_id},
            'res2': {'chain_id': res2_chain_id, 'res_id': res2_id}
        }

    if dyad_bp_info is None:
        return {'error': 'Could not identify a dyad base pair.'}

    dyad_res1_atoms = interacting_dna[(interacting_dna.chain_id == dyad_bp_info['res1']['chain_id']) & (interacting_dna.res_id == dyad_bp_info['res1']['res_id'])]
    dyad_res2_atoms = interacting_dna[(interacting_dna.chain_id == dyad_bp_info['res2']['chain_id']) & (interacting_dna.res_id == dyad_bp_info['res2']['res_id'])]
    logging.info(f"Identified dyad BP: {dyad_res1_atoms[0].chain_id}:{dyad_res1_atoms[0].res_id} - {dyad_res2_atoms[0].chain_id}:{dyad_res2_atoms[0].res_id}")

    # --- New plane-defining BP logic ---
    all_bps = struc.base_pairs(dna_structure)
    bp_partner_map = {}
    for i, j in all_bps:
        res1_info = (dna_structure.chain_id[i], dna_structure.res_id[i])
        res2_info = (dna_structure.chain_id[j], dna_structure.res_id[j])
        bp_partner_map[res1_info] = res2_info
        bp_partner_map[res2_info] = res1_info

    def find_closest_bp(ref_chain, target_resid, window_size, bp_map):
        best_residue_info = None
        min_dist = float('inf')
        for offset in range(-window_size, window_size + 1):
            current_resid = target_resid + offset
            current_res_info = (ref_chain, current_resid)
            if current_res_info in bp_map:
                dist = abs(offset)
                if dist < min_dist:
                    min_dist = dist
                    best_residue_info = current_res_info
        if best_residue_info:
            return best_residue_info, bp_map[best_residue_info]
        else:
            return None, None

    dyad_refs = [
        (dyad_res1_atoms[0].chain_id, dyad_res1_atoms[0].res_id),
        (dyad_res2_atoms[0].chain_id, dyad_res2_atoms[0].res_id)
    ]
    
    plane_bp_found = False
    last_error_info = ""
    search_window = 5 # Search +/- 5 residues around the target

    for ref_chain, ref_resid in dyad_refs:
        target_neg_id = ref_resid - 20
        target_pos_id = ref_resid + 20
        
        plane_neg_res1_info, plane_neg_res2_info = find_closest_bp(ref_chain, target_neg_id, search_window, bp_partner_map)
        plane_pos_res1_info, plane_pos_res2_info = find_closest_bp(ref_chain, target_pos_id, search_window, bp_partner_map)

        if plane_neg_res1_info and plane_pos_res1_info:
            plane_bp_found = True
            break 
        else:
            last_error_info = f"Could not find valid BPs near {ref_chain}:{target_neg_id} or {ref_chain}:{target_pos_id}"
    
    if not plane_bp_found:
        return {'error': f'Could not find plane-defining base pairs. {last_error_info}'}

    plane_res_neg1_atoms = dna_structure[(dna_structure.chain_id == plane_neg_res1_info[0]) & (dna_structure.res_id == plane_neg_res1_info[1])]
    plane_res_neg2_atoms = dna_structure[(dna_structure.chain_id == plane_neg_res2_info[0]) & (dna_structure.res_id == plane_neg_res2_info[1])]
    plane_res_pos1_atoms = dna_structure[(dna_structure.chain_id == plane_pos_res1_info[0]) & (dna_structure.res_id == plane_pos_res1_info[1])]
    plane_res_pos2_atoms = dna_structure[(dna_structure.chain_id == plane_pos_res2_info[0]) & (dna_structure.res_id == plane_pos_res2_info[1])]

    if any(a.array_length() == 0 for a in [plane_res_neg1_atoms, plane_res_neg2_atoms, plane_res_pos1_atoms, plane_res_pos2_atoms]):
        return {'error': 'Could not find all atom records for plane-defining base pairs at dyad +/- 20 bp.'}

    logging.info("Successfully found plane-defining residues.")

    dna_segments = []
    half_len = (NCP_DNA_LENGTH - 1) // 2
    dna_segments.append((dyad_res1_atoms[0].chain_id, dyad_res1_atoms[0].res_id - half_len, dyad_res1_atoms[0].res_id + half_len))
    dna_segments.append((dyad_res2_atoms[0].chain_id, dyad_res2_atoms[0].res_id - half_len, dyad_res2_atoms[0].res_id + half_len))

    return {
        'central_bp': (dyad_res1_atoms[0], dyad_res2_atoms[0]),
        'plane_bp_neg': (plane_res_neg1_atoms[0], plane_res_neg2_atoms[0]),
        'plane_bp_pos': (plane_res_pos1_atoms[0], plane_res_pos2_atoms[0]),
        'dna_segments': dna_segments,
    }

def generate_config_file(pdb_id, ncp_octamers, histone_map, structure, output_file, args):
    
    ncp_data = []
    for octamer_atoms in ncp_octamers:
        features = find_ncp_features(octamer_atoms, histone_map, structure, args)
        if features and 'error' not in features:
            res1, res2 = features['central_bp']
            # Sort key is the minimum of the two (chain_id, res_id) tuples from the central BP
            sort_key = min((res1.chain_id, res1.res_id), (res2.chain_id, res2.res_id))
            ncp_data.append({'sort_key': sort_key, 'octamer': octamer_atoms, 'features': features})
        else:
            # Handle cases where features could not be identified
            # Give it a sort key that places it at the end
            ncp_data.append({'sort_key': ('~', float('inf')), 'octamer': octamer_atoms, 'features': features})

    # Sort NCPs based on the defined key
    ncp_data.sort(key=lambda x: x['sort_key'])

    with open(output_file, 'w') as f:
        f.write(f"# NCP Configuration File for {pdb_id}\n")
        f.write(f"# Generated using biotite for base-pair analysis.\n")
        f.write(f"# NCPs are sorted by their position on the DNA strand.\n")
        
        for i, data in enumerate(ncp_data):
            octamer_atoms = data['octamer']
            features = data['features']
            
            f.write(f"\nNCP: {i + 1}\n")
            histone_ids = sorted(np.unique(octamer_atoms.chain_id))
            f.write(f"HISTONE_CHAINS: {' '.join(histone_ids)}\n")
            map_items = [f"{hid}:{histone_map.get(hid, '?')}" for hid in histone_ids]
            f.write(f"HISTONE_MAP: {' '.join(map_items)}\n")

            if features and 'error' not in features:
                def res_to_str(res): return f"{res.chain_id}:{res.res_id}"
                f.write(f"CENTRAL_BP: {res_to_str(features['central_bp'][0])} {res_to_str(features['central_bp'][1])}\n")
                f.write(f"PLANE_BP_NEG: {res_to_str(features['plane_bp_neg'][0])} {res_to_str(features['plane_bp_neg'][1])}\n")
                f.write(f"PLANE_BP_POS: {res_to_str(features['plane_bp_pos'][0])} {res_to_str(features['plane_bp_pos'][1])}\n")
                for chain, start, end in features['dna_segments']:
                    f.write(f"DNA_SEGMENT: {chain} {start} {end}\n")
            elif features and 'error' in features:
                f.write(f"# Could not identify features for this NCP. Reason: {features['error']}\n")
            else:
                f.write("# Could not identify features for this NCP for an unknown reason.\n")


def main():
    parser = argparse.ArgumentParser(description="Identify NCPs and features for analysis using biotite.")
    parser.add_argument("pdb_file", help="Path to the input PDB file.")
    parser.add_argument("-f", "--histone-fasta", default="histones.fa", 
                       help="Path to FASTA file of canonical histone sequences (default: histones.fa)")
    parser.add_argument("-o", "--output-config", default=None, 
                       help="Output config file. Defaults to <pdb_filename>_config.txt")
    parser.add_argument("--strict", action="store_true", 
                       help="Use stricter criteria for base pairing (not used with biotite)." )
    
    args = parser.parse_args()

    pdb_id = os.path.splitext(os.path.basename(args.pdb_file))[0]
    log_file = f"{pdb_id}_id.log"
    logging.basicConfig(filename=log_file, level=logging.INFO, 
                        format='%(asctime)s - %(levelname)s - %(message)s', filemode='w')

    logging.info(f"Starting NCP identification for {args.pdb_file}")

    if not os.path.exists(args.pdb_file): 
        logging.error(f"PDB file not found: {args.pdb_file}")
        return print(f"Error: PDB file not found: {args.pdb_file}")
    if not os.path.exists(args.histone_fasta): 
        logging.error(f"FASTA file not found: {args.histone_fasta}")
        return print(f"Error: FASTA file not found: {args.histone_fasta}")

    db_path = os.path.splitext(args.histone_fasta)[0]
    output_config = args.output_config if args.output_config else f"{pdb_id}_config.txt"

    if not os.path.exists(f"{db_path}.phr"):
        logging.info(f"Creating BLAST database from {args.histone_fasta}...")
        try:
            subprocess.run(["makeblastdb", "-in", args.histone_fasta, "-dbtype", "prot", "-out", db_path, "-parse_seqids"], check=True, capture_output=True)
        except Exception as e: 
            logging.error(f"Error creating BLAST db: {e}")
            return print(f"Error creating BLAST db: {e}")

    logging.info(f"Loading structure with biotite from {args.pdb_file}...")
    try:
        # Explicitly use the PDB file parser to avoid issues with non-standard extensions
        pdb_f = PDBFile.read(args.pdb_file)
        # Get only the first model to handle multi-model PDBs
        structure = pdb_f.get_structure(model=1)
    except Exception as e:
        logging.error(f"Could not load structure with biotite: {e}")
        return print(f"Error: Could not load PDB file with biotite.")

    logging.info("Identifying histones...")
    sequences = get_chain_sequences(structure)
    histone_map = identify_histones_by_blast(sequences, db_path)
    if not histone_map: 
        logging.error("Could not identify any histones.")
        return print("Could not identify any histones.")

    logging.info("Clustering histone octamers...")
    ncp_octamers = cluster_chains_into_ncps(structure, histone_map)
    if not ncp_octamers:
        logging.error("Could not cluster any histone octamers.")
        return print("Could not cluster any histone octamers.")

    logging.info(f"Identifying features and writing config file to {output_config}...")
    generate_config_file(pdb_id, ncp_octamers, histone_map, structure, output_config, args)
    print(f"\nIdentification complete. Config written to {output_config}. See {log_file} for details.")

if __name__ == "__main__":
    main()
