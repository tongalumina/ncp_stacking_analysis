#!/usr/bin/env python
"""
ncp_identification.py

Identifies histone octamers (NCPs) and key DNA features from a PDB file.
Implements a robust, symmetry-based method for finding the central dyad axis.
"""

import argparse
import os
import subprocess
import numpy as np
from Bio.PDB import PDBParser, NeighborSearch
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Data.IUPACData import protein_letters_3to1_extended
from collections import defaultdict
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning

warnings.filterwarnings("ignore", category=PDBConstructionWarning)

# --- Constants ---
DNA_RESIDUES = ['DA', 'DC', 'DG', 'DT']
INTERACTION_DISTANCE = 10.0
NCP_DNA_LENGTH = 147

# --- Helper Functions ---
def get_center_of_mass(entity_or_atom_list):
    if hasattr(entity_or_atom_list, 'get_atoms'):
        atom_list = list(entity_or_atom_list.get_atoms())
    else:
        atom_list = entity_or_atom_list
    if not atom_list: return None
    coords = [atom.get_coord() for atom in atom_list]
    return np.mean(coords, axis=0)

def find_partner_for_residue(target_res, partner_strand_residues):
    COMPLEMENTS = {'DA': 'DT', 'DT': 'DA', 'DC': 'DG', 'DG': 'DC'}
    target_name = target_res.get_resname().strip()
    complement_name = COMPLEMENTS.get(target_name)
    if not complement_name: return None
    best_partner, min_dist = None, float('inf')
    target_com = get_center_of_mass(target_res)
    for res in partner_strand_residues:
        if res.get_parent().id != target_res.get_parent().id and res.get_resname().strip() == complement_name:
            dist = np.linalg.norm(target_com - get_center_of_mass(res))
            if dist < min_dist:
                min_dist, best_partner = dist, res
    if min_dist < 15.0: return best_partner
    return None

# --- Core Logic ---
def get_chain_sequences(structure):
    sequences = {}
    for chain in structure.get_chains():
        seq = ""
        for residue in chain:
            if residue.get_resname().strip() in DNA_RESIDUES or 'CA' not in residue: continue
            res_name = residue.get_resname().strip().title()
            one_letter_code = protein_letters_3to1_extended.get(res_name, 'X')
            seq += one_letter_code
        if seq: sequences[chain.id] = SeqRecord(Seq(seq), id=chain.id, description="")
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
        print(f"BLAST error: {e}")
        return None
    finally:
        if os.path.exists(temp_fasta_query): os.remove(temp_fasta_query)
    return histone_map

def cluster_chains_into_ncps(structure, histone_map):
    histone_chains = {c.id: c for c in structure.get_chains() if c.id in histone_map}
    if not histone_chains: return []
    chain_coms = {cid: get_center_of_mass(c) for cid, c in histone_chains.items()}
    chains_by_type = defaultdict(list)
    for cid, h_type in histone_map.items(): chains_by_type[h_type].append(histone_chains[cid])
    ncp_octamers, assigned_ids = [], set()
    if not chains_by_type['H3']: return []
    for h3 in chains_by_type['H3']:
        if h3.id in assigned_ids: continue
        unassigned = [(cid, com) for cid, com in chain_coms.items() if cid not in assigned_ids and cid != h3.id]
        unassigned.sort(key=lambda x: np.linalg.norm(chain_coms[h3.id] - x[1]))
        potential_ids = [h3.id] + [cid for cid, com in unassigned[:7]]
        stoichiometry = defaultdict(int)
        for cid in potential_ids: stoichiometry[histone_map[cid]] += 1
        if all(stoichiometry.get(ht, 0) == 2 for ht in ['H2A', 'H2B', 'H3', 'H4']):
            ncp_octamers.append([histone_chains[cid] for cid in potential_ids])
            assigned_ids.update(potential_ids)
    ncp_octamers.sort(key=lambda o: get_center_of_mass([a for c in o for a in c.get_atoms()])[2])
    return ncp_octamers

def find_ncp_features(octamer, histone_map, structure):
    # 1. Find all interacting BPs for the current octamer
    histone_atoms = [a for c in octamer for a in c.get_atoms()]
    all_dna_residues = [r for r in structure.get_residues() if r.get_resname() in DNA_RESIDUES]
    ns = NeighborSearch(histone_atoms)
    interacting_dna_residues = {r for r in all_dna_residues if ns.search(get_center_of_mass(r), INTERACTION_DISTANCE, level='R')}
    if not interacting_dna_residues: return {'error': 'No interacting DNA residues found.'}

    strand_map = defaultdict(list)
    for res in interacting_dna_residues: strand_map[res.get_parent().id].append(res)
    if len(strand_map) < 2: return {'error': 'Could not identify two distinct DNA strands.'}
    main_strand_ids = sorted(strand_map.keys(), key=lambda k: len(strand_map[k]), reverse=True)[:2]
    strand1_res, strand2_res = strand_map[main_strand_ids[0]], strand_map[main_strand_ids[1]]
    all_bps = { (r1, p) for r1 in strand1_res if (p := find_partner_for_residue(r1, strand2_res)) }
    if not all_bps: return {'error': 'Could not form any base pairs.'}

    # 2. Define symmetry axes from global properties
    h3_chains = [c for c in octamer if histone_map.get(c.id) == 'H3']
    h3_coms = [get_center_of_mass(c) for c in h3_chains]
    v_h3_h3 = h3_coms[0] - h3_coms[1]

    dna_coords = np.array([get_center_of_mass(bp[0]) for bp in all_bps])
    _, _, vh = np.linalg.svd(dna_coords - np.mean(dna_coords, axis=0))
    v_superhelix = vh[2,:]

    dyad_axis = np.cross(v_h3_h3, v_superhelix)
    dyad_axis /= np.linalg.norm(dyad_axis)

    # 3. Find the dyad BP (closest to the calculated dyad axis)
    gho_com = get_center_of_mass(histone_atoms)
    min_dist, dyad_bp = float('inf'), None
    for bp in all_bps:
        bp_com = get_center_of_mass(list(bp[0].get_atoms()) + list(bp[1].get_atoms()))
        vec_to_bp = bp_com - gho_com
        dist_to_axis = np.linalg.norm(vec_to_bp - np.dot(vec_to_bp, dyad_axis) * dyad_axis)
        if dist_to_axis < min_dist:
            min_dist, dyad_bp = dist_to_axis, bp
    if not dyad_bp: return {'error': 'Could not find a dyad BP close to the symmetry axis.'}

    # 4. Find plane-defining BPs (+/- 20 from the dyad)
    center_res = next((r for r in dyad_bp if r.get_parent().id == main_strand_ids[0]), dyad_bp[0])
    center_idx = center_res.id[1]
    parent_chain = center_res.get_parent()
    
    # Get the full chain objects for robust partner searching
    chain1_obj = structure[0][main_strand_ids[0]]
    chain2_obj = structure[0][main_strand_ids[1]]
    partner_strand_obj = chain2_obj if parent_chain.id == main_strand_ids[0] else chain1_obj

    plane_res_neg = next((r for r in parent_chain if r.id[1] == center_idx - 20), None)
    plane_res_pos = next((r for r in parent_chain if r.id[1] == center_idx + 20), None)
    if not plane_res_neg or not plane_res_pos: return {'error': 'Could not find residues +/- 20 from center.'}

    plane_partner_neg = find_partner_for_residue(plane_res_neg, partner_strand_obj)
    plane_partner_pos = find_partner_for_residue(plane_res_pos, partner_strand_obj)
    if not plane_partner_neg or not plane_partner_pos: return {'error': 'Could not find partners for plane-defining residues.'}

    # 5. Define the 147bp DNA segment around the dyad
    dna_segments = []
    half_len = (NCP_DNA_LENGTH - 1) // 2
    for res in dyad_bp:
        dna_segments.append((res.get_parent().id, res.id[1] - half_len, res.id[1] + half_len))

    return {
        'central_bp': dyad_bp,
        'plane_bp_neg': (plane_res_neg, plane_partner_neg),
        'plane_bp_pos': (plane_res_pos, plane_partner_pos),
        'dna_segments': dna_segments
    }

def generate_config_file(pdb_id, ncp_octamers, histone_map, structure, output_file):
    with open(output_file, 'w') as f:
        f.write(f"# NCP Configuration File for {pdb_id}\n")
        for i, octamer in enumerate(ncp_octamers):
            f.write(f"\nNCP: {i + 1}\n")
            histone_ids = sorted([c.id for c in octamer])
            f.write(f"HISTONE_CHAINS: {' '.join(histone_ids)}\n")
            map_items = [f"{hid}:{histone_map.get(hid, '?')}" for hid in histone_ids]
            f.write(f"HISTONE_MAP: {' '.join(map_items)}\n")

            features = find_ncp_features(octamer, histone_map, structure)
            if features and 'error' not in features:
                def bp_to_str(bp): return f"{bp[0].get_parent().id}:{bp[0].id[1]} {bp[1].get_parent().id}:{bp[1].id[1]}"
                f.write(f"CENTRAL_BP: {bp_to_str(features['central_bp'])}\n")
                f.write(f"PLANE_BP_NEG: {bp_to_str(features['plane_bp_neg'])}\n")
                f.write(f"PLANE_BP_POS: {bp_to_str(features['plane_bp_pos'])}\n")
                for chain, start, end in features['dna_segments']:
                    f.write(f"DNA_SEGMENT: {chain} {start} {end}\n")
            elif features and 'error' in features:
                f.write(f"# Could not identify features for this NCP. Reason: {features['error']}\n")
            else:
                f.write("# Could not identify features for this NCP for an unknown reason.\n")

def main():
    parser = argparse.ArgumentParser(description="Identify NCPs and features for analysis.")
    parser.add_argument("pdb_file", help="Path to the input PDB file.")
    parser.add_argument("--histone-fasta", required=True, help="Path to FASTA file of canonical histone sequences.")
    parser.add_argument("--output-config", default=None, help="Output config file. Defaults to <pdb_filename>_config.txt")
    args = parser.parse_args()

    if not os.path.exists(args.pdb_file): return print(f"Error: PDB file not found: {args.pdb_file}")
    if not os.path.exists(args.histone_fasta): return print(f"Error: FASTA file not found: {args.histone_fasta}")

    pdb_id = os.path.splitext(os.path.basename(args.pdb_file))[0]
    db_path = os.path.splitext(args.histone_fasta)[0]
    output_config = args.output_config if args.output_config else f"{pdb_id}_config.txt"

    if not os.path.exists(f"{db_path}.phr"):
        print(f"Creating BLAST database from {args.histone_fasta}...")
        try:
            subprocess.run(["makeblastdb", "-in", args.histone_fasta, "-dbtype", "prot", "-out", db_path, "-parse_seqids"], check=True, capture_output=True)
        except Exception as e: return print(f"Error creating BLAST db: {e}")

    print(f"Loading structure from {args.pdb_file}...")
    structure = PDBParser().get_structure(pdb_id, args.pdb_file)
    
    print("Identifying histones...")
    sequences = get_chain_sequences(structure)
    histone_map = identify_histones_by_blast(sequences, db_path)
    if not histone_map: return print("Could not identify any histones.")

    print("Clustering histone octamers...")
    ncp_octamers = cluster_chains_into_ncps(structure, histone_map)
    if not ncp_octamers: return print("Could not cluster any histone octamers.")

    print(f"Identifying features and writing config file...")
    generate_config_file(pdb_id, ncp_octamers, histone_map, structure, output_config)
    print(f"\nIdentification complete. Config written to {output_config}")

if __name__ == "__main__":
    main()
