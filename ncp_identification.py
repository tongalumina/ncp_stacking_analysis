#!/usr/bin/env python
"""
ncp_identification.py

Identifies histone octamers (NCPs) and key DNA features from a PDB file.

Workflow:
1.  Reads a PDB structure.
2.  Extracts protein sequences and identifies histone chains (H2A, H2B, H3, H4)
    by running a BLAST search against a user-provided histone sequence database.
3.  Clusters the identified histones into octamers based on proximity.
4.  For each NCP, identifies the central DNA base pair that defines the local
    symmetry axis. This is done by finding the midpoint between the two
    C-alpha atoms of Arginine 49 (R49) on the two H3 histones, and then
    finding the closest DNA base pair to this point.
5.  Identifies the 147bp DNA segments wrapped around each octamer.
6.  Writes the results to a configuration file for downstream analysis.

Usage:
python ncp_identification.py <PDB_FILE> --histone-fasta <FASTA_FILE> --output-config <CONFIG_FILE>
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
INTERACTION_DISTANCE = 10.0  # Angstroms

def get_chain_sequences(structure):
    """Extracts protein sequences from a PDB structure."""
    sequences = {}
    for chain in structure.get_chains():
        seq = ""
        # Reset residue list for each chain
        residues = list(chain.get_residues())
        for residue in residues:
            # Skip non-protein or DNA residues
            if residue.get_resname().strip() in DNA_RESIDUES or 'CA' not in residue:
                continue
            res_name = residue.get_resname().strip().title()
            one_letter_code = protein_letters_3to1_extended.get(res_name, 'X')
            seq += one_letter_code
        if seq:
            sequences[chain.id] = SeqRecord(Seq(seq), id=chain.id, description="")
    return sequences

def identify_histones_by_blast(sequences, db_path):
    """Identifies histone type for each chain using BLASTP."""
    histone_map = {}
    # Create a temporary file for all sequences
    temp_fasta_query = "temp_query.fasta"
    SeqIO.write(sequences.values(), temp_fasta_query, "fasta")

    # Run BLASTP
    cmd = [
        "blastp", "-query", temp_fasta_query, "-db", db_path,
        "-outfmt", "6 qseqid sseqid pident", "-max_target_seqs", "1"
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"Error running BLAST: {e}")
        print("Please ensure NCBI BLAST+ is installed and in your system's PATH.")
        os.remove(temp_fasta_query)
        return None

    os.remove(temp_fasta_query)

    # Process BLAST results
    for line in result.stdout.strip().split('\n'):
        if not line:
            continue
        parts = line.split()
        chain_id, subject_id, pident = parts[0], parts[1], float(parts[2])
        histone_type = subject_id.split('_')[0] # Assuming DB format like 'H3_human'
        if pident > 50:  # Confidence threshold
            histone_map[chain_id] = histone_type

    return histone_map

def get_center_of_mass(entity_or_atom_list):
    """Calculates the center of mass of a PDB entity or a list of atoms."""
    if hasattr(entity_or_atom_list, 'get_atoms'): # It's an entity
        atom_list = list(entity_or_atom_list.get_atoms())
    else: # Assume it's a list of atoms
        atom_list = entity_or_atom_list
    
    if not atom_list:
        return None
    coords = [atom.get_coord() for atom in atom_list]
    return np.mean(coords, axis=0)

def cluster_chains_into_ncps(structure, histone_map):
    """Clusters histone chains into octamers based on proximity."""
    histone_chains = {chain.id: chain for chain in structure.get_chains() if chain.id in histone_map}
    if not histone_chains:
        return []

    # Pre-calculate COMs for all histone chains
    chain_coms = {cid: get_center_of_mass(c) for cid, c in histone_chains.items()}
    
    # Group chains by histone type
    chains_by_type = defaultdict(list)
    for chain_id, h_type in histone_map.items():
        chains_by_type[h_type].append(histone_chains[chain_id])

    ncp_histone_octamers = []
    assigned_chain_ids = set()

    # Use H3 as a seed for clustering
    if not chains_by_type['H3']:
        return []

    for h3_chain in chains_by_type['H3']:
        if h3_chain.id in assigned_chain_ids:
            continue

        h3_com = chain_coms[h3_chain.id]
        
        # Find the 7 closest histone chains
        unassigned_chains = [(cid, com) for cid, com in chain_coms.items() if cid not in assigned_chain_ids and cid != h3_chain.id]
        unassigned_chains.sort(key=lambda x: np.linalg.norm(h3_com - x[1]))
        
        closest_7_ids = [cid for cid, com in unassigned_chains[:7]]
        potential_octamer_ids = [h3_chain.id] + closest_7_ids
        
        # Check stoichiometry
        stoichiometry = defaultdict(int)
        for cid in potential_octamer_ids:
            stoichiometry[histone_map[cid]] += 1
            
        is_octamer = (stoichiometry['H3'] == 2 and stoichiometry['H4'] == 2 and
                      stoichiometry['H2A'] == 2 and stoichiometry['H2B'] == 2)

        if is_octamer:
            octamer_chains = [histone_chains[cid] for cid in potential_octamer_ids]
            ncp_histone_octamers.append(octamer_chains)
            assigned_chain_ids.update(potential_octamer_ids)

    # Sort NCPs by Z-coordinate of their COM for consistent ordering
    ncp_histone_octamers.sort(key=lambda octamer: get_center_of_mass([atom for chain in octamer for atom in chain.get_atoms()]).tolist()[2])
    return ncp_histone_octamers

def find_watson_crick_partner(target_res, partner_strand_residues):
    """    Finds the best partner for a residue.
    1. Tries to find the best complementary partner based on H-bond atom distance.
    2. If no complement is found, falls back to the geometrically closest residue.
    Returns a tuple: (partner_residue, distance, is_complementary_flag).
    """
    COMPLEMENTS = {'DA': 'DT', 'DT': 'DA', 'DC': 'DG', 'DG': 'DC'}
    PURINES = ['DA', 'DG']

    target_res_name = target_res.get_resname().strip()
    complement_name = COMPLEMENTS.get(target_res_name)

    best_partner = None
    min_dist = float('inf')
    is_complementary = False

    # 1. First, search for the best COMPLEMENTARY partner
    if complement_name:
        try:
            target_atom = target_res['N1'] if target_res_name in PURINES else target_res['N3']
            for partner_res in partner_strand_residues:
                if partner_res.get_resname().strip() == complement_name:
                    try:
                        partner_atom = partner_res['N1'] if partner_res.get_resname().strip() in PURINES else partner_res['N3']
                        dist = target_atom - partner_atom
                        if dist < min_dist:
                            min_dist = dist
                            best_partner = partner_res
                    except KeyError:
                        continue
            if best_partner:
                is_complementary = True
        except KeyError:
            pass # Target atom not found, will proceed to geometric search

    # 2. If no complementary partner found, fall back to GEOMETRIC closest
    if not best_partner:
        target_com = get_center_of_mass(target_res)
        for partner_res in partner_strand_residues:
            dist = np.linalg.norm(target_com - get_center_of_mass(partner_res))
            if dist < min_dist:
                min_dist = dist
                best_partner = partner_res
        is_complementary = False # It's a geometric match

    return best_partner, min_dist, is_complementary

def find_central_bp_and_dna_segment(structure, octamer, histone_map):
    """
    Identifies the central base pair and the 147bp DNA segment for an NCP.
    """
    # 1. Find H3 R49 anchor point
    h3_chains = [c for c in octamer if histone_map[c.id] == 'H3']
    if len(h3_chains) != 2:
        return None, None, None

    try:
        r49_ca_1 = h3_chains[0][49]['CA'].get_coord()
        r49_ca_2 = h3_chains[1][49]['CA'].get_coord()
    except KeyError:
        print(f"Warning: H3 Arginine 49 not found in chains {h3_chains[0].id} or {h3_chains[1].id}. Skipping NCP.")
        return None, None, None
        
    h3_r49_midpoint = (r49_ca_1 + r49_ca_2) / 2.0

    # 2. Find interacting DNA chains
    histone_atoms = [atom for chain in octamer for atom in chain.get_atoms()]
    all_dna_residues = [res for res in structure.get_residues() if res.get_resname() in DNA_RESIDUES]
    
    ns = NeighborSearch(histone_atoms)
    interacting_dna_residues = {res for res in all_dna_residues if ns.search(get_center_of_mass(res), INTERACTION_DISTANCE, level='R')}

    if not interacting_dna_residues:
        return None, None, None

    # 3. Find the closest (pseudo) base pair to the anchor point
    min_dist = float('inf')
    central_bp_geometric = None
    
    dna_by_chain = defaultdict(list)
    for res in interacting_dna_residues:
        dna_by_chain[res.get_parent().id].append(res)
    
    if len(dna_by_chain) < 2: return None, None, None
    
    main_dna_chain_ids = sorted(dna_by_chain.keys(), key=lambda k: len(dna_by_chain[k]), reverse=True)[:2]
    strand1_res = sorted(dna_by_chain[main_dna_chain_ids[0]], key=lambda r: r.id[1])
    strand2_res = sorted(dna_by_chain[main_dna_chain_ids[1]], key=lambda r: r.id[1])

    for res1 in strand1_res:
        partner = min(strand2_res, key=lambda res2: np.linalg.norm(get_center_of_mass(res1) - get_center_of_mass(res2)))
        bp_com = (get_center_of_mass(res1) + get_center_of_mass(partner)) / 2.0
        dist = np.linalg.norm(bp_com - h3_r49_midpoint)
        if dist < min_dist:
            min_dist = dist
            central_bp_geometric = (res1, partner)

    if not central_bp_geometric:
        return None, None, None

    # 4. Find the actual Watson-Crick partners (FIX 1: search entire chain)
    res_a, res_b = central_bp_geometric
    full_strand1_residues = list(res_a.get_parent().get_residues())
    full_strand2_residues = list(res_b.get_parent().get_residues())

    partner_of_a = find_watson_crick_partner(res_a, full_strand2_residues)
    partner_of_b = find_watson_crick_partner(res_b, full_strand1_residues)
    actual_pairing_info = (res_a, partner_of_a, res_b, partner_of_b)

    # 5. Define the 147bp DNA segment (reverted to user-preferred logic)
    dna_segments = []
    for central_res in central_bp_geometric:
        chain_id = central_res.get_parent().id
        center_res_num = central_res.id[1]
        start_res_num = center_res_num - 73
        end_res_num = center_res_num + 73
        dna_segments.append({'chain': chain_id, 'start': start_res_num, 'end': end_res_num})
        
    central_bp_info = [
        {'chain': res.get_parent().id, 'resid': res.id[1]} for res in central_bp_geometric
    ]

    return central_bp_info, dna_segments, actual_pairing_info

def generate_config_file(pdb_id, ncp_octamers, histone_map, structure, output_file):
    """Generates the final configuration file."""
    with open(output_file, 'w') as f:
        f.write(f"# NCP Configuration File for {pdb_id}\n")
        f.write("# Generated by ncp_identification.py\n")
        f.write("# Defines histone octamers, DNA segments, and the central base pair for symmetry axis calculation.\n\n")

        for i, octamer in enumerate(ncp_octamers):
            ncp_id = i + 1
            f.write(f"NCP: {ncp_id}\n")
            
            histone_ids = sorted([c.id for c in octamer])
            f.write(f"HISTONE_CHAINS: {' '.join(histone_ids)}\n")

            type_info = [f"{hid}({histone_map.get(hid, '?')})" for hid in histone_ids]
            f.write(f"# Histone types: {' '.join(type_info)}\n")

            # Add the machine-readable histone map
            map_items = [f"{hid}:{histone_map.get(hid, '?')}" for hid in histone_ids]
            f.write(f"HISTONE_MAP: {' '.join(map_items)}\n")
            
            central_bp, dna_segments, actual_pairing = find_central_bp_and_dna_segment(structure, octamer, histone_map)
            
            if central_bp and dna_segments:
                bp_str = ' '.join([f"{bp['chain']}:{bp['resid']}" for bp in central_bp])
                f.write(f"CENTRAL_BP: {bp_str}\n")

                # Add detailed comment about actual pairing
                res_a, (partner_a, dist_a, is_comp_a), res_b, (partner_b, dist_b, is_comp_b) = actual_pairing
                
                def get_partner_note(res, partner, dist, is_comp):
                    res_str = f"{res.get_parent().id}:{res.id[1]}"
                    if not partner:
                        return f"# Note: Actual partner of {res_str} is None"
                    
                    partner_str = f"{partner.get_parent().id}:{partner.id[1]}"
                    dist_str = f"{dist:.2f}"
                    
                    note = "complementary"
                    if not is_comp:
                        note = "geometric closest, MISMATCH"
                    elif is_comp and dist > 3.5:
                        note += ", distant"
                        
                    return f"# Note: Partner of {res_str} is {partner_str} (dist: {dist_str} Å, {note})"

                f.write(get_partner_note(res_a, partner_a, dist_a, is_comp_a) + '\n')
                f.write(get_partner_note(res_b, partner_b, dist_b, is_comp_b) + '\n')

                for seg in dna_segments:
                    f.write(f"DNA_SEGMENT: {seg['chain']} {seg['start']} {seg['end']}\n")
            else:
                f.write("# Could not identify DNA for this NCP.\n")
            
            f.write("\n")
    print(f"Configuration file generated: {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Identify NCPs and generate a configuration file for analysis.")
    parser.add_argument("pdb_file", help="Path to the input PDB file.")
    parser.add_argument("-f", "--histone-fasta", required=True, help="Path to FASTA file of canonical histone sequences (for BLAST db).")
    parser.add_argument("-o", "--output-config", default=None, help="Name of the output configuration file. Defaults to <pdb_filename>_config.txt")
    args = parser.parse_args()

    if not os.path.exists(args.pdb_file):
        print(f"Error: PDB file not found at {args.pdb_file}"); return
    if not os.path.exists(args.histone_fasta):
        print(f"Error: Histone FASTA file not found at {args.histone_fasta}"); return

    pdb_id = os.path.splitext(os.path.basename(args.pdb_file))[0]
    db_path = os.path.splitext(args.histone_fasta)[0]

    # Determine output config file name
    output_config = args.output_config if args.output_config else f"{pdb_id}_config.txt"

    # 1. Create BLAST database if it doesn't exist
    if not os.path.exists(f"{db_path}.phr"):
        print(f"Creating BLAST database from {args.histone_fasta}...")
        cmd_makedb = ["makeblastdb", "-in", args.histone_fasta, "-dbtype", "prot", "-out", db_path, "-parse_seqids"]
        try:
            subprocess.run(cmd_makedb, check=True, capture_output=True)
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            print(f"Error creating BLAST database: {e}")
            print("Please ensure NCBI BLAST+ is installed and in your system's PATH.")
            return

    # 2. Load Structure and Identify Histones
    print(f"Loading structure from {args.pdb_file}...")
    parser = PDBParser()
    structure = parser.get_structure(pdb_id, args.pdb_file)
    
    print("Extracting chain sequences...")
    sequences = get_chain_sequences(structure)
    
    print("Identifying histones via BLAST...")
    histone_map = identify_histones_by_blast(sequences, db_path)
    if not histone_map:
        print("Could not identify any histones. Aborting."); return
    print(f"Identified {len(histone_map)} histone chains.")

    # 3. Cluster NCPs
    print("Clustering chains into histone octamers...")
    ncp_octamers = cluster_chains_into_ncps(structure, histone_map)
    if not ncp_octamers:
        print("Could not cluster any complete histone octamers. Aborting."); return
    print(f"Found {len(ncp_octamers)} NCP(s).")

    # 4. Generate Configuration File
    print("Identifying central base pairs and DNA segments...")
    generate_config_file(pdb_id, ncp_octamers, histone_map, structure, output_config)
    print("\nIdentification complete.")


if __name__ == "__main__":
    main()
