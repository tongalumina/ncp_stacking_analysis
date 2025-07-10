#!/usr/bin/env python

"""
ncp_analysis.py

This script performs a systematic analysis of nucleosome core particle (NCP)
stacking geometry, based on the methods described in the paper:
"A systematic analysis of nucleosome core particle and nucleosome-nucleosome
stacking structure" by Korolev, Lyubartsev, & Nordenskiöld (2018).
Scientific Reports, 8(1), 1543. DOI: 10.1038/s41598-018-19875-0.

This script takes a PDB ID as input, downloads the structure, identifies histone
and DNA chains using BLAST, clusters them into NCPs, and performs a pairwise
geometric analysis on all adjacent NCPs in the stack.

It also generates a PyMOL script (.pml) to visualize the NCPs, their coordinate
systems, and the calculated geometric parameters for each pair.

Requires:
- BioPython
- NumPy
- NCBI BLAST+ (blastp and makeblastdb must be in the system's PATH)
"""

import argparse
import os
import subprocess
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Data.IUPACData import protein_letters_3to1, protein_letters_3to1_extended
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import warnings
import urllib.request
from collections import defaultdict

warnings.filterwarnings("ignore", category=PDBConstructionWarning)

# --- Constants ---
DNA_RESIDUES = ['DA', 'DC', 'DG', 'DT']

# --- Main Functions ---

def fetch_pdb(pdb_id):
    """Downloads a PDB file from the RCSB PDB database."""
    # Prioritize the biological assembly file (.pdb1)
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb1"
    file_path = f"{pdb_id}.pdb1"
    if not os.path.exists(file_path):
        print(f"Downloading biological assembly {pdb_id}.pdb1...")
        try:
            urllib.request.urlretrieve(url, file_path)
        except urllib.error.HTTPError:
            print(f"Biological assembly not found. Trying standard PDB {pdb_id}.pdb...")
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            file_path = f"{pdb_id}.pdb"
            if not os.path.exists(file_path):
                 urllib.request.urlretrieve(url, file_path)
    print("Download complete.")
    return file_path

def get_chain_sequences(structure):
    """Extracts protein sequences from a PDB structure."""
    sequences = {}
    for model in structure:
        for chain in model:
            seq = ""
            for residue in chain:
                # Skip DNA residues and non-amino acid residues
                if residue.get_resname() in DNA_RESIDUES or not residue.has_id('CA'): # Check for Alpha Carbon for amino acids
                    continue
                
                res_name = residue.get_resname().strip().title()
                one_letter_code = protein_letters_3to1_extended.get(res_name, 'X')
                if one_letter_code == 'X':
                    print(f"Warning: Unknown residue type '{res_name}' in chain {chain.id}. Skipping.")
                    continue
                seq += one_letter_code
            if seq:
                sequences[chain.id] = SeqRecord(Seq(seq), id=chain.id)
                print(f"Extracted sequence for chain {chain.id}: {seq[:50]}...") # Print first 50 chars
    return sequences

def identify_histones_by_blast(sequences, db_path):
    """Identifies histone type for each chain using BLASTp."""
    histone_map = {}
    for chain_id, seq_record in sequences.items():
        # Write sequence to a temporary fasta file
        temp_fasta = f"temp_{chain_id}.fasta"
        SeqIO.write(seq_record, temp_fasta, "fasta")

        # Run blastp
        cmd = [
            "blastp", "-query", temp_fasta, "-db", db_path,
            "-outfmt", "6 qseqid sseqid pident", "-max_target_seqs", "1"
        ]
        print(f"Running BLASTp command: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        os.remove(temp_fasta)

        if result.returncode != 0:
            print(f"BLASTp error for chain {chain_id}: {result.stderr}")
            continue

        if not result.stdout:
            print(f"No BLASTp hits for chain {chain_id}.")
            continue

        print(f"BLASTp output for chain {chain_id}:\n{result.stdout.strip()}")

        parts = result.stdout.strip().split()
        if len(parts) < 3:
            print(f"Unexpected BLASTp output format for chain {chain_id}: {result.stdout.strip()}")
            continue

        subject_id = parts[1]
        pident = float(parts[2])
        
        # Extract histone type directly from subject_id (e.g., H3 from H3)
        histone_type = subject_id
        
        # Normalize histone types (e.g., H3.1 -> H3, H2B1A -> H2B)
        if histone_type.startswith("H3"): histone_type = "H3"
        elif histone_type.startswith("H4"): histone_type = "H4"
        elif histone_type.startswith("H2A"): histone_type = "H2A"
        elif histone_type.startswith("H2B"): histone_type = "H2B"

        if pident > 50: # Confidence threshold
            histone_map[chain_id] = histone_type
            print(f"Identified Chain {chain_id} as {histone_type} (pident: {pident:.2f}%)")
        else:
            print(f"Chain {chain_id} identified as {subject_id} but percent identity ({pident:.2f}%) is below threshold (50%).")
    return histone_map

def get_center_of_mass(atom_list):
    """Calculates the center of mass of a list of atoms."""
    coords = [atom.get_coord() for atom in atom_list]
    return np.mean(coords, axis=0)

def get_plane_normal(atom_list):
    """
    Calculates the normal to the best-fit plane for a list of atoms using SVD.
    """
    coords = np.array([atom.get_coord() for atom in atom_list])
    centroid = np.mean(coords, axis=0)
    coords_centered = coords - centroid
    _, _, vh = np.linalg.svd(coords_centered)
    normal = vh[2, :]
    return normal / np.linalg.norm(normal)

def cluster_chains_into_ncps(structure, histone_map):
    """Clusters histone chains into potential NCP histone octamers based on proximity and stoichiometry."""
    all_chains = {chain.id: chain for model in structure for chain in model}
    
    histone_chains_by_type = defaultdict(list)
    for chain_id, chain_obj in all_chains.items():
        if chain_id in histone_map:
            histone_chains_by_type[histone_map[chain_id]].append(chain_obj)

    ncp_histone_octamers = []
    assigned_histone_chains_ids = set()

    # Iterate until no more histone octamers can be formed
    while True:
        found_octamer_in_iteration = False
        
        # Try to find a starting H3 chain that hasn't been assigned yet
        for h3_chain in histone_chains_by_type["H3"]:
            if h3_chain.id in assigned_histone_chains_ids:
                continue

            potential_octamer_chains = {h3_chain.id: h3_chain}
            
            # Greedily add other unassigned histones close to this H3
            h3_com = get_center_of_mass(list(h3_chain.get_atoms()))
            
            candidate_histones = []
            for h_type, chains in histone_chains_by_type.items():
                for chain_obj in chains:
                    if chain_obj.id not in assigned_histone_chains_ids and chain_obj.id != h3_chain.id:
                        candidate_histones.append((chain_obj, get_center_of_mass(list(chain_obj.get_atoms()))))
            
            # Sort candidate histones by distance to the starting H3
            candidate_histones.sort(key=lambda x: np.linalg.norm(h3_com - x[1]))

            for other_h_chain, other_h_com in candidate_histones:
                if len(potential_octamer_chains) == 8: # Octamer formed
                    break
                potential_octamer_chains[other_h_chain.id] = other_h_chain

            # Check if a full octamer (8 histones) is formed
            if len(potential_octamer_chains) != 8:
                continue # Not a full octamer, try next H3

            # Verify stoichiometry (2x H3, 2x H4, 2x H2A, 2x H2B)
            stoichiometry = defaultdict(int)
            for h_chain_id in potential_octamer_chains:
                stoichiometry[histone_map[h_chain_id]] += 1
            
            if not (stoichiometry["H3"] == 2 and stoichiometry["H4"] == 2 and 
                    stoichiometry["H2A"] == 2 and stoichiometry["H2B"] == 2):
                continue # Incorrect stoichiometry, try next H3

            # If we reached here, we found a valid histone octamer
            ncp_histone_octamers.append(list(potential_octamer_chains.values()))
            
            # Mark chains as assigned
            for h_chain_id in potential_octamer_chains:
                assigned_histone_chains_ids.add(h_chain_id)
            
            found_octamer_in_iteration = True
            break # Found an octamer, break to restart search for next octamer
        
        if not found_octamer_in_iteration:
            break # No more octamers found in this iteration

    # Sort octamers along the z-axis (assuming a typical stacking orientation)
    ncp_histone_octamers.sort(key=lambda octamer: get_center_of_mass([atom for chain in octamer for atom in chain.get_atoms()])[2])
    return ncp_histone_octamers

def calculate_ncp_basis(histone_octamer, all_dna_chains):
    """Calculates the geometric basis for a single NCP, associating DNA residues."""
    gho_atoms = []
    for histone_chain in histone_octamer:
        for res in list(histone_chain.get_residues())[15:-15]: # Approximate globular part
            gho_atoms.extend(res.get_atoms())

    if not gho_atoms:
        print(f"Debug: gho_atoms is empty for NCP. Histone chains: {[c.id for c in histone_octamer]}")
        return None

    gho_com = get_center_of_mass(gho_atoms)

    # Find DNA residues associated with this histone octamer
    dna_residues_for_ncp = []
    for dna_chain in all_dna_chains:
        for residue in dna_chain:
            # Only consider DNA residues (DA, DC, DG, DT)
            if residue.get_resname() in DNA_RESIDUES:
                # Calculate distance from residue to octamer COM
                res_com = get_center_of_mass(list(residue.get_atoms()))
                if np.linalg.norm(gho_com - res_com) < 70: # Cutoff for DNA residue association
                    dna_residues_for_ncp.append(residue)
    
    if not dna_residues_for_ncp:
        print(f"Debug: No DNA residues found for this NCP. Histone octamer COM: {gho_com}")
        return None

    # Extract phosphates and central DNA atoms from associated DNA residues
    dna_phosphates = []
    dna_atoms_for_axis = []
    for residue in dna_residues_for_ncp:
        if 'P' in residue:
            dna_phosphates.append(residue['P'])
        # Get central DNA atoms for symmetry axis calculation (adjust residue ID range as needed)
        # This part might need further refinement based on actual PDB numbering
        if 70 < residue.id[1] < 78: # Example range, adjust based on typical wrapped DNA length
            dna_atoms_for_axis.extend(residue.get_atoms())

    if not dna_phosphates:
        print(f"Debug: dna_phosphates is empty for NCP. DNA residues found: {len(dna_residues_for_ncp)}")
        return None
    if not dna_atoms_for_axis:
        print(f"Debug: dna_atoms_for_axis is empty for NCP. DNA residues found: {len(dna_residues_for_ncp)}")
        return None

    plane_normal = get_plane_normal(dna_phosphates)
    central_bp_com = get_center_of_mass(dna_atoms_for_axis)
    symmetry_axis = central_bp_com - gho_com
    symmetry_axis /= np.linalg.norm(symmetry_axis)

    if np.dot(symmetry_axis, plane_normal) < 0:
        plane_normal *= -1

    return {"com": gho_com, "axis": symmetry_axis, "normal": plane_normal, "dna_residues": dna_residues_for_ncp}

def calculate_stacking_parameters(basis1, basis2):
    """Calculates the 7 NCP-NCP stacking parameters."""
    # This function remains the same as before
    params = {}
    c1, ax1, n1 = basis1["com"], basis1["axis"], basis1["normal"]
    c2, ax2, n2 = basis2["com"], basis2["axis"], basis2["normal"]
    c1_c2_vec = c2 - c1
    params["Distance"] = np.linalg.norm(c1_c2_vec)
    params["Rise"] = np.dot(c1_c2_vec, n1)
    proj_c2_on_plane1 = c2 - params["Rise"] * n1
    shift_vec = proj_c2_on_plane1 - c1
    params["Shift"] = np.linalg.norm(shift_vec)
    if params["Shift"] > 1e-6:
        ref_vec = ax1 - np.dot(ax1, n1) * n1
        ref_vec /= np.linalg.norm(ref_vec)
        cos_phi = np.dot(shift_vec, ref_vec) / params["Shift"]
        sign = np.sign(np.dot(np.cross(ref_vec, shift_vec), n1))
        params["Shift Orientation"] = np.rad2deg(np.arccos(np.clip(cos_phi, -1.0, 1.0))) * sign
    else:
        params["Shift Orientation"] = 0.0
    proj_ax2_on_plane1 = ax2 - np.dot(ax2, n1) * n1
    proj_ax2_norm = np.linalg.norm(proj_ax2_on_plane1)
    if proj_ax2_norm > 1e-6:
        proj_ax2_on_plane1 /= proj_ax2_norm
        ref_vec = ax1 - np.dot(ax1, n1) * n1
        ref_vec /= np.linalg.norm(ref_vec)
        cos_delta = np.dot(proj_ax2_on_plane1, ref_vec)
        sign = np.sign(np.dot(np.cross(ref_vec, proj_ax2_on_plane1), n1))
        params["Symmetry Axes Orientation"] = np.rad2deg(np.arccos(np.clip(cos_delta, -1.0, 1.0))) * sign
    else:
        params["Symmetry Axes Orientation"] = 0.0
    cos_tilt = np.dot(n1, n2)
    params["Tilt"] = np.rad2deg(np.arccos(np.clip(cos_tilt, -1.0, 1.0)))
    proj_n2_on_plane1 = n2 - np.dot(n2, n1) * n1
    proj_n2_norm = np.linalg.norm(proj_n2_on_plane1)
    if proj_n2_norm > 1e-6:
        proj_n2_on_plane1 /= proj_n2_norm
        ref_vec = ax1 - np.dot(ax1, n1) * n1
        ref_vec /= np.linalg.norm(ref_vec)
        cos_tilt_dir = np.dot(proj_n2_on_plane1, ref_vec)
        sign = np.sign(np.dot(np.cross(ref_vec, proj_n2_on_plane1), n1))
        params["Tilt Direction"] = np.rad2deg(np.arccos(np.clip(cos_tilt_dir, -1.0, 1.0))) * sign
    else:
        params["Tilt Direction"] = 0.0
    return params

def generate_pymol_script_polynucleosome(pdb_id, ncps, all_params):
    """Generates a PyMOL script for visualizing the entire nucleosome stack."""
    pml_file = f"visualize_{pdb_id}_stack.pml"
    
    script_parts = [
        f"load {pdb_id}_stack.pdb",
        "bg_color white",
        "hide everything",
        "show cartoon",
        "color gray80, all"
    ]

    colors = ["palecyan", "lightorange", "lightmagenta", "lightgreen"]
    for i, ncp in enumerate(ncps):
        ncp_name = f"ncp{i+1}"
        # Get chain IDs for histones
        histone_chain_ids = [c.id for c in ncp["histones"]]
        # Get chain IDs for DNA residues (from the original chains)
        dna_chain_ids = list(set([res.get_parent().id for res in ncp["dna_residues"]]))
        
        all_ncp_chain_ids = histone_chain_ids + dna_chain_ids
        script_parts.append(f"select {ncp_name}, chain {'+'.join(all_ncp_chain_ids)}")
        script_parts.append(f"color {colors[i % len(colors)]}, {ncp_name}")

    for i, params in enumerate(all_params):
        pair_name = f"pair_{i+1}_{i+2}"
        script_parts.append(f"print('\n--- Analysis for NCP{i+1}-NCP{i+2} ---')")
        for key, value in params.items():
            unit = "A" if key in ["Distance", "Rise", "Shift"] else "deg"
            script_parts.append(f"print('{key:<25}: {value:8.2f} {unit}')")

    script_parts.append("zoom center")
    
    with open(pml_file, "w") as f:
        f.write("\n".join(script_parts))
    print(f"\nGenerated PyMOL script: {pml_file}")

class NCPSelect(Select):
    def __init__(self, chains_to_keep):
        self.chains_to_keep = chains_to_keep
    def accept_chain(self, chain):
        return chain.id in self.chains_to_keep

def main():
    parser = argparse.ArgumentParser(description="Analyze polynucleosome stacking from a PDB file.")
    parser.add_argument("pdb_id", help="The PDB ID of the structure (e.g., 8Y3C).")
    parser.add_argument("histone_db", help="Path to the BLAST database of histone sequences.")
    args = parser.parse_args()

    # 1. Fetch PDB
    pdb_file = fetch_pdb(args.pdb_id.upper())
    if not pdb_file:
        return

    # 2. Parse Structure and Identify Histones
    p = PDBParser()
    structure = p.get_structure(args.pdb_id, pdb_file)
    sequences = get_chain_sequences(structure)
    histone_map = identify_histones_by_blast(sequences, args.histone_db)

    if not histone_map:
        print("Error: Could not identify any histone chains. Aborting.")
        return

    # 3. Cluster chains into NCPs (histone octamers only at this stage)
    print("\nClustering histone octamers...")
    ncp_histone_octamers = cluster_chains_into_ncps(structure, histone_map)
    print(f"Found {len(ncp_histone_octamers)} potential histone octamers.")

    if len(ncp_histone_octamers) < 2:
        print("Error: Less than two histone octamers were found. Cannot perform stacking analysis.")
        return

    # 4. Associate DNA and calculate basis for each NCP
    ncps = [] # This will store the full NCP objects (histones + DNA)
    all_dna_chains = [chain for model in structure for chain in model if any(res.get_resname() in DNA_RESIDUES for res in chain)]

    for i, octamer in enumerate(ncp_histone_octamers):
        print(f"Calculating basis for NCP{i+1} (associating DNA)...")
        # Pass all_dna_chains to calculate_ncp_basis for DNA residue association
        basis = calculate_ncp_basis(octamer, all_dna_chains)
        if basis:
            # Add the full NCP object (histones + associated DNA residues) to the ncps list
            ncps.append({"histones": octamer, "dna_residues": basis["dna_residues"], "com": basis["com"], "axis": basis["axis"], "normal": basis["normal"]})
        else:
            print(f"Warning: Could not calculate basis for NCP{i+1}. Skipping.")

    if len(ncps) < 2:
        print("Error: Less than two complete NCPs (histones + DNA) were found. Cannot perform stacking analysis.")
        return

    # 5. Perform pairwise analysis
    all_params = []
    ncp_bases = [{
        "com": ncp["com"],
        "axis": ncp["axis"],
        "normal": ncp["normal"]
    } for ncp in ncps] # Extract basis from full NCP objects

    for i in range(len(ncp_bases) - 1):
        print(f"\n--- Analyzing stack between NCP{i+1} and NCP{i+2} ---")
        params = calculate_stacking_parameters(ncp_bases[i], ncp_bases[i+1])
        all_params.append(params)
        for key, value in params.items():
            unit = "A" if key in ["Distance", "Rise", "Shift"] else "deg"
            print(f"{key:<25}: {value:8.2f} {unit}")

    # 6. Save a clean PDB of the stack and generate PyMOL script
    # The NCPSelect class needs to be updated to handle DNA residues, not just chains
    class NCPSelectForPDB(Select):
        def __init__(self, residues_to_keep):
            self.residues_to_keep = residues_to_keep
            self.residue_ids = set([(res.get_parent().id, res.get_id()) for res in residues_to_keep])

        def accept_residue(self, residue):
            return (residue.get_parent().id, residue.get_id()) in self.residue_ids

    all_ncp_residues = []
    for ncp in ncps:
        for h_chain in ncp["histones"]:
            all_ncp_residues.extend(list(h_chain.get_residues()))
        all_ncp_residues.extend(ncp["dna_residues"])

    io = PDBIO()
    io.set_structure(structure)
    stack_pdb_file = f"{args.pdb_id}_stack.pdb"
    io.save(stack_pdb_file, NCPSelectForPDB(all_ncp_residues))
    print(f"\nSaved clustered NCP stack to {stack_pdb_file}")

    generate_pymol_script_polynucleosome(args.pdb_id, ncps, all_params)

if __name__ == "__main__":
    main()

def calculate_ncp_basis(histone_octamer, all_dna_chains):
    """Calculates the geometric basis for a single NCP, associating DNA residues."""
    gho_atoms = []
    for histone_chain in histone_octamer:
        for res in list(histone_chain.get_residues())[15:-15]: # Approximate globular part
            gho_atoms.extend(res.get_atoms())

    if not gho_atoms:
        print(f"Debug: gho_atoms is empty for NCP. Histone chains: {[c.id for c in histone_octamer]}")
        return None

    gho_com = get_center_of_mass(gho_atoms)

    # Find DNA residues associated with this histone octamer
    dna_residues_for_ncp = []
    for dna_chain in all_dna_chains:
        for residue in dna_chain:
            # Only consider DNA residues (DA, DC, DG, DT)
            if residue.get_resname() in DNA_RESIDUES:
                # Calculate distance from residue to octamer COM
                res_com = get_center_of_mass(list(residue.get_atoms()))
                if np.linalg.norm(gho_com - res_com) < 70: # Cutoff for DNA residue association
                    dna_residues_for_ncp.append(residue)
    
    if not dna_residues_for_ncp:
        print(f"Debug: No DNA residues found for this NCP. Histone octamer COM: {gho_com}")
        return None

    # Extract phosphates and central DNA atoms from associated DNA residues
    dna_phosphates = []
    dna_atoms_for_axis = []
    for residue in dna_residues_for_ncp:
        if 'P' in residue:
            dna_phosphates.append(residue['P'])
        # Get central DNA atoms for symmetry axis calculation (adjust residue ID range as needed)
        # This part might need further refinement based on actual PDB numbering
        if 70 < residue.id[1] < 78: # Example range, adjust based on typical wrapped DNA length
            dna_atoms_for_axis.extend(residue.get_atoms())

    if not dna_phosphates:
        print(f"Debug: dna_phosphates is empty for NCP. DNA residues found: {len(dna_residues_for_ncp)}")
        return None
    if not dna_atoms_for_axis:
        print(f"Debug: dna_atoms_for_axis is empty for NCP. DNA residues found: {len(dna_residues_for_ncp)}")
        return None

    plane_normal = get_plane_normal(dna_phosphates)
    central_bp_com = get_center_of_mass(dna_atoms_for_axis)
    symmetry_axis = central_bp_com - gho_com
    symmetry_axis /= np.linalg.norm(symmetry_axis)

    if np.dot(symmetry_axis, plane_normal) < 0:
        plane_normal *= -1

    return {"com": gho_com, "axis": symmetry_axis, "normal": plane_normal, "dna_residues": dna_residues_for_ncp}

def calculate_stacking_parameters(basis1, basis2):
    """Calculates the 7 NCP-NCP stacking parameters."""
    # This function remains the same as before
    params = {}
    c1, ax1, n1 = basis1["com"], basis1["axis"], basis1["normal"]
    c2, ax2, n2 = basis2["com"], basis2["axis"], basis2["normal"]
    c1_c2_vec = c2 - c1
    params["Distance"] = np.linalg.norm(c1_c2_vec)
    params["Rise"] = np.dot(c1_c2_vec, n1)
    proj_c2_on_plane1 = c2 - params["Rise"] * n1
    shift_vec = proj_c2_on_plane1 - c1
    params["Shift"] = np.linalg.norm(shift_vec)
    if params["Shift"] > 1e-6:
        ref_vec = ax1 - np.dot(ax1, n1) * n1
        ref_vec /= np.linalg.norm(ref_vec)
        cos_phi = np.dot(shift_vec, ref_vec) / params["Shift"]
        sign = np.sign(np.dot(np.cross(ref_vec, shift_vec), n1))
        params["Shift Orientation"] = np.rad2deg(np.arccos(np.clip(cos_phi, -1.0, 1.0))) * sign
    else:
        params["Shift Orientation"] = 0.0
    proj_ax2_on_plane1 = ax2 - np.dot(ax2, n1) * n1
    proj_ax2_norm = np.linalg.norm(proj_ax2_on_plane1)
    if proj_ax2_norm > 1e-6:
        proj_ax2_on_plane1 /= proj_ax2_norm
        ref_vec = ax1 - np.dot(ax1, n1) * n1
        ref_vec /= np.linalg.norm(ref_vec)
        cos_delta = np.dot(proj_ax2_on_plane1, ref_vec)
        sign = np.sign(np.dot(np.cross(ref_vec, proj_ax2_on_plane1), n1))
        params["Symmetry Axes Orientation"] = np.rad2deg(np.arccos(np.clip(cos_delta, -1.0, 1.0))) * sign
    else:
        params["Symmetry Axes Orientation"] = 0.0
    cos_tilt = np.dot(n1, n2)
    params["Tilt"] = np.rad2deg(np.arccos(np.clip(cos_tilt, -1.0, 1.0)))
    proj_n2_on_plane1 = n2 - np.dot(n2, n1) * n1
    proj_n2_norm = np.linalg.norm(proj_n2_on_plane1)
    if proj_n2_norm > 1e-6:
        proj_n2_on_plane1 /= proj_n2_norm
        ref_vec = ax1 - np.dot(ax1, n1) * n1
        ref_vec /= np.linalg.norm(ref_vec)
        cos_tilt_dir = np.dot(proj_n2_on_plane1, ref_vec)
        sign = np.sign(np.dot(np.cross(ref_vec, proj_n2_on_plane1), n1))
        params["Tilt Direction"] = np.rad2deg(np.arccos(np.clip(cos_tilt_dir, -1.0, 1.0))) * sign
    else:
        params["Tilt Direction"] = 0.0
    return params

def generate_pymol_script_polynucleosome(pdb_id, ncps, all_params):
    """Generates a PyMOL script for visualizing the entire nucleosome stack."""
    pml_file = f"visualize_{pdb_id}_stack.pml"
    
    script_parts = [
        f"load {pdb_id}_stack.pdb",
        "bg_color white",
        "hide everything",
        "show cartoon",
        "color gray80, all"
    ]

    colors = ["palecyan", "lightorange", "lightmagenta", "lightgreen"]
    for i, ncp in enumerate(ncps):
        ncp_name = f"ncp{i+1}"
        # Get chain IDs for histones
        histone_chain_ids = [c.id for c in ncp["histones"]]
        # Get chain IDs for DNA residues (from the original chains)
        dna_chain_ids = list(set([res.get_parent().id for res in ncp["dna_residues"]]))
        
        all_ncp_chain_ids = histone_chain_ids + dna_chain_ids
        script_parts.append(f"select {ncp_name}, chain {'+'.join(all_ncp_chain_ids)}")
        script_parts.append(f"color {colors[i % len(colors)]}, {ncp_name}")

    for i, params in enumerate(all_params):
        pair_name = f"pair_{i+1}_{i+2}"
        script_parts.append(f"print('\n--- Analysis for NCP{i+1}-NCP{i+2} ---')")
        for key, value in params.items():
            unit = "A" if key in ["Distance", "Rise", "Shift"] else "deg"
            script_parts.append(f"print('{key:<25}: {value:8.2f} {unit}')")

    script_parts.append("zoom center")
    
    with open(pml_file, "w") as f:
        f.write("\n".join(script_parts))
    print(f"\nGenerated PyMOL script: {pml_file}")


class NCPSelect(Select):
    def __init__(self, chains_to_keep):
        self.chains_to_keep = chains_to_keep
    def accept_chain(self, chain):
        return chain.id in self.chains_to_keep

def main():
    parser = argparse.ArgumentParser(description="Analyze polynucleosome stacking from a PDB file.")
    parser.add_argument("pdb_id", help="The PDB ID of the structure (e.g., 8Y3C).")
    parser.add_argument("histone_db", help="Path to the BLAST database of histone sequences.")
    args = parser.parse_args()

    # 1. Fetch PDB
    pdb_file = fetch_pdb(args.pdb_id.upper())
    if not pdb_file:
        return

    # 2. Parse Structure and Identify Histones
    p = PDBParser()
    structure = p.get_structure(args.pdb_id, pdb_file)
    sequences = get_chain_sequences(structure)
    histone_map = identify_histones_by_blast(sequences, args.histone_db)

    if not histone_map:
        print("Error: Could not identify any histone chains. Aborting.")
        return

    # 3. Cluster chains into NCPs (histone octamers only at this stage)
    print("\nClustering histone octamers...")
    ncp_histone_octamers = cluster_chains_into_ncps(structure, histone_map)
    print(f"Found {len(ncp_histone_octamers)} potential histone octamers.")

    if len(ncp_histone_octamers) < 2:
        print("Error: Less than two histone octamers were found. Cannot perform stacking analysis.")
        return

    # 4. Associate DNA and calculate basis for each NCP
    ncps = [] # This will store the full NCP objects (histones + DNA)
    all_dna_chains = [chain for model in structure for chain in model if any(res.get_resname() in DNA_RESIDUES for res in chain)]

    for i, octamer in enumerate(ncp_histone_octamers):
        print(f"Calculating basis for NCP{i+1} (associating DNA)...")
        # Pass all_dna_chains to calculate_ncp_basis for DNA residue association
        basis = calculate_ncp_basis(octamer, all_dna_chains)
        if basis:
            # Add the full NCP object (histones + associated DNA residues) to the ncps list
            ncps.append({"histones": octamer, "dna_residues": basis["dna_residues"], "com": basis["com"], "axis": basis["axis"], "normal": basis["normal"]})
        else:
            print(f"Warning: Could not calculate basis for NCP{i+1}. Skipping.")

    if len(ncps) < 2:
        print("Error: Less than two complete NCPs (histones + DNA) were found. Cannot perform stacking analysis.")
        return

    # 5. Perform pairwise analysis
    all_params = []
    ncp_bases = [{
        "com": ncp["com"],
        "axis": ncp["axis"],
        "normal": ncp["normal"]
    } for ncp in ncps] # Extract basis from full NCP objects

    for i in range(len(ncp_bases) - 1):
        print(f"\n--- Analyzing stack between NCP{i+1} and NCP{i+2} ---")
        params = calculate_stacking_parameters(ncp_bases[i], ncp_bases[i+1])
        all_params.append(params)
        for key, value in params.items():
            unit = "A" if key in ["Distance", "Rise", "Shift"] else "deg"
            print(f"{key:<25}: {value:8.2f} {unit}")

    # 6. Save a clean PDB of the stack and generate PyMOL script
    # The NCPSelect class needs to be updated to handle DNA residues, not just chains
    class NCPSelectForPDB(Select):
        def __init__(self, residues_to_keep):
            self.residues_to_keep = residues_to_keep
            self.residue_ids = set([(res.get_parent().id, res.get_id()) for res in residues_to_keep])

        def accept_residue(self, residue):
            return (residue.get_parent().id, residue.get_id()) in self.residue_ids

    all_ncp_residues = []
    for ncp in ncps:
        for h_chain in ncp["histones"]:
            all_ncp_residues.extend(list(h_chain.get_residues()))
        all_ncp_residues.extend(ncp["dna_residues"])

    io = PDBIO()
    io.set_structure(structure)
    stack_pdb_file = f"{args.pdb_id}_stack.pdb"
    io.save(stack_pdb_file, NCPSelectForPDB(all_ncp_residues))
    print(f"\nSaved clustered NCP stack to {stack_pdb_file}")

    generate_pymol_script_polynucleosome(args.pdb_id, ncps, all_params)

if __name__ == "__main__":
    main()