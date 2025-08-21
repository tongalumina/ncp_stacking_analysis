#!/usr/bin/env python

import sys
from Bio.PDB import PDBParser, PDBIO

# Check if a PDB file name is provided as a command-line argument
if len(sys.argv) != 2:
    print("Usage: python script.py input.pdb")
    sys.exit(1)

# Get the input PDB file name from command line
input_pdb = sys.argv[1]

# Generate output file name with '_merged' suffix
output_pdb = input_pdb.replace('.pdb', '_merged.pdb')

# Load the PDB file
parser = PDBParser(QUIET=True)
structure = parser.get_structure("structure", input_pdb)

# Create a new structure for the merged model
new_structure = structure.copy()
new_structure[0].id = 0  # Ensure only one model

# Define chain mapping for MODEL 2 (adjust as needed for your PDB)
chain_mapping = {
    "A": "Q", "B": "R", "C": "S", "D": "T",
    "E": "U", "F": "V", "G": "W", "H": "X",
    "a": "q", "b": "r", "c": "s", "d": "t",
    "e": "u", "f": "v", "g": "w", "h": "x",
    "I": "Y", "J": "Z"
}

# Rename chains in MODEL 2 (assuming MODEL 2 is structure[1])
for model in structure:
    if model.id == 1:  # MODEL 2
        for chain in model:
            if chain.id in chain_mapping:
                chain.id = chain_mapping[chain.id]
        # Move chains from MODEL 2 to MODEL 1
        for chain in model:
            new_structure[0].add(chain.copy())

# Remove MODEL 2
if 1 in new_structure:
    del new_structure[1]

# Save the new PDB file
io = PDBIO()
io.set_structure(new_structure)
io.save(output_pdb)
print(f"Merged PDB saved as {output_pdb}")
