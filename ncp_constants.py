"""
ncp_constants.py

This file contains shared constants used by the NCP analysis scripts.
"""

# --- DNA and Interaction Constants ---
DNA_RESIDUES = ['DA', 'DC', 'DG', 'DT']
"""A list of standard DNA residue names."""

INTERACTION_DISTANCE = 10.0
"""The distance (in Angstroms) to define histone-DNA interaction."""

NCP_DNA_LENGTH = 147
"""The canonical length of DNA (in base pairs) wrapped in a nucleosome core particle."""

# --- Histone Core Definitions ---
HISTONE_CORE_RANGES = {
    'H2A': (21, 116),
    'H2B': (34, 121),
    'H3': (44, 135),
    'H4': (24, 102),
}
"""Defines the structured core region for each histone type based on residue IDs,
used to calculate a more stable center of mass by excluding flexible tails."""
