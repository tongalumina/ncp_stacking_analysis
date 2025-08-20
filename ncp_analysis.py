#!/usr/bin/env python
"""
ncp_analysis.py

Calculates NCP-NCP stacking parameters and generates a visualization script.
This is a revised version that incorporates more robust visualization methods.
"""

import argparse
import os
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select
from collections import defaultdict
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning

warnings.filterwarnings("ignore", category=PDBConstructionWarning)

# --- Class for PDB cleaning ---
class NCPSelect(Select):
    """Selects only the residues that are part of an NCP definition."""
    def __init__(self, ncp_configs, histone_map, core_ranges):
        self.residues_to_keep = set()
        for ncp_def in ncp_configs:
            # Add histone core residues
            for chain_id in ncp_def['histone_chains']:
                h_type = histone_map.get(chain_id)
                if h_type and h_type in core_ranges:
                    start, end = core_ranges[h_type]
                    for res_id in range(start, end + 1):
                        self.residues_to_keep.add((chain_id, res_id))
            
            # Add DNA segment residues
            for chain_id, start, end in ncp_def['dna_segments']:
                for res_id in range(start, end + 1):
                    self.residues_to_keep.add((chain_id, res_id))

    def accept_residue(self, residue):
        return (residue.get_parent().id, residue.id[1]) in self.residues_to_keep

# --- Constants ---
HISTONE_CORE_RANGES = {
    'H2A': (21, 116),
    'H2B': (34, 121),
    'H3': (44, 135),
    'H4': (24, 102),
}

# --- Configuration Parsing ---
def parse_config_file(config_file):
    """Parse the configuration file to extract NCP definitions."""
    ncps_config = []
    with open(config_file, 'r') as f:
        current_ncp = None
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or ':' not in line:
                continue
            
            key, value = line.split(':', 1)
            key, value = key.strip(), value.strip()

            if key == 'NCP':
                if current_ncp:
                    ncps_config.append(current_ncp)
                current_ncp = {'id': int(value), 'histone_chains': [], 'dna_segments': [], 'central_bp': [], 'histone_map': {}}
            elif key == 'HISTONE_CHAINS' and current_ncp:
                current_ncp['histone_chains'] = value.split()
            elif key == 'HISTONE_MAP' and current_ncp:
                current_ncp['histone_map'] = dict(pair.split(':') for pair in value.split())
            elif key == 'CENTRAL_BP' and current_ncp:
                for part in value.split():
                    chain_id, resid = part.split(':')
                    current_ncp['central_bp'].append((chain_id, int(resid)))
            elif key == 'DNA_SEGMENT' and current_ncp:
                parts = value.split()
                current_ncp['dna_segments'].append((parts[0], int(parts[1]), int(parts[2])))

    if current_ncp:
        ncps_config.append(current_ncp)
    return ncps_config

# --- Geometric Calculations ---
def get_center_of_mass(atom_list):
    """Calculates the center of mass for a list of atoms."""
    if not atom_list: return None
    coords = [atom.get_coord() for atom in atom_list]
    return np.mean(coords, axis=0)

def get_plane_normal(coords):
    """Calculates the normal vector to a plane defined by a set of coordinates."""
    centroid = np.mean(coords, axis=0)
    coords_centered = coords - centroid
    _, _, vh = np.linalg.svd(coords_centered)
    normal = vh[2, :]
    return normal / np.linalg.norm(normal)

def calculate_ncp_basis(ncp_id, histone_chains, dna_residues, central_bp_residues, histone_map):
    """Calculates the COM, symmetry axis, and plane normal for a single NCP."""
    # 1. Calculate Globular Histone Octamer (gHO) COM
    gho_atoms = []
    for chain in histone_chains:
        h_type = histone_map.get(chain.id)
        if not h_type or h_type not in HISTONE_CORE_RANGES: continue
        start, end = HISTONE_CORE_RANGES[h_type]
        for res in chain:
            if start <= res.id[1] <= end:
                gho_atoms.extend(a for a in res.get_atoms() if a.element != 'H')
    if not gho_atoms: return None
    gho_com = get_center_of_mass(gho_atoms)

    # 2. Calculate Symmetry Axis
    if not central_bp_residues: return None
    central_bp_com = get_center_of_mass([a for r in central_bp_residues for a in r.get_atoms() if a.element != 'H'])
    if gho_com is None or central_bp_com is None: return None
    
    symmetry_axis = central_bp_com - gho_com
    if np.linalg.norm(symmetry_axis) < 1e-6: return None
    symmetry_axis /= np.linalg.norm(symmetry_axis)

    # 3. Calculate NCP Plane Normal
    plane_defining_points = [a.get_coord() for r in dna_residues for a in r.get_atoms() if a.name == 'P']
    if len(plane_defining_points) < 3: return None
    plane_normal = get_plane_normal(plane_defining_points)

    # Ensure normal vector points in a direction consistent with the symmetry axis
    if np.dot(symmetry_axis, plane_normal) < 0:
        plane_normal *= -1

    return {"com": gho_com, "axis": symmetry_axis, "normal": plane_normal}

def calculate_stacking_parameters(basis1, basis2):
    """Calculates the 6 stacking parameters based on two NCP bases."""
    c1, ax1, n1 = basis1["com"], basis1["axis"], basis1["normal"]
    c2, ax2, n2 = basis2["com"], basis2["axis"], basis2["normal"]
    c1_c2_vec = c2 - c1
    params = {}
    params["Distance"] = np.linalg.norm(c1_c2_vec)
    params["Rise"] = np.dot(c1_c2_vec, n1)
    proj_c2_on_plane1 = c2 - params["Rise"] * n1
    shift_vec = proj_c2_on_plane1 - c1
    params["Shift"] = np.linalg.norm(shift_vec)
    
    # Define a reference vector in plane 1 for angle calculations
    ref_vec = ax1 - np.dot(ax1, n1) * n1
    if np.linalg.norm(ref_vec) < 1e-6: ref_vec = np.cross(ax1, n1)
    ref_vec /= np.linalg.norm(ref_vec)
    
    # Shift Orientation (phi)
    if params["Shift"] > 1e-6:
        unit_shift_vec = shift_vec / params["Shift"]
        cos_phi = np.clip(np.dot(unit_shift_vec, ref_vec), -1.0, 1.0)
        angle_phi = np.rad2deg(np.arccos(cos_phi))
        sign_phi = np.sign(np.dot(np.cross(ref_vec, unit_shift_vec), n1))
        params["Shift Orientation"] = angle_phi * sign_phi
    else: params["Shift Orientation"] = 0.0

    # Symmetry Axes Orientation (delta)
    proj_ax2_on_plane1 = ax2 - np.dot(ax2, n1) * n1
    proj_ax2_norm = np.linalg.norm(proj_ax2_on_plane1)
    if proj_ax2_norm > 1e-6:
        unit_proj_ax2 = proj_ax2_on_plane1 / proj_ax2_norm
        cos_delta = np.clip(np.dot(unit_proj_ax2, ref_vec), -1.0, 1.0)
        angle_delta = np.rad2deg(np.arccos(cos_delta))
        sign_delta = np.sign(np.dot(np.cross(ref_vec, unit_proj_ax2), n1))
        params["Symmetry Axes Orientation"] = angle_delta * sign_delta
    else: params["Symmetry Axes Orientation"] = 0.0

    # Tilt
    cos_tilt = np.clip(np.dot(n1, n2), -1.0, 1.0)
    params["Tilt"] = np.rad2deg(np.arccos(cos_tilt))

    # Tilt Direction
    proj_n2_on_plane1 = n2 - np.dot(n2, n1) * n1
    proj_n2_norm = np.linalg.norm(proj_n2_on_plane1)
    if proj_n2_norm > 1e-6:
        unit_proj_n2 = proj_n2_on_plane1 / proj_n2_norm
        cos_tilt_dir = np.clip(np.dot(unit_proj_n2, ref_vec), -1.0, 1.0)
        angle_tilt_dir = np.rad2deg(np.arccos(cos_tilt_dir))
        sign_tilt_dir = np.sign(np.dot(np.cross(ref_vec, unit_proj_n2), n1))
        params["Tilt Direction"] = angle_tilt_dir * sign_tilt_dir
    else: params["Tilt Direction"] = 0.0
    return params

# --- Visualization ---
def generate_pymol_script(output_prefix, ncp_configs, ncp_bases, all_params, histone_map):
    """Generates a PyMOL script for visualizing the analysis results."""
    pml_file = f"{output_prefix}_visualization.pml"
    stack_pdb_file = f"{output_prefix}_stack.pdb"
    
    def format_vec(v): return f"[{v[0]:.3f}, {v[1]:.3f}, {v[2]:.3f}]"

    cgo_objects = []
    colors = ["palecyan", "lightorange", "lightmagenta", "palegreen"]
    
    # --- CGO objects for Basis Vectors and Planes ---
    for i, ncp in enumerate(ncp_bases):
        ncp_id = ncp['id']
        com, axis, normal = ncp['com'], ncp['axis'], ncp['normal']
        color = colors[i % len(colors)]

        # Axis and Normal Arrows
        cgo_objects.extend([f"cgo_arrow {format_vec(com)}, {format_vec(com + axis * 30)}, radius=0.4, color=red, name=axis{ncp_id}",
                            f"cgo_arrow {format_vec(com)}, {format_vec(com + normal * 30)}, radius=0.4, color=blue, name=normal{ncp_id}"])

        # Plane (as a semi-transparent disk)
        cgo_objects.append("BEGIN, TRIANGLE_FAN")
        cgo_objects.append(f"COLOR, {color}")
        cgo_objects.append(f"ALPHA, 0.4") # Set transparency
        cgo_objects.append(f"VERTEX, {com[0]:.3f}, {com[1]:.3f}, {com[2]:.3f}")
        v_ref = np.cross(normal, axis) / np.linalg.norm(np.cross(normal, axis))
        for k in range(37):
            angle = k * 10 * np.pi / 180
            pt = com + 40 * (np.cos(angle) * v_ref + np.sin(angle) * np.cross(normal, v_ref))
            cgo_objects.append(f"VERTEX, {pt[0]:.3f}, {pt[1]:.3f}, {pt[2]:.3f}")
        cgo_objects.append("END")

    # --- CGO objects for Stacking Parameter Labels ---
    if all_params:
        params = all_params[0]
        ncp1, ncp2 = ncp_bases[0], ncp_bases[1]
        label_pos = (ncp1['com'] + ncp2['com']) / 2.0
        y_offset = 0
        for key, value in params.items():
            unit = " A" if key in ["Distance", "Rise", "Shift"] else " deg"
            label_text = f'{key}: {value:.1f}{unit}'
            pos = label_pos + np.array([25, y_offset, 0]) # Offset labels for clarity
            cgo_objects.append(f"pseudoatom label_{key.replace(' ','_')}, pos={format_vec(pos)}")
            cgo_objects.append(f'label label_{key.replace(' ','_')}, "{label_text}" ')
            y_offset -= 6

    # --- Construct final PyMOL script ---
    script = [f"load {stack_pdb_file}", "bg_color white", "hide everything", "show cartoon", "set cartoon_transparency, 0.5"]
    
    # Add selections for each NCP
    for i, ncp_def in enumerate(ncp_configs):
        ncp_name = f"ncp{ncp_def['id']}"
        h_selectors = [f"(chain {c} and resi {HISTONE_CORE_RANGES[histone_map[c]][0]}-{HISTONE_CORE_RANGES[histone_map[c]][1]})" for c in ncp_def['histone_chains']]
        d_selectors = [f"(chain {seg[0]} and resi {seg[1]}-{seg[2]})" for seg in ncp_def['dna_segments']]
        script.append(f"select {ncp_name}, {' or '.join(h_selectors + d_selectors)}")
        script.append(f"color {colors[i % len(colors)]}, {ncp_name}")
    
    # Load all CGOs at once
    script.append(f"load_cgo([{', '.join(cgo_objects)}], viz, 1)")
    script.append("set label_color, black, label_*")
    script.append("set label_size, -0.8")
    script.append("hide labels, label_*") # Hide pseudoatom dots
    script.append("zoom center")

    with open(pml_file, 'w') as f:
        f.write("\n".join(script))
    print(f"Generated PyMOL script: {pml_file}")

# --- Main Execution ---
def main():
    """Main function to run the analysis pipeline."""
    # 1. Argument Parsing
    parser = argparse.ArgumentParser(description="Calculate NCP stacking parameters and generate visualization.")
    parser.add_argument("pdb_file", help="Path to the input PDB file.")
    parser.add_argument("--config", required=True, help="Path to the NCP configuration file.")
    parser.add_argument("--output-prefix", default=None, help="Prefix for all output files. Defaults to PDB filename.")
    args = parser.parse_args()

    if not os.path.exists(args.pdb_file): return print(f"Error: PDB file not found at {args.pdb_file}")
    if not os.path.exists(args.config): return print(f"Error: Config file not found at {args.config}")

    pdb_id = os.path.splitext(os.path.basename(args.pdb_file))[0]
    output_prefix = args.output_prefix if args.output_prefix else pdb_id

    # 2. Configuration and Structure Loading
    print("Parsing configuration file...")
    ncps_config = parse_config_file(args.config)
    if not ncps_config: return print("Error: No NCPs found in config file.")

    print(f"Loading structure from {args.pdb_file}...")
    p = PDBParser()
    structure = p.get_structure(pdb_id, args.pdb_file)
    all_chains_map = {c.id: c for c in structure.get_chains()}

    histone_map = {k: v for ncp in ncps_config for k, v in ncp.get('histone_map', {}).items()}

    # 3. Per-NCP Basis Calculation
    ncp_bases = []
    for ncp_def in ncps_config:
        print(f"Processing NCP {ncp_def['id']}...")
        histone_chains = [all_chains_map[cid] for cid in ncp_def['histone_chains'] if cid in all_chains_map]
        dna_residues = [res for seg in ncp_def['dna_segments'] for res in all_chains_map.get(seg[0], []) if seg[1] <= res.id[1] <= seg[2]]
        central_bp_residues = [all_chains_map[bp[0]][bp[1]] for bp in ncp_def['central_bp'] if bp[0] in all_chains_map and bp[1] in all_chains_map[bp[0]]]

        if len(histone_chains) != 8 or not dna_residues or not central_bp_residues:
            print(f"Warning: Incomplete data for NCP {ncp_def['id']}. Skipping.")
            continue

        basis = calculate_ncp_basis(ncp_def['id'], histone_chains, dna_residues, central_bp_residues, histone_map)
        if basis: ncp_bases.append({'id': ncp_def['id'], **basis})

    if len(ncp_bases) < 2: return print("Analysis requires at least two valid NCPs. Aborting.")

    # 4. Stacking Parameter Calculation
    ncp_bases.sort(key=lambda x: x['id'])
    all_params = [calculate_stacking_parameters(ncp_bases[i], ncp_bases[i+1]) for i in range(len(ncp_bases) - 1)]

    # 5. Save Numerical Results
    results_file = f"{output_prefix}_results.txt"
    print(f"\n--- Saving Numerical Analysis to {results_file} ---")
    with open(results_file, 'w') as f:
        f.write(f"NCP Stacking Analysis Results for {pdb_id}\n")
        for i, params in enumerate(all_params):
            f.write(f"\n--- NCP {ncp_bases[i]['id']} - NCP {ncp_bases[i+1]['id']} ---\n")
            for key, value in params.items():
                unit = "A" if key in ["Distance", "Rise", "Shift"] else "deg"
                f.write(f"{key:<25}: {value:8.2f} {unit}\n")
    print(f"Analysis results saved to: {results_file}")

    # 6. Create Visualization Files
    print("\n--- Creating Visualization Files ---")
    print("Creating a cleaned PDB file for visualization...")
    io = PDBIO()
    io.set_structure(structure)
    stack_pdb_file = f"{output_prefix}_stack.pdb"
    io.save(stack_pdb_file, NCPSelect(ncps_config, histone_map, HISTONE_CORE_RANGES))
    print(f"Saved clustered NCP stack to {stack_pdb_file}")

    print("Generating PyMOL script...")
    generate_pymol_script(output_prefix, ncps_config, ncp_bases, all_params, histone_map)
    print("\nProcess complete.")

if __name__ == "__main__":
    main()