# NCP Stacking Analysis

This project provides a set of Python scripts for the analysis of Nucleosome Core Particle (NCP) stacking in PDB structures. It is a reimplementation of the methodology described in the paper by Korolev, Lyubartsev, & Nordenskiöld (2018).

**Reference:**
Korolev, N., Lyubartsev, A.P. & Nordenskiöld, L. A systematic analysis of nucleosome core particle and nucleosome-nucleosome stacking structure. *Sci Rep* **8**, 1543 (2018). https://doi.org/10.1038/s41598-018-19875-0

## Methodology

The analysis follows a two-step workflow:

1.  **Identification:** The first script, `ncp_identification.py`, processes a PDB file to identify all NCPs. It uses BLAST to determine the identity of histone chains, clusters them into octamers, and locates key DNA features (like the central dyad) for each NCP. The output is a configuration file (`*_config.txt`) that defines the components of each identified NCP.

2.  **Analysis:** The second script, `ncp_analysis.py`, takes the PDB file and the generated configuration file as input. It calculates the geometric stacking parameters (Distance, Rise, Shift, Tilt, etc.) between adjacent NCPs. The script outputs a detailed results file (`*_results.txt`) and a PyMOL visualization script (`*_visualization.pml`) to inspect the stacked structures and their coordinate systems.

## Scripts

-   `ncp_identification.py`: Identifies NCPs in a PDB file and generates a configuration file. It requires a local BLAST database of canonical histone sequences.
-   `ncp_analysis.py`: Calculates NCP stacking parameters using the PDB file and the configuration file from the identification step.
-   `ncp_constants.py`: A file containing shared constants (e.g., histone core residue ranges) used by both scripts.
-   `pdb_merge_2models.py`: A utility script to merge multiple models in a PDB file into a single model for analysis (e.g., for PDB ID 1ZBB).

## Requirements

-   Python 3.x
-   NCBI BLAST+ (must be available in the system's PATH)
-   Python packages as listed in `requirements.txt`. Install them using:
    ```bash
    pip install -r requirements.txt
    ```

## Usage

### Step 1: Prepare the BLAST database

If you haven't already, create a BLAST database from the provided `histones.fa` file.

```bash
makeblastdb -in histones.fa -dbtype prot
```

### Step 2: Run NCP Identification

Run the identification script on your PDB file.

**Command:**
```bash
python ncp_identification.py <PDB_FILE> -f <HISTONE_FASTA> [-o <OUTPUT_CONFIG_FILE>]
```

**Example:**
```bash
python ncp_identification.py 1zbb_ba.pdb -f histones.fa -o 1zbb_ba_config.txt
```
This will generate `1zbb_ba_config.txt` and a log file `1zbb_ba_id.log`.

### Step 3: Run NCP Stacking Analysis

Run the analysis script using the original PDB file and the newly generated config file.

**Command:**
```bash
python ncp_analysis.py <PDB_FILE> -c <CONFIG_FILE> [-o <OUTPUT_PREFIX>]
```

**Example:**
```bash
python ncp_analysis.py 1zbb_ba.pdb -c 1zbb_ba_config.txt -o 1zbb_ba
```
This will generate three files:
-   `1zbb_ba_results.txt`: The calculated stacking parameters.
-   `1zbb_ba_stack.pdb`: A cleaned PDB file containing only the identified NCPs.
-   `1zbb_ba_visualization.pml`: The PyMOL script for visualization.
