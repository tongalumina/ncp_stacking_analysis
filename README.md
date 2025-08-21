## NCP Stacking Analysis

This is a Python implementation of the nucleosome core particle (NCP) stacking
analysis using the NCP-centred coordinate system described in the following paper:

Korolev, N., Lyubartsev, A.P. & Nordenskiöld, L. A systematic analysis of
nucleosome core particle and nucleosome-nucleosome stacking structure. Sci Rep
8, 1543 (2018).

https://doi.org/10.1038/s41598-018-19875-0

The scripts were mainly generated using Gemini CLI 2.5 pro (CLI version
0.1.22), with extensive validation using PDB: 1ZBB.


Two main scripts are implemented:
- _ncp_identification.py_ will use `blastp` to find out the identities of the
  polypeptide chains, find the dyad base pair, set proper coordinate system,
  and generate a config file for inspection.
- _ncp_analysis.py_ will carry out all the calulation, the parameters are saved
  in a result file, it will also generate a PyMOL script for visualization.

The result file contains:
Distance, Rise, Shift, Tilt, Shift Orientation, Symmetry Axes Orientation, Tilt Direction

Additionally, a helper script is implemented:
- _pdb_merge_2models.py_ merges the two models in a biological assembly in 1ZBB into one model for analysis.

### Requirements
The scripts have been tested on Ubuntu/WSL2.

- `ncbi-blast+` package
- Python package requirement are defined in the _requirements.txt_ file

### Files
- 1zbb_ba.pdb, merged two models in PDB 1ZBB into one, remapping chain ID of model 2
- histones.fa, histone sequence FASTA file for blastp to find out the identity of chains.

