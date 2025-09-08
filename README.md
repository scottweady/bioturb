# Bioturb

This repository contains scripts for processing simulation data generated with the software package CellSimulator (TODO: url).

Data processing scripts for computing autocorrelation functions and mean-square displacement are located in ./scripts and should be called as

    python scripts/correlation.py {nu} {L}
    python scripts/msd.py {nu} {L}

These will load data, stored in .vtp format, from the directory ./data/raw and write into ./data/processed. The vtp files are located in

    .data/raw/v1.0_nu{nu}_L{L}_A5/result/result{t0}-{tf}

Useful functions are located in scripts/sbt.py.

## Dependencies

This project requires the python (3.10.10) libraries

-finufft
-natsort
-numpy
-pyyaml
-vtk

It is recommeded to install these in a virtual environment. Be warned, the vtk module is large (~300MB).

This repository is associated with the manuscript "Correlations, mean-field limits, and transition to the concentrated regime in motile particle suspensions," B. Palmer et al. (https://arxiv.org/abs/2505.18299).
