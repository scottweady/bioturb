# Bioturb

This repository contains scripts for generating and processing simulation data associated with the manuscript "Correlations, mean-field limits, and transition to the concentrated regime in motile particle suspensions," B. Palmer, S. Weady et al. (https://arxiv.org/abs/2505.18299). The source code for the discrete particle simulations can be found in the directory `./CellSimulatorSBT`, which also includes compilation instructions. The source code for the continuum simulations is located in `./bingham-cpp`, which also includes compilation instructions. The scripts for data processing and analysis are located in `./analysis`.

## Data analysis

Data processing scripts for computing autocorrelation functions and mean-square displacement from the discrete particle data can be called as

```
python analysis/correlation.py {nu} {L} {overwrite (1, optional)}
python analysis/msd.py {nu} {L} {overwrite (1, optional)}
```

These use i/o and processing functions defined in `./analysis/sbt.py`. Given an input volume fraction nu and box size L, the analysis scripts will search the directories specificed in `./analysis/config.py` for corresponding `.vtp` files, and write into the target directory specified in the same `config.py` file. An example with nu = 0.0625 and L = 25 is located in `.data/examples`.

## Dependencies

The data processing codes require the python libraries `finufft, natsort, numpy, pyyaml, vtk`, and have been tested with python 3.11.11. It is recommeded to install these in a virtual environment. The vtk module is large and takes some time to load. Detailed version information can be found in `requirements.txt`. 