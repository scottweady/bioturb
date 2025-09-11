# Bioturb

This repository contains scripts for processing simulation data generated with the software package CellSimulator.

Data processing scripts for computing autocorrelation functions and mean-square displacement are located in `./analysis` and should be called as

```
python analysis/correlation.py {nu} {L} {overwrite (0 or 1, optional)}
python analysis/msd.py {nu} {L} {overwrite (0 or 1, optional)}
```

These will load data, stored in `.vtp` format, from the directory `./data/raw` and write into `./data/processed`. The vtp files are located in

```
.data/raw/v1.0_nu{nu}_L{L}_A5/result/result{t0}-{tf}
```

Useful functions are located in `analysis/sbt.py`.

## Dependencies

This project requires the python libraries `finufft, natsort, numpy, pyyaml, vtk`. Detailed version information can be found in `requirements.txt`.

It is recommeded to install these in a virtual environment. The vtk module is large and takes some time to load.

## References

This repository is associated with the manuscript "Correlations, mean-field limits, and transition to the concentrated regime in motile particle suspensions," B. Palmer et al. (https://arxiv.org/abs/2505.18299).
