# mesopelagic-airgun-response

Code for reproducing statistical analysis and figures in a study of mesopelagic dive responses to seismic airgun exposure, using AUV-mounted split-beam echosounder tracking.

## Contents

- `fig1.m`–`fig5.m` — one standalone script per manuscript figure. `fig4.m`/`fig5.m` also print the statistical test results reported in the text.
- `utils/` — small helper functions shared across scripts (axis styling, font handling, ship-track distance calculation).
- `data/` — input data (see below; not included in this repository).

## Requirements

MATLAB R2019b or later (for `tiledlayout`/`nexttile`), plus:
- Statistics and Machine Learning Toolbox (`kruskalwallis`, `ranksum`, `boxplot`, `prctile`, `ksdensity`, `normcdf`)
- Signal Processing Toolbox (`spectrogram`, `hann`)
- Image Processing Toolbox (`imresize`, `imgaussfilt`)

## Data

Input data is not tracked in this repository (several files exceed GitHub's size limits) and is instead archived separately on Zenodo: https://doi.org/10.5281/zenodo.21684608.

Download the archive and place it as a `data/` folder alongside these scripts, preserving the subfolder structure (`data/auv/`, `data/ship/`, `data/airgun/`, plus the two bathymetry `.nc` files at the top level of `data/`).

## Citation

If you use this code, please cite the associated paper:

[Citation]
