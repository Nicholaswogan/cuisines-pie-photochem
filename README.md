# cuisines-pie-photochem

This repository hosts code for the CUISINES Photochemical model Inter-comparison for Exoplanet (PIE) Science. The goal of the project is to compare photochemical models. This repository is for the [photochem](https://github.com/Nicholaswogan/photochem) model.

## Installation

```
conda create -n pie -c conda-forge photochem=0.4.2 matplotlib jupyter
```

## Units of photochem output files

The units of the atmosphere files that photochem outputs are as follows
- `alt` is km
- `press` is bars
- `den` is molecules/cm^3
- `temp` is K
- `eddy` is cm^2/s
- Columns labeled by species are volume mixing ratios.
- Species with "aer" at the end of their name are particles (e.g. H2SO4aer). The column gives the effective particle volume mixing ratio.
- Species with "_r" at the end of their name are particle radii in cm (e.g. H2SO4aer_r)