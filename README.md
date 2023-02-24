# cuisines-pie-photochem

This repository hosts code for the CUISINES Photochemical model Inter-comparison for Exoplanet (PIE) Science. The goal of the project is to compare photochemical models. This repository is for the [photochem](https://github.com/Nicholaswogan/photochem) model.

## Installation

To install photochem, create a dedicated `conda` environment and execute the commands below. You will need a Fortran (`gfortran`) and C compiler ('clang' or GNU `gcc`).

```
# Create conda environment
conda create -n photochem -c conda-forge python=3.10 numpy scipy pyyaml scikit-build cython jupyter cantera matplotlib numba h5py

# activate
conda activate photochem

# photochem
git clone --recursive --branch=dev https://github.com/Nicholaswogan/photochem
cd photochem
git checkout 6c241b213542b3cc28457ec038d457c0a2705927
python -m pip install --no-deps --no-build-isolation . -v
cd ..
rm -rf photochem

# clima
git clone --recursive --branch=dev https://github.com/Nicholaswogan/clima
cd clima
git checkout dd51fed4bdbf7d8080bed8d2231a79ba25ce59c7
python -m pip install --no-deps --no-build-isolation . -v
cd ..
rm -rf clima
```