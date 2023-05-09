# cuisines-pie-photochem

This repository hosts code for the CUISINES Photochemical model Inter-comparison for Exoplanet (PIE) Science. The goal of the project is to compare photochemical models. This repository is for the [photochem](https://github.com/Nicholaswogan/photochem) model.

## Installation

To install photochem, create a dedicated `conda` environment and execute the commands below. You will need a Fortran (`gfortran`) and C compiler ('clang' or GNU `gcc`).

```
# Create conda environment
conda create -n pie -c conda-forge python=3.10 numpy scipy pyyaml scikit-build cython jupyter cantera matplotlib numba h5py ninja

# activate
conda activate pie

# photochem
git clone --recursive --branch=dev https://github.com/Nicholaswogan/photochem
cd photochem
git checkout 03623d29519079927c250763c73b6d8c02b0bec1
python -m pip install --no-deps --no-build-isolation . -v
cd ..
rm -rf photochem

# clima
git clone --recursive --branch=dev https://github.com/Nicholaswogan/clima
cd clima
git checkout c2abeefd048973008c6fd6962a21371b41a7a74b
python -m pip install --no-deps --no-build-isolation . -v
cd ..
rm -rf clima
```