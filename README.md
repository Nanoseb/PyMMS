# PyMMS
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4428181.svg)](https://doi.org/10.5281/zenodo.4428181)

This python module generates Manufactured Solutions for RANS equations using Sympy. It also directly generates a compilable Fortran module file containing the analytical source terms to be used in a CFD code.

The following turbulence models have been implemented:
- Spalart-Allmaras (1992) (both Standard and "noft2" variant)
- One equation Eddy viscosity model from Menter (1997)

## Installation
Can be installed and used like any python package, for example:
```
pip3 install https://github.com/Nanoseb/PyMMS/archive/master.zip
```

## Usage
See examples folder:
- `basic_example.py`: a very basic test case to show how PyMMS can be used
- `3D_RANS_unsteady.py`: implementation of the Manufactured Solution defined in [Eca et al. (2012)](https://doi.org/10.1080/10618562.2012.717617)


# How to cite?

This code can be cited with:
```
@software{lemaire_sebastien_2021_4428181,
  author       = {Lemaire, Sébastien},
  title        = {PyMMS},
  month        = jan,
  year         = 2021,
  publisher    = {Zenodo},
  version      = {v1.0},
  doi          = {10.5281/zenodo.4428181},
  url          = {https://doi.org/10.5281/zenodo.4428181}
}
```
