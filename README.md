Small Molecule Conformer Generation
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/SM_ConfGen/workflows/CI/badge.svg)](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/SM_ConfGen/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/Small Molecule Conformer Generation/branch/main/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/Small Molecule Conformer Generation/branch/main)


Small molecule conformer generation using TREMD

### Installation
'''
git clone git@github.com:ajfriedman22/SM_ConfGen.git
cd SM_ConfGen
pip install .
'''

### Useage
N = Number of temperature replicates specified. If the default temperatures are selected N=27.
'''
mpirun -np N prep_ConfGen -y params.yaml
run_ConfGen -y params.yaml
mpirun -np N analyze_ConfGen -y params.yaml
'''

### Copyright

Copyright (c) 2024, Anika Friedman


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
