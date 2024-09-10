Small Molecule Conformer Generation
==============================

Small molecule conformer generation using TREMD

### Installation
'''pre
git clone git@github.com:ajfriedman22/SM_ConfGen.git
mamba env create -f env.yaml
mamba activate sm_confgen
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
