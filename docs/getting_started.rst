1. Getting Started
===============
:code:`sm_confgen` is a Python package providing a method of generating ligand 
conformers from Temperature Replic Exchange Molecular Dynamics (TREMD) 
simulations in GROMACS.

2. Installation
===============
2.1. Requirements
-----------------
Before installing :code:`sm_confgen`, one should have a working version of `GROMACS`_. Please refer to the 
GROMACS documentation for full installation instructions. All other dependencies are either included in the 
provided conda environment or will be installed automatically by pip.

.. _`GROMACS`: https://manual.gromacs.org/current/install-guide/index.html 

2.2. Installation from source
-----------------------------
One can also install :code:`sm_confgen` from the source code, which is available in our
`GitHub repository`_. Specifically, one can execute the following commands:
::

    git clone https://github.com/ajfriedman22/SM_ConfGen
    cd sm_confgen/
    mamba env create -n env.yaml
    mamba activate sm_confgen
    pip install .

If you would like to install the package in the editable mode, simply append the last command with the flag :code:`-e`
so that changes you make in the source code will take effect without re-installation of the package. This is particularly
useful if you would like to contribute to the development of the package. (Pull requests and issues are always welcome!)