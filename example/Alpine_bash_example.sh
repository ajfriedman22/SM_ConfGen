#!/bin/bash

#SBATCH -N 1
#SBATCH -p amilan
#SBATCH --ntasks=4
#SBATCH --job-name=name
#SBATCH -t 24:00:00
#SBATCH --output=name.out

module purge
ml anaconda
conda activate sm_confgen

mpiexec -np 2 prep_ConfGen -y params_setup.yaml

ml gcc/11.2.0
ml openmpi/4.1.1
ml gromacs/2022.4

run_ConfGen -y params.yaml
analyze_ConfGen -y params.yaml
