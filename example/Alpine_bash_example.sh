#!/bin/bash

#SBATCH -N 1
#SBATCH -p amilan
#SBATCH --ntasks=54
#SBATCH --job-name=name
#SBATCH -t 24:00:00
#SBATCH --output=name.out

ml anaconda
conda activate sm_confgen

mpirun -np 27 prep_ConfGen -y params_setup.yaml
run_ConfGen -y params.yaml
analyze_ConfGen -y params.yaml
