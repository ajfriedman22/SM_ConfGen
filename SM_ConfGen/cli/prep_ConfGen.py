import time
import copy
import shutil
import traceback
import argparse
import numpy as np
from mpi4py import MPI
from datetime import datetime
import sys
import os

from SM_ConfGen.utils import utils
from SM_ConfGen.sm_confgen import SM_REMD

def initialize(args):
    parser = argparse.ArgumentParser(
        description='This code runs a REXEE simulation given necessary inputs.')
    parser.add_argument('-y',
                        '--yaml',
                        type=str,
                        default='params.yaml',
                        help='The input YAML file that contains REXEE parameters. (Default: params.yaml)')
    parser.add_argument('-o',
                        '--output',
                        type=str,
                        default='ConfGen_log.txt',
                        help='The output file for logging how replicas interact with each other. \
                            (Default: ConfGen_log.txt)')
    parser.add_argument('-m',
                        '--maxwarn',
                        type=int,
                        default=0,
                        help='The maximum number of warnings in parameter specification to be ignored.')
    args_parse = parser.parse_args(args)

    return args_parse

def main():
    t1 = time.time()
    args = initialize(sys.argv[1:])

    sys.stdout = utils.Logger(logfile=args.output)
    sys.stderr = utils.Logger(logfile=args.output)

    # Step 1: Set up MPI rank and instantiate ReplicaExchangeEE to set up REXEE parameters
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()  # Note that this is a GLOBAL variable

    if rank == 0:
        print(f'Current time: {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}')
        print(f'Command line: {" ".join(sys.argv)}\n')

    TREMD = SM_REMD(args.yaml)

    if rank == 0:
        # Print out simulation parameters
        TREMD.print_params()

        # Print out warnings and fail if needed
        for i in TREMD.warnings:
            print(f'\n{i}\n')

        if len(TREMD.warnings) > args.maxwarn:
            print(f"The execution failed due to warning(s) about parameter spcificaiton. Check the warnings, or consider setting maxwarn in the input YAML file if you find them harmless.")  # noqa: E501, F541
            comm.Abort(101)
    
    #Step 2: Parameterize System (if files not available)
    if rank == 0:
        if not os.path.exists('prep/conf.gro') or not os.path.exists('prep/topol.top'):
            if not os.path.exists('prep'):
                os.mkdir('prep')
            TREMD.parameterize_system()
    
    #Step 3: Solvate and Neutralize System (if files not available)
    if rank == 0:
        if not os.path.exists('prep/ions.gro'):
            TREMD.solvate_system()

    #Step 4: Perform Energy Minimization (if files not available)
    if rank == 0:
        if not os.path.exists('prep/min.gro'):
            # 4-1 Run Energy Minimization
            TREMD.run_grompp('EM')
            TREMD.run_mdrun('EM')
        if not os.path.exists('prep/min.gro'): 
            # 4-2 Check that Energy minimization ran properly
            raise Exception('Energy minimization did not finish!')
    # Synchronyzation point for all ranks
    comm.barrier()

    #Step 5: Preform NVT equilibration (if files not available)
    if not os.path.exists('prep/rep_0/nvt.gro') and not os.path.exists(f'prep/rep_{TREMD.n_rep-1}/nvt.gro'): 
        # 5-1. Set up input files for all simulations
        if rank < TREMD.n_rep:
            os.chdir('prep')
            if not os.path.exists(f'rep_{rank}/nvt.mdp'):
                if not os.path.exists(f'rep_{rank}'):
                    os.mkdir(f'rep_{rank}')
                MDP = TREMD.initialize_MDP(rank, True, TREMD.mdp[1])
                MDP.write(f'rep_{rank}/nvt.mdp', skipempty=True)
            os.chdir('../')

        # Synchronyzation point for all ranks
        comm.barrier()

        # 5-2. Run the first set of simulations
        TREMD.run_grompp('NVT')
        TREMD.run_mdrun('NVT')
    
    # Synchronyzation point for all ranks
    comm.barrier()

    #Step 6: Preform NPT equilibration (if files not available)
    if not os.path.exists('prep/rep_0/npt.gro') and not os.path.exists(f'prep/rep_{TREMD.n_rep-1}/npt.gro'):
        # 6-1. Set up input files for all simulations
        if rank < TREMD.n_rep:
            os.chdir('prep')
            if not os.path.exists(f'rep_{rank}/npt.mdp'):
                MDP = TREMD.initialize_MDP(rank, True, TREMD.mdp[2])
                MDP.write(f'rep_{rank}/npt.mdp', skipempty=True)
            os.chdir('../')
        # Synchronyzation point for all ranks
        comm.barrier()

        # 6-2. Run the first set of simulations
        TREMD.run_grompp('NPT')
        TREMD.run_mdrun('NPT')
    # Synchronyzation point for all ranks
    comm.barrier()

    #Step 7: Run TREMD (check for chechpoint)
    check_pt = False 
    if not os.path.exists(f'rep_{rank}/md.cpt'): #Use checkpoint if present for all replicas
        check_pt = False

    # Synchronyzation point for all ranks
    comm.barrier()

    if not check_pt:
        if rank < TREMD.n_rep:
            if not os.path.exists(f'rep_{rank}'):
                os.mkdir(f'rep_{rank}')
            MDP = TREMD.initialize_MDP(rank, True, TREMD.mdp[3])
            MDP.write(f'rep_{rank}/md.mdp', skipempty=True)
        TREMD.run_grompp('MD')
        # Synchronyzation point for all ranks
        comm.barrier()