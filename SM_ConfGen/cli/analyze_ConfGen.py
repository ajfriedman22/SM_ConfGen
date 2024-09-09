import time
import argparse
from mpi4py import MPI
from datetime import datetime
import sys
import os

from SM_ConfGen.utils import utils
from SM_ConfGen.sm_confgen import SM_REMD

def initialize(args):
    parser = argparse.ArgumentParser(
        description='This code analyzes the output from TREMD to generate conformers given necessary inputs.')
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

   
    TREMD = SM_REMD(args.yaml)
    
    # Step 2: Create analysis directory if not present
    if not os.path.exists('analysis'):
        os.mkdir('analysis')

    # Step 1: Process trajectory for lowest temperature replicate
    if not os.path.exists('analysis/center.xtc') or not os.path.exists('analysis/md.gro'):
        TREMD.process_traj()

    #Step 2: Compute Dihedral Peaks
    if not os.path.exists('analysis/dihe_ind_max.csv'):
        TREMD.compute_dihedral_peaks()

    #Step 3: Determine sampled conformations and cluster
    TREMD.clust_dihedrals()
    