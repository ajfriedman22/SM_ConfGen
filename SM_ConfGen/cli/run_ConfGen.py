import os
from subprocess import Popen
from SM_ConfGen.sm_confgen import SM_REMD
import argparse
import time
import sys

def f(cmd):
    os.system(cmd) 

def initialize(args):
    parser = argparse.ArgumentParser(
        description='This code runs a runs a set of TREMD simulations given necessary inputs.')
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

    TREMD = SM_REMD(args.yaml)

    arguments = ['mpirun', '-np', str(TREMD.n_rep), TREMD.gmx_executable, 'mdrun', '-deffnm', 'md', '-multidir']

    #Add multiple directories for replica exchange run
    for i in range(TREMD.n_rep):
        arguments.append(f'rep_{i}')
        
    arguments.append('-replex')
    arguments.append(str(TREMD.replex_rate))

    if TREMD.runtime_args is not None:
        # Turn the dictionary into a list with the keys alternating with values
        add_args = [elem for pair in TREMD.runtime_args.items() for elem in pair]
        arguments.extend(add_args)
    
    restart = True
    for r in range(TREMD.n_rep):
        if not os.path.exists(f'rep_{r}/md.cpt'): #Use checkpoint if present for all replicas
            restart = False

    if restart:
        arguments.append(['-cpi', 'md'])

    #Run Simulation
    process = Popen(arguments)
    process.wait()
