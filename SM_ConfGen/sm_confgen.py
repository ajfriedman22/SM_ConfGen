"""Provide the primary functions."""
import os
import sys
import copy
import yaml
import subprocess
from SM_ConfGen.utils import utils
from SM_ConfGen.utils import analysis_func as af
from SM_ConfGen.utils.exceptions import ParameterError
from SM_ConfGen.utils import gmx_parser
from mpi4py import MPI
from openff.toolkit import Molecule, ForceField
from openff.units import unit
import numpy as np
from openff.interchange import Interchange
from os.path import dirname, join as joinpath
import mdtraj as md
import pandas as pd
DATADIR = joinpath(dirname(__file__), 'data')

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def canvas(with_attribution=True):
    """
    Placeholder function to show example docstring (NumPy format).

    Replace this function and doc string for your own project.

    Parameters
    ----------
    with_attribution : bool, Optional, default: True
        Set whether or not to display who the quote is from.

    Returns
    -------
    quote : str
        Compiled string including quote and optional attribution.
    """

    quote = "The code is but a canvas to our imagination."
    if with_attribution:
        quote += "\n\t- Adapted from Henry David Thoreau"
    return quote

class SM_REMD:
    """
    This class provides a variety of functions useful for setting up and running Temperature 
    Replica Exchange Molecular Dynamics Simulations for the purpose of generating small molecule conformers.
    Upon instantiation, all parameters in the YAML file will be assigned to an attribute in the class. 
    In addition to these variables, below is a list of attributes of the class. (All the the attributes 
    are assigned by :obj:`set_params` unless otherwise noted.)

    :ivar gmx_path: The absolute path of the GROMACS exectuable.
    :ivar gmx_version: The version of the GROMACS executable.
    :ivar yaml: The input YAML file used to instantiate the class. Assigned by the :code:`__init__` function.
    :ivar warnings: Warnings about parameter specification in either YAML or MDP files.
    :ivar template: The instance of the :obj:`MDP` class based on the template MDP file.
    :ivar dt: The simulation timestep in ps.
    :ivar temp_range: An array of temperatures (Kelvin) to run simulations.
    :ivar n_rep: The number of parrallel replica exchange simulations

    """
    

    def __init__(self, yaml_file, analysis=False):
        self.yaml = yaml_file
        self.set_params(analysis)

    def set_params(self, analysis):
        """
        Sets up or reads in the user-defined parameters from a yaml file and an MDP template.
        This function is called to instantiate the class in the :code:`__init__` function of
        class. Specifically, it does the following:

          1. Sets up constants.
          2. Reads in parameters from a YAML file.
          3. Handles YAML parameters.
          4. Checks if the parameters in the YAML file are well-defined.
          5. Reads in parameters from the MDP template.

        After instantiation, the class instance will have attributes corresponding to
        each of the parameters specified in the YAML file.

        :param yaml_file: The file name of the YAML file for specifying the parameters for REXEE.
        :type yaml_file: str
        :param analysis: Whether the instantiation of the class is for data analysis of REXEE simulations.
            The default is :code:`False`
        :type analysis: bool

        :raises ParameterError:
              - If a required parameter is not specified in the YAML file.
              - If a specified parameter is not recognizable.
              - If a specified option is not available for a parameter.
              - If the data type or range (e.g., positive or negative) of a parameter is not correct.
              - If an invalid MDP file is detected.
        """

        self.warnings = []  # Store warnings, if any.

        # Step 0: Set up constants

        # Step 1: Read in parameters from the YAML file.
        with open(self.yaml) as f:
            params = yaml.load(f, Loader=yaml.FullLoader)

        for attr in params:
            setattr(self, attr, params[attr])

        # Step 2: Handle the compulsory YAML parameters
        required_args = [
            "gmx_executable",
            "input_structure",
            "output_name"
        ]

        for i in required_args:
            if hasattr(self, i) is False or getattr(self, i) is None:
                raise ParameterError(
                    f"Required parameter '{i}' not specified in {self.yaml}."
                )  # noqa: F405

        # Step 3: Handle the optional YAML parameters
        # Key: Optional argument; Value: Default value
        optional_args = {
            "mdp": [f'{DATADIR}/default_mdp/em.mdp', f'{DATADIR}/default_mdp/nvt.mdp', f'{DATADIR}/default_mdp/npt.mdp', f'{DATADIR}/default_mdp/TREMD.mdp'],
            "temp_range": np.array([300, 306, 311.08, 316.23, 321.45, 326.74, 332.10, 337.54, 343.05, 348.63, 354.29, 360.02, 365.84, 371.73, 377.71, 383.76, 389.90, 396.12, 402.43, 408.83, 415.31, 421.88, 428.55, 435.30, 442.14, 449.07, 450],dtype=float),
            "nst_sim": None,
            "dt": None,
            "replex_rate": 1000,
            "grompp_args": None,
            "runtime_args": None,
            "water_model": 'TIP3P',
            "conformer_threshold": 1,
            "peak_threshold": 0.0005,
            "min_probability": 0.001, 
        }

        for i in optional_args:
            if hasattr(self, i) is False or getattr(self, i) is None:
                setattr(self, i, optional_args[i])

        # all_args: Arguments that can be specified in the YAML file.
        all_args = required_args + list(optional_args.keys())
        for i in params:
            if i not in all_args:
                self.warnings.append(f'Warning: Parameter "{i}" specified in the input YAML file is not recognizable.')

        # Step 4: Check if the parameters in the YAML file are well-defined
        #Check string inputs
        params_str = ['gmx_executable', 'input_structure', 'output_name']
        for i in params_str:
            if type(getattr(self, i)) != str:
                raise ParameterError(f"The parameter '{i}' should be a string.")
        
        #Check that input files exsist and have appropriate file extensions
        if isinstance(self.mdp, str):
            if not self.mdp.lower().endswith(('.mdp')):#Add file extension
                self.mdp = self.mdp + '.mdp'
            if not os.path.exists(self.mdp):
                raise ParameterError(f'{self.mdp} file does not exsist')
            self.warnings.append(f'Using Default MDP tempates for energy minimization and equilibration. Using {self.mdp} as template for production run only')
            self.mdp = ['data/default_mdp/em.mdp', 'data/default_mdp/nvt.mdp', 'data/default_mdp/npt.mdp', self.mdp]
        elif isinstance(self.mdp, list) and len(self.mdp) != 4:
            raise ParameterError('Input MDP file must either be for production simulation only or for energy minimization, NVT equilibration, NPT equilibration, and production simulations')
        elif isinstance(self.mdp, list) and len(self.mdp) == 4:
            file_w_extension = []
            for file in self.mdp:
                if not file.lower().endswith(('.mdp')):#Add file extension
                    file = file + '.mdp'
                file_w_extension.append(file)
                if not os.path.exists(file):
                    raise ParameterError(f'{file} file does not exsist') 
            self.mdp = file_w_extension
        else:
            raise ParameterError('Input MDP file must either be list with len(mdp) == 4 or str')

        if not self.input_structure.lower().endswith(('.pdb', '.sdf')):
            if os.path.exists(self.input_structure + '.pdb'):
                self.input_structure = self.input_structure + '.pdb'
            elif os.path.exists(self.input_structure + '.sdf'):
                self.input_structure = self.input_structure + '.sdf'
            else:
                raise ParameterError(f'{self.input_structure} file does not exsist')
        elif not os.path.exists(self.input_structure):
            raise ParameterError(f'{self.input_structure} file does not exsist')
        
        #Check that temperature range are integers are 
        if not all(isinstance(x, float) for x in self.temp_range) and not all(isinstance(x, int) for x in self.temp_range):
            raise ParameterError('The input temperature range must be of type int or float')
        if self.temp_range[-1] > 500:
            self.warnings.append(f'WARNING: Maximum temperature of {self.temp_range[-1]} is > 500 K which may cause issues solvent')
        self.n_rep = len(self.temp_range)

        #Check simulation parameters
        int_params = ['nst_sim', 'dt', 'replex_rate', 'conformer_threshold']
        for i in int_params:
            if type(getattr(self, i)) != int and getattr(self, i) != None:
                raise ParameterError(f"The parameter '{i}' should be an integer not {type(getattr(self, i))}.")
        
        float_params = ['peak_threshold', 'min_probability']
        for i in float_params:
            if type(getattr(self, i)) != float and getattr(self, i) != None:
                raise ParameterError(f"The parameter '{i}' should be a float not {type(getattr(self, i))}.")
            
        #Check that option input args are formatted correctly
        params_dict = ['grompp_args', 'runtime_args']
        for i in params_dict:
            if getattr(self, i) is not None and not isinstance(getattr(self, i), dict):
                raise ParameterError(f"The parameter '{i}' should be a dictionary.")
        
        #Check that water model is suported
        if self.water_model != 'TIP3P':
            raise ParameterError(f'Water model {self.water_model} is not supported only TIP3P is an acceptable input')
        
        # Step 5: Read in and check parameters from the MDP template
        self.template = gmx_parser.MDP(self.mdp[-1])
        if self.dt == None:
            self.dt = self.template["dt"]  # ps

        if self.nst_sim is None:
            self.nst_sim = self.template["nsteps"]
        
        if self.dt > 0.004:
            self.warnings.append(f'WARNING: Time step dt={self.dt} is greater than the reccommended 4fs. This may cause system instability.')
        
        # Step 6. Check the executables
        if analysis is False:
            self._check_gmx_executable()

    def _check_gmx_executable(self):
        """
        Checks if the GROMACS executable can be used and gets its absolute path and version.
        """
        try:
            result = subprocess.run(['which', self.gmx_executable], capture_output=True, text=True, check=True)
            self.gmx_path = result.stdout.strip()  # this can be exactly the same as self.gmx_executable

            version_output = subprocess.run([self.gmx_path, "-version"], capture_output=True, text=True, check=True)
            for line in version_output.stdout.splitlines():
                if "GROMACS version" in line:
                    self.gmx_version = line.split()[-1]
                    break
        except subprocess.CalledProcessError:
            print(f"{self.gmx_executable} is not available on this system.")
        except Exception as e:
            print(f"An error occurred:\n{e}")

    def print_params(self, params_analysis=False):
        """
        Prints important parameters related to the TREMD simulation.

        Parameters
        ----------
        params_analysis : bool, optional
            If True, additional parameters related to data analysis will be printed. Default is False.
        """
        print("Important parameters of ConfGen")
        print("============================")
        print(f"Python version: {sys.version}")
        print(f"GROMACS executable: {self.gmx_path}")  # we print the full path here
        print(f'Simulation Input Structure: {self.input_structure}')
        print(f'Simulation MDP Files: {self.mdp}')
        print(f"Number of replicas: {self.n_rep}")
        print(f"Temperatures for Replicas: {self.temp_range}")
        print(f"Length of each replica: {self.dt * self.nst_sim / 1000} ns")
        print(f"Additional grompp arguments: {self.grompp_args}")
        print(f"Additional runtime arguments: {self.runtime_args}")

        if params_analysis is True: #### FILL IN LATTER ######
            print()
    
    def initialize_MDP(self, idx:int , equil:bool, mdp_file:str):
        """
        Initializes the MDP object for all temperatures

        Parameters
        ----------
        idx : int
            Index of the simulation whose MDP parameters need to be initialized.
        equil : bool
            Is this equilibration or production simulation
        mdp_file : str
            The file path for the mdp template
        Returns
        -------
        MDP : :obj:`.gmx_parser.MDP` obj
            An updated object of :obj:`.gmx_parser.MDP` that can be used to write MDP files.
        """
        mdp = gmx_parser.MDP(mdp_file)
        MDP = copy.deepcopy(mdp)
        if not equil:
            MDP["nsteps"] = self.nst_sim
            MDP["dt"] = self.dt
        MDP["ref_t"] = self.temp_range[idx]

        return MDP
    
    def parameterize_system(self):
        """
        Prepares the system files for all replicate simulations

        """
        # Create interchange object
        mol = Molecule.from_file(self.input_structure)
        sage = ForceField("openff-2.0.0.offxml")
        cubic_box = unit.Quantity(30 * np.eye(3), unit.angstrom)
        interchange = Interchange.from_smirnoff(topology=[mol], force_field=sage, box=cubic_box)

        # Export topology and coordinates
        interchange.to_gro('prep/conf.gro')
        interchange.to_top('prep/init-topol.top')

        #Get water parameters depending on the water model
        if self.water_model == 'TIP3P':
            file_path_water = f'{DATADIR}/tip3p_ions.itp'

        # Edit topology to add necessary elements
        init_topology = open('prep/init-topol.top', 'r').readlines()
        topology = open('prep/topol.top', 'w')
        sect = None
        for line in init_topology:
            str_line = line.split(' ')
            while '' in str_line: str_line.remove('')
            if sect == 'atomtypes' and len(str_line) == 1: # Add water and ion atom types from data
                topology.write('Na          11      22.99    0.0000  A   2.43928e-01  3.65846e-02\nCl          17      35.45    0.0000  A   4.47766e-01  1.48913e-01\nOW           8      16.00    0.0000  A   3.15061e-01  6.36386e-01\nHW           1       1.008   0.0000  A   0.00000e+00  0.00000e+00\n')
                topology.write(line)
                sect = None
            elif sect == 'dihedrals' and len(str_line) == 1: # Add position restraint section to topology
                topology.write('; Include Position restraint file\n#ifdef POSRES\n#include "posre.itp"\n#endif\n')
                topology.write(line)
                sect = None
            elif sect == 'excl' and len(str_line) == 1:
                topology.write(f'#include "{file_path_water}"\n') # Add water and ion molecule types from data
                topology.write(line)
                sect = None
            else:
                topology.write(line)
            
            # Determine section
            if 'atomtypes' in line:
                sect = 'atomtypes'
            elif 'dihedrals' in line:
                sect = 'dihedrals'
            elif 'exclusions' in line:
                sect = 'excl'        
        topology.close()

        #Create position restraints file
        create_posre = [self.gmx_executable, 'genrestr', '-f', 'prep/conf.gro', '-o', 'prep/posre.itp']
        returncode, stdout, stderr = utils.run_gmx_cmd(create_posre, '1\n')
        if returncode != 0:
            print(f'Error (return code: {returncode}):\n{stderr}')

    def solvate_system(self):
        """
        Solvate and neutralize charge for the system

        """
        #Create Box
        create_box = [self.gmx_executable, 'editconf', '-f', 'prep/conf.gro', '-o', 'prep/box.gro', '-bt', 'dodecahedron', '-d', '1']
        returncode, stdout, stderr = utils.run_gmx_cmd(create_box)
        if returncode != 0:
            print(f'Error (return code: {returncode}):\n{stderr}')

        #Solvate Box
        solv_box = [self.gmx_executable, 'solvate', '-cp', 'prep/box.gro', '-p', 'prep/topol.top', '-o', 'prep/solv.gro', '-cs']
        returncode, stdout, stderr = utils.run_gmx_cmd(solv_box)
        if returncode != 0:
            print(f'Error (return code: {returncode}):\n{stderr}')

        #Add neuralizing ions
        neutral_box = [self.gmx_executable, 'grompp', '-f', f'{DATADIR}/default_mdp/ions.mdp', '-c', 'prep/solv.gro', '-p', 'prep/topol.top', '-o', 'prep/ions.tpr']
        # Add additional arguments if any
        if self.grompp_args is not None:
            # Turn the dictionary into a list with the keys alternating with values
            add_args = [elem for pair in self.grompp_args.items() for elem in pair]
            neutral_box.extend(add_args)
        returncode, stdout, stderr = utils.run_gmx_cmd(neutral_box)
        if returncode != 0:
            print(f'Error (return code: {returncode}):\n{stderr}')
        neutral_run = [self.gmx_executable, 'genion', '-s', 'prep/ions.tpr', '-o', 'prep/ions.gro', '-p', 'prep/topol.top', '-pname', 'NA', '-nname', 'CL', '-conc', '0.15', '-neutral']
        returncode, stdout, stderr = utils.run_gmx_cmd(neutral_run, '4\n')
        if returncode != 0:
            print(f'Error (return code: {returncode}):\n{stderr}')   

        #Check that files generated properly
        if not os.path.exists('prep/ions.gro'):
            raise Exception('Error is System Solvation!')

    def run_grompp(self, step):
        """
        Prepares TPR files for the simulation ensemble using the GROMACS :code:`grompp` command.

        Parameters
        ----------
        n : int
            The iteration index (starting from 0).
        step : str
            Whether this is an energy minimization, NVT equilibration, NPT equilibration or production run
        """
        args_list = []
        top = 'prep/topol.top'
        if step == 'EM':
            arguments = [self.gmx_executable, 'grompp']

            # Input Files
            mdp = self.mdp[0]
            gro = 'prep/ions.gro'
            output_base = 'prep/min'

            # Add input file arugments
            arguments.extend(['-f', mdp, '-c', gro, '-p', top])

            # Add output file arguments
            arguments.extend([
                "-o", f"{output_base}.tpr",
                "-po", f"{output_base}_mdout.mdp"
            ])

            # Add additional arguments if any
            if self.grompp_args is not None:
                # Turn the dictionary into a list with the keys alternating with values
                add_args = [elem for pair in self.grompp_args.items() for elem in pair]
                arguments.extend(add_args)
            
            returncode, stdout, stderr = utils.run_gmx_cmd(arguments)
            if returncode != 0:
                print(f'Error on rank {rank} (return code: {returncode}):\n{stderr}')
        else:
            for i in range(self.n_rep):
                arguments = [self.gmx_executable, 'grompp']

                # Input files
                if step == 'NVT':
                    mdp = f'prep/rep_{i}/nvt.mdp'
                    gro = 'prep/min.gro'
                    restr = 'prep/min.gro'
                    cpt = None
                    output_base = f'prep/rep_{i}/nvt'
                elif step == 'NPT':
                    mdp = f'prep/rep_{i}/npt.mdp'
                    gro = f'prep/rep_{i}/nvt.gro'
                    restr = f'prep/rep_{i}/nvt.gro'
                    cpt = f'prep/rep_{i}/nvt.cpt'
                    output_base = f'prep/rep_{i}/npt'
                else:
                    mdp = f'rep_{i}/md.mdp'
                    gro = f'prep/rep_{i}/npt.gro'
                    restr = None
                    cpt = f'prep/rep_{i}/npt.cpt'
                    output_base = f'rep_{i}/md'

                # Add input file arguments
                arguments.extend(['-f', mdp, '-c', gro, '-p', top])
                if restr != None:
                    arguments.extend(['-r', restr])
                if cpt != None:
                    arguments.extend(['-t', cpt])

                # Add output file arguments
                arguments.extend([
                    "-o", f"{output_base}.tpr",
                    "-po", f"{output_base}_mdout.mdp"
                ])

                # Add additional arguments if any
                if self.grompp_args is not None:
                    # Turn the dictionary into a list with the keys alternating with values
                    add_args = [elem for pair in self.grompp_args.items() for elem in pair]
                    arguments.extend(add_args)

                args_list.append(arguments)

            # Run the GROMACS grompp commands in parallel
            returncode = None  # Initialize as None for all ranks (necessary for the case when -np > n_rep, which is rare)
            if rank == 0:
                print('Generating TPR files ...')
            if rank < self.n_rep:
                returncode, stdout, stderr = utils.run_gmx_cmd(args_list[rank])
                if returncode != 0:
                    print(f'Error on rank {rank} (return code: {returncode}):\n{stderr}')

            # gather return codes at rank 0
            code_list = comm.gather(returncode, root=0)

            if rank == 0:
                # Filter out None values which represent ranks that did not execute the command
                code_list = [code for code in code_list if code is not None]
                if code_list != [0] * self.n_rep:
                    MPI.COMM_WORLD.Abort(1)   # Doesn't matter what non-zero returncode we put here as the code from GROMACS will be printed before this point anyway.  # noqa: E501

    def run_mdrun(self, step):
        """
        Executes GROMACS mdrun commands in parallel.

        Parameters
        ----------
        step : str
            Whether this is an energy minimization, NVT equilibration, NPT equilibration or production run
        """
        # We will change the working directory so the mdrun command should be the same for all replicas.
        arguments = [self.gmx_executable, 'mdrun']

        # Add input file arguments
        arguments.extend(['-deffnm'])
        if step == 'EM':
            arguments.extend(['min'])
        elif step == 'NVT':
            arguments.extend(['nvt'])
        else:
            arguments.extend(['npt'])

        if self.runtime_args is not None:
            # Turn the dictionary into a list with the keys alternating with values
            add_args = [elem for pair in self.runtime_args.items() for elem in pair]
            arguments.extend(add_args)

        os.chdir('prep')
        if step == 'EM':
            returncode, stdout, stderr = utils.run_gmx_cmd(arguments)
            if returncode != 0:
                print(f'Error (return code: {returncode}):\n{stderr}')
        else:
            # Run the GROMACS mdrun commands in parallel
            returncode = None  # Initialize as None for all ranks (necessary for the case when -np > n_rep, which is rare)
            if rank == 0:
                print(f'Running {step} equilibration simulations ...')
            if rank < self.n_rep:
                os.chdir(f'rep_{rank}/')
                returncode, stdout, stderr = utils.run_gmx_cmd(arguments)
                if returncode != 0:
                    print(f'Error on rank {rank} (return code: {returncode}):\n{stderr}')
                os.chdir('../../')

            # gather return codes at rank 0
            code_list = comm.gather(returncode, root=0)

            if rank == 0:
                # Filter out None values which represent ranks that did not execute the command
                code_list = [code for code in code_list if code is not None]
                if code_list != [0] * self.n_rep:
                    MPI.COMM_WORLD.Abort(1)   # Doesn't matter what non-zero returncode we put here as the code from GROMACS will be printed before this point anyway.  # noqa: E501
        os.chdir('../')

    def process_traj(self):
        """
        Center trajectory for analysis

        """
        if rank == 0:
            # Check that trajectory is present to process
            if not os.path.exists(f'rep_{rank}/md.xtc'):
                raise Exception(f'Error: Trajectory rep_{rank}/md.xtc does not exsist')
            if not os.path.exists(f'rep_{rank}/md.gro'):
                raise Exception(f'Error: File rep_{rank}/md.gro does not exsist')

        if rank == 0:
            # Center trajectory
            arguments = [self.gmx_executable, 'trjconv', '-s', 'rep_0/md.tpr', '-f', 'rep_0/md.xtc', '-pbc', 'cluster', '-center', '-o', 'analysis/center.xtc']
            returncode, stdout, stderr = utils.run_gmx_cmd(arguments, prompt_input='1\n1\n1\n')
            if returncode != 0:
                print(f'Error (return code: {returncode}):\n{stderr}')
        
            # Create GRO File with just the ligand
            arguments = [self.gmx_executable, 'trjconv', '-s', 'rep_0/md.gro', '-f', 'rep_0/md.gro', '-o', 'analysis/md.gro']
            returncode, stdout, stderr = utils.run_gmx_cmd(arguments, prompt_input='1\n')
            if returncode != 0:
                print(f'Error (return code: {returncode}):\n{stderr}')
        
    def compute_dihedral_peaks(self):
        """
        Compute the value for all unique ligand dihedrals to determine which are multi-modal
        """
        # Step 1: Get indices for all dihedrals in ligand
        if rank == 0:
            dihe_name, dihe_ind  = af.define_dihedral('prep/topol.top')

        # Step 2: Compute the dihedral angles
        if rank == 0:
            traj = md.load('analysis/center.xtc', top='analysis/md.gro')
            dihedral = md.compute_dihedrals(traj, indices=dihe_ind)
            #Convert to degree
            dihedral = dihedral*(180/np.pi)
        
        # Step 3: Determine which peaks are multimodal and save that info
        if rank == 0 and not os.path.exists('analysis/dihedrals'):
            os.mkdir('analysis/dihedrals')
        
        if rank == 0:
            df_max = pd.DataFrame(columns=['Dihedral Name', 'Dihderal Atom Index 1', 'Dihderal Atom Index 2', 'Dihderal Atom Index 3', 'Dihderal Atom Index 4', 'Maxima'])
            for i in range(len(dihe_name)):
                maxima, dihe_dist = af.deter_multimodal(dihedral, i, dihe_name, self.peak_threshold, self.min_probability)
                af.plot_torsion(dihe_dist, maxima, dihe_name[i])
                #If multiple peaks add to dataframe
                if len(maxima) > 1:
                    df = pd.DataFrame({'Dihedral Name': dihe_name[i], 'Dihderal Atom Index 1': dihe_ind[i][0], 'Dihderal Atom Index 2': dihe_ind[i][1], 'Dihderal Atom Index 3': dihe_ind[i][2], 'Dihderal Atom Index 4': dihe_ind[i][3], 'Maxima': maxima})
                    df_max = pd.concat([df_max, df])
            df_max.to_csv('analysis/dihe_ind_max.csv')

    def clust_dihedrals(self):
        """
        Determine which dihedral combinations are sampled and cluster conformers

        """
        from itertools import product

        if rank == 0:
            # Step 1: Load dihedral peaks and determine possible options
            dihe_name, dihe_ind, max_values, peak_options = af.input_torsion()
        
            # Step 2: Compute dihedral angles for ligand
            traj = md.load('analysis/center.xtc', top='analysis/md.gro')
            dihedral = md.compute_dihedrals(traj, indices=dihe_ind)
            dihedral = dihedral*(180/np.pi) #Convert to degree

            # Step 3: Determine which dihderal peak is being sampled per frame
            num_dihe = len(dihe_name)
            dihe_peak_sampled = np.zeros((traj.n_frames, num_dihe))
            for t in range(traj.n_frames):
                for i in range(num_dihe):
                    max_value_i = np.array(max_values[i], dtype=float)
                    value = dihedral[t,i]
                    dihe_peak_sampled[t,i] = af.find_nearest(max_value_i, value)

            # Step 4: Clasify conformers from possible combinations
            #Name possible dihedral conformations
            conf = list(product(*peak_options))
            conformer = np.zeros((len(conf), len(conf[0])))
            for c in range(len(conf)):
                conf_c = conf[c]
                conformer[c,:] = conf_c
        
            #Classify dihedrals into conformations
            count = np.zeros(len(conformer))
            frame_select = []
            for t in range(traj.n_frames):
                find_conf = (dihe_peak_sampled[t,:] == conformer).all(axis=1)
                conf_idx = find_conf.nonzero()[0][0]
                if count[conf_idx] == 0:
                    frame_select.append(t)    
                count[conf_idx] += 1
            per = 100*(count/traj.n_frames)

            # Step 5: Print conformer angle combinations, percent ligand is in conformation, and frame in which the ligand is in that conformation
            conformer_list, traj_confs, per_non_zero = af.process_confs(traj, frame_select, per, f'{self.output_name}_dihe', self.conformer_threshold)

            #Cluster conformers
            frames_dihe_clust, per_dihe_clust, group = af.clust_conf(traj_confs, per_non_zero, f'{self.output_name}_dihe_clust')
            cluster_list, traj_clust_confs, clust_per = af.process_confs(traj_confs, frames_dihe_clust, per_dihe_clust, f'{self.output_name}_dihe_clust', self.conformer_threshold, 'Cluster')
            df = pd.DataFrame({'Cluster ID': cluster_list, 'Grouped Confs': group})
            df.to_csv(f'analysis/{self.output_name}_dihe_clust_def.csv')

if __name__ == "__main__":
    # Do something if this file is invoked on its own
    print(canvas())
