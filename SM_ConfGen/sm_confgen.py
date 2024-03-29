"""Provide the primary functions."""
import os
import sys
import copy
import yaml
import random
import warnings
import subprocess
from SM_ConfGen.utils import utils
from SM_ConfGen.utils.exception import ParameterError
from SM_ConfGen.utils import gmx_parser
from mpi4py import MPI

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
            "mdp",
        ]

        for i in required_args:
            if hasattr(self, i) is False or getattr(self, i) is None:
                raise ParameterError(
                    f"Required parameter '{i}' not specified in {self.yaml}."
                )  # noqa: F405

        # Step 3: Handle the optional YAML parameters
        # Key: Optional argument; Value: Default value
        optional_args = {
            "mdp": ['data/default_mdp/em.mdp', 'data/default_mdp/nvt.mdp', 'data/default_mdp/npt.mdp', 'data/default_mdp/TREMD.mdp'],
            "temp_range": [300, 306, 311.08, 316.23, 321.45, 326.74, 332.10, 337.54, 343.05, 348.63, 354.29, 360.02, 365.84, 371.73, 377.71, 383.76, 389.90, 396.12, 402.43, 408.83, 415.31, 421.88, 428.55, 435.30, 442.14, 449.07, 450],
            "nst_sim": None,
            "dt": None,
            "replex_rate": 1000,
            "grompp_args": None,
            "runtime_args": None,
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
        params_str = ['gmx_executible', 'input_file']
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

        if not self.input_file.lower().endswith(('.pdb', '.sdf')):
            if os.path.exists(self.input_file + '.pdb'):
                self.input_file = self.input_file + '.pdb'
            elif os.path.exists(self.input_file + '.sdf'):
                self.input_file = self.input_file + '.sdf'
            else:
                raise ParameterError(f'{self.input_file} file does not exsist')
        elif not os.path.exists(self.input_file):
            raise ParameterError(f'{self.input_file} file does not exsist')
        
        #Check that temperature range are integers are 
        if not all(isinstance(x, float) for x in self.temp_range) and not all(isinstance(x, int) for x in self.temp_range):
            raise ParameterError('The input temperature range must be of type int or float')
        if self.temp_range[-1] > 500:
            self.warnings.append(f'WARNING: Maximum temperature of {self.temp_range[-1]} is > 500 K which may cause issues solvent')
        self.n_rep = len(self.temp_range)

        #Check simulation parameters
        float_or_int_params = ['nst_sim', 'dt', 'replex_rate']
        for i in float_or_int_params:
            if type(getattr(self, i)) != int and type(getattr(self, i)) != float and type(getattr(self, i)) != None:
                raise ParameterError(f"The parameter '{i}' should be an integer or float.")
        
        #Check that option input args are formatted correctly
        params_dict = ['grompp_args', 'runtime_args']
        for i in params_dict:
            if getattr(self, i) is not None and not isinstance(getattr(self, i), dict):
                raise ParameterError(f"The parameter '{i}' should be a dictionary.")
        
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
        print("Important parameters of EXEE")
        print("============================")
        print(f"Python version: {sys.version}")
        print(f"GROMACS executable: {self.gmx_path}")  # we print the full path here
        print(f"GROMACS version: {self.gmx_version}")
        #print(f"SM ConfGen version: {SM_ConfGen.__version__}")
        print(f'Simulation inputs: {self.input_file}, {self.mdp}')
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
    
    def parameterize_system(self, ):
        """
        Prepares the system files for all replicate simulations

        Parameters
        ----------

        """

    def solvate_system(self):
        """
        Solvate and neutralize charge for the system

        """
        #Create Box
        create_box = []
        returncode, stdout, stderr = utils.run_gmx_cmd(create_box)
        if returncode != 0:
            print(f'Error (return code: {returncode}):\n{stderr}')

        #Solvate Box
        solv_box = []
        returncode, stdout, stderr = utils.run_gmx_cmd(solv_box)
        if returncode != 0:
            print(f'Error (return code: {returncode}):\n{stderr}')
        
        #Add neuralizing ions
        neutral_box = []
        # Add additional arguments if any
        if self.grompp_args is not None:
            # Turn the dictionary into a list with the keys alternating with values
            add_args = [elem for pair in self.grompp_args.items() for elem in pair]
            neutral_box.extend(add_args)
        returncode, stdout, stderr = utils.run_gmx_cmd(neutral_box)
        if returncode != 0:
            print(f'Error (return code: {returncode}):\n{stderr}')

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
            gro = 'prep/conf.gro'
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
            for i in range(self.n_sim):
                arguments = [self.gmx_executable, 'grompp']

                # Input files
                if step == 'NVT':
                    mdp = f'prep/rep_{i}/nvt.mdp'
                    gro = 'prep/min.gro'
                    cpt = None
                    output_base = f'prep/rep_{i}/nvt'
                elif step == 'NPT':
                    mdp = f'prep/rep_{i}/npt.mdp'
                    gro = f'prep/rep_{i}/nvt.gro'
                    cpt = f'prep/rep_{i}/nvt.cpt'
                    output_base = f'prep/rep_{i}/npt'
                else:
                    mdp = f'rep_{i}/md.mdp'
                    gro = f'prep/rep_{i}/npt.gro'
                    cpt = f'prep/rep_{i}/npt.cpt'
                    output_base = f'rep_{i}/md'

                # Add input file arguments
                arguments.extend(['-f', mdp, '-c', gro, '-p', top])
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
            returncode = None  # Initialize as None for all ranks (necessary for the case when -np > n_sim, which is rare)
            if rank == 0:
                print('Generating TPR files ...')
            if rank < self.n_sim:
                returncode, stdout, stderr = utils.run_gmx_cmd(args_list[rank])
                if returncode != 0:
                    print(f'Error on rank {rank} (return code: {returncode}):\n{stderr}')

            # gather return codes at rank 0
            code_list = comm.gather(returncode, root=0)

            if rank == 0:
                # Filter out None values which represent ranks that did not execute the command
                code_list = [code for code in code_list if code is not None]
                if code_list != [0] * self.n_sim:
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
            arguments.extend('min')
        elif step == 'NVT':
            arguments.extend('nvt')
        else:
            arguments.extend('npt')

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
            returncode = None  # Initialize as None for all ranks (necessary for the case when -np > n_sim, which is rare)
            if rank == 0:
                print(f'Running {step} equilibration simulations ...')
            if rank < self.n_sim:
                os.chdir(f'rep_{rank}/')
                returncode, stdout, stderr = utils.run_gmx_cmd(arguments)
                if returncode != 0:
                    print(f'Error on rank {rank} (return code: {returncode}):\n{stderr}')
                if self.rm_cpt is True:
                    # if the simulation went wrong, there would be no checkpoint file
                    try:
                        os.remove('state.cpt')
                    except Exception:
                        print('\n--------------------------------------------------------------------------\n')
                        MPI.COMM_WORLD.Abort(1)
                os.chdir('../../')

            # gather return codes at rank 0
            code_list = comm.gather(returncode, root=0)

            if rank == 0:
                # Filter out None values which represent ranks that did not execute the command
                code_list = [code for code in code_list if code is not None]
                if code_list != [0] * self.n_sim:
                    MPI.COMM_WORLD.Abort(1)   # Doesn't matter what non-zero returncode we put here as the code from GROMACS will be printed before this point anyway.  # noqa: E501
        os.chdir('../')
    
    def run_TREMD(self, restart:bool):
        """
        Run temperature replica exchange molecular dynamics simulations

        Parameters
        ----------
        restart : bool
            Is this simulation restarting from checkpoint?
        """
        arguments = ['mpirun', '-np', self.n_rep, self.gmx_executable, 'mdrun', '-deffnm', 'md', '-multidir']

        #Add multiple directories for replica exchange run
        for i in range(self.n_rep):
            arguments.append(f'rep_{i}')
        
        arguments.append(['-replex', self.replex_rate])

        if self.runtime_args is not None:
            # Turn the dictionary into a list with the keys alternating with values
            add_args = [elem for pair in self.runtime_args.items() for elem in pair]
            arguments.extend(add_args)
        
        if restart:
            arguments.append(['-cpi', 'md'])

        #Run Simulation
        returncode, stdout, stderr = utils.run_gmx_cmd(arguments)
        if returncode != 0:
                print(f'Error (return code: {returncode}):\n{stderr}')

if __name__ == "__main__":
    # Do something if this file is invoked on its own
    print(canvas())
