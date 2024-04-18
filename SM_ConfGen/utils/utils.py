import subprocess
import sys
from mpi4py import MPI

class Logger:
    """
    Redirects the STDOUT and STDERR to a specified output file while preserving them on screen.

    Parameters
    ----------
    logfile : str
        Name of the output file to write the logged messages.

    Attributes
    ----------
    terminal : file object
        The file object that represents the original STDOUT (i.e., the screen).
    log : file object
        The file object that represents the logfile where messages will be written.
    """

    def __init__(self, logfile):
        """
        Initializes a Logger instance.

        Parameters
        ----------
        logfile : str
            Name of the output file to write the logged messages.
        """
        self.terminal = sys.stdout
        self.log = open(logfile, "a")

    def write(self, message):
        """
        Writes the given message to both the STDOUT and the logfile.

        Parameters
        ----------
        message : str
            The message to be written to STDOUT and logfile.
        """
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        """
        This method is needed for Python 3 compatibility. This handles the flush command by doing nothing.
        You might want to specify some extra behavior here.
        """
        # self.terminal.log()
        pass

def run_gmx_cmd(arguments, prompt_input=None):
    """
    Runs a GROMACS command as a subprocess

    Parameters
    ----------
    arguments : list
        A list of arguments that compose of the GROMACS command to run, e.g.
        :code:`['gmx', 'mdrun', '-deffnm', 'sys']`.
    prompt_input : str or None
        The input to be passed to the GROMACS command when it prompts for input.

    Returns
    -------
    return_code : int
        The exit code of the GROMACS command. Any number other than 0 indicates an error.
    stdout : str or None
        The STDOUT of the process.
    stderr: str or None
        The STDERR or the process.

    """
    try:
        result = subprocess.run(arguments, capture_output=True, text=True, input=prompt_input, check=True)
        return_code, stdout, stderr = result.returncode, result.stdout, None
    except subprocess.CalledProcessError as e:
        return_code, stdout, stderr = e.returncode, None, e.stderr

    return return_code, stdout, stderr

def _autoconvert(s):
    """
    Converts input to a numerical type if possible. Used for the MDP parser.
    Modified from `utilities.py in GromacsWrapper <https://github.com/Becksteinlab/GromacsWrapper>`_.
    Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>

    Parameters
    ----------
    s : str or any
        The input value to be converted to a numerical type if possible. If :code:`s` is not a string,
        it is returned as is.

    Returns
    -------
    numerical : int, float, numpy.ndarray, or any
        The converted numerical value. If :code:`s` can be converted to a single numerical value,
        that value is returned as an :code:`int` or :code:`float`. If :code:`s` can be converted to
        multiple numerical values, a :code:`numpy.ndarray` containing those values is returned.
        If :code:`s` cannot be converted to a numerical value, :code:`s` is returned as is.

    Raises
    ------
    ValueError
        If :code:`s` cannot be converted to a numerical value.
    """
    if type(s) is not str:
        return s
    for converter in int, float, str:  # try them in increasing order of lenience
        try:
            s = [converter(i) for i in s.split()]
            if len(s) == 1:
                return s[0]
            else:
                return s
                """
                if len(s) != 0 and type(s[0]) == str:
                    # For the case like pull_coord1_dim = Y Y Y
                    return s
                else:
                    return np.array(s)
                """
        except (ValueError, AttributeError):
            pass
    raise ValueError("Failed to autoconvert {0!r}".format(s))
