import subprocess

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
