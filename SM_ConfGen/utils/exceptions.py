class ParameterError(Exception):
    """Error raised when detecting improperly specified parameters in the YAML file."""


class ParseError(Exception):
    """Error raised when parsing of a file failed. Modified from GromacsWrapper."""