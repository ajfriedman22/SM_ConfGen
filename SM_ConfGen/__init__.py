"""Small molecule conformer generation using TREMD"""

# Add imports here
from .sm_confgen import *


from ._version import __version__

from . import _version
__version__ = _version.get_versions()['version']
