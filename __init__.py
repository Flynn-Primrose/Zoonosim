'''
Initialize Zoonosim by importing all modules.

Convention is to use "import zoonosim as zn", and then to use all functions and
classes directly, e.g. zn.Sim() rather than zn.sim.Sim()
'''

# Check that requirements are met and set options
from . import requirements
from .config import options

from .version import __version__, __versiondate__, __license__
if config.options.verbose:
    print(__license__)