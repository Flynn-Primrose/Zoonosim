'''
Initialize Zoonosim by importing all modules.

Convention is to use "import zoonosim as zn", and then to use all functions and
classes directly, e.g. zn.Sim() rather than zn.sim.Sim()
'''

# Check that requirements are met and set options
from . import requirements
from . import defaults as znd

from .version import __version__, __versiondate__, __license__
if znd.verbose:
    print(__license__)