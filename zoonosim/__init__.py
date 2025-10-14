'''
Initialize zoonosim by importing all modules.

Convention is to use "import zoonosim as zn", and then to use all functions and
classes directly, e.g. zn.Sim() rather than zn.sim.Sim()
'''

# Check that requirements are met and set options
from . import requirements
from . import options as zno

from .version import __version__, __versiondate__, __license__
if zno.options.verbose:
    print(__license__)


from .defaults import *
from .misc import *
from .utils import *
from .parameters import *
from .base import *
from .rosters.rosters import *
from .rosters.humans import *
from .rosters.flocks import *
from .rosters.barns import *
from .rosters.waters import *
from .interventions import *
from .analysis import *
from .immunity import *
from .population import *
from .testing import *
from .sim import *
from .run import *
