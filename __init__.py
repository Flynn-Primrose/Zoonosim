'''
Initialize Zoonosim by importing all modules
'''

from .version import __version__, __versiondate__, __license__
if settings.options.verbose:
    print(__license__)