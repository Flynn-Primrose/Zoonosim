# Miscellaneous functions used exclusively for the base module

import sciris as sc

from .. import version as znv
from .. import misc as znm

__all__ = ['set_metadata']

def set_metadata(obj, **kwargs):
    ''' Set standard metadata for an object '''
    obj.created = kwargs.get('created', sc.now())
    obj.version = kwargs.get('version', znv.__version__) 
    obj.git_info = kwargs.get('git_info', znm.git_info())
    return