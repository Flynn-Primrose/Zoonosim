'''
Defines the Subroster class. A subroster is a list of agents of a single type as well as all values tracked about those agents (e.g. infection status, symptomatic status, age, etc.). 
'''

import sciris as sc
import numpy as np

from ..base import BaseRoster
from .. import utils as znu

__all__= ['Subroster']

class Subroster(BaseRoster):
    '''
    A class to handle all the boilerplate for subrosters -- note that as with the
    BaseSim vs Sim classes, everything interesting happens in the derived classes,
    whereas this class exists to handle the less interesting implementation details.
    '''



    def validate(self, roster_pars=None, die=True, verbose=False):
        '''
        Perform validation on the roster object.

        Args:
            agents_pars (dict): dictionary of parameters from the agents object to ensure they match the current roster object
            die (bool): whether to raise an exception if validation fails
            verbose (bool): detail to print
        '''

        # Check that parameters match
        if roster_pars is not None:
            mismatches = {}
            keys = ['n_variants'] # These are the keys that must match the agents object
            for key in keys:
                roster_v = roster_pars.get(key)
                subroster_v = self.pars.get(key)
                if roster_v is not None and subroster_v is not None:
                    if roster_v != subroster_v:
                        mismatches[key] = sc.objdict(roster=roster_v, subroster=subroster_v)
            if len(mismatches):
                errormsg = 'Validation failed due to the following mismatches between the agent object and the subroster parameters:\n'
                for k,v in mismatches.items():
                    errormsg += f'  {k}: roster={v.roster}, subroster={v.subroster}'
                raise ValueError(errormsg)

        # Check that the length of each array is consistent
        expected_len = len(self)
        expected_variants = self.pars['n_variants']
        for key in self.keys():
            if self[key].ndim == 1:
                actual_len = len(self[key])
            else: # If it's 2D, variants need to be checked separately
                actual_variants, actual_len = self[key].shape
                if actual_variants != expected_variants:
                    if verbose:
                        print(f'Resizing "{key}" from {actual_variants} to {expected_variants}')
                    self._resize_arrays(keys=key, new_size=(expected_variants, expected_len))
            if actual_len != expected_len: # pragma: no cover
                if die:
                    errormsg = f'Length of key "{key}" did not match population size ({actual_len} vs. {expected_len})'
                    raise IndexError(errormsg)
                else:
                    if verbose:
                        print(f'Resizing "{key}" from {actual_len} to {expected_len}')
                    self._resize_arrays(keys=key)

        return
    
    def check_inds(self, current, date, filter_inds=None):
        ''' Return indices for which the current state is false and which are assigned a date on or before the current date

        Args:
            current (array): list of boolean values that represent a current state
            date (array): list that contains either a date or a Nan
        '''
        if filter_inds is None:
            not_current = znu.false(current)
        else:
            not_current = znu.ifalsei(current, filter_inds)
        has_date = znu.idefinedi(date, not_current)
        inds     = znu.itrue(self.t >= date[has_date], has_date)
        return inds









