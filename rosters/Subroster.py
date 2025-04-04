'''
Defines the Subroster class. A subroster is a list of agents of a single type as well as all values tracked about those agents (e.g. infection status, symptomatic status, age, etc.). 
'''

import sciris as sc
import numpy as np

from ..base import BaseRoster

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
                    errormsg += f'  {k}: agents={v.roster}, roster={v.subroster}'
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


    def __add__(self, roster2):
        ''' Combine two subroster arrays '''
        newroster = sc.dcp(self)
        keys = list(self.keys())
        for key in keys:
            nrval = newroster[key]
            r2val = roster2[key]
            if nrval.ndim == 1:
                newroster.set(key, np.concatenate([nrval, r2val], axis=0), die=False) # Allow size mismatch
            elif nrval.ndim == 2:
                newroster.set(key, np.concatenate([nrval, r2val], axis=1), die=False)
            else:
                errormsg = f'Not sure how to combine arrays of {nrval.ndim} dimensions for {key}'
                raise NotImplementedError(errormsg)

        # Validate
        newroster.pars['pop_size'] += roster2.pars['pop_size']
        newroster.validate()

        # Reassign UIDs so they're unique
        newroster.set('uid', np.arange(len(newroster))) # This is going to be a problem since each roster only holds a subset of the agents. TODO:revisit

        return newroster


    def __radd__(self, roster2):
        ''' Allows sum() to work correctly '''
        if not roster2: return self
        else:           return self.__add__(roster2)






