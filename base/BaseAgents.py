'''
Base object for the Agents object. The Agents object will itself contain a list of all agents as well as all state information that is common to all agent types. 
It will also include a pars object with the parameters needed to generate the population as well as a specific roster for each agent type. Each type specific roster 
will contain a list of all agents of that type as well as all state information that is specific to that data type.
'''


import sciris as sc
import numpy as np
import pandas as pd

from . import FlexPretty
from .. import defaults as znd

class BaseAgents(FlexPretty):
    '''
    A class to handle all the boilerplate for the Agent object -- note that as with the
    BaseSim vs Sim classes, everything interesting happens in the derived classes,
    whereas this class exists to handle the less interesting implementation details.
    '''

    def set_pars(self, pars=None):
        '''
        Re-link the parameters stored in the people object to the sim containing it,
        and perform some basic validation.
        '''
        orig_pars = self.__dict__.get('pars') # Get the current parameters using dict's get method
        if pars is None:
            if orig_pars is not None: # If it has existing parameters, use them
                pars = orig_pars
            else:
                pars = {}
        elif sc.isnumber(pars): # Interpret as the number of farms in the model (needed to generate population)
            pars = {'n_farms':pars} # Ensure it's a dictionary

        # Copy from old parameters to new parameters
        if isinstance(orig_pars, dict):
            for k,v in orig_pars.items():
                if k not in pars:
                    pars[k] = v

        # Do minimal validation -- needed here since pop_size should be converted to an int when first set
        if 'n_farms' not in pars:
            errormsg = f'The parameter "n_farms" must be included in a population; keys supplied were:\n{sc.newlinejoin(pars.keys())}'
            raise sc.KeyNotFoundError(errormsg)
        pars['n_farms'] = int(pars['n_farms'])
        pars.setdefault('n_variants', 1)
        # pars.setdefault('location', None) # Currently there is no functionality related to different locations.
        self.pars = pars # Actually store the pars
        return
    
    def validate(self, sim_pars=None, die=True, verbose=False):
        '''
        Perform validation on the Agents object.

        Args:
            sim_pars (dict): dictionary of parameters from the sim to ensure they match the current People object
            die (bool): whether to raise an exception if validation fails
            verbose (bool): detail to print
        '''

        # Check that parameters match
        if sim_pars is not None:
            mismatches = {}
            keys = ['n_farms', 'n_variants'] # These are the keys used in generating the population
            for key in keys:
                sim_v = sim_pars.get(key)
                agents_v = self.pars.get(key)
                if sim_v is not None and agents_v is not None:
                    if sim_v != agents_v:
                        mismatches[key] = sc.objdict(sim=sim_v, agents = agents_v)
            if len(mismatches):
                errormsg = 'Validation failed due to the following mismatches between the sim and the agent parameters:\n'
                for k,v in mismatches.items():
                    errormsg += f'  {k}: sim={v.sim}, agents={v.agents}'
                raise ValueError(errormsg)

        # Check that the keys match
        contact_layer_keys = set(self.contacts.keys())
        layer_keys    = set(self.layer_keys())
        if contact_layer_keys != layer_keys:
            errormsg = f'Parameters layers {layer_keys} are not consistent with contact layers {contact_layer_keys}'
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

        # Check that the layers are valid
        for layer in self.contacts.values():
            layer.validate()

        return
