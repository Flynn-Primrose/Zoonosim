'''
Defines the BaseRoster class. A roster is a list of agents of a single type as well as all values tracked about those agents (e.g. infection status, symptomatic status, age, etc.). 
'''

import sciris as sc
import numpy as np
import pandas as pd

from . import FlexPretty

from .. import defaults as znd

__all__= ['BaseRoster']

class BaseRoster(FlexPretty):
    '''
    A class to handle all the boilerplate for rosters -- note that as with the
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

        # Copy from old parameters to new parameters
        if isinstance(orig_pars, dict):
            for k,v in orig_pars.items():
                if k not in pars:
                    pars[k] = v

        # Do minimal validation -- needed here since pop_size should be converted to an int when first set
        if 'agent_type' not in pars:
            errormsg = f'The parameter "agent_type" most be included in every roster object; keys supplied were: \n{sc.newlinejoin(pars.keys())}'
        elif 'pop_size' not in pars:
            errormsg = f'The parameter "pop_size" must be included in a population; keys supplied were:\n{sc.newlinejoin(pars.keys())}'
            raise sc.KeyNotFoundError(errormsg)
        pars['pop_size'] = int(pars['pop_size'])
        pars.setdefault('n_variants', 1)
        self.pars = pars # Actually store the pars
        return


    def validate(self, agents_pars=None, die=True, verbose=False):
        '''
        Perform validation on the roster object.

        Args:
            agents_pars (dict): dictionary of parameters from the agents object to ensure they match the current roster object
            die (bool): whether to raise an exception if validation fails
            verbose (bool): detail to print
        '''

        # Check that parameters match
        if agents_pars is not None:
            mismatches = {}
            keys = ['n_variants'] # These are the keys that must match the agents object
            for key in keys:
                agents_v = agents_pars.get(key)
                rstr_v = self.pars.get(key)
                if agents_v is not None and rstr_v is not None:
                    if agents_v != rstr_v:
                        mismatches[key] = sc.objdict(agents=agents_v, roster=rstr_v)
            if len(mismatches):
                errormsg = 'Validation failed due to the following mismatches between the agent object and the roster parameters:\n'
                for k,v in mismatches.items():
                    errormsg += f'  {k}: agents={v.agents}, roster={v.roster}'
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


    def _resize_arrays(self, new_size=None, keys=None):
        ''' Resize arrays if any mismatches are found '''

        # Handle None or tuple input (representing variants and pop_size)
        if new_size is None:
            new_size = len(self)
        pop_size = new_size if not isinstance(new_size, tuple) else new_size[1]
        self.pars['pop_size'] = pop_size

        # Reset sizes
        if keys is None:
            keys = self.keys()
        keys = sc.promotetolist(keys)
        for key in keys:
            self[key].resize(new_size, refcheck=False) # Don't worry about cross-references to the arrays

        return


    def lock(self):
        ''' Lock the people object to prevent keys from being added '''
        self._lock = True
        return


    def unlock(self):
        ''' Unlock the people object to allow keys to be added '''
        self._lock = False
        return


    def __getitem__(self, key):
        ''' Allow people['attr'] instead of getattr(people, 'attr')
        '''
        try:
            return self.__dict__[key]
        except: # pragma: no cover
            errormsg = f'Key "{key}" is not a valid attribute of the {self.pars['agent_type']} roster'
            raise AttributeError(errormsg)


    def __setitem__(self, key, value):
        ''' Ditto '''
        if self._lock and key not in self.__dict__: # pragma: no cover
            errormsg = f'Key "{key}" is not a current attribute of {self.pars['agent_type']} roster, and the roster object is locked; see roster.unlock()'
            raise AttributeError(errormsg)
        self.__dict__[key] = value
        return


    def __len__(self):
        ''' This is just a scalar, but validate() and _resize_arrays() make sure it's right '''
        return int(self.pars['pop_size'])


    def __iter__(self):
        ''' Iterate over agents '''
        for i in range(len(self)):
            yield self[i]


    def __add__(self, roster2):
        ''' Combine two roster arrays '''
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


    def _brief(self):
        '''
        Return a one-line description of the roster -- used internally and by repr();
        see roster.brief() for the user version.
        '''
        try:
            layerstr = ', '.join([str(k) for k in self.layer_keys()])
            string   = f'{self.pars['agent_type']}(n={len(self):0n}; layers: {layerstr})'
        except Exception as E: # pragma: no cover
            string = sc.objectid(self)
            string += f'Warning, multisim appears to be malformed:\n{str(E)}'
        return string


    def summarize(self, output=False):
        ''' Print a summary of the roster -- same as brief '''
        return self.brief(output=output)


    def set(self, key, value, die=True):
        ''' Ensure sizes and dtypes match '''
        current = self[key]
        value = np.array(value, dtype=self._dtypes[key]) # Ensure it's the right type
        if die and len(value) != len(current): # pragma: no cover
            errormsg = f'Length of new array does not match current ({len(value)} vs. {len(current)})'
            raise IndexError(errormsg)
        self[key] = value
        return


    def get(self, key):
        ''' Convenience method -- key can be string or list of strings '''
        if isinstance(key, str):
            return self[key]
        elif isinstance(key, list):
            arr = np.zeros((len(self), len(key)))
            for k,ky in enumerate(key):
                arr[:,k] = self[ky]
            return arr


    def true(self, key):
        ''' Return indices matching the condition '''
        return self[key].nonzero()[0]


    def false(self, key):
        ''' Return indices not matching the condition '''
        return (~self[key]).nonzero()[0]


    def defined(self, key):
        ''' Return indices of agents who are not-nan '''
        return (~np.isnan(self[key])).nonzero()[0]


    def undefined(self, key):
        ''' Return indices of agents who are nan '''
        return np.isnan(self[key]).nonzero()[0]


    def count(self, key):
        ''' Count the number of agents for a given key '''
        return np.count_nonzero(self[key])


    def count_by_variant(self, key, variant):
        ''' Count the number of agents for a given key '''
        return np.count_nonzero(self[key][variant,:])

    def r_count(self, key, i_start, i_end):
        ''' Count the number of agents for a given key in the half-open
        interval [i_start, i_end), for the given range of agent id's. '''
        return np.count_nonzero(self[key][i_start:i_end])


    def r_count_by_variant(self, key, variant, i_start, i_end):
        ''' Count the number of agents for a given key. In a range, as above. '''
        return np.count_nonzero(self[key][variant,i_start:i_end])


    def count_not(self, key):
        ''' Count the number of agents who do not have a property for a given key '''
        return len(self[key]) - self.count(key)


    def keys(self):
        ''' Returns keys for all properties of the roster object '''
        return self.meta.all_states[:]


    def person_keys(self):
        ''' Returns keys specific to an agent (e.g., their age) '''
        return self.meta.specific[:]


    def state_keys(self):
        ''' Returns keys for different states of a person (e.g., symptomatic) '''
        return self.meta.states[:]


    def date_keys(self):
        ''' Returns keys for different event dates (e.g., date a person became symptomatic) '''
        return self.meta.dates[:]


    def dur_keys(self):
        ''' Returns keys for different durations (e.g., the duration from exposed to infectious) '''
        return self.meta.durs[:]


    def layer_keys(self):
        ''' Get the available contact keys -- try contacts first, then beta_layer '''
        try:
            keys = list(self.contacts.keys())
        except: # If not fully initialized
            try:
                keys = list(self.pars['beta_layer'].keys())
            except:  # pragma: no cover # If not even partially initialized
                keys = []
        return keys


    def indices(self):
        ''' The indices of each agent array '''
        return np.arange(len(self))


    def to_df(self):
        ''' Convert to a Pandas dataframe '''
        df = pd.DataFrame.from_dict({key:self[key] for key in self.keys()})
        return df


    def to_arr(self):
        ''' Return as numpy array '''
        arr = np.empty((len(self), len(self.keys())), dtype=znd.default_float)
        for k,key in enumerate(self.keys()):
            if key == 'uid':
                arr[:,k] = np.arange(len(self))
            else:
                arr[:,k] = self[key]
        return arr


    def to_list(self):
        ''' Return all agents as a list '''
        return list(self)


    def from_list(self, agents, resize=True):
        ''' Convert a list of agents back into a Roster object '''

        # Handle population size
        pop_size = len(agents)
        if resize:
            self._resize_arrays(new_size=pop_size)

        # Iterate over people -- slow!
        for a,agent in enumerate(agents):
            for key in self.keys():
                self[key][a] = getattr(agent, key)

        return