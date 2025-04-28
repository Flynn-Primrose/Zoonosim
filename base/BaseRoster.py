'''
Defines the BaseRoster class,  which is a base class for managing a roster of agents.
'''

import sciris as sc
import numpy as np
import pandas as pd

from . import FlexPretty

from .. import defaults as znd

__all__ = ['BaseRoster']

class BaseRoster(FlexPretty):
    '''
    Base class for rosters, which are collections of agents.
    '''

    def set_pars(self, pars=None):
        '''
        Re-link the parameters stored in the roster object to the object containing it.
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

        # Do minimal validation -- needed here since pop_size must exist should be converted to an int when first set
        if 'pop_size' not in pars:
            errormsg = f'The parameter "pop_size" must be included in a population; keys supplied were:\n{sc.newlinejoin(pars.keys())}'
            raise sc.KeyNotFoundError(errormsg)
        pars['pop_size'] = int(pars['pop_size'])
        pars.setdefault('n_variants', 1)
        self.pars = pars # Actually store the pars
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
        ''' Lock the roster object to prevent keys from being added '''
        self._lock = True
        return


    def unlock(self):
        ''' Unlock the roster object to allow keys to be added '''
        self._lock = False
        return


    def __getitem__(self, key):
        ''' Allow roster['attr'] instead of getattr(roster, 'attr')
        '''
        try:
            return self.__dict__[key]
        except: # pragma: no cover
            errormsg = f'Key "{key}" is not a valid attribute of the this roster'
            raise AttributeError(errormsg)


    def __setitem__(self, key, value):
        ''' Ditto '''
        if self._lock and key not in self.__dict__: # pragma: no cover
            errormsg = f'Key "{key}" is not a current attribute of this roster, and the roster object is locked; see roster.unlock()'
            raise AttributeError(errormsg)
        self.__dict__[key] = value
        return


    def __len__(self):
        ''' This is just a scalar, but validate() and _resize_arrays() make sure it's right '''
        #return int(self.pars['pop_size'])
        return len(self.uid) # NOTE: this is probably slower but it's necessary because self.pars['pop_size'] will always equal the size of the total population.


    def __iter__(self):
        ''' Iterate over agents '''
        for i in range(len(self)):
            yield self[i]

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


    def agent_keys(self):
        ''' Returns keys specific to an agent (e.g., their age) '''
        return self.meta.agent[:]


    def state_keys(self):
        ''' Returns keys for different states of an agent (e.g., symptomatic) '''
        return self.meta.states[:]


    def date_keys(self):
        ''' Returns keys for different event dates (e.g., date a person became symptomatic) '''
        return self.meta.dates[:]


    def dur_keys(self):
        ''' Returns keys for different durations (e.g., the duration from exposed to infectious) '''
        return self.meta.durs[:]


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
    
    def __add__(self, roster2):
        ''' Combine two BaseRoster arrays '''
        newroster = sc.dcp(self)
        keys = list(self.keys())
        for key in keys:
            nrval = newroster[key]
            r2val = roster2[key]

            if isinstance(nrval, BaseRoster) or isinstance(r2val, BaseRoster):
                if isinstance(nrval, BaseRoster) and isinstance(r2val, BaseRoster):
                    newroster.set(key, nrval + r2val, die=False) # NOTE: I'm not confident this will work
                else: 
                    errormsg = f'Cannot add a roster to a non-roster object: {key}'
                    raise NotImplementedError(errormsg)
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