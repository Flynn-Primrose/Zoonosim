'''
Defines the BaseRoster class. A roster is a list of agents of a single species as well as all values tracked about those agents (e.g. infection status, symptomatic status, age, etc.). 
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
        elif sc.isnumber(pars): # Interpret as a population size
            pars = {'pop_size':pars} # Ensure it's a dictionary

        # Copy from old parameters to new parameters
        if isinstance(orig_pars, dict):
            for k,v in orig_pars.items():
                if k not in pars:
                    pars[k] = v

        # Do minimal validation -- needed here since pop_size should be converted to an int when first set
        if 'pop_size' not in pars:
            errormsg = f'The parameter "pop_size" must be included in a population; keys supplied were:\n{sc.newlinejoin(pars.keys())}'
            raise sc.KeyNotFoundError(errormsg)
        pars['pop_size'] = int(pars['pop_size'])
        pars.setdefault('n_variants', 1)
        pars.setdefault('location', None)
        self.pars = pars # Actually store the pars
        return


    def validate(self, sim_pars=None, die=True, verbose=False):
        '''
        Perform validation on the People object.

        Args:
            sim_pars (dict): dictionary of parameters from the sim to ensure they match the current People object
            die (bool): whether to raise an exception if validation fails
            verbose (bool): detail to print
        '''

        # Check that parameters match
        if sim_pars is not None:
            mismatches = {}
            keys = ['pop_size', 'pop_type', 'location', 'pop_infected', 'frac_susceptible', 'n_variants'] # These are the keys used in generating the population
            for key in keys:
                sim_v = sim_pars.get(key)
                ppl_v = self.pars.get(key)
                if sim_v is not None and ppl_v is not None:
                    if sim_v != ppl_v:
                        mismatches[key] = sc.objdict(sim=sim_v, people=ppl_v)
            if len(mismatches):
                errormsg = 'Validation failed due to the following mismatches between the sim and the people parameters:\n'
                for k,v in mismatches.items():
                    errormsg += f'  {k}: sim={v.sim}, people={v.people}'
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
            If the key is an integer, alias `people.person()` to return a `Person` instance
        '''
        try:
            return self.__dict__[key]
        except: # pragma: no cover
            if isinstance(key, int):
                return self.person(key)
            else:
                errormsg = f'Key "{key}" is not a valid attribute of people'
                raise AttributeError(errormsg)


    def __setitem__(self, key, value):
        ''' Ditto '''
        if self._lock and key not in self.__dict__: # pragma: no cover
            errormsg = f'Key "{key}" is not a current attribute of people, and the people object is locked; see people.unlock()'
            raise AttributeError(errormsg)
        self.__dict__[key] = value
        return


    def __len__(self):
        ''' This is just a scalar, but validate() and _resize_arrays() make sure it's right '''
        return int(self.pars['pop_size'])


    def __iter__(self):
        ''' Iterate over people '''
        for i in range(len(self)):
            yield self[i]


    def __add__(self, people2):
        ''' Combine two people arrays '''
        newpeople = sc.dcp(self)
        keys = list(self.keys())
        for key in keys:
            npval = newpeople[key]
            p2val = people2[key]
            if npval.ndim == 1:
                newpeople.set(key, np.concatenate([npval, p2val], axis=0), die=False) # Allow size mismatch
            elif npval.ndim == 2:
                newpeople.set(key, np.concatenate([npval, p2val], axis=1), die=False)
            else:
                errormsg = f'Not sure how to combine arrays of {npval.ndim} dimensions for {key}'
                raise NotImplementedError(errormsg)

        # Validate
        newpeople.pars['pop_size'] += people2.pars['pop_size']
        newpeople.validate()

        # Reassign UIDs so they're unique
        newpeople.set('uid', np.arange(len(newpeople)))

        return newpeople


    def __radd__(self, people2):
        ''' Allows sum() to work correctly '''
        if not people2: return self
        else:           return self.__add__(people2)


    def _brief(self):
        '''
        Return a one-line description of the people -- used internally and by repr();
        see people.brief() for the user version.
        '''
        try:
            layerstr = ', '.join([str(k) for k in self.layer_keys()])
            string   = f'People(n={len(self):0n}; layers: {layerstr})'
        except Exception as E: # pragma: no cover
            string = sc.objectid(self)
            string += f'Warning, multisim appears to be malformed:\n{str(E)}'
        return string


    def summarize(self, output=False):
        ''' Print a summary of the people -- same as brief '''
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
        ''' Return indices of people who are not-nan '''
        return (~np.isnan(self[key])).nonzero()[0]


    def undefined(self, key):
        ''' Return indices of people who are nan '''
        return np.isnan(self[key]).nonzero()[0]


    def count(self, key):
        ''' Count the number of people for a given key '''
        return np.count_nonzero(self[key])


    def count_by_variant(self, key, variant):
        ''' Count the number of people for a given key '''
        return np.count_nonzero(self[key][variant,:])

    def r_count(self, key, i_start, i_end):
        ''' Count the number of people for a given key in the half-open
        interval [i_start, i_end), for the given range of people id's. '''
        return np.count_nonzero(self[key][i_start:i_end])


    def r_count_by_variant(self, key, variant, i_start, i_end):
        ''' Count the number of people for a given key. In a range, as above. '''
        return np.count_nonzero(self[key][variant,i_start:i_end])


    def count_not(self, key):
        ''' Count the number of people who do not have a property for a given key '''
        return len(self[key]) - self.count(key)


    def keys(self):
        ''' Returns keys for all properties of the people object '''
        return self.meta.all_states[:]


    def person_keys(self):
        ''' Returns keys specific to a person (e.g., their age) '''
        return self.meta.person[:]


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
        ''' The indices of each people array '''
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


    def person(self, ind):
        ''' Method to create person from the people '''
        p = Person()
        for key in self.meta.all_states:
            data = self[key]
            if data.ndim == 1:
                val = data[ind]
            elif data.ndim == 2:
                val = data[:,ind]
            else:
                errormsg = f'Cannot extract data from {key}: unexpected dimensionality ({data.ndim})'
                raise ValueError(errormsg)
            setattr(p, key, val)

        contacts = {}
        for lkey, layer in self.contacts.items():
            contacts[lkey] = layer.find_contacts(ind)
        p.contacts = contacts

        return p


    def to_list(self):
        ''' Return all people as a list '''
        return list(self)


    def from_list(self, people, resize=True):
        ''' Convert a list of people back into a People object '''

        # Handle population size
        pop_size = len(people)
        if resize:
            self._resize_arrays(new_size=pop_size)

        # Iterate over people -- slow!
        for p,person in enumerate(people):
            for key in self.keys():
                self[key][p] = getattr(person, key)

        return


    def to_graph(self): # pragma: no cover
        '''
        Convert all people to a networkx MultiDiGraph, including all properties of
        the people (nodes) and contacts (edges).

        **Example**::

            import covasim as cv
            import networkx as nx
            sim = cv.Sim(pop_size=50, pop_type='hybrid', contacts=dict(h=3, s=10, w=10, c=5)).run()
            G = sim.people.to_graph()
            nodes = G.nodes(data=True)
            edges = G.edges(keys=True)
            node_colors = [n['age'] for i,n in nodes]
            layer_map = dict(h='#37b', s='#e11', w='#4a4', c='#a49')
            edge_colors = [layer_map[G[i][j][k]['layer']] for i,j,k in edges]
            edge_weights = [G[i][j][k]['beta']*5 for i,j,k in edges]
            nx.draw(G, node_color=node_colors, edge_color=edge_colors, width=edge_weights, alpha=0.5)
        '''
        import networkx as nx

        # Copy data from people into graph
        G = self.contacts.to_graph()
        for key in self.keys():
            data = {k:v for k,v in enumerate(self[key])}
            nx.set_node_attributes(G, data, name=key)

        # Include global layer weights
        for u,v,k in G.edges(keys=True):
            edge = G[u][v][k]
            edge['beta'] *= self.pars['beta_layer'][edge['layer']]

        return G


    def save(self, filename=None, force=False, **kwargs):
        '''
        Save to disk as a gzipped pickle.

        Note: by default this function raises an exception if trying to save a
        run or partially run People object, since the changes that happen during
        a run are usually irreversible.

        Args:
            filename (str or None): the name or path of the file to save to; if None, uses stored
            force (bool): whether to allow saving even of a run or partially-run People object
            kwargs: passed to ``sc.makefilepath()``

        Returns:
            filename (str): the validated absolute path to the saved file

        **Example**::

            sim = cv.Sim()
            sim.initialize()
            sim.people.save() # Saves to a .ppl file
        '''

        # Check if we're trying to save an already run People object
        if self.t > 0 and not force:
            errormsg = f'''
The People object has already been run (t = {self.t}), which is usually not the
correct state to save it in since it cannot be re-initialized. If this is intentional,
use sim.people.save(force=True). Otherwise, the correct approach is:

    sim = cv.Sim(...)
    sim.initialize() # Create the people object but do not run
    sim.people.save() # Save people immediately after initialization
    sim.run() # The People object is
'''
            raise RuntimeError(errormsg)

        # Handle the filename
        if filename is None:
            filename = 'covasim.ppl'
        filename = sc.makefilepath(filename=filename, **kwargs)
        self.filename = filename # Store the actual saved filename
        cvm.save(filename=filename, obj=self)

        return filename


    @staticmethod
    def load(filename, *args, **kwargs):
        '''
        Load from disk from a gzipped pickle.

        Args:
            filename (str): the name or path of the file to load from
            args (list): passed to ``cv.load()``
            kwargs (dict): passed to ``cv.load()``

        Returns:
            people (People): the loaded people object

        **Example**::

            people = cv.people.load('my-people.ppl')
        '''
        people = cvm.load(filename, *args, **kwargs)
        if not isinstance(people, BasePeople): # pragma: no cover
            errormsg = f'Cannot load object of {type(people)} as a People object'
            raise TypeError(errormsg)
        return people



    def init_contacts(self, reset=False):
        ''' Initialize the contacts dataframe with the correct columns and data types '''

        # Create the contacts dictionary
        contacts = Contacts(layer_keys=self.layer_keys())

        if self.contacts is None or reset: # Reset all
            self.contacts = contacts
        else: # Only replace specified keys
            for key,layer in contacts.items():
                self.contacts[key] = layer
        return


    def add_contacts(self, contacts, lkey=None, beta=None):
        '''
        Add new contacts to the array. See also contacts.add_layer().
        '''

        # Validate the supplied contacts
        if isinstance(contacts, (Contacts, dict)): # If it's a Contacts object or a dict, we can use it directly
            if isinstance(contacts, dict) and lkey is not None: # Edge case: a dict for a single layer has been provided
                new_contacts = {}
                new_contacts[lkey] = contacts
            else:
                if 'p1' in contacts: # Avoid the mistake of passing a single layer
                    errormsg = 'To supply a single layer as a dict, you must supply an lkey as well'
                    raise ValueError(errormsg)
                new_contacts = contacts # Main use case
        elif isinstance(contacts, Layer):
            if lkey is None: # If no layer key is supplied, use the first layer
                lkey = self.layer_keys()[0]
            new_contacts = {}
            new_contacts[lkey] = contacts
        elif isinstance(contacts, list): # Assume it's a list of contacts by person, not an edgelist
            new_contacts = self.make_edgelist(contacts) # Assume contains key info
        else: # pragma: no cover
            errormsg = f'Cannot understand contacts of type {type(contacts)}; expecting dataframe, array, or dict'
            raise TypeError(errormsg)

        # Ensure the columns are right and add values if supplied
        for lkey, new_layer in new_contacts.items():
            n = len(new_layer['p1'])
            if 'beta' not in new_layer.keys() or len(new_layer['beta']) != n:
                if beta is None:
                    beta = 1.0
                beta = znd.default_float(beta)
                new_layer['beta'] = np.ones(n, dtype=znd.default_float)*beta

            # Create the layer if it doesn't yet exist
            if lkey not in self.contacts:
                self.contacts[lkey] = Layer(label=lkey)

            # Actually include them, and update properties if supplied
            for col in self.contacts[lkey].keys(): # Loop over the supplied columns
                self.contacts[lkey][col] = np.concatenate([self.contacts[lkey][col], new_layer[col]])
            self.contacts[lkey].validate()

        return


    def make_edgelist(self, contacts):
        '''
        Parse a list of people with a list of contacts per person and turn it
        into an edge list.
        '''

        # Handle layer keys
        lkeys = self.layer_keys()
        if len(contacts):
            contact_keys = contacts[0].keys() # Pull out the keys of this contact list
            lkeys += [key for key in contact_keys if key not in lkeys] # Extend the layer keys

        # Initialize the new contacts
        new_contacts = Contacts(layer_keys=lkeys)
        for lkey in lkeys:
            new_contacts[lkey]['p1']    = [] # Person 1 of the contact pair
            new_contacts[lkey]['p2']    = [] # Person 2 of the contact pair

        # Populate the new contacts
        for p,cdict in enumerate(contacts):
            for lkey,p_contacts in cdict.items():
                n = len(p_contacts) # Number of contacts
                new_contacts[lkey]['p1'].extend([p]*n) # e.g. [4, 4, 4, 4]
                new_contacts[lkey]['p2'].extend(p_contacts) # e.g. [243, 4538, 7,19]

        # Turn into a dataframe
        for lkey in lkeys:
            new_layer = Layer(label=lkey)
            for ckey,value in new_contacts[lkey].items():
                new_layer[ckey] = np.array(value, dtype=new_layer.meta[ckey])
            new_contacts[lkey] = new_layer

        return new_contacts


    @staticmethod
    def remove_duplicates(df):
        ''' Sort the dataframe and remove duplicates -- note, not extensively tested '''
        p1 = df[['p1', 'p2']].values.min(1) # Reassign p1 to be the lower-valued of the two contacts
        p2 = df[['p1', 'p2']].values.max(1) # Reassign p2 to be the higher-valued of the two contacts
        df['p1'] = p1
        df['p2'] = p2
        df.sort_values(['p1', 'p2'], inplace=True) # Sort by p1, then by p2
        df.drop_duplicates(['p1', 'p2'], inplace=True) # Remove duplicates
        df = df[df['p1'] != df['p2']] # Remove self connections
        df.reset_index(inplace=True, drop=True)
        return df