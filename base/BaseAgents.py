'''
Base object for the Agents object. The Agents object will itself contain a list of all agents as well as all state information that is common to all agent types. 
It will also include a pars object with the parameters needed to generate the population as well as a specific roster for each agent type. Each type specific roster 
will contain a list of all agents of that type as well as all state information that is specific to that data type.
'''


import sciris as sc
import numpy as np
import pandas as pd

from . import FlexPretty
from . import Contacts
from . import Layer

from .. import defaults as znd
from .. import misc as znm

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
                errormsg = 'Validation failed due to the following mismatches between the sim and the agents parameters:\n'
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
    def to_graph(self): # pragma: no cover
        '''
        Convert all agents to a networkx MultiDiGraph, including all properties of
        the agents (nodes) and contacts (edges).

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
        run or partially run Agents object, since the changes that happen during
        a run are usually irreversible.

        Args:
            filename (str or None): the name or path of the file to save to; if None, uses stored
            force (bool): whether to allow saving even of a run or partially-run People object
            kwargs: passed to ``sc.makefilepath()``

        Returns:
            filename (str): the validated absolute path to the saved file

        **Example**::

            sim = zn.Sim()
            sim.initialize()
            sim.agents.save() # Saves to a .agnts file
        '''

        # Check if we're trying to save an already run PeopleAgents object
        if self.t > 0 and not force:
            errormsg = f'''
The Agents object has already been run (t = {self.t}), which is usually not the
correct state to save it in since it cannot be re-initialized. If this is intentional,
use sim.people.save(force=True). Otherwise, the correct approach is:

    sim = zn.Sim(...)
    sim.initialize() # Create the Agents object but do not run
    sim.agents.save() # Save Agents immediately after initialization
    sim.run() #
'''
            raise RuntimeError(errormsg)

        # Handle the filename
        if filename is None:
            filename = 'zoonosim.agnts'
        filename = sc.makefilepath(filename=filename, **kwargs)
        self.filename = filename # Store the actual saved filename
        znm.save(filename=filename, obj=self)

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
        agents = znm.load(filename, *args, **kwargs)
        if not isinstance(agents, BaseAgents): # pragma: no cover
            errormsg = f'Cannot load object of {type(agents)} as an Agents object'
            raise TypeError(errormsg)
        return agents



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
