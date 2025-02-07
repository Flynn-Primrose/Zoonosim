import sciris as sc
import numpy as np
import pandas as pd

from . import FlexDict
from .. import defaults as znd # TODO: is there a better way to organize default values?
from .. import utils as znu # TODO: is there a better way to organize utility functions?

__all__ = ['Layer']

class Layer(FlexDict):
    '''
    A small class holding a single layer of contact edges (connections) between agents.

    The input is typically three arrays: agent 1 of the connection, agent 2 of
    the connection, and the weight of the connection. Connections are undirected;
    each agent is both a source and sink.

    This class is usually not invoked directly by the user, but instead is called
    as part of the population creation.

    Args:
        p1 (array): an array of N connections, representing people on one side of the connection
        p2 (array): an array of people on the other side of the connection
        beta (array): an array of weights for each connection
        label (str): the name of the layer (optional)
        kwargs (dict): other keys copied directly into the layer

    Note that all arguments (except for label) must be arrays of the same length,
    although not all have to be supplied at the time of creation (they must all
    be the same at the time of initialization, though, or else validation will fail).

    **Examples**::

        # Generate an average of 10 contacts for 1000 people
        n = 10_000
        n_people = 1000
        p1 = np.random.randint(n_people, size=n)
        p2 = np.random.randint(n_people, size=n)
        beta = np.ones(n)
        layer = zn.Layer(p1=p1, p2=p2, beta=beta, label='rand')
        layer = zn.Layer(dict(p1=p1, p2=p2, beta=beta), label='rand') # Alternate method

        # Convert one layer to another with extra columns
        index = np.arange(n)
        self_conn = p1 == p2
        layer2 = zn.Layer(**layer, index=index, self_conn=self_conn, label=layer.label)
    '''

    def __init__(self, *args, label=None, **kwargs):
        self.meta = {
            'p1':    znd.default_int,   # Person 1
            'p2':    znd.default_int,   # Person 2
            'beta':  znd.default_float, # Default transmissibility for this contact type
        }
        self.basekey = 'p1' # Assign a base key for calculating lengths and performing other operations
        self.label = label

        # Handle args
        kwargs = sc.mergedicts(*args, kwargs)

        # Initialize the keys of the layers
        for key,dtype in self.meta.items():
            self[key] = np.empty((0,), dtype=dtype)

        # Set data, if provided
        for key,value in kwargs.items():
            self[key] = np.array(value, dtype=self.meta.get(key))

        # Set beta if not provided
        key = 'beta'
        if key not in kwargs.keys():
            self[key] = np.ones(len(self), dtype=self.meta[key])

        return


    def __len__(self):
        try:
            return len(self[self.basekey])
        except: # pragma: no cover
            return 0


    def __repr__(self):
        ''' Convert to a dataframe for printing '''
        namestr = self.__class__.__name__
        labelstr = f'"{self.label}"' if self.label else '<no label>'
        keys_str = ', '.join(self.keys())
        output = f'{namestr}({labelstr}, {keys_str})\n' # e.g. Layer("h", p1, p2, beta)
        output += self.to_df().__repr__()
        return output


    def __contains__(self, item):
        """
        Check if a person is present in a layer

        Args:
            item: Person index

        Returns: True if person index appears in any interactions

        """
        return (item in self['p1']) or (item in self['p2'])

    @property
    def members(self):
        """
        Return sorted array of all members
        """
        return np.unique([self['p1'], self['p2']])


    def meta_keys(self):
        ''' Return the keys for the layer's meta information -- i.e., p1, p2, beta '''
        return self.meta.keys()


    def validate(self, force=True):
        '''
        Check the integrity of the layer: right types, right lengths.

        If dtype is incorrect, try to convert automatically; if length is incorrect,
        do not.
        '''
        n = len(self[self.basekey])
        for key,dtype in self.meta.items():
            if dtype:
                actual = self[key].dtype
                expected = dtype
                if actual != expected:
                    self[key] = np.array(self[key], dtype=expected) # Probably harmless, so try to convert to correct type
            actual_n = len(self[key])
            if n != actual_n:
                errormsg = f'Expecting length {n} for layer key "{key}"; got {actual_n}' # We can't fix length mismatches
                raise TypeError(errormsg)
        return


    def pop_inds(self, inds):
        '''
        "Pop" the specified indices from the edgelist and return them as a dict.
        Returns in the right format to be used with layer.append().

        Args:
            inds (int, array, slice): the indices to be removed
        '''
        output = {}
        for key in self.meta_keys():
            output[key] = self[key][inds] # Copy to the output object
            self[key] = np.delete(self[key], inds) # Remove from the original
        return output


    def append(self, contacts):
        '''
        Append contacts to the current layer.

        Args:
            contacts (dict): a dictionary of arrays with keys p1,p2,beta, as returned from layer.pop_inds()
        '''
        for key in self.keys():
            new_arr = contacts[key]
            n_curr = len(self[key]) # Current number of contacts
            n_new = len(new_arr) # New contacts to add
            n_total = n_curr + n_new # New size
            self[key] = np.resize(self[key], n_total) # Resize to make room, preserving dtype
            self[key][n_curr:] = new_arr # Copy contacts into the layer
        return


    def to_df(self):
        ''' Convert to dataframe '''
        df = pd.DataFrame.from_dict(self)
        return df


    def from_df(self, df, keys=None):
        ''' Convert from a dataframe '''
        if keys is None:
            keys = self.meta_keys()
        for key in keys:
            self[key] = df[key].to_numpy()
        return self


    def to_graph(self): # pragma: no cover
        '''
        Convert to a networkx DiGraph

        **Example**::

            import networkx as nx
            sim = zn.Sim(pop_size=20, pop_type='hybrid').run()
            G = sim.people.contacts['h'].to_graph()
            nx.draw(G)
        '''
        import networkx as nx
        data = [np.array(self[k], dtype=dtype).tolist() for k,dtype in [('p1', int), ('p2', int), ('beta', float)]]
        G = nx.DiGraph()
        G.add_weighted_edges_from(zip(*data), weight='beta')
        nx.set_edge_attributes(G, self.label, name='layer')
        return G


    def find_contacts(self, inds, as_array=True):
        """
        Find all contacts of the specified people

        For some purposes (e.g. contact tracing) it's necessary to find all of the contacts
        associated with a subset of the people in this layer. Since contacts are bidirectional
        it's necessary to check both P1 and P2 for the target indices. The return type is a Set
        so that there is no duplication of indices (otherwise if the Layer has explicit
        symmetric interactions, they could appear multiple times). This is also for performance so
        that the calling code doesn't need to perform its own unique() operation. Note that
        this cannot be used for cases where multiple connections count differently than a single
        infection, e.g. exposure risk.

        Args:
            inds (array): indices of people whose contacts to return
            as_array (bool): if true, return as sorted array (otherwise, return as unsorted set)

        Returns:
            contact_inds (array): a set of indices for pairing partners

        Example: If there were a layer with
        - P1 = [1,2,3,4]
        - P2 = [2,3,1,4]
        Then find_contacts([1,3]) would return {1,2,3}
        """

        # Check types
        if not isinstance(inds, np.ndarray):
            inds = sc.promotetoarray(inds)
        if inds.dtype != np.int64:  # pragma: no cover # This is int64 since indices often come from cv.true(), which returns int64
            inds = np.array(inds, dtype=np.int64)

        # Find the contacts
        contact_inds = znu.find_contacts(self['p1'], self['p2'], inds)
        if as_array:
            contact_inds = np.fromiter(contact_inds, dtype=znd.default_int)
            contact_inds.sort()  # Sorting ensures that the results are reproducible for a given seed as well as being identical to previous versions of Covasim

        return contact_inds


    def update(self, people, frac=1.0):
        '''
        Regenerate contacts on each timestep.

        This method gets called if the layer appears in ``sim.pars['dynam_layer']``.
        The Layer implements the update procedure so that derived classes can customize
        the update e.g. implementing over-dispersion/other distributions, random
        clusters, etc.

        Typically, this method also takes in the ``people`` object so that the
        update can depend on person attributes that may change over time (e.g.
        changing contacts for people that are severe/critical).

        Args:
            people (People): the Covasim People object, which is usually used to make new contacts
            frac (float): the fraction of contacts to update on each timestep
        '''
        # Choose how many contacts to make
        pop_size   = len(people) # Total number of people
        n_contacts = len(self) # Total number of contacts
        n_new = int(np.round(n_contacts*frac)) # Since these get looped over in both directions later
        inds = znu.choose(n_contacts, n_new)

        # Create the contacts, not skipping self-connections
        self['p1'][inds]   = np.array(znu.choose_r(max_n=pop_size, n=n_new), dtype=znd.default_int) # Choose with replacement
        self['p2'][inds]   = np.array(znu.choose_r(max_n=pop_size, n=n_new), dtype=znd.default_int)
        self['beta'][inds] = np.ones(n_new, dtype=znd.default_float)
        return