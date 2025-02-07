import sciris as sc

from . import FlexDict
from . import Layer

__all__ = ['Contacts']

class Contacts(FlexDict): 
    '''
    A simple (for now) class for storing different contact layers.

    Args:
        data (dict): a dictionary that looks like a Contacts object
        layer_keys (list): if provided, create an empty Contacts object with these layers
        kwargs (dict): additional layer(s), merged with data
    '''
    def __init__(self, data=None, layer_keys=None, **kwargs):
        data = sc.mergedicts(data, kwargs)
        if layer_keys is not None:
            for lkey in layer_keys:
                self[lkey] = Layer(label=lkey)
        if data:
            for lkey,layer_data in data.items():
                self[lkey] = Layer(**layer_data)
        return

    def __repr__(self):
        ''' Use slightly customized repr'''
        keys_str = ', '.join([str(k) for k in self.keys()])
        output = f'Contacts({keys_str})\n'
        for key in self.keys():
            output += f'\n"{key}": '
            output += self[key].__repr__() + '\n'
        return output


    def __len__(self):
        ''' The length of the contacts is the length of all the layers '''
        output = 0
        for key in self.keys():
            try:
                output += len(self[key])
            except: # pragma: no cover
                pass
        return output


    def add_layer(self, **kwargs):
        '''
        Small method to add one or more layers to the contacts. Layers should
        be provided as keyword arguments.

        **Example**::

            hospitals_layer = zn.Layer(label='hosp')
            sim.people.contacts.add_layer(hospitals=hospitals_layer)
        '''
        for lkey,layer in kwargs.items():
            if not isinstance(layer, Layer):
                try:
                    layer = Layer(layer, label=lkey)
                except Exception as E:
                    exc = type(E)
                    errormsg = f'Could not parse {type(layer)} as layer: must be Layer or dict'
                    raise exc(errormsg) from E
            layer.validate()
            self[lkey] = layer
        return


    def pop_layer(self, *args):
        '''
        Remove the layer(s) from the contacts.

        **Example**::

            sim.people.contacts.pop_layer('hospitals')

        Note: while included here for convenience, this operation is equivalent
        to simply popping the key from the contacts dictionary.
        '''
        for lkey in args:
            self.pop(lkey)
        return


    def to_graph(self): # pragma: no cover
        '''
        Convert all layers to a networkx MultiDiGraph

        **Example**::

            import networkx as nx
            sim = cv.Sim(pop_size=50, pop_type='hybrid').run()
            G = sim.people.contacts.to_graph()
            nx.draw(G)
        '''
        import networkx as nx
        H = nx.MultiDiGraph()
        for lkey,layer in self.items():
            G = layer.to_graph()
            H = nx.compose(H, nx.MultiDiGraph(G))
        return H