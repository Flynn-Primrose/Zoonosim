'''
Contacts via shared water sources
'''

import sciris as sc
from ..base import LayerGroup

class WaterSource(LayerGroup):

    def __init__(self, wsid=None, **kwargs):
        '''
        Class constructor for an empty water source

        Args:
            **wsid (int): Unique water source identifier
            **member_uids (np.array) : ids of group members
        '''
        super(WaterSource, self).__init__(wsid = wsid, **kwargs)
        self.validate()
        return

    def validate(self):
        """
        Check that information supplied to make a water source is valid and update
        to the correct type if necessary.
        """
        super().validate(layer_str='water source')
        return
    
def get_water_source(pop, wsid):
    """
    Return household with id: wsid.

    Args:
        pop (contacts.Pop) : population
        wsid (int)   : water source id number

    Returns:
        contacts.FarmWaterSource: A populated water source.
    """
    if not isinstance(wsid, int):
        raise TypeError(f"wsid must be an int. Instead supplied fid with type: {type(wsid)}.")
    if len(pop.water_sources) <= wsid:
        raise IndexError(f"FWater source id (wsid): {wsid} out of range. There are {len(pop.water_sources)} water sources stored in this object.")
    return pop.water_sources[wsid]
    
def add_water_source(pop, water_source):
    """
    Add a water source to the list of water sources.

    Args:
        pop (contacts.Pop)             : population
        water_source (contacts.WaterSource) : Water source with at minimum the wsid, member_uids.
    """
    if not isinstance(water_source, WaterSource):
        raise ValueError('water_source is not a contacts.WaterSource object.')

    # ensure wsid to match the index in the list
    if water_source['wsid'] != len(pop.water_sources):
        water_source['wsid'] = len(pop.water_sources)
    pop.water_sources.append(water_source)
    pop.n_water_sources = len(pop.water_sources)
    return
    
def initialize_empty_water_sources(pop, n_water_sources=None):
    """
    Array of empty water sources.

    Args:
        pop (contacts.Pop)       : population
        n_water_sources (int) : the number of water sources to initialize
    """
    if n_water_sources is not None and isinstance(n_water_sources, int):
        pop.n_water_sources = n_water_sources
    else:
        pop.n_water_sources = 0

    pop.water_sources = [WaterSource() for nws in range(pop.n_water_sources)]
    return

def populate_water_sources(pop, water_sources):
    """
    Populate all of the water sources. Store each water source at the index corresponding to it's wsid.

    Args:
        pop (contacts.Pop)      : population
        water_sources (list) : list of lists where each sublist represents a farm and contains the ids of the farm members
    """
    # initialize an empty set of water sources
    # if previously you had 10 water sources and now you want to repopulate with
    # this method and only supply 5 water sources, this method will overwrite the list to produce only 5 water sources
    initialize_empty_water_sources(pop, len(water_sources))

    # now populate water sources
    for nh, hh in enumerate(water_sources):
        kwargs = dict(wsid=nh,
                      member_uids=hh,
                      )
        water_source = WaterSource()
        water_source.set_layer_group(**kwargs)
        pop.water_sources[water_source['wsid']] = sc.dcp(water_source)

    return