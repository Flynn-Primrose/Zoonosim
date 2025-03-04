'''
Contacts via shared farm
'''
import sciris as sc
from ..base import LayerGroup

class Farm(LayerGroup):

    def __init__(self, fid=None, **kwargs):
        '''
        Class constructor for an empty farm

        Args:
            **fid (int): Unique farm identifier
            **member_uids (np.array) : ids of group members
        '''
        super(Farm, self).__init__(fid = fid, **kwargs)
        self.validate()
        return

    def validate(self):
        """
        Check that information supplied to make a household is valid and update
        to the correct type if necessary.
        """
        super().validate(layer_str='farm')
        return
    
def get_farm(pop, fid):
    """
    Return equipment with id: fid.

    Args:
        pop (contacts.Pop) : population
        fid (int)   : farm id number

    Returns:
        contacts.Farm: A populated farm.
    """
    if not isinstance(fid, int):
        raise TypeError(f"fid must be an int. Instead supplied fid with type: {type(fid)}.")
    if len(pop.farms) <= fid:
        raise IndexError(f"Farm id (fid): {fid} out of range. There are {len(pop.farms)} farms stored in this object.")
    return pop.farms[fid]
    
def add_farm(pop, farm):
    """
    Add a farm to the list of farms.

    Args:
        pop (contacts.Pop)             : population
        farm (contacts.Farm) : Farm with at minimum the fid, member_uids.
    """
    if not isinstance(farm, Farm):
        raise ValueError('farm is not a contacts.Farm object.')

    # ensure fid to match the index in the list
    if farm['fid'] != len(pop.farms):
        farm['fid'] = len(pop.farms)
    pop.farms.append(farm)
    pop.n_farms = len(pop.farms)
    return
    
def initialize_empty_farms(pop, n_farms=None):
    """
    Array of empty farms.

    Args:
        pop (contacts.Pop)       : population
        n_farms (int) : the number of farms to initialize
    """
    if n_farms is not None and isinstance(n_farms, int):
        pop.n_farms = n_farms
    else:
        pop.n_farms = 0

    pop.farms = [Farm() for nf in range(pop.n_farms)]
    return

def populate_farms(pop, farms):
    """
    Populate all of the farms. Store each farm at the index corresponding to it's fid.

    Args:
        pop (contacts.Pop)      : population
        farms (list) : list of lists where each sublist represents a farm and contains the ids of the farm members
        age_by_uid (dict) : dictionary mapping each person's id to their age
    """
    # initialize an empty set of farms
    # if previously you had 10 farms and now you want to repopulate with
    # this method and only supply 5 farms, this method will overwrite the list to produce only 5 farms
    initialize_empty_farms(pop, len(farms))

    # now populate households
    for nh, hh in enumerate(farms):
        kwargs = dict(fid=nh,
                      member_uids=hh,
                      )
        farm = Farm()
        farm.set_layer_group(**kwargs)
        pop.farms[farm['fid']] = sc.dcp(farm)

    return