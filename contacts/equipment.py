'''
Contact via shared equipment
'''
import sciris as sc
from ..base import LayerGroup

class Equipment(LayerGroup):

    def __init__(self, eid=None, **kwargs):
        '''
        Class constructor for an empty piece of equipment

        Args:
            **eid (int): Unique equipment identifier
            **member_uids (np.array) : ids of group members
        '''
        super(Equipment, self).__init__(eid = eid, **kwargs)
        self.validate()
        return

    def validate(self):
        """
        Check that information supplied to make a piece of equipment is valid and update
        to the correct type if necessary.
        """
        super().validate(layer_str='equipment')
        return
    
def get_Equipment(pop, eid):
    """
    Return equipment with id: eid.

    Args:
        pop (contacts.Pop) : population
        eid (int)   : equipment id number

    Returns:
        contacts.Equipment: A populated farm.
    """
    if not isinstance(eid, int):
        raise TypeError(f"eid must be an int. Instead supplied eid with type: {type(eid)}.")
    if len(pop.equipment) <= eid:
        raise IndexError(f"FEquipment id (eid): {eid} out of range. There are {len(pop.equipment)} pieces of equipment stored in this object.")
    return pop.equipment[eid]
    
def add_equipment(pop, equipment):
    """
    Add a piece of equipment to the list of equipment.

    Args:
        pop (contacts.Pop)             : population
        equipment (contacts.Equipment) : piece of equipment with at minimum the eid, member_uids.
    """
    if not isinstance(equipment, Equipment):
        raise ValueError('farm is not a contacts.Farm object.')

    # ensure fid to match the index in the list
    if equipment['eid'] != len(pop.equipment):
        equipment['eid'] = len(pop.equipment)
    pop.equipment.append(equipment)
    pop.n_equipment = len(pop.equipment)
    return
    
def initialize_empty_equipment(pop, n_equipment=None):
    """
    Array of empty equipment.

    Args:
        pop (contacts.Pop)       : population
        n_equipment (int) : the number of pieces of equipment to initialize
    """
    if n_equipment is not None and isinstance(n_equipment, int):
        pop.n_equipment = n_equipment
    else:
        pop.n_equipment = 0

    pop.equipment = [Equipment() for ne in range(pop.n_equipment)]
    return

def populate_equipment(pop, equipment):
    """
    Populate all of the equipment. Store each piece at the index corresponding to it's eid.

    Args:
        pop (contacts.Pop)      : population
        equipment (list) : list of lists where each sublist represents a farm and contains the ids of the farm members
    """
    # initialize an empty set of equipment
    initialize_empty_equipment(pop, len(equipment))

    # now populate households
    for nh, hh in enumerate(equipment):
        kwargs = dict(eid=nh,
                      member_uids=hh,
                      )
        farm = Equipment()
        farm.set_layer_group(**kwargs)
        pop.equipment[equipment['eid']] = sc.dcp(equipment)

    return