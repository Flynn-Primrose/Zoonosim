

import numpy as np
import sciris as sc
from ..config import Defaults as znd

__all__ = ['LayerGroup']

class LayerGroup(dict):
    """
    A generic class for individual setting group and some methods to operate on each.

    Args:
        kwargs (dict) : data dictionary for the setting group

    Notes:
        Settings currently supported include : 
    """

    def __init__(self, **kwargs):
        """
        Class constructor for an base empty setting group.

        Args:
            **member_uids (np.array) : ids of group members
        """
        # set up default values
        default_kwargs = znd.default_layer_info
        kwargs = sc.mergedicts(default_kwargs, kwargs)
        self.update(kwargs)
        self.validate()

        return

    def set_layer_group(self, **kwargs):
        """Set layer group values."""
        for key, value in kwargs.items():
            self[key] = value
        self.validate()

        return

    def __len__(self):
        """Return the length as the number of members in the layer group"""
        return len(self['member_uids'])

    def validate(self, layer_str=''):
        """
        Check that information supplied to make a household is valid and update
        to the correct type if necessary.
        """
        for key in self.keys():
            if key in ['member_uids']:
                try:
                    self[key] = sc.promotetoarray(self[key], dtype=int)
                except:
                    errmsg = f"Could not convert key {key} to an np.array() with type int. This key only takes arrays with int values."
                    raise TypeError(errmsg)
            else:
                if not isinstance(self[key], (int, np.int32, np.int64)):
                    if self[key] is not None:
                        errmsg = f"error: Expected type int or None for {layer_str} key {key}. Instead the type of this value is {type(self[key])}."
                        raise TypeError(errmsg)

        return

    def member_ages(self, age_by_uid, subgroup_member_uids=None):
        """
        Return the ages of members in the layer group given the pop object.

        Args:
            age_by_uid (np.ndarray) : mapping of age to uid
            subgroup_member_uids (np.ndarray, list) : subgroup of uids to return ages for

        Returns:
            nd.ndarray : ages of members in group or subgroup
        """
        if len(age_by_uid) == 0:
            print("age_by_uid is empty. Returning an empty array for member_ages.")
            return np.array([])

        if subgroup_member_uids is None:
            return np.array(age_by_uid[self['member_uids']])
        else:
            subgroup_member_uids = sc.tolist(subgroup_member_uids)
            return np.array(age_by_uid[subgroup_member_uids])