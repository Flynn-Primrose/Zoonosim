"""
Defines default settings for generating a population.
"""
import numpy as np
import sciris as sc



default_human_pop = 100
default_poultry_pop = 10

default_layer_info = dict(
    member_uids=np.array([], dtype=int),
    )

# available globally if needed or via defaults.py --- stores information about

settings = sc.objdict()
settings.max_age = 101


def reset_settings_by_key(key, value):
    """
    Reset a key in the globally available settings dictionary with a new value.

    Returns:
        None
    """
    settings[key] = value


def reset_settings(new_config):
    """
    Reset multiple keys in the globally available settings dictionary based on a new
    dictionary of values.

    Args:
        new_config (dict) : a dictionary with new values mapped to keys

    Returns:
        None.
    """
    for key, value in new_config.items():
        reset_settings_by_key(key, value)
