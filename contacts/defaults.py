"""
Defaults for synthpops files and data types.

NOTE: The contents of this file should eventually be migrated to the config submodule for consistency.
"""
import numpy as np
import sciris as sc
import os


default_pop_size = 200

# specify default valid probability distributions - users can easily supply
# their own list if interested in other properties
valid_probability_distributions = [
    'population_age_distributions',
    'household_size_distribution',
    # 'school_size_distribution',
]

default_data = {
    'defaults': {
        'country_location': 'Canada',
        'state_location': 'Quebec',
        'location': 'Montreal',
        'sheet_name': 'Canada',
        'nbrackets': 20,
    }
}

default_layer_info = dict(
    member_uids=np.array([], dtype=int),
    )

def default_datadir_path(): # TODO: Direct this to zoonosim/data/contact_data
    """Return the path to the internal data folder."""
    thisdir = os.path.dirname(os.path.abspath(__file__)) 
    return os.path.join(thisdir, 'data')

# available globally if needed or via defaults.py --- stores information about
# location information to search for data when unavailable for the location
# specified or the parent locations specified
settings = sc.objdict()

settings.thisdir = os.path.dirname(os.path.abspath(__file__))
settings.localdatadir = default_datadir_path()
settings.datadir = settings.localdatadir

settings.relative_path = []


settings.max_age = 101
settings.nbrackets = 20
settings.valid_nbracket_ranges = [16, 18, 20]

settings.country_location = None
settings.state_location = None
settings.location = None
settings.sheet_name = None

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


def reset_default_settings():
    """
    By default, synthpops uses default settings for Seattle, Washington, USA.
    After having changed the values in the settings dictionary, calling this
    method can easily reset the settings dictionary to the values for Montreal, 
    Quebec, Canada.

    Returns:
        None
    """
    reset_settings(default_data['defaults']) # pragma: no cover