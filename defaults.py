'''
Defines default values for Zoonosim.
'''
import sciris as sc
import numpy as np
import numba as nb
import os

# Specify all externally visible functions this file defines -- other things are available as e.g. zn.defaults.default_int
__all__ = ['default_float', 'default_int', 'nbfloat', 'nbint', 'verbose', 'get_default_colors']

result_float = np.float64 # Always use float64 for results, for simplicity

# Uncomment the following lines to use float32 and int32 instead of float64 and int64, I'm always using float64 and int64 for simplicity.
#default_float = np.float32
#default_int   = np.int32
#nbfloat       = nb.float32
#nbint         = nb.int32

default_float = np.float64
default_int   = np.int64
nbfloat       = nb.float64
nbint         = nb.int64

verbose = 0.1
warnings = 'warn' # must be one of 'warn', 'print', or 'error'

numba_parallel = str(os.getenv('ZOONOSIM_NUMBA_PARALLEL', 'none'))
numba_cache = bool(int(os.getenv('ZOONOSIM_NUMBA_CACHE', 1)))

default_pop_pars = {
    'avg_barns_per_farm': 1.0,
    'avg_humans_per_barn': 1.0,
    'avg_barn_occupancy': 1.0, # Probability that a barn is occupied by a flock at any given time
    'avg_flock_size': 1.0,
    'avg_water_per_farm': 1.0,
}

def get_default_colors():
    '''
    Specify plot colors -- used in sim.py.

    NB, includes duplicates since stocks and flows are named differently.

    Colour palette: 
    '''
    c = sc.objdict()
    # TODO: add colours as needed
    c.default = '#1f77b4' # Default colour
    return c

# Parameters that can vary by variant
variant_pars = [
    'rel_beta',
    'rel_symp_prob',
    'rel_severe_prob',
    'rel_death_prob',
]



# States to be used for results.
# TODO: modify to accommodate multiple agent types
result_stocks = {
    'susceptible': 'Number susceptible',
    'exposed':     'Number exposed',
    'infectious':  'Number infectious',
    'symptomatic': 'Number symptomatic',
    'severe':      'Number of severe cases',
    'recovered':   'Number recovered',
    'dead':        'Number dead',
    'diagnosed':   'Number of confirmed cases',
    'known_dead':  'Number of confirmed deaths',
    'quarantined': 'Number in quarantine',
    'vaccinated':  'Number of people vaccinated',
}

result_stocks_by_variant = {
    'exposed_by_variant':    'Number exposed by variant',
    'infectious_by_variant': 'Number infectious by variant',
}

# The types of result that are counted as flows -- used in sim.py; value is the label suffix
result_flows = {
    'infections':   'infections',
    'reinfections': 'reinfections',
    'infectious':   'infectious',
    'symptomatic':  'symptomatic cases',
    'severe':       'severe cases',
    'recoveries':   'recoveries',
    'deaths':       'deaths',
    'tests':        'tests',
    'diagnoses':    'diagnoses',
    'known_deaths': 'known deaths',
    'quarantined':  'quarantines started',
    'doses':        'vaccine doses',
    'vaccinated':   'vaccinated people'
}

result_flows_by_variant = {
    'infections_by_variant':  'infections by variant',
    'symptomatic_by_variant': 'symptomatic by variant',
    'severe_by_variant':      'severe by variant',
    'infectious_by_variant':  'infectious by variant',
}

# Define new and cumulative flows
new_result_flows = [f'new_{key}' for key in result_flows.keys()]
cum_result_flows = [f'cum_{key}' for key in result_flows.keys()]
new_result_flows_by_variant = [f'new_{key}' for key in result_flows_by_variant.keys()]
cum_result_flows_by_variant = [f'cum_{key}' for key in result_flows_by_variant.keys()]