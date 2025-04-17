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

default_str = str # Use str for strings, for simplicity
default_float = np.float64
default_int   = np.int64
nbfloat       = nb.float64
nbint         = nb.int64

verbose = 0.1
warnings = 'warn' # must be one of 'warn', 'print', or 'error'

numba_parallel = str(os.getenv('ZOONOSIM_NUMBA_PARALLEL', 'none'))
numba_cache = bool(int(os.getenv('ZOONOSIM_NUMBA_CACHE', 1)))

default_pop_pars = {
    'avg_barns_per_farm': 3.0,
    'avg_humans_per_barn': 2.0,
    'avg_flock_size': 1000.0,
    'avg_water_per_farm': 1.0,
}

default_flock_breeds = ['breeder', 'layer', 'broiler'] # Breeds of flocks present
default_flock_breed_freqs = [0.1, 0.2, 0.7] # frequency of the different breed types



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


# tracked stocks and flows -- used in sim.py; value is the label suffix

human_stocks = {
    'susceptible_humans': 'Number of susceptible humans',
    'exposed_humans':    'Number of exposed humans',
    'infectious_humans':   'Number of infectious humans',
    'symptomatic_humans':  'Number of symptomatic humans',
    'severe_humans':       'Number of severe humans',
    'recovered_humans':    'Number of recovered humans',
    'dead_humans':         'Number of dead humans',
}

flock_stocks = {
    'susceptible_flocks': 'Number of susceptible flocks',
    'exposed_flocks':     'Number of exposed flocks',
    'infectious_flocks':  'Number of infectious flocks',
    'symptomatic_flocks': 'Number of symptomatic flocks',
    'severe_flocks':      'Number of severe flocks',
}

barn_stocks = {
    'contaminated_barns': 'Number of contaminated barns',
}

water_stocks = {
    'contaminated_water': 'Number of contaminated waters',
}

all_stocks = {**human_stocks, **flock_stocks, **barn_stocks, **water_stocks}

human_stocks_by_variant = {
    'exposed_by_variant_humans':    'Number exposed humans by variant',
    'infectious_by_variant_humans': 'Number infectious humans by variant',
    'symptomatic_by_variant_humans': 'Number symptomatic humans by variant',
    'severe_by_variant_humans':      'Number severe humans by variant',
}

flock_stocks_by_variant = {
    'exposed_by_variant_flocks':     'Number exposed flocks by variant',
    'infectious_by_variant_flocks':  'Number infectious flocks by variant',
    'symptomatic_by_variant_flocks': 'Number symptomatic flocks by variant',
}


barn_stocks_by_variant = {
    'contaminated_by_variant_barns': 'Number contaminated barns by variant',
}

water_stocks_by_variant = {
    'contaminated_by_variant_water': 'Number contaminated water by variant',
}

all_stocks_by_variant = {**human_stocks_by_variant, **flock_stocks_by_variant, **barn_stocks_by_variant, **water_stocks_by_variant}


# The types of result that are counted as flows -- used in sim.py; value is the label suffix
human_flows = {
    'human_reinfections': 'reinfections',
    'human_infections':   'infections',
    'human_infectious':   'infectious',
    'human_symptomatic':  'symptomatic cases',
    'human_severe':       'severe cases',
    'human_recoveries':   'recoveries',
    'human_deaths':       'deaths',
    'human_tests':        'tests',
    'human_diagnoses':    'diagnoses',
    'human_known_deaths': 'known deaths',
    'human_quarantined':  'quarantines started',
    'human_doses':        'vaccine doses',
    'human_vaccinated':   'vaccinated people'
}

flock_flows = {
    'flock_infections':   'infections',
    'flock_infectious':   'infectious',
    'flock_symptomatic':  'symptomatic cases',
}

barn_flows = {
    'barn_contaminated': 'contaminated',
}

water_flows = {
    'water_contaminated': 'contaminated',
}

all_flows = {**human_flows, **flock_flows, **barn_flows, **water_flows}

human_flows_by_variant = {
    'human_infections_by_variant':  'infections by variant',
    'human_symptomatic_by_variant': 'symptomatic by variant',
    'human_severe_by_variant':      'severe by variant',
    'human_infectious_by_variant':  'infectious by variant',
}

flock_flows_by_variant = {
    'flock_infections_by_variant':  'infections by variant',
    'flock_infectious_by_variant':  'infectious by variant',
    'flock_symptomatic_by_variant': 'symptomatic by variant',
}

barn_flows_by_variant = {
    'barn_contaminated_by_variant': 'contaminated by variant',
}

water_flows_by_variant = {
    'water_contaminated_by_variant': 'contaminated by variant',
}

all_flows_by_variant = {**human_flows_by_variant, **flock_flows_by_variant, **barn_flows_by_variant, **water_flows_by_variant}

# Define new and cumulative flows
new_result_flows = [f'new_{key}' for key in all_flows.keys()]
cum_result_flows = [f'cum_{key}' for key in all_flows.keys()]
new_result_flows_by_variant = [f'new_{key}' for key in all_flows_by_variant.keys()]
cum_result_flows_by_variant = [f'cum_{key}' for key in all_flows_by_variant.keys()]