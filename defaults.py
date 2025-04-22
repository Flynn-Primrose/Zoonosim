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
    c.default = '#000000' # Default colour (black)
    # Stocks
    c.susceptible_humans = '#00FF00' # Green
    c.exposed_humans = '#FFFF00' # Yellow
    c.infectious_humans = '#FF0000' # Red
    c.symptomatic_humans = '#FF8000' # Orange
    c.severe_humans = '#800000' # Dark red
    c.recovered_humans = '#00FFFF' # Cyan
    c.dead_humans = '#0000FF' # Blue
    c.susceptible_flocks = '#00FF00' # Green
    c.exposed_flocks = '#FFFF00' # Yellow
    c.infectious_flocks = '#FF0000' # Red
    c.symptomatic_flocks = '#FF8000' # Orange
    c.severe_flocks = '#800000' # Dark red
    c.contaminated_barns = '#00FFFF' # Cyan
    c.contaminated_water = '#0000FF' # Blue
    # Stocks by variant
    c.exposed_by_variant_humans = '#FFFF00' # Yellow
    c.infectious_by_variant_humans = '#FF0000' # Red
    c.symptomatic_by_variant_humans = '#FF8000' # Orange
    c.severe_by_variant_humans = '#800000' # Dark red
    c.exposed_by_variant_flocks = '#FFFF00' # Yellow
    c.infectious_by_variant_flocks = '#FF0000' # Red
    c.symptomatic_by_variant_flocks = '#FF8000' # Orange
    c.severe_by_variant_flocks = '#800000' # Dark red
    c.contaminated_by_variant_barns = '#00FFFF' # Cyan
    c.contaminated_by_variant_water = '#0000FF' # Blue
    # Flows
    c.human_reinfections = '#FF00FF' # Magenta
    c.human_infections = '#FF00FF' # Magenta
    c.human_infectious = '#FF00FF' # Magenta
    c.human_symptomatic = '#FF00FF' # Magenta
    c.human_severe = '#FF00FF' # Magenta
    c.human_recoveries = '#FF00FF' # Magenta
    c.human_deaths = '#FF00FF' # Magenta
    c.human_tests = '#FF00FF' # Magenta
    c.human_diagnoses = '#FF00FF' # Magenta
    c.human_known_deaths = '#FF00FF' # Magenta
    c.human_quarantined = '#FF00FF' # Magenta
    c.human_doses = '#FF00FF' # Magenta
    c.human_vaccinated = '#FF00FF' # Magenta
    c.flock_infections = '#FF00FF' # Magenta
    c.flock_infectious = '#FF00FF' # Magenta
    c.flock_symptomatic = '#FF00FF' # Magenta
    c.barn_contaminated = '#FF00FF' # Magenta
    c.water_contaminated = '#FF00FF' # Magenta
    # Flows by variant
    c.human_infections_by_variant = '#FF00FF' # Magenta
    c.human_symptomatic_by_variant = '#FF00FF' # Magenta
    c.human_severe_by_variant = '#FF00FF' # Magenta
    c.human_infectious_by_variant = '#FF00FF' # Magenta
    c.flock_infections_by_variant = '#FF00FF' # Magenta
    c.flock_infectious_by_variant = '#FF00FF' # Magenta
    c.flock_symptomatic_by_variant = '#FF00FF' # Magenta
    c.barn_contaminated_by_variant = '#FF00FF' # Magenta
    c.water_contaminated_by_variant = '#FF00FF' # Magenta
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
new_human_flows = [f'new_{key}' for key in human_flows.keys()]
cum_human_flows = [f'cum_{key}' for key in human_flows.keys()]

flock_flows = {
    'flock_infections':   'infections',
    'flock_infectious':   'infectious',
    'flock_symptomatic':  'symptomatic cases',
}
new_flock_flows = [f'new_{key}' for key in flock_flows.keys()]
cum_flock_flows = [f'cum_{key}' for key in flock_flows.keys()]

barn_flows = {
    'barn_contaminated': 'contaminated',
}
new_barn_flows = [f'new_{key}' for key in barn_flows.keys()]
cum_barn_flows = [f'cum_{key}' for key in barn_flows.keys()]

water_flows = {
    'water_contaminated': 'contaminated',
}
new_water_flows = [f'new_{key}' for key in water_flows.keys()]
cum_water_flows = [f'cum_{key}' for key in water_flows.keys()]

all_flows = {**human_flows, **flock_flows, **barn_flows, **water_flows}

human_flows_by_variant = {
    'human_infections_by_variant':  'infections by variant',
    'human_symptomatic_by_variant': 'symptomatic by variant',
    'human_severe_by_variant':      'severe by variant',
    'human_infectious_by_variant':  'infectious by variant',
}
new_human_flows_by_variant = [f'new_{key}' for key in human_flows_by_variant.keys()]
cum_human_flows_by_variant = [f'cum_{key}' for key in human_flows_by_variant.keys()]

flock_flows_by_variant = {
    'flock_infections_by_variant':  'infections by variant',
    'flock_infectious_by_variant':  'infectious by variant',
    'flock_symptomatic_by_variant': 'symptomatic by variant',
}
new_flock_flows_by_variant = [f'new_{key}' for key in flock_flows_by_variant.keys()]
cum_flock_flows_by_variant = [f'cum_{key}' for key in flock_flows_by_variant.keys()]

barn_flows_by_variant = {
    'barn_contaminated_by_variant': 'contaminated by variant',
}
new_barn_flows_by_variant = [f'new_{key}' for key in barn_flows_by_variant.keys()]
cum_barn_flows_by_variant = [f'cum_{key}' for key in barn_flows_by_variant.keys()]

water_flows_by_variant = {
    'water_contaminated_by_variant': 'contaminated by variant',
}
new_water_flows_by_variant = [f'new_{key}' for key in water_flows_by_variant.keys()]
cum_water_flows_by_variant = [f'cum_{key}' for key in water_flows_by_variant.keys()]

all_flows_by_variant = {**human_flows_by_variant, **flock_flows_by_variant, **barn_flows_by_variant, **water_flows_by_variant}

# Define new and cumulative flows
new_result_flows = [f'new_{key}' for key in all_flows.keys()]
cum_result_flows = [f'cum_{key}' for key in all_flows.keys()]
new_result_flows_by_variant = [f'new_{key}' for key in all_flows_by_variant.keys()]
cum_result_flows_by_variant = [f'cum_{key}' for key in all_flows_by_variant.keys()]