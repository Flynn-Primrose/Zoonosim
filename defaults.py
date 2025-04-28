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
nbbool        = nb.bool_

verbose = 0.1
warnings = 'warn' # must be one of 'warn', 'print', or 'error'

safe_opts = [1, '1', 'safe'] 
full_opts = [2, '2', 'full'] 

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
    c.susceptible = '#00FF00' # Green
    c.exposed = '#FFFF00' # Yellow
    c.infectious = '#FF0000' # Red
    c.symptomatic = '#FF8000' # Orange
    c.severe = '#800000' # Dark red
    c.recovered = '#00FFFF' # Cyan
    c.dead = '#0000FF' # Blue
    c.contaminated = '#00FFFF' # Cyan

    # Stocks by variant
    c.exposed_by_variant = '#FFFF00' # Yellow
    c.infectious_by_variant = '#FF0000' # Red
    c.symptomatic_by_variant = '#FF8000' # Orange
    c.severe_by_variant = '#800000' # Dark red
    c.contaminated_by_variant_barns = '#00FFFF' # Cyan
    # Flows
    c.reinfections = '#FF00FF' # Magenta
    c.infections = '#FF00FF' # Magenta
    c.infectious = '#FF00FF' # Magenta
    c.symptomatic = '#FF00FF' # Magenta
    c.severe = '#FF00FF' # Magenta
    c.recoveries = '#FF00FF' # Magenta
    c.deaths = '#FF00FF' # Magenta
    c.tests = '#FF00FF' # Magenta
    c.diagnoses = '#FF00FF' # Magenta
    c.known_deaths = '#FF00FF' # Magenta
    c.quarantined = '#FF00FF' # Magenta
    c.doses = '#FF00FF' # Magenta
    c.vaccinated = '#FF00FF' # Magenta
    # Flows by variant
    c.infections_by_variant = '#FF00FF' # Magenta
    c.symptomatic_by_variant = '#FF00FF' # Magenta
    c.severe_by_variant = '#FF00FF' # Magenta
    c.infectious_by_variant = '#FF00FF' # Magenta
    c.infections_by_variant = '#FF00FF' # Magenta
    c.contaminated_by_variant = '#FF00FF' # Magenta
    return c

# Parameters that can vary by variant
variant_pars = [
    'rel_beta',
    'rel_symp_prob',
    'rel_severe_prob',
    'rel_death_prob',
]

default_human_prognoses = dict(
    age_cutoffs   = np.array([0,       10,      20,      30,      40,      50,      60,      70,      80,      90,]),     # Age cutoffs (lower limits)
    sus_ORs       = np.array([0.34,    0.67,    1.00,    1.00,    1.00,    1.00,    1.24,    1.47,    1.47,    1.47]),    # Odds ratios for relative susceptibility -- from Zhang et al., https://science.sciencemag.org/content/early/2020/05/04/science.abb8001; 10-20 and 60-70 bins are the average across the ORs
    trans_ORs     = np.array([1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00]),    # Odds ratios for relative transmissibility -- no evidence of differences
    comorbidities = np.array([1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00]),    # Comorbidities by age -- set to 1 by default since already included in disease progression rates
    symp_probs    = np.array([0.50,    0.55,    0.60,    0.65,    0.70,    0.75,    0.80,    0.85,    0.90,    0.90]),    # Overall probability of developing symptoms (based on https://www.medrxiv.org/content/10.1101/2020.03.24.20043018v1.full.pdf, scaled for overall symptomaticity)
    severe_probs  = np.array([0.00050, 0.00165, 0.00720, 0.02080, 0.03430, 0.07650, 0.13280, 0.20655, 0.24570, 0.24570]), # Overall probability of developing severe symptoms (derived from Table 1 of https://www.imperial.ac.uk/media/imperial-college/medicine/mrc-gida/2020-03-16-COVID19-Report-9.pdf)
    death_probs   = np.array([0.00002, 0.00002, 0.00010, 0.00032, 0.00098, 0.00265, 0.00766, 0.02439, 0.08292, 0.16190]), # Overall probability of dying -- from O'Driscoll et al., https://www.nature.com/articles/s41586-020-2918-0; last data point from Brazeau et al., https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/covid-19/report-34-ifr/
        )

default_flock_prognoses = dict(
    breeds = np.array(['breeder', 'broiler', 'layer']),
    sus_ORs = np.array([1.00, 1.00, 1.00]),
    trans_ORs = np.array([1.00, 1.00, 1.00]),
    symp_probs = np.array([1.00, 1.00, 1.00])
)

default_barn_prognoses = dict(
    sus_ORs = 1.00,
    trans_ORs = 1.00
)

default_water_prognoses = dict(
    sus_ORs = 1.00,
    trans_ORs = 1.00
)


# tracked stocks and flows -- used in sim.py; value is the label suffix

human_stocks = {
    'susceptible': 'Number of susceptible humans',
    'exposed':    'Number of exposed humans',
    'infectious':   'Number of infectious humans',
    'symptomatic':  'Number of symptomatic humans',
    'severe':       'Number of severe humans',
    'recovered':    'Number of recovered humans',
    'dead':         'Number of dead humans',
}

flock_stocks = {
    'susceptible': 'Number of susceptible flocks',
    'exposed':     'Number of exposed flocks',
    'infectious':  'Number of infectious flocks',
    'symptomatic': 'Number of symptomatic flocks',
}

barn_stocks = {
    'contaminated': 'Number of contaminated barns',
}

water_stocks = {
    'contaminated': 'Number of contaminated waters',
}

#all_stocks = {**human_stocks, **flock_stocks, **barn_stocks, **water_stocks}

human_stocks_by_variant = {
    'exposed_by_variant':    'Number exposed humans by variant',
    'infectious_by_variant': 'Number infectious humans by variant',
    'symptomatic_by_variant': 'Number symptomatic humans by variant',
    'severe_by_variant':      'Number severe humans by variant',
}

flock_stocks_by_variant = {
    'exposed_by_variant':     'Number exposed flocks by variant',
    'infectious_by_variant':  'Number infectious flocks by variant',
    'symptomatic_by_variant': 'Number symptomatic flocks by variant',
}


barn_stocks_by_variant = {
    'contaminated_by_variant': 'Number contaminated barns by variant',
}

water_stocks_by_variant = {
    'contaminated_by_variant': 'Number contaminated water by variant',
}

#all_stocks_by_variant = {**human_stocks_by_variant, **flock_stocks_by_variant, **barn_stocks_by_variant, **water_stocks_by_variant}


# The types of result that are counted as flows -- used in sim.py; value is the label suffix
human_flows = {
    'reinfections': 'human reinfections',
    'infections':   'human infections',
    'infectious':   'infectious humans',
    'symptomatic':  'symptomatic human cases',
    'severe':       'severe human cases',
    'recoveries':   'human recoveries',
    'deaths':       'human deaths',
    'tests':        'human tests',
    'diagnoses':    'human diagnoses',
    'known_deaths': 'known human deaths',
    'quarantined':  'human quarantines started',
    'doses':        'human vaccine doses',
    'vaccinated':   'vaccinated humans'
}
new_human_flows = [f'new_{key}' for key in human_flows.keys()]
cum_human_flows = [f'cum_{key}' for key in human_flows.keys()]

flock_flows = {
    'infections':   'flock infections',
    'infectious':   'infectious flocks',
    'symptomatic':  'symptomatic flock cases',
}
new_flock_flows = [f'new_{key}' for key in flock_flows.keys()]
cum_flock_flows = [f'cum_{key}' for key in flock_flows.keys()]

barn_flows = {
    'contaminated': 'contaminated barns',
}
new_barn_flows = [f'new_{key}' for key in barn_flows.keys()]
cum_barn_flows = [f'cum_{key}' for key in barn_flows.keys()]

water_flows = {
    'contaminated': 'contaminated waterbodies',
}
new_water_flows = [f'new_{key}' for key in water_flows.keys()]
cum_water_flows = [f'cum_{key}' for key in water_flows.keys()]

#all_flows = {**human_flows, **flock_flows, **barn_flows, **water_flows}

human_flows_by_variant = {
    'infections_by_variant':  'human infections by variant',
    'symptomatic_by_variant': 'symptomatic humans by variant',
    'severe_by_variant':      'severe humans by variant',
    'infectious_by_variant':  'infectious humans by variant',
}
new_human_flows_by_variant = [f'new_{key}' for key in human_flows_by_variant.keys()]
cum_human_flows_by_variant = [f'cum_{key}' for key in human_flows_by_variant.keys()]

flock_flows_by_variant = {
    'infections_by_variant':  'flock infections by variant',
    'infectious_by_variant':  'infectious flocks by variant',
    'symptomatic_by_variant': 'symptomatic flocks by variant',
}
new_flock_flows_by_variant = [f'new_{key}' for key in flock_flows_by_variant.keys()]
cum_flock_flows_by_variant = [f'cum_{key}' for key in flock_flows_by_variant.keys()]

barn_flows_by_variant = {
    'contaminated_by_variant': 'contaminated barns by variant',
}
new_barn_flows_by_variant = [f'new_{key}' for key in barn_flows_by_variant.keys()]
cum_barn_flows_by_variant = [f'cum_{key}' for key in barn_flows_by_variant.keys()]

water_flows_by_variant = {
    'contaminated_by_variant': 'contaminated waterbodies by variant',
}
new_water_flows_by_variant = [f'new_{key}' for key in water_flows_by_variant.keys()]
cum_water_flows_by_variant = [f'cum_{key}' for key in water_flows_by_variant.keys()]

#all_flows_by_variant = {**human_flows_by_variant, **flock_flows_by_variant, **barn_flows_by_variant, **water_flows_by_variant}

# Define new and cumulative flows
#new_result_flows = [f'new_{key}' for key in all_flows.keys()]
#cum_result_flows = [f'cum_{key}' for key in all_flows.keys()]
#new_result_flows_by_variant = [f'new_{key}' for key in all_flows_by_variant.keys()]
#cum_result_flows_by_variant = [f'cum_{key}' for key in all_flows_by_variant.keys()]