'''
Defines default values for Zoonosim.
'''
import sciris as sc
import numpy as np
import numba as nb
import pylab as pl
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


show = True # Show figures by default
close = False # Close figures by default
returnfig = True # Return figures by default
dpi = pl.rcParams['figure.dpi'] # Default DPI for figures
font = pl.rcParams['font.family'] # Default font family for figures
fontsize = pl.rcParams['font.size'] # Default font size for figures
backend = pl.get_backend() # Default backend for figures
sep = ',' # Default thousands separator for text output
precision = 32 # Default arithmetic precision for Numba -- 32-bit by default for efficiency

# Define the 'overview plots', i.e. the most useful set of plots to explore different aspects of a simulation
overview_plots = [
    'new_human_infections',
    'new_flock_infections',
    'new_barn_contaminated',
    'new_water_contaminated',
]

overview_variant_plots = [
    'new_human_infections_by_variant',
    'new_flock_infections_by_variant',
    'new_barn_contaminated_by_variant',
    'new_water_contaminated_by_variant',
]

def get_default_plots(which='default', kind='sim', sim=None):
    '''
    Specify which quantities to plot; used in sim.py.

    Args:
        which (str): either 'default' or 'overview'
    '''
    which = str(which).lower() # To make comparisons easier

    # Check that kind makes sense
    sim_kind   = 'sim'
    scens_kind = 'scens'
    kindmap = {
        None:      sim_kind,
        'sim':     sim_kind,
        'default': sim_kind,
        'msim':    scens_kind,
        'scen':    scens_kind,
        'scens':   scens_kind,
    }
    if kind not in kindmap.keys():
        errormsg = f'Expecting "sim" or "scens", not "{kind}"'
        raise ValueError(errormsg)
    else:
        is_sim = kindmap[kind] == sim_kind

    # Default plots -- different for sims and scenarios
    if which in ['none', 'default']:

        if is_sim:
            plots = sc.odict({
                'Total counts': [
                    'cum_human_infections',
                    'n_human_infectious',
                    'cum_flock_infections',
                    'n_flock_infectious',
                    'cum_barn_contaminated',
                    'n_barn_contaminated',
                    'cum_water_contaminated',
                    'n_water_contaminated',
                ],
                'Daily counts': [
                    'new_human_infections',
                    'new_flock_infections',
                    'new_barn_contaminated',
                    'new_water_contaminated',
                ],
                'Health outcomes': [
                    'cum_human_infectious',
                    'cum_human_symptomatic',
                    'cum_human_severe',
                    'cum_human_deaths',
                ],
            })

        else: # pragma: no cover
            plots = sc.odict({
                'Cumulative infections': [
                    'cum_human_infections',
                    'cum_flock_infections',
                    'cum_barn_contaminated',
                    'cum_water_contaminated',
                ],
                'New infections per day': [
                    'new_human_infections',
                    'new_flock_infections',
                    'new_barn_contaminated',
                    'new_water_contaminated',
                ],
                'Cumulative deaths': [
                    'cum_human_deaths',
                ],
            })

    # Show an overview
    elif which == 'overview': # pragma: no cover
        plots = sc.dcp(overview_plots)

    # Plot absolutely everything
    elif which == 'all': # pragma: no cover
        plots = sim.result_keys('all')

    # Show an overview plus variants
    elif 'overview' in which and 'variant' in which: # pragma: no cover
        plots = sc.dcp(overview_plots) + sc.dcp(overview_variant_plots)

    # Show default but with variants
    elif which.startswith('variant'): # pragma: no cover
        if is_sim:
            plots = sc.odict({
                'Cumulative infections by variant': [
                    'cum_human_infections_by_variant',
                    'cum_flock_infections_by_variant',
                    'cum_barn_contaminated_by_variant',
                    'cum_water_contaminated_by_variant',
                ],
                'New infections by variant': [
                    'new_human_infections_by_variant',
                    'new_flock_infections_by_variant',
                    'new_barn_contaminated_by_variant',
                    'new_water_contaminated_by_variant',
                ],
                'Health outcomes': [
                    'cum_human_infections',
                    'cum_human_severe',
                    'cum_human_deaths',
                ],
            })

        else: # pragma: no cover
            plots = sc.odict({
                    'Cumulative infections by variant': [
                        'cum_human_infections_by_variant',
                        'cum_flock_infections_by_variant',
                        'cum_barn_contaminated_by_variant',
                        'cum_water_contaminated_by_variant',
                    ],
                    'New infections by variant': [
                        'new_human_infections_by_variant',
                        'new_flock_infectious_by_variant',
                        'new_barn_contaminated_by_variant',
                        'new_water_contaminated_by_variant',
                    ],
                    'New diagnoses': [
                        'new_human_diagnoses',
                    ],
                    'Cumulative deaths': [
                        'cum_human_deaths',
                    ],
            })

    # Plot agent type specific plots
    elif which == 'human': # pragma: no cover
        plots = [
            'n_human_susceptible',
            'n_human_exposed',
            'n_human_infectious',
            'n_human_symptomatic',
            'n_human_severe',
            'n_human_recovered',
            'n_human_dead',
        ],

    else: # pragma: no cover
        errormsg = f'The choice which="{which}" is not supported: choices are "default", "overview", "all", "variant", "overview-variant", or "seir", along with any result key (see sim.results_keys(\'all\') for options)'
        raise ValueError(errormsg)

    return plots



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


default_pop_pars = {
    'avg_barns_per_farm': 3.0,
    'avg_humans_per_barn': 2.0,
    'avg_flock_size': 1000.0,
    'avg_water_per_farm': 1.0,
}

default_flock_breeds = ['duck', 'layer', 'broiler'] # Breeds of flocks present
default_flock_breed_freqs = [0.1, 0.2, 0.7] # frequency of the different breed types
default_suspicious_mortality_rate = 0.002 # mortality rate that triggers a flock to be considered suspicious
default_suspicious_symptomatic_rate = 0.005 # symptomatic rate that triggers a flock to be considered suspicious
default_suspicious_consumption_rate = 1.5 # water rate that triggers a flock to be considered suspicious (L/bird/day)

# Parameters that can vary by variant
variant_pars = [
    'rel_beta',
    'rel_symp_prob',
    'rel_severe_prob',
    'rel_death_prob',
]

default_human_prognoses = dict(
    age_cutoffs   = np.array([0,       10,      20,      30,      40,      50,      60,      70,      80,      90,]),     # Age cutoffs (lower limits)
    sus_ORs       = np.array([0.34,    0.67,    1.00,    1.00,    1.00,    1.00,    1.24,    1.47,    1.47,    1.47]),    # Odds ratios for relative susceptibility 
    trans_ORs     = np.array([1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00]),    # Odds ratios for relative transmissibility
    comorbidities = np.array([1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00]),    # Comorbidities by age -- set to 1 by default since already included in disease progression rates
    symp_probs    = np.array([0.50,    0.55,    0.60,    0.65,    0.70,    0.75,    0.80,    0.85,    0.90,    0.90]),    # Overall probability of developing symptoms 
    severe_probs  = np.array([0.01,    0.01,    0.05,    0.05,    0.05,    0.07,    0.09,    0.1,     0.2,     0.2]), # Overall probability of developing severe symptoms
    death_probs   = np.array([0.001,   0.001,   0.005,   0.007,   0.009,   0.01,    0.01,   0.05,     0.07,    0.09]), # Overall probability of dying
        )

default_flock_prognoses = dict(
    breeds = np.array(['duck', 'broiler', 'layer'], dtype=default_str),
    sus_ORs = np.array([1.00, 1.00, 1.00]),
    trans_ORs = np.array([1.00, 1.00, 1.00]),
    baseline_symptomatic_rate = np.array([0.001, 0.001, 0.001]),
    mean_symptomatic_rate_increase = np.array([0.005, 0.005, 0.005]),
    baseline_mortality_rate = np.array([0.001, 0.001, 0.001]),
    mean_mortality_rate_increase = np.array([0.002, 0.002, 0.002]),
    baseline_water_rate = np.array([1.00, 1.00, 1.00]),
    mean_water_rate_increase = np.array([1.00, 1.00, 1.00]),
)

default_barn_prognoses = dict(
    sus_ORs = 1.00,
    trans_ORs = 1.00
)

default_water_prognoses = dict(
    sus_ORs = 1.00,
    trans_ORs = 1.00
)

default_production_cycle = dict(
    breeds = np.array(['duck', 'broiler', 'layer'], dtype=default_str),
    cycle_dur = [dict(dist = 'normal_pos', par1 = 600, par2 = 50),
                 dict(dist = 'normal_pos', par1 = 45, par2 = 5),
                 dict(dist = 'normal_pos', par1 = 150, par2=25)],
    flock_size = [dict(dist = 'normal_pos', par1 = 1000, par2 = 100),
                dict(dist = 'normal_pos', par1 = 20000, par2 = 1000),
                dict(dist = 'normal_pos', par1 = 10000, par2 = 500)]
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