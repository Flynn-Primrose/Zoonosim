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
default_bool  = np.bool_
nbfloat       = nb.float64
nbint         = nb.int64
nbbool        = nb.bool_

verbose = 0.1
warnings = 'warn' # must be one of 'warn', 'print', or 'error'

safe_opts = [1, '1', 'safe'] 
full_opts = [2, '2', 'full'] 

numba_parallel = str(os.getenv('ZOONOSIM_NUMBA_PARALLEL', 'none'))
numba_cache = bool(int(os.getenv('ZOONOSIM_NUMBA_CACHE', 1)))

# Define the 'overview plots', i.e. the most useful set of plots to explore different aspects of a simulation
overview_plots = [
    'new_human_infections',
    'new_flock_infectious',
    'new_barn_contaminated',
    'new_water_contaminated',
]

overview_variant_plots = [
    'new_human_infections_by_variant',
    'new_flock_infections_by_variant',
    'new_barn_contaminated_by_variant',
    'new_water_contaminated_by_variant',
]

human_plots = sc.odict({
    'Total counts': [
        'n_human_infectious',
        'n_human_exposed',
        'n_human_symptomatic',
        'n_human_severe',
        'n_human_dead',
    ],
    'Cumulative counts': [
        'cum_human_infections',
        'cum_human_exposed',
        'cum_human_symptomatic',
        'cum_human_severe',
        'cum_human_dead',
    ],
    'Daily counts': [
        'new_human_infections',
        'new_human_exposed',
        'new_human_symptomatic',
        'new_human_severe',
        'new_human_dead',
    ],
})

flock_plots = sc.odict({
    'Total counts': [
        'n_flock_infectious',
        'n_flock_exposed',
        'n_flock_suspected',
        'n_flock_quarantined',
    ],
    'Cumulative counts': [
        'cum_flock_infectious',
        'cum_flock_exposed',
        'cum_flock_suspected',
        'cum_flock_quarantined',
    ],
    'Daily counts': [
        'new_flock_infectious',
        'new_flock_exposed',
        'new_flock_suspected',
        'new_flock_quarantined',
    ],
})

barn_plots = sc.odict({
    'Total counts': [
        'n_barn_contaminated',
        'n_barn_uncontaminated',
        'n_barn_composting',
        'n_barn_cleaning',
    ],
    'Cumulative counts': [
        'cum_barn_contaminated',
        'cum_barn_uncontaminated',
        'cum_barn_composting',
        'cum_barn_cleaning',
    ],
    'Daily counts': [
        'new_barn_contaminated',
        'new_barn_uncontaminated',
        'new_barn_composting',
        'new_barn_cleaning',
    ],
})

water_plots = sc.odict({
    'Total counts': [
        'n_water_contaminated',
        'n_water_uncontaminated',
    ],
    'Cumulative counts': [
        'cum_water_contaminated',
        'cum_water_uncontaminated',
    ],
    'Daily counts': [
        'new_water_contaminated',
        'new_water_uncontaminated',
    ],
})

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
                    'n_human_infectious',
                    'n_flock_infectious',
                    'n_barn_contaminated',
                    'n_water_contaminated',
                ],
                'Cumulative counts': [
                    'cum_human_infections',
                    'cum_flock_infectious',
                    'cum_barn_contaminated',
                    'cum_water_contaminated',
                ],
                'Daily counts': [
                    'new_human_infections',
                    'new_flock_infectious',
                    'new_barn_contaminated',
                    'new_water_contaminated',
                ],
                'Health outcomes': [
                    'cum_human_infectious',
                    'cum_human_symptomatic',
                    'cum_human_severe',
                    'cum_human_dead',
                ],
            })

        else: # pragma: no cover
            plots = sc.odict({
                'Cumulative infections': [
                    'cum_human_infections',
                    'cum_flock_infectious',
                    'cum_barn_contaminated',
                    'cum_water_contaminated',
                ],
                'New infections per day': [
                    'new_human_infections',
                    'new_flock_infectious',
                    'new_barn_contaminated',
                    'new_water_contaminated',
                ],
                'Cumulative deaths': [
                    'cum_human_dead',
                ],
            })

    # Show an overview
    elif which == 'overview': # pragma: no cover
        plots = sc.dcp(overview_plots)

    # Plot absolutely everything
    elif which == 'all': # pragma: no cover
        plots = sim.result_keys('all')

    elif which == 'human': # pragma: no cover
        plots = sc.dcp(human_plots)

    elif which == 'flock': # pragma: no cover
        plots = sc.dcp(flock_plots)

    elif which == 'barn': # pragma: no cover
        plots = sc.dcp(barn_plots)

    elif which == 'water': # pragma: no cover
        plots = sc.dcp(water_plots)

    # Show an overview plus variants
    elif 'overview' in which and 'variant' in which: # pragma: no cover
        plots = sc.dcp(overview_plots) + sc.dcp(overview_variant_plots)

    # Show default but with variants
    elif which.startswith('variant'): # pragma: no cover
        if is_sim:
            plots = sc.odict({
                'Cumulative infections by variant': [
                    'cum_human_infections_by_variant',
                    'cum_flock_infectious_by_variant',
                    'cum_barn_contaminated_by_variant',
                    'cum_water_contaminated_by_variant',
                ],
                'New infections by variant': [
                    'new_human_infections_by_variant',
                    'new_flock_infectious_by_variant',
                    'new_barn_contaminated_by_variant',
                    'new_water_contaminated_by_variant',
                ],
                'Health outcomes': [
                    'cum_human_infections',
                    'cum_human_severe',
                    'cum_human_dead',
                ],
            })

        else: # pragma: no cover
            plots = sc.odict({
                    'Cumulative infections by variant': [
                        'cum_human_infections_by_variant',
                        'cum_flock_infectious_by_variant',
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

    else: # pragma: no cover
        errormsg = f'The choice which="{which}" is not supported: choices are "default", "overview", "all", "variant", "overview-variant", or "seir", along with any result key (see sim.results_keys(\'all\') for options)'
        raise ValueError(errormsg)

    return plots



def get_default_colors(agent_type):
    '''
    Specify plot colors -- used in sim.py.

    Args:
        agent_type (str): the type of agent to get colors for; one of 'human', 'flock', 'barn', or 'water'

    Colour palette: 
    '''
    c = sc.objdict()
    match agent_type.lower():
        case 'human':
            c.susceptible = "#00FF00FF" # Green
            c.exposed = "#00FF00A6" # 
            c.symptomatic = "#00FF0067" # 
            c.infectious = "#00FF001F" # 
            c.infections = '#00FF002F' # 
            c.reinfections = '#00FF001F' # 
            c.severe = "#00FF9D" # 
            c.recovered = "#00CCFF" # 
            c.dead = "#00FF0000" # 
            c.known_dead = '#00FF0000' # 
            c.diagnosed = '#00FF0000' # 
            c.tested = "#00FF0000" # 
            c.quarantined = "#00FF0000" # 
            c.vaccinated = "#00FF0000" # 
            c.doses = "#00FF0000" # 
            c.exposed_by_variant = '#00FF00A6' # 
            c.infectious_by_variant = '#00FF001F' # 
            c.infections_by_variant = '#00FF002F' # 
            c.symptomatic_by_variant = '#00FF0067' # 
            c.severe_by_variant = '#00FF9D' # 
        case 'flock':
            c.susceptible = '#FFFF00' # 
            c.exposed = "#FFFF00BC" # Yellow
            c.infectious = "#FFFF0090" # 
            c.suspected = "#FFD000" # 
            c.quarantined = "#C8FF00" # 
            c.exposed_by_variant = '#FFFF00BC' # 
            c.infectious_by_variant = '#FFFF0090' # 
        case 'barn':
            c.uncontaminated = "#33FF00" # 
            c.contaminated = "#FF0000FF" # Red
            c.composting = "#D18A07BA"
            c.cleaning = "#1D20ACDA"
            c.contaminated_by_variant = "#00FFFF3A" # Red
        case 'water':
            c.uncontaminated = "#0004FF" # 
            c.contaminated = "#0004FF84" # 
            c.contaminated_by_variant = "#0004FF1F" # 
        case 'default':
            c.default = '#000000' # Black
        case _:
            errormsg = f'Unknown agent type "{agent_type}" for get_default_colors()'
            raise ValueError(errormsg)
    return c


default_pop_pars = {
    'avg_barns_per_farm': 3.0,
    'avg_humans_per_barn': 2.0,
    'avg_flock_size': 1000.0,
    'avg_water_per_farm': 1.0,
}

default_flock_breeds = ['duck', 'layer', 'broiler'] # Breeds of flocks present
default_flock_breed_freqs = [0.1, 0.2, 0.7] # frequency of the different breed types
default_suspicious_mortality_rate = 0.005 # mortality rate that triggers a flock to be considered suspicious
default_suspicious_symptomatic_rate = 0.005 # symptomatic rate that triggers a flock to be considered suspicious
default_suspicious_consumption_rate = 2 # water rate that triggers a flock to be considered suspicious (L/bird/day)

# Parameters that can vary by variant
variant_pars = [
    'rel_beta',
    'rel_symp_prob',
    'rel_severe_prob',
    'rel_death_prob',
]

default_human_prognoses = dict(
    age_cutoffs   = np.array([0,       10,      20,      30,      40,      50,      60,      70,      80,      90,]),     # Age cutoffs (lower limits)
    sus_ORs       = np.array([0.25,    0.50,    0.75,    1.0,     1.25,    1.50,    1.75,    2.0,     2.25,    2.50]),    # Odds ratios for relative susceptibility 
    trans_ORs     = np.array([0.10,    0.10,    0.10,    0.10,    0.10,    0.10,    0.10,    0.10,    0.10,    0.10]),    # Odds ratios for relative transmissibility
    comorbidities = np.array([1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00]),    # Comorbidities by age -- set to 1 by default since already included in disease progression rates
    symp_probs    = np.array([0.50,    0.55,    0.60,    0.65,    0.70,    0.75,    0.80,    0.85,    0.90,    0.90]),    # Overall probability of developing symptoms 
    severe_probs  = np.array([0.01,    0.01,    0.05,    0.05,    0.05,    0.07,    0.09,    0.1,     0.2,     0.2]), # Overall probability of developing severe symptoms
    death_probs   = np.array([0.001,   0.001,   0.005,   0.007,   0.009,   0.01,    0.01,   0.05,     0.07,    0.09]), # Overall probability of dying
        )

default_flock_prognoses = dict(
    breeds = np.array(['duck', 'broiler', 'layer'], dtype=default_str),
    sus_ORs = np.array([2.00, 1.00, 0.75]),
    trans_ORs = np.array([1.00, 1.00, 1.00]),
    baseline_symptomatic_rate = np.array([0.001, 0.001, 0.001]),
    mean_symptomatic_rate_increase = np.array([0.001, 0.0005, 0.0001]),
    baseline_mortality_rate = np.array([0.001, 0.001, 0.001]),
    mean_mortality_rate_increase = np.array([0.005, 0.002, 0.002]),
    baseline_water_rate = np.array([1.00, 1.00, 1.00]),
    mean_water_rate_increase = np.array([1.50, 1.00, 0.75]),
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
    'suspected':   'Number of suspected flocks',
    'infectious':  'Number of infectious flocks',
    'quarantined': 'Number of quarantined flocks',
}

barn_stocks = {
    'composting': 'Number of barns that have finished composting',
    'cleaning': 'Number of cleaned barns',
    'contaminated': 'Number of contaminated barns',
    'uncontaminated': 'Number of uncontaminated barns',
}

water_stocks = {
    'contaminated': 'Number of contaminated waters',
    'uncontaminated': 'Number of uncontaminated waters',
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
    'exposed':    'human exposures',
    'reinfections': 'human reinfections',
    'infections':   'human infections',
    'infectious':   'infectious humans',
    'symptomatic':  'symptomatic human cases',
    'severe':       'severe human cases',
    'recovered':   'human recoveries',
    'dead':       'human deaths',
    'tested':        'human tests',
    'diagnosed':    'human diagnoses',
    'known_dead': 'known human deaths',
    'quarantined':  'human quarantines started',
    'doses':        'human vaccine doses',
    'vaccinated':   'vaccinated humans'
}
new_human_flows = [f'new_{key}' for key in human_flows.keys()]
cum_human_flows = [f'cum_{key}' for key in human_flows.keys()]

flock_flows = {
    'exposed':     'exposed flocks',
    'infectious':   'infectious flocks',
    'suspected' : 'flocks suspected of being infectious',
    'quarantined': 'flocks quarantined',
}
new_flock_flows = [f'new_{key}' for key in flock_flows.keys()]
cum_flock_flows = [f'cum_{key}' for key in flock_flows.keys()]

barn_flows = {
    'contaminated': 'contaminated barns',
    'uncontaminated': 'uncontaminated barns',
    'composting' : 'barns finished composting',
    'cleaning' : 'cleaned barns'
}
new_barn_flows = [f'new_{key}' for key in barn_flows.keys()]
cum_barn_flows = [f'cum_{key}' for key in barn_flows.keys()]

water_flows = {
    'contaminated': 'contaminated waterbodies',
    'uncontaminated' : 'uncontaminated waterbodies'
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
    'infectious_by_variant':  'infectious flocks by variant',
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

