'''
Defines default values for Zoonosim.
'''
import sciris as sc
import numpy as np
import numba as nb

# Specify all externally visible functions this file defines -- other things are available as e.g. zn.defaults.default_int
__all__ = ['default_float', 'default_int', 'nbfloat', 'nbint', 'get_default_colors']

result_float = np.float64 # Always use float64 for results, for simplicity

# Uncomment the following lines to use float32 and int32 instead of float64 and int64, I'm always using float64 and int64 for simplicity.
#default_float = np.float32
#default_int   = np.int32
#nbfloat       = nb.float32
#nbint         = nb.int32

default_str = object # 
default_float = np.float64
default_int   = np.int64
default_bool  = np.bool_
nbfloat       = nb.float64
nbint         = nb.int64
nbbool        = nb.bool_


safe_opts = [1, '1', 'safe'] 
full_opts = [2, '2', 'full'] 

# Define the 'overview plots', i.e. the most useful set of plots to explore different aspects of a simulation
overview_plots = [
    'cum_human_infections',
    'cum_flock_infectious',
    'cum_barn_contaminated',
    'cum_water_contaminated',
]

overview_variant_plots = sc.odict({
    'Human Infections by Variant' : ['new_human_infections_by_variant'],
    'Flock Infections by Variant' : ['new_flock_infections_by_variant'],
    'Contaminated Barns by Variant': ['new_barn_contaminated_by_variant'],
    'Contaminated Water by Variant' : ['new_water_contaminated_by_variant'],
})

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
        'cum_barn_composted',
        'cum_barn_cleaned',
    ],
    'Daily counts': [
        'new_barn_contaminated',
        'new_barn_uncontaminated',
        'new_barn_composted',
        'new_barn_cleaned',
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

testing_plots = sc.odict({
    'new counts': [
        'new_RAT_tests',
        'new_PCR_tests',
        'new_diagnoses_custom',
    ],
    'cumulative counts': [
        'cum_RAT_tests',
        'cum_PCR_tests',
        'cum_diagnoses_custom',
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

    elif which == 'testing': # pragma: no cover
        plots = sc.dcp(testing_plots)

    # Show an overview plus variants
    elif 'overview' in which and 'variant' in which: # pragma: no cover
        plots = sc.dcp(overview_plots) + sc.dcp(overview_variant_plots)

    # Show default but with variants
    elif which.startswith('variant'): # pragma: no cover
        if is_sim:
            plots = sc.odict({
                'Human infections by variant': [
                    #'cum_human_infections_by_variant',
                    'new_human_infections_by_variant',  
                ],
                'Flock infections by variant': [
                    #'cum_flock_infectious_by_variant',
                    'new_flock_infectious_by_variant',  
                ],
                'Barn contaminations by variant': [
                    #'cum_barn_contaminated_by_variant',
                    'new_barn_contaminated_by_variant',
                ],
                'Water contaminations by variant' : [
                    #'cum_water_contaminated_by_variant',
                    'new_water_contaminated_by_variant',
                ]
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
            # Spare colors are '#fb8072', '#80b1d3'
            c.susceptible = '#1f77b4' # 
            c.exposed = '#aec7e8' # 
            c.symptomatic = '#ff7f0e' # 
            c.infectious = '#ffbb78' # 
            c.infections = '#2ca02c' # 
            c.reinfections = '#98df8a' # 
            c.severe = '#d62728' # 
            c.recovered = '#ff9896' # 
            c.dead = '#9467bd' # 
            c.known_dead = '#c5b0d5' # 
            c.pop_nabs = '#8c564b' #
            c.pop_protection = '#c49c94' #
            c.pop_symp_protection = '#e377c2' #
            c.diagnosed = '#f7b6d2' # 
            c.tested = '#7f7f7f' # 
            c.quarantined = '#c7c7c7' # 
            c.vaccinated = '#bcbd22' # 
            c.doses = '#dbdb8d' # 
            c.exposed_by_variant = '#17becf' # 
            c.infectious_by_variant = '#9edae5' # 
            c.infections_by_variant = '#8dd3c7' # 
            c.symptomatic_by_variant = '#ffffb3' # 
            c.severe_by_variant = '#bebada' # 
        case 'flock':
            c.susceptible = '#66c2a5' # 
            c.exposed = '#fc8d62' # 
            c.infectious = '#8da0cb' # 
            c.suspected = '#e78ac3' # 
            c.quarantined = '#a6d854' # 
            c.exposed_by_variant = '#ffd92f' # 
            c.infectious_by_variant = '#e5c494' # 
        case 'barn': 
            c.uncontaminated = '#1b9e77' # 
            c.contaminated = '#d95f02' # 
            c.composting = '#7570b3'
            c.composted = '#e7298a'
            c.cleaning = '#66a61e'
            c.cleaned = '#e6ab02'
            c.contaminated_by_variant = '#a6761d' # 
        case 'water':
            c.uncontaminated = '#fbb4ae' # 
            c.contaminated = '#b3cde3' # 
            c.contaminated_by_variant = '#ccebc5' # 
        case 'misc':
            c.misc1 = '#e41a1c'
            c.misc2 = '#377eb8'
            c.misc3 = '#4daf4a'
            c.misc4 = '#984ea3'
            c.misc5 = '#ff7f00'
            c.misc6 = '#ffff33'
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
default_suspicious_mortality_rate = 0.003 # mortality rate that triggers a flock to be considered suspicious
default_suspicious_symptomatic_rate = 0.001 # symptomatic rate that triggers a flock to be considered suspicious
default_suspicious_consumption_rate = 1.33 # water rate that triggers a flock to be considered suspicious (L/bird/day)

# Parameters that can vary by variant
# TODO: Refactor to include sub pars for each species type.
variant_pars = [
    'human',
    'flock',
    'barn',
    'water'
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

# ILI default parameters
default_bkg_ILI = 0.02 # percentage of the population that is infected with ILI at any given time
default_Avian_to_ILI = False # whether to allow humans with avian influenza to be infected with ILI


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
    'composting': 'Number of barns that are being composted',
    'cleaning': 'Number of barns that are being cleaned',
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

result_flows_smartwatches = {'alerted': 'smartwatch alerts',
                             'alerts_tp': 'true positive smartwatch alerts',
                             'alerts_fn': 'false negative smartwatch alerts',
                             'alerts_tn': 'true negative smartwatch alerts',
                             'alerts_fp': 'false positives smartwatch alerts',
                             'Q_w_i': 'smartwatch users incorrectly quarantined',
                             'Q_w': 'smartwatch users quarantined'}

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
    'composted' : 'barns finished composting',
    'cleaned' : 'cleaned barns'
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
    'exposed_by_variant':     'exposed flocks by variant',
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

