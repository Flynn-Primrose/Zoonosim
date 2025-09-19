import Zoonosim as zn
import numpy as np

new_param = {
    'n_farms': 100,
    'n_days': 1000,
    'beta': {
        'human': 1.0,
        'flock': 1.0,
        'barn': 1.0,
        'water': 1.0
    },
    'wild': {
        'human': {
            'rel_beta': 1.0,
            'rel_symp_prob': 1.0,
            'rel_sev_prob': 1.0,
            'rel_death_prob': 1.0,
            'rel_asymp_fact': 1.0
        },
        'flock': {
            'rel_beta': 1.0,
            'rel_symp_prob': 1.0,
            'rel_sev_prob': 1.0,
            'rel_death_prob': 1.0,
        },
        'barn': {
            'rel_beta': 1.0,
            'rel_dur_contamination': 1.0
        },
        'water': {
            'rel_beta': 1.0,
            'rel_dur_contamination': 1.0
        }
    },
    'transmission_pars': {
        'human': {
            'beta_dist': dict(dist = 'neg_binomial', par1 = 1.0, par2 = 5.0, step = 0.01),
            'viral_dist': None, # this should work since we never actually use this distribution
            'enable_vl': True,
            'viral_levels': dict(min_vl = 0.5, max_vl = 2.0)
        },
        'flock': {
            'beta_dist': dict(dist = 'neg_binomial', par1 = 1.0, par2 = 5.0, step = 0.01),
        },
        'barn': {
            'beta_dist': dict(dist = 'neg_binomial', par1 = 1.0, par2 = 5.0, step = 0.01),
        },
        'water': {
            'beta_dist': dict(dist = 'neg_binomial', par1 = 1.0, par2 = 5.0, step = 0.01),
        }
    },
    'immunity_pars': {
        'human': {
        'use_waning': True, # Whether or not to use waning immunity; set to True for humans by default
        'nab_init':dict(dist='normal', par1=0, par2=2),  
        'nab_decay':dict(form='nab_growth_decay', growth_time=21, decay_rate1=np.log(2) / 50, decay_time1=150, decay_rate2=np.log(2) / 250, decay_time2=365),
        'nab_kin': None, # Constructed during sim initialization using the nab_decay parameters
        'nab_boost': 2.0,
        'nab_eff':dict(alpha_inf=1.00, alpha_inf_diff=2.0, beta_inf=1.0, alpha_symp_inf=-1.0, beta_symp_inf=0.1, alpha_sev_symp=-0.01, beta_sev_symp=0.1), # Parameters to map nabs to efficacy
        'rel_imm_symp':dict(asymp=1.0, mild=1.0, severe=1.0), # Relative immunity from natural infection varies by symptoms. Assumption.
        'immunity':None, # Matrix of immunity and cross-immunity factors, set by init_immunity() in immunity.py
        'trans_redux':1.0 # Reduction in transmission for breakthrough infection
        },
        'flock': {
            'use_waning': False, # Whether or not to use waning immunity; set to False for flocks by default
        },
        'barn': {
            'use_waning': False, # Whether or not to use waning immunity; set to False for barns by default
        },
        'water': {
            'use_waning': False, # Whether or not to use waning immunity; set to False for water by default
        }
    },
    'dur': {
        'human': {
            # Duration: disease progression
            'exp2inf': dict(dist='lognormal_int', par1=14.0, par2=2.0), # Duration from exposed to infectious
            'inf2sym': dict(dist='lognormal_int', par1=7.0, par2=1.0), # Duration from infectious to symptomatic
            'sym2sev': dict(dist='lognormal_int', par1=7.0, par2=1.0), # Duration from symptomatic to severe symptoms

            # Duration: Recovery
            'asym2rec': dict(dist='lognormal_int', par1=7.0,  par2=1.0), # Duration for asymptomatic people to recover
            'mild2rec': dict(dist='lognormal_int', par1=7.0,  par2=1.0), # Duration for people with mild symptoms to recover
            'sev2rec': dict(dist='lognormal_int', par1=14.0, par2=2.0), # Duration for people with severe symptoms to recover
            'sev2die': dict(dist='lognormal_int', par1=21.0, par2=3.0), # Duration from critical symptoms to death
        },
        'flock': {
            'exp2inf': dict(dist='lognormal_int', par1=7.0, par2=1.0), # Duration from exposed to infectious. 
            'inf2out': dict(dist='lognormal_int', par1=7.0, par2=1.0), # Duration from infectious to recovery/removal. 
            'susp2res': dict(dist='lognormal_int', par1=7.0, par2=1.0), # Duration from suspicion to a definitive test result. 
        },
        'barn': {
            'contamination': dict(dist='lognormal_int', par1=3.0, par2=1.5), # Duration of contamination. 
            'composting': dict(dist='lognormal_int', par1=7.0, par2=2.0), # Duration of composting. 
            'cleaning': dict(dist='lognormal_int', par1=7.0, par2=2.0), # Duration of cleaning process. 
        },
        'water': {
            'contamination': dict(dist='lognormal_int', par1=3.0, par2=1.5), # Duration of contamination. 
        }
    },
    'prognosis': {
        'human': {
            'age_cutoffs'   : np.array([0,       10,      20,      30,      40,      50,      60,      70,      80,      90,]),     # Age cutoffs (lower limits)
            'sus_ORs'       : np.array([0.25,    0.50,    0.75,    1.0,     1.25,    1.50,    1.75,    2.0,     2.25,    2.50]),    # Odds ratios for relative susceptibility 
            'trans_ORs'     : np.array([0.50,    0.50,    0.50,    0.50,    0.50,    0.50,    0.50,    0.50,    0.50,    0.50]),    # Odds ratios for relative transmissibility
            'comorbidities' : np.array([1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00]),    # Comorbidities by age -- set to 1 by default since already included in disease progression rates
            'symp_probs'    : np.array([0.60,    0.60,    0.70,    0.70,    0.80,    0.80,    0.90,    0.90,    0.96,    0.96]),    # Overall probability of developing symptoms 
            'severe_probs'  : np.array([0.1,     0.1,     0.2,     0.3,     0.4,     0.5,     0.6,     0.7,     0.8,     0.9]), # Overall probability of developing severe symptoms
            'death_probs'   : np.array([0.1,     0.1,     0.2,     0.3,     0.4,     0.5,     0.6,     0.7,     0.8,     0.9]), # Overall probability of dying
        },
        'flock': {
            'breeds'                           : np.array(['duck', 'broiler', 'layer'], dtype=zn.default_str),
            'sus_ORs'                          : np.array([2.00, 1.00, 0.75]),
            'trans_ORs'                        : np.array([1.00, 1.00, 1.00]),
            'baseline_symptomatic_rate'        : np.array([0.001, 0.001, 0.001]),
            'mean_symptomatic_rate_increase'   : np.array([0.003, 0.003, 0.003]),
            'baseline_mortality_rate'          : np.array([0.001, 0.001, 0.001]),
            'mean_mortality_rate_increase'     : np.array([0.003, 0.003, 0.003]),
            'baseline_water_rate'              : np.array([1.00, 1.00, 1.00]),
            'mean_water_rate_increase'         : np.array([1.50, 1.00, 0.75]),
        },
        'barn': {
            'sus_ORs'   : 1.00,
            'trans_ORs' : 1.00
        },
        'water': {
            'sus_ORs'   : 1.00,
            'trans_ORs' : 1.00
        }
    },
    'production_cycle': {
            'breeds'      : np.array(['duck', 'broiler', 'layer'], dtype=zn.default_str),
            'cycle_dur'   : [dict(dist = 'normal_pos', par1 = 175, par2 = 50),
                         dict(dist = 'normal_pos', par1 = 90, par2 = 25),
                         dict(dist = 'normal_pos', par1 = 175, par2=50)],
            'flock_size'  : [dict(dist = 'normal_pos', par1 = 10000, par2 = 1000),
                            dict(dist = 'normal_pos', par1 = 50000, par2 = 5000),
                            dict(dist = 'normal_pos', par1 = 25000, par2 = 1000)]
    }

}

sim = zn.Sim(new_param)
sim.export_pars('parameter_change.json')