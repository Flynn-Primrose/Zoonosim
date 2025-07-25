'''
Set the parameters for Zoonosim.
'''

import numpy as np
import sciris as sc
from . import defaults as znd
from .Options import options as zno

__all__ = ['make_pars']

def make_pars(set_prognoses = False, version = None, **kwargs):

    '''
    Create the parameters for the simulation. Typically, this function is used
    internally rather than called by the user; e.g. typical use would be to do
    sim = zn.Sim() and then inspect sim.pars, rather than calling this function
    directly.

    Args:
        set_prognoses (bool): whether or not to create prognoses (else, added when the population is created)
        prog_by_age   (bool): whether or not to use age-based severity, mortality etc.
        kwargs        (dict): any additional kwargs are interpreted as parameter names
        version       (str):  if supplied, use parameters from this Zoonosim version

    Returns:
        pars (dict): the parameters of the simulation
    '''

    pars = {}

    # Population pars
    pars['agent_types'] = ['human', 'flock', 'barn', 'water'] # Every type of agent in the model
    pars['n_farms'] = 25 # Number of farms in the simulation. This is used to generate the rest of the population
    pars['pop_size'] = None # The total number of agents in the simulation. This should be equal to the sum of the population sizes of all agent types.
    pars['pop_size_by_type'] = {}
    

    for type in pars['agent_types']:
        pars['pop_size_by_type'][type] = None # This will be set after the population has been created

    pars['pop_scale'] = 1.0 # Scale factor for the population size. We don't use this functionality so this will always be 1.0. It is included for compatibility with modules inherited from Covasim/Pathosim.
    pars['rescale'] = False # We will never use rescaling, so this is always False. It is included for compatibility with modules inherited from Covasim/Pathosim.

    pars['initial_conditions'] = {
        'human': 0, # Number of initial humans exposed
        'flock': 0, # Number of initial flocks exposed
        'barn': 0, # Number of initial contaminated barns
        'water': 0 # Number of initial contaminated water
    }

    # Simulation parameters
    pars['start_day']  = '2025-01-01' # Start day of the simulation
    pars['end_day']    = None         # End day of the simulation
    pars['n_days']     = 150           # Number of days to run, if end_day isn't specified
    pars['rand_seed']  = None            # Random seed, if None, don't reset
    pars['verbose']    = zno.verbose  # Whether or not to display information during the run -- options are 0 (silent), 0.1 (some; default), 1 (default), 2 (everything)


    # Testing par dictionaries. Note: This is used to toggle symptoms and bkg ILI, since we only need symptoms and bkg ILI when we need testing. 
    pars['enable_testobjs'] = False 
    pars['testobjs'] = None  # This naming convention is possibly ambiguous - do NOT supply actually test objects (do that for pars['testing'])


    pars['enable_smartwatches'] = False
    pars['smartwatch_pars'] = {
        'mean_fpr'                  :   0.08, # mean false positive rate
        'use_variable_fpr'          : True, # Whether to use a variable false positive rate
        'day_i'                     : np.arange(-21, 22, 1), #
        'loc'                       : 3.25, # Day of max probability of alert, relative to the day of symptom onset.
        'alpha'                     : 1, # Scales the probability of receiving an alert
        'usage_rate'                : 1, # Out of people who have smartwatches, the amount who use download the alerting app and stick with it.
        'compliance_rate'           :   0.05, # probability of quarantining if a smartwatch detects symptoms
        'participation_rate'        :   0.3,  # proportion of the population that has a smartwatch
    }


    # Network parameters, generally initialized after the population has been constructed
    #pars['contacts']        = None  # The number of contacts per layer; set by reset_layer_pars() below
    pars['dynam_layer']     = None  # Which layers are dynamic; set by reset_layer_pars() below
    pars['beta_layer']      = None  # Transmissibility per layer; set by reset_layer_pars() below
    pars['quar_factor']     = None  # Quarantine multiplier on transmissibility and susceptibility; set by reset_layer_pars() below
    pars['quar_period']     = 14  # Number of days to quarantine for. Assumption based on standard policies


    # Parameters that control settings and defaults for multi-variant runs
    pars['n_variants'] = 1 # The number of variants circulating in the population

    pars['beta'] = {} # The transmissibility of the disease for each agent type.
    pars['beta']['human'] = 0.3 # The transmissibility of the disease for humans. This is a dummy variable!
    pars['beta']['flock'] = 0.7 # The transmissibility of the disease for flocks. This is a dummy variable!
    pars['beta']['barn'] = 0.2 # The transmissibility of the disease for barns. This is a dummy variable!
    pars['beta']['water'] = 0.2 # The transmissibility of the disease for water. This is a dummy variable!

    # Default variant parameters
    # These are all dummy values
    pars['wild'] = dict(
        human = dict(
            rel_beta = 0.33,
            rel_symp_prob = 0.33,
            rel_severe_prob = 0.25,
            rel_death_prob = 0.01,
            rel_asymp_fact = 0.5

        ),
        flock = dict(
            rel_beta = 1.0,
            rel_symp_prob = 0.75,
            rel_severe_prob = 1.0,
            rel_death_prob = 1.0
        ),
        barn = dict(
            rel_beta = 1.0,
            rel_dur_contamination = 1.0
        ),
        water = dict(
            rel_beta = 1.0,
            rel_dur_contamination = 1.0
        )
    )


    # Basic disease transmission parameters
    pars['transmission_pars'] = {}

    pars['transmission_pars']['human'] = {
        'beta_dist': dict(dist='neg_binomial', par1=1.0, par2=0.45, step=0.01), # Distribution to draw individual level transmissibility
        'viral_dist':dict(frac_time=0.3, load_ratio=2, high_cap=4), # The time varying viral load (transmissibility)
        'enable_vl':True, # Specifies whether we should use the updated viral load calculation; False = use native calculation
        'viral_levels':dict(min_vl=0.75, max_vl=2) # Specifies the range within which viral load should be scaled so it can contribute to relative transmissibility
    }

    pars['transmission_pars']['flock'] = {
        'beta_dist': dict(dist='neg_binomial', par1=1.0, par2=0.45, step=0.01), # NOTE: Dummy variables
    }

    pars['transmission_pars']['barn'] = {
        'beta_dist': dict(dist='neg_binomial', par1=1.0, par2=0.45, step=0.01), # NOTE: Dummy variables
    }

    pars['transmission_pars']['water'] = {
        'beta_dist': dict(dist='neg_binomial', par1=1.0, par2=0.45, step=0.01), # NOTE: Dummy variables
    }

    # Number of imported cases per day for each agent type
    pars['n_imports']  = {
        'human': 0,
        'flock': 0,
        'barn' : 0.1,
        'water': 0.1,
    } 




    # Parameters used to calculate immunity
    pars['immunity_pars'] = {}

    pars['immunity_pars']['human'] = {
        'use_waning': True, # Whether or not to use waning immunity; set to True for humans by default
        'nab_init':dict(dist='normal', par1=0, par2=2),  # Parameters for the distribution of the initial level of log2(nab) following natural infection, taken from fig1b of https://doi.org/10.1101/2021.03.09.21252641
        'nab_decay':dict(form='nab_growth_decay', growth_time=21, decay_rate1=np.log(2) / 50, decay_time1=150, decay_rate2=np.log(2) / 250, decay_time2=365), # NOTE: I have no idea where this comes from
        'nab_kin': None, # Constructed during sim initialization using the nab_decay parameters
        'nab_boost': 1.5, # Multiplicative factor applied to a person's nab levels if they get reinfected. No data on this, assumption.
        'nab_eff':dict(alpha_inf=1.08, alpha_inf_diff=1.812, beta_inf=0.967, alpha_symp_inf=-0.739, beta_symp_inf=0.038, alpha_sev_symp=-0.014, beta_sev_symp=0.079), # Parameters to map nabs to efficacy
        'rel_imm_symp':dict(asymp=0.85, mild=1, severe=1.5), # Relative immunity from natural infection varies by symptoms. Assumption.
        'immunity':None, # Matrix of immunity and cross-immunity factors, set by init_immunity() in immunity.py
        'trans_redux':0.59 # Reduction in transmission for breakthrough infection
    }
    pars['immunity_pars']['flock'] = {'use_waning': False} # No waning immunity for flock agents
    pars['immunity_pars']['barn'] = {'use_waning': False} # No waning immunity for barn agents
    pars['immunity_pars']['water'] = {'use_waning': False} # No waning immunity for water agents




    # Duration Parameters

    pars['dur'] = {}

    pars['dur']['human'] = {

        # Duration: disease progression
        'exp2inf': dict(dist='lognormal_int', par1=3.0, par2=1.5), # Duration from exposed to infectious
        'inf2sym': dict(dist='lognormal_int', par1=1.5, par2=0.5), # Duration from infectious to symptomatic
        'sym2sev': dict(dist='lognormal_int', par1=5.0, par2=2.0), # Duration from symptomatic to severe symptoms

        # Duration: Recovery
        'asym2rec': dict(dist='lognormal_int', par1=8.0,  par2=2.0), # Duration for asymptomatic people to recover
        'mild2rec': dict(dist='lognormal_int', par1=8.0,  par2=2.0), # Duration for people with mild symptoms to recover
        'sev2rec': dict(dist='lognormal_int', par1=14.0, par2=6.0), # Duration for people with severe symptoms to recover
        'sev2die': dict(dist='lognormal_int', par1=10.0, par2=5.0), # Duration from critical symptoms to death, 18.8 days total
    }

    pars['dur']['flock'] = {
        'exp2inf': dict(dist='lognormal_int', par1=1.5, par2=0.5), # Duration from exposed to infectious. NOTE: This data is just a guess, and should be replaced with real data
        'inf2out': dict(dist='lognormal_int', par1=2.0, par2=1.0), # Duration from infectious to recovery/removal. NOTE: This data is just a guess, and should be replaced with real data
        'susp2res': dict(dist='lognormal_int', par1=5.0, par2=1.0), # Duration from suspicion to a definitive test result. NOTE: This data is just a guess, and should be replaced with real data
    }
    pars['dur']['barn'] = {
        'contamination': dict(dist='lognormal_int', par1=14, par2=5.0), # Duration of contamination. NOTE: This data is just a guess, and should be replaced with real data
        'composting': dict(dist='lognormal_int', par1=7.0, par2=1.5), # Duration of composting. NOTE: This data is just a guess, and should be replaced with real data
        'cleaning': dict(dist='lognormal_int', par1=7.0, par2=1.5), # Duration of cleaning process. NOTE: This data is just a guess, and should be replaced with real data
    }
    pars['dur']['water'] = {
        'contamination': dict(dist='lognormal_int', par1=14, par2=5.0), # Duration of contamination. NOTE: This data is just a guess, and should be replaced with real data
    }




    # Prognoses
    pars['prognoses'] = {}
    pars['prognoses']['human'] = relative_human_prognoses(znd.default_human_prognoses)
    pars['prognoses']['flock'] = relative_flock_prognoses(znd.default_flock_prognoses)
    pars['prognoses']['barn'] = relative_barn_prognoses(znd.default_barn_prognoses)
    pars['prognoses']['water'] = relative_water_prognoses(znd.default_water_prognoses)

    pars['production_cycle'] = znd.default_production_cycle
    
    # Background ILI parameters
    pars['bkg_ILI'] = znd.default_bkg_ILI # Weekly background ILI rates, used to infect human agents with ILI
    pars['Avian_to_ILI'] = znd.default_Avian_to_ILI # should we allow humans infected with avian influenza to also be infected with ILI? This is a boolean value, default is False

        # Events and interventions
    pars['interventions'] = []   # The interventions present in this simulation; populated by the user
    pars['surveillance'] = []    # The surveillance systems present in this simualtion; populated by the user
    pars['testing'] = []         # The testing systems present in this simualtion; populated by the user. These can be populated externally, or internally by supplying a parameter dictionary. 
    pars['analyzers']     = []   # Custom analysis functions; populated by the user
    pars['timelimit']     = None # Time limit for the simulation (seconds)
    pars['stopping_func'] = None # A function to call to stop the sim partway through

    # Health system parameters
    pars['n_beds_hosp']    = None # The number of hospital (adult acute care) beds available for severely ill patients (default is no constraint)
    pars['no_hosp_factor'] = 2.0  # Multiplier for how much more likely severely ill people are to become critical if no hospital beds are available

        # Handle vaccine and variant parameters
    pars['vaccine_pars'] = {} # Vaccines that are being used; populated during initialization
    pars['vaccine_map']  = {} #Reverse mapping from number to vaccine key
    pars['variants']     = [] # Additional variants of the virus; populated by the user, see immunity.py
    pars['variant_map']  = {0:'wild'} # Reverse mapping from number to variant key
    pars['variant_pars'] = dict(wild={}) # Populated just below
    for sp in znd.variant_pars:
        if sp in pars['wild'].keys():
            pars['variant_pars']['wild'][sp] = pars['wild'][sp]

    # Update with any supplied parameter values and generate things that need to be generated
    pars.update(kwargs)
    # reset_layer_pars(pars)
    #if set_prognoses: # If not set here, gets set when the population is initialized
        # pars['prognoses'] = get_prognoses(pars['prog_by_age'], version=version) # Default to age-specific prognoses NOTE: this will have to change to accommodate multiple agent types

    return pars

# Define which parameters need to be specified as a dictionary by layer -- define here so it's available at the module level for sim.py
layer_pars = ['beta_layer', 'dynam_layer', 'quar_factor']

def reset_layer_pars(pars, layer_keys=None, force=False):
    '''
    Helper function to set layer-specific parameters. If layer keys are not provided,
    then set them based on the population type. This function is not usually called
    directly by the user, although it can sometimes be used to fix layer key mismatches
    (i.e. if the contact layers in the population do not match the parameters). More
    commonly, however, mismatches need to be fixed explicitly.

    NOTE:   fb = Flock-Barn
            bw = Barn-Water
            fw = Flock-Water
            hb = Human-Barn
            hf = Human-Flock
            hh = Human-Human

    Args:
        pars (dict): the parameters dictionary
        layer_keys (list): the layer keys of the population, if available
        force (bool): reset the parameters even if they already exist
    '''

    layer_defaults = dict(
        beta_layer = dict(fb=1.0, bw=1.0, fw = 1.0, hb = 1.0, hf = 1.0, hh = 1.0), # Transmissibility per layer -- set to one by default
        dynam_layer = dict(fb=0.0, bw=0.0, fw = 0.0, hb = 0.0, hf = 0.0, hh = 0.0), # Dynamic layer -- set to zero by default
        quar_factor = dict(fb=1.0, bw=1.0, fw = 1.0, hb = 1.0, hf = 1.0, hh = 1.0), # Quarantine factor -- set to one by default
    )

    default_val = 1.0 # Default value for parameters that are not specified
    default_layer_keys = list(layer_defaults['beta_layer'].keys()) # Get the default layer keys from the first parameter

        # Actually set the parameters
    for pkey in layer_pars:
        par = {} # Initialize this parameter

        # If forcing, we overwrite any existing parameter values
        if force:
            par_dict = layer_defaults[pkey] # Just use defaults
        else:
            par_dict = sc.mergedicts(layer_defaults[pkey], pars.get(pkey, None)) # Use user-supplied parameters if available, else default

        # Figure out what the layer keys for this parameter are (may be different between parameters)
        if layer_keys:
            par_layer_keys = layer_keys # Use supplied layer keys
        else:
            par_layer_keys = list(sc.odict.fromkeys(default_layer_keys + list(par_dict.keys())))  # If not supplied, use the defaults, plus any extra from the par_dict; adapted from https://www.askpython.com/python/remove-duplicate-elements-from-list-python

        # Construct this parameter, layer by layer
        for lkey in par_layer_keys: # Loop over layers
            par[lkey] = par_dict.get(lkey, default_val) # Get the value for this layer if available, else use the default for random
        pars[pkey] = par # Save this parameter to the dictionary

    return 

def relative_human_prognoses(prognoses):
    '''
    Convenience function to revert absolute human prognoses into relative (conditional)
    ones. Internally, Zoonosim uses relative prognoses.
    '''
    out = sc.dcp(prognoses)
    out['death_probs']  /= out['severe_probs']   # Conditional probability of dying, given critical symptoms
    out['severe_probs'] /= out['symp_probs']   # Conditional probability of symptoms becoming severe, given symptomatic
    return out

def relative_flock_prognoses(prognoses):
        '''
    Convenience function to revert absolute flock prognoses into relative (conditional)
    ones. Internally, Zoonosim uses relative prognoses.
    '''
        out = sc.dcp(prognoses)
        return(out)

def relative_barn_prognoses(prognoses):
        '''
    Convenience function to revert absolute flock prognoses into relative (conditional)
    ones. Internally, Zoonosim uses relative prognoses.
    '''
        out = sc.dcp(prognoses)
        return(out)

def relative_water_prognoses(prognoses):
        '''
    Convenience function to revert absolute flock prognoses into relative (conditional)
    ones. Internally, Zoonosim uses relative prognoses.
    '''
        out = sc.dcp(prognoses)
        return(out)

def absolute_human_prognoses(prognoses):
    '''
    Convenience function to revert relative (conditional) prognoses into absolute
    ones. Used to convert internally used relative prognoses into more readable
    absolute ones.

    **Example**::

        sim = cv.Sim()
        abs_progs = cv.parameters.absolute_prognoses(sim['prognoses'])
    '''
    out = sc.dcp(prognoses)
    out['severe_probs'] *= out['symp_probs']   # Absolute probability of severe symptoms
    out['death_probs']  *= out['severe_probs']   # Absolute probability of dying
    return out


def get_variant_choices():
    '''
    Define valid pre-defined variant names
    '''
    # List of choices currently available: new ones can be added to the list along with their aliases
    choices = {
        'HPAI':  ['HPAI', 'hpai'], # Highly pathogenic avian influenza
        'LPAI': ['LPAI', 'lpai'], # Low pathogenic avian influenza

    }
    mapping = {name:key for key,synonyms in choices.items() for name in synonyms} # Flip from key:value to value:key
    return choices, mapping


def get_vaccine_choices():
    '''
    Define valid pre-defined vaccine names
    '''
    # List of choices currently available: new ones can be added to the list along with their aliases
    choices = {
        'default': ['default', None],
    }
    mapping = {name:key for key,synonyms in choices.items() for name in synonyms} # Flip from key:value to value:key
    return choices, mapping


def _get_from_pars(pars, default=False, key=None, defaultkey='default'):
    ''' Helper function to get the right output from vaccine and variant functions '''

    # If a string was provided, interpret it as a key and swap
    if isinstance(default, str):
        key, default = default, key

    # Handle output
    if key is not None:
        try:
            return pars[key]
        except Exception as E:
            errormsg = f'Key "{key}" not found; choices are: {sc.strjoin(pars.keys())}'
            raise sc.KeyNotFoundError(errormsg) from E
    elif default:
        return pars[defaultkey]
    else:
        return pars


def get_variant_pars(default=False, variant=None):
    '''
    Define the default parameters for the different variants
    '''
    pars = dict(
        

        HPAI = dict(
            human = dict(
                rel_beta        = 1.0, # Default values
                #rel_gamma      = 1.0, # Default values
                rel_symp_prob   = 1.0, # Default values
                rel_severe_prob = 1.0, # Default values
                rel_death_prob  = 1.0, # Default values
                rel_asymp_fact  = 1.0
            ),
            flock = dict(
                rel_beta        = 1.0, # Default values
                rel_gamma      = 1.0, # Default values
                rel_symp_prob   = 1.0, # Default values
                rel_death_prob  = 1.0, # Default values
            ),
            barn = dict(
                rel_beta        = 1.0, # Default values
                # rel_gamma      = 1.0, # Default values, not used for barns
                rel_dur_contamination = 1.0, # Default values
            ),
            water = dict(
                rel_beta        = 1.0, # Default values
                # rel_gamma      = 1.0, # Default values, not used for water
                rel_dur_contamination = 1.0, # Default values
            ),
        ),

        LPAI = dict(
            human = dict(
                rel_beta        = 1.0, # guessed values
                # rel_gamma      = 1.0, # guessed value, not used for humans
                rel_symp_prob   = 0.1, # guess but LPAI is less severe than HPAI
                rel_severe_prob = 0.25, # guess but LPAI is less severe than HPAI
                rel_death_prob  = 0.25, # guess but LPAI is less severe than HPAI
                rel_asymp_fact  = 0.5
            ),
            flock = dict(
                rel_beta        = 1.0, # guessed values
                rel_gamma       = 0.5, # guessed values
                rel_symp_prob   = 0.25, # guess but LPAI is less severe than HPAI
                rel_death_prob  = 0.05, # guess but LPAI is less severe than HPAI
            ),
            barn = dict(
                rel_beta        = 1.0, # guessed values
                # rel_gamma      = 1.0, # guessed values, not used for barns
                rel_dur_contamination = 0.5, # guess but LPAI is less severe than HPAI
            ),
            water = dict(
                rel_beta        = 1.0, # guessed values
                # rel_gamma      = 1.0, # guessed values, not used for water
                rel_dur_contamination = 0.5, # guess but LPAI is less severe than HPAI
            ),
        ),
    )

    return _get_from_pars(pars, default, key=variant, defaultkey='wild')


def get_cross_immunity(default=False, variant=None):
    '''
    Get the cross immunity between each variant in a sim
    '''
    pars = dict(

        HPAI = dict(
            HPAI  = 1.0, # Default for own-immunity
            LPAI  = 0.5, # Guess for cross-immunity

        ),

        LPAI = dict(
            HPAI  = 0.5, # Guess for cross-immunity
            LPAI  = 1.0, # Default for own-immunity
        ),
    )

    return _get_from_pars(pars, default, key=variant, defaultkey='wild')


def get_vaccine_variant_pars(default=False, vaccine=None):
    '''
    Define the effectiveness of each vaccine against each variant
    '''
    pars = dict(

        default = dict(
            HPAI  = 1.0, # Default values
            LPAI  = 1.0, # Default values
        ),
    )

    return _get_from_pars(pars, default=default, key=vaccine)


def get_vaccine_dose_pars(default=False, vaccine=None):
    '''
    Define the parameters for each vaccine
    '''

    pars = dict(

        default = dict(
            nab_init  = dict(dist='normal', par1=0, par2=2), # Initial distribution of NAbs
            nab_boost = 2, # Factor by which a dose increases existing NABs
            doses     = 1, # Number of doses for this vaccine
            interval  = None, # Interval between doses
        ),
    )

    return _get_from_pars(pars, default, key=vaccine)