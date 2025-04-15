'''
Set the parameters for Zoonosim.
'''

import numpy as np
import sciris as sc
from . import defaults as znd

__all__ = ['make_pars']

def make_pars(set_prognoses = False, prog_by_age = False, version = None, **kwargs):

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
    pars['n_farms'] = 10 # Number of farms in the simulation. This is used to generate the rest of the population
    pars['pop_size'] = None # The total number of agents in the simulation. This should be equal to the sum of the population sizes of all agent types.
    pars['pop_size_by_type'] = {}

    for type in pars['agent_types']:
        pars['pop_size_by_type'][type] = None # This will be set after the population has been created

    # Simulation parameters
    pars['start_day']  = '2020-03-01' # Start day of the simulation
    pars['end_day']    = None         # End day of the simulation
    pars['n_days']     = 60           # Number of days to run, if end_day isn't specified
    pars['rand_seed']  = 1            # Random seed, if None, don't reset
    pars['verbose']    = znd.verbose  # Whether or not to display information during the run -- options are 0 (silent), 0.1 (some; default), 1 (default), 2 (everything)

        # Rescaling parameters - I'm not sure if we're going to be using these or not.
    # pars['pop_scale']         = 1    # Factor by which to scale the population -- e.g. pop_scale=10 with pop_size=100e3 means a population of 1 million
    # pars['scaled_pop']        = None # The total scaled population, i.e. the number of agents times the scale factor
    # pars['rescale']           = True # Enable dynamic rescaling of the population -- starts with pop_scale=1 and scales up dynamically as the epidemic grows
    # pars['rescale_threshold'] = 0.05 # Fraction susceptible population that will trigger rescaling if rescaling
    # pars['rescale_factor']    = 1.2  # Factor by which the population is rescaled on each step
    # pars['frac_susceptible']  = 1.0  # What proportion of the population is susceptible to infection

    # Network parameters, generally initialized after the population has been constructed
    pars['contacts']        = None  # The number of contacts per layer; set by reset_layer_pars() below
    pars['dynam_layer']     = None  # Which layers are dynamic; set by reset_layer_pars() below
    pars['beta_layer']      = None  # Transmissibility per layer; set by reset_layer_pars() below

    # Basic disease transmission parameters
    pars['transmission_pars'] = {}

    pars['transmission_pars']['human'] = {
        'beta_dist': dict(dist='neg_binomial', par1=1.0, par2=0.45, step=0.01), # Distribution to draw individual level transmissibility; dispersion from https://www.researchsquare.com/article/rs-29548/v1
        'viral_dist':dict(frac_time=0.3, load_ratio=2, high_cap=4), # The time varying viral load (transmissibility); estimated from Lescure 2020, Lancet, https://doi.org/10.1016/S1473-3099(20)30200-0
        'asymp_factor': 1.0,  # Multiply beta by this factor for asymptomatic cases; no statistically significant difference in transmissibility: https://www.sciencedirect.com/science/article/pii/S1201971220302502
        'enable_vl':False, # Specifies whether we should use the updated viral load calculation; False = use native calculation NOTE: Im note sure how to set this yet
        'viral_levels':dict(min_vl=0.75, max_vl=2) # Specifies the range within which viral load should be scaled so it can contribute to relative transmissibility
    }

    pars['transmission_pars']['flock'] = {

    }

    pars['transmission_pars']['barn'] = {
        
    }

    pars['transmission_pars']['water'] = {
        
    }

    # Parameters that control settings and defaults for multi-variant runs
    pars['n_imports']  = {} 
    for type in pars['agent_types']:# Average daily number of imported cases for the given agent type (actual number is drawn from Poisson distribution)
        pars['n_imports'][type] = 0
    pars['n_variants'] = 1 # The number of variants circulating in the population

    # Parameters used to calculate immunity
    pars['immunity_pars'] = {}

    pars['immunity_pars']['human'] = {
        'use_waning':True, # Whether to use dynamically calculated immunity NOTE: I think this will always be true in our case
        'nab_init':dict(dist='normal', par1=0, par2=2),  # Parameters for the distribution of the initial level of log2(nab) following natural infection, taken from fig1b of https://doi.org/10.1101/2021.03.09.21252641
        'nab_decay':dict(form='nab_growth_decay', growth_time=21, decay_rate1=np.log(2) / 50, decay_time1=150, decay_rate2=np.log(2) / 250, decay_time2=365), # NOTE: I have no idea where this comes from
        'nab_kin': None, # Constructed during sim initialization using the nab_decay parameters
        'nab_boost': 1.5, # Multiplicative factor applied to a person's nab levels if they get reinfected. No data on this, assumption.
        'nab_eff':dict(alpha_inf=1.08, alpha_inf_diff=1.812, beta_inf=0.967, alpha_symp_inf=-0.739, beta_symp_inf=0.038, alpha_sev_symp=-0.014, beta_sev_symp=0.079), # Parameters to map nabs to efficacy
        'rel_imm_symp':dict(asymp=0.85, mild=1, severe=1.5), # Relative immunity from natural infection varies by symptoms. Assumption.
        'immunity':None, # Matrix of immunity and cross-immunity factors, set by init_immunity() in immunity.py
        'trans_redux':0.59 # Reduction in transmission for breakthrough infections, https://www.medrxiv.org/content/10.1101/2021.07.13.21260393v
    }
    pars['immunity_pars']['poultry'] = None
    pars['immunity_pars']['barn'] = None
    pars['immunity_pars']['water'] = None


    # Variant-specific disease transmission parameters. By default, these are set up for a single variant, but can all be modified for multiple variants
    pars['rel_beta']        = 1.0 # Relative transmissibility varies by variant

    # Duration Parameters

    pars['dur'] = {}

    pars['dur']['human'] = {

        # Duration: disease progression
        'exp2inf': dict(dist='lognormal_int', par1=4.5, par2=1.5), # Duration from exposed to infectious; see Lauer et al., https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7081172/, appendix table S2, subtracting inf2sym duration
        'inf2sym': dict(dist='lognormal_int', par1=1.1, par2=0.9), # Duration from infectious to symptomatic; see Linton et al., https://doi.org/10.3390/jcm9020538, from Table 2, 5.6 day incubation period - 4.5 day exp2inf from Lauer et al.
        'sym2sev': dict(dist='lognormal_int', par1=6.6, par2=4.9), # Duration from symptomatic to severe symptoms; see Linton et al., https://doi.org/10.3390/jcm9020538, from Table 2, 6.6 day onset to hospital admission (deceased); see also Wang et al., https://jamanetwork.com/journals/jama/fullarticle/2761044, 7 days (Table 1)

        # Duration: Recovery
        'asym2rec': dict(dist='lognormal_int', par1=8.0,  par2=2.0), # Duration for asymptomatic people to recover; see Wölfel et al., https://www.nature.com/articles/s41586-020-2196-x
        'mild2rec': dict(dist='lognormal_int', par1=8.0,  par2=2.0), # Duration for people with mild symptoms to recover; see Wölfel et al., https://www.nature.com/articles/s41586-020-2196-x
        'sev2rec': dict(dist='lognormal_int', par1=18.1, par2=6.3), # Duration for people with severe symptoms to recover, 24.7 days total; see Verity et al., https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30243-7/fulltext; 18.1 days = 24.7 onset-to-recovery - 6.6 sym2sev; 6.3 = 0.35 coefficient of variation * 18.1; see also https://doi.org/10.1017/S0950268820001259 (22 days) and https://doi.org/10.3390/ijerph17207560 (3-10 days)
        'sev2die': dict(dist='lognormal_int', par1=10.7, par2=4.8), # Duration from critical symptoms to death, 18.8 days total; see Verity et al., https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30243-7/fulltext; 10.7 = 18.8 onset-to-death - 6.6 sym2sev - 1.5 sev2crit; 4.8 = 0.45 coefficient of variation * 10.7

        # Severity:

        'rel_symp_prob':1.0,  # Scale factor for proportion of symptomatic cases
        'rel_severe_prob':1.0,  # Scale factor for proportion of symptomatic cases that become severe
        'rel_crit_prob'  :1.0,  # Scale factor for proportion of severe cases that become critical
        'rel_death_prob' :1.0,  # Scale factor for proportion of critical cases that result in death
        'prog_by_age'    :prog_by_age, # Whether to set disease progression based on the person's age NOTE: I'm not sure if we will end up using this or not
        'prognoses'    :None # The actual arrays of prognoses by age; this is populated later
    }

    pars['dur']['flock'] = None
    pars['dur']['barn'] = None
    pars['dur']['water'] = None

    # Efficacy of non-pharmaceutical interventions (NPIs)
    pars['NPIs'] = {}

    pars['NPIs']['human'] = {
        'quar_factor': None, # Quarantine multiplier on transmissibility and susceptibility; set by reset_layer_pars() below
        'quar_period': 14 # Number of days to quarantine for; assumption based on standard policies
    }

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
        if sp in pars.keys():
            pars['variant_pars']['wild'][sp] = pars[sp]

    # Update with any supplied parameter values and generate things that need to be generated
    pars.update(kwargs)
    # reset_layer_pars(pars)
    #if set_prognoses: # If not set here, gets set when the population is initialized
        # pars['prognoses'] = get_prognoses(pars['prog_by_age'], version=version) # Default to age-specific prognoses NOTE: this will have to change to accommodate multiple agent types

    #     # If version is specified, load old parameters
    # if version is not None:
    #     version_pars = znm.get_version_pars(version, verbose=pars['verbose']) 
    #     if sc.compareversions(version, '<3.0.0'): # Waning was introduced in 3.0, so is always false before
    #         version_pars['use_waning'] = False
    #     for key in pars.keys(): # Only loop over keys that have been populated
    #         if key in version_pars: # Only replace keys that exist in the old version
    #             pars[key] = version_pars[key]
    return pars

def get_prognoses(agent_type, version=None):
    '''
    Return the default parameter values for prognoses

    The prognoses probabilities are conditional given the previous disease state.

    Args:
        agent_type (string): type of agent for which we need to retrieve a prognoses

    Returns:
        prog_pars (dict): the dictionary of prognoses probabilities
    '''
    
    prognoses = {}    

    prognoses['human'] = dict(
        age_cutoffs   = np.array([0,       10,      20,      30,      40,      50,      60,      70,      80,      90,]),     # Age cutoffs (lower limits)
        sus_ORs       = np.array([0.34,    0.67,    1.00,    1.00,    1.00,    1.00,    1.24,    1.47,    1.47,    1.47]),    # Odds ratios for relative susceptibility -- from Zhang et al., https://science.sciencemag.org/content/early/2020/05/04/science.abb8001; 10-20 and 60-70 bins are the average across the ORs
        trans_ORs     = np.array([1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00]),    # Odds ratios for relative transmissibility -- no evidence of differences
        comorbidities = np.array([1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00]),    # Comorbidities by age -- set to 1 by default since already included in disease progression rates
        symp_probs    = np.array([0.50,    0.55,    0.60,    0.65,    0.70,    0.75,    0.80,    0.85,    0.90,    0.90]),    # Overall probability of developing symptoms (based on https://www.medrxiv.org/content/10.1101/2020.03.24.20043018v1.full.pdf, scaled for overall symptomaticity)
        severe_probs  = np.array([0.00050, 0.00165, 0.00720, 0.02080, 0.03430, 0.07650, 0.13280, 0.20655, 0.24570, 0.24570]), # Overall probability of developing severe symptoms (derived from Table 1 of https://www.imperial.ac.uk/media/imperial-college/medicine/mrc-gida/2020-03-16-COVID19-Report-9.pdf)
        death_probs   = np.array([0.00002, 0.00002, 0.00010, 0.00032, 0.00098, 0.00265, 0.00766, 0.02439, 0.08292, 0.16190]), # Overall probability of dying -- from O'Driscoll et al., https://www.nature.com/articles/s41586-020-2918-0; last data point from Brazeau et al., https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/covid-19/report-34-ifr/
            )
    prognoses['human'] = relative_prognoses(prognoses['human']) # Convert to conditional probabilities

    prognoses['poultry'] = dict( #TODO: develop poultry prognosis
            
        )

    prognoses['barn'] = dict( #TODO: develop barn prognosis

        )

    prognoses['water'] = dict( #TODO: develop water prognosis

        )
    return prognoses


    # # If version is specified, load old parameters
    # if by_age and version is not None:
    #     version_prognoses = znm.get_version_pars(version, verbose=False)['prognoses']
    #     for key in version_prognoses.keys(): # Only loop over keys that have been populated
    #         if key in version_prognoses: # Only replace keys that exist in the old version
    #             prognoses[key] = np.array(version_prognoses[key])

    # Check that lengths match
    # expected_len = len(prognoses['age_cutoffs'])
    # for key,val in prognoses.items():
    #     this_len = len(prognoses[key])
    #     if this_len != expected_len: # pragma: no cover
    #         errormsg = f'Lengths mismatch in prognoses: {expected_len} age bins specified, but key "{key}" has {this_len} entries'
    #         raise ValueError(errormsg)

    return prognoses

def relative_prognoses(prognoses):
    '''
    Convenience function to revert absolute prognoses into relative (conditional)
    ones. Internally, Covasim uses relative prognoses.
    '''
    out = sc.dcp(prognoses)
    out['death_probs']  /= out['severe_probs']   # Conditional probability of dying, given critical symptoms
    out['severe_probs'] /= out['symp_probs']   # Conditional probability of symptoms becoming severe, given symptomatic
    return out


def absolute_prognoses(prognoses):
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
            rel_beta        = 1.0, # Default values
            rel_symp_prob   = 1.0, # Default values
            rel_severe_prob = 1.0, # Default values
            rel_death_prob  = 1.0, # Default values
        ),

        LPAI = dict(
            rel_beta        = 1.0, # guessed values
            rel_symp_prob   = 0.1, # guess but LPAI is less severe than HPAI
            rel_severe_prob = 0.25, # guess but LPAI is less severe than HPAI
            rel_death_prob  = 0.25, # guess but LPAI is less severe than HPAI
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