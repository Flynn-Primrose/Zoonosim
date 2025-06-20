
import numpy as np
import sciris as sc
from collections import defaultdict
from .. import version as znv
from .. import utils as znu
from .. import defaults as znd
from .Subroster import Subroster
#from .. import base as znb

__all__ = ['Humans']

class HumanMeta(sc.prettyobj):
    ''' Defines all keys used for human agents '''

    def __init__(self):

        # Set the properties of a person
        self.agent = [
            'uid',              # Int
            'age',              # Float
            'sex',              # Float
            'symp_prob',        # Float
            'severe_prob',      # Float
            'death_prob',       # Float
            'rel_trans',        # Float
            'rel_sus',          # Float
            'viral_load',       # Float
            'rescaled_vl',      # Float
            'n_infections',     # Int
            'n_breakthroughs',  # Int
            'cons_days_in_quar',    # Int
            'cons_days_neg_rat',    # Int
        ]

        # Set the states that a person can be in: these are all booleans per person -- used in people.py
        self.states = [
            'susceptible',
            'exposed',
            'infectious',
            'symptomatic',
            'severe',
            'tested',
            'diagnosed',
            'recovered',
            'known_dead',
            'dead',
            'quarantined',
            'vaccinated',
        ]

        # Variant states -- these are ints
        self.variant_states = [
            'exposed_variant',
            'infectious_variant',
            'recovered_variant',
        ]

        # Variant states -- these are ints, by variant
        self.by_variant_states = [
            'symptomatic_by_variant',
            'exposed_by_variant',
            'infectious_by_variant',
            'infections_by_variant',
            'severe_by_variant'
        ]

        # Immune states, by variant
        self.imm_states = [
            'sus_imm',  # Float, by variant
            'symp_imm', # Float, by variant
            'sev_imm',  # Float, by variant
        ]

        # Neutralizing antibody states
        self.nab_states = [
            'peak_nab',    # Float, peak neutralization titre relative to convalescent plasma
            'nab',         # Float, current neutralization titre relative to convalescent plasma
            't_nab_event', # Int, time since nab-conferring event
        ]

        # Additional vaccination states
        self.vacc_states = [
            'doses',          # Number of doses given per person
            'vaccine_source', # index of vaccine that individual received
        ]

        # Set the dates various events took place: these are floats per agent
        self.dates = [f'date_{state}' for state in self.states] # Convert each state into a date
        self.dates.append('date_pos_test') # Store the date when a person tested which will come back positive
        self.dates.append('date_end_quarantine') # Store the date when a person comes out of quarantine
        # Duration of different states: these are floats per human.
        self.durs = [
            'dur_exp2inf',
            'dur_inf2sym',
            'dur_sym2sev',
            'dur_sev2crit',
            'dur_disease',
        ]

        # Timing of control points for viral load
        self.ctrl_points = [
            'x_p_inf',
            'y_p_inf',
            'x_p1',
            'y_p1',
            'x_p2',
            'y_p2',
            'x_p3',
            'y_p3',
        ]
        self.all_recordable_states = self.agent + self.states + self.variant_states + self.nab_states + self.vacc_states + self.dates + self.durs + self.ctrl_points
        self.all_states = self.agent + self.states + self.variant_states + self.by_variant_states + self.imm_states + self.nab_states + self.vacc_states + self.dates + self.durs + self.ctrl_points

        # Validate
        self.state_types = ['agent', 'states', 'variant_states', 'by_variant_states', 'imm_states',
                            'nab_states', 'vacc_states', 'dates', 'durs', 'ctrl_points', 'all_states']
        for state_type in self.state_types:
            states = getattr(self, state_type)
            n_states        = len(states)
            n_unique_states = len(set(states))
            if n_states != n_unique_states: # pragma: no cover
                errormsg = f'In {state_type}, only {n_unique_states} of {n_states} state names are unique'
                raise ValueError(errormsg)

        return
    

class Humans(Subroster):
    '''
    A class to perform all the operations on the human agents -- usually not invoked directly.

    Note that this class handles the mechanics of updating the actual humans, while
    ``cv.BaseRoster`` takes care of housekeeping (saving, loading, exporting, etc.).
    Please see the BaseRoster class for additional methods.

    Args:
        pars (dict): the sim parameters, e.g. sim.pars
        strict (bool): whether or not to only create keys that are already in self.meta.agent; otherwise, let any key be set
        kwargs (dict): the actual data, e.g. from a popdict, being specified

    **Examples**::
    '''

    def __init__(self, pars, strict=True, **kwargs):

        # Handle pars and population size
        self.set_pars(pars)
        self.version = znv.__version__ # Store version info

        # Other initialization
        self.t = 0 # Keep current simulation time
        self._lock = False # Prevent further modification of keys
        self.meta = HumanMeta() # Store list of keys and dtypes
        self.infection_log = [] # Record of infections - keys for ['source','target','date','layer']
        
        pop_size = self.pars['pop_size_by_type']['human']

        # Set person properties -- all floats except for UID
        for key in self.meta.agent:
            if key == 'uid':
                self[key] = np.zeros(pop_size, dtype=znd.default_int) 
            elif key in ['n_infections', 'n_breakthroughs']:
                self[key] = np.zeros(pop_size, dtype=znd.default_int)
            elif key in ['viral_load']:
                self[key] = np.zeros(pop_size, dtype=znd.default_float)
            elif key in ['rescaled_vl']:  # for tracking purposes
                self[key] = np.zeros(pop_size, dtype=znd.default_float)
            else:
                self[key] = np.full(pop_size, np.nan, dtype=znd.default_float)

        # Set health states -- only susceptible is true by default -- booleans except exposed by variant which should return the variant that ind is exposed to
        for key in self.meta.states:
            val = (key in ['susceptible']) # Default value is True for susceptible, False otherwise
            self[key] = np.full(pop_size, val, dtype=bool)

        # Set variant states, which store info about which variant a person is exposed to
        for key in self.meta.variant_states:
            self[key] = np.full(pop_size, np.nan, dtype=znd.default_float)
        for key in self.meta.by_variant_states:
            self[key] = np.full((self.pars['n_variants'], pop_size), False, dtype=bool)

        # Set immunity and antibody states
        for key in self.meta.imm_states:  # Everyone starts out with no immunity
            self[key] = np.zeros((self.pars['n_variants'], pop_size), dtype=znd.default_float)
        for key in self.meta.nab_states:  # Everyone starts out with no antibodies
            dtype = znd.default_int if key == 't_nab_event' else znd.default_float
            self[key] = np.zeros(pop_size, dtype=dtype)
        for key in self.meta.vacc_states:
            self[key] = np.zeros(pop_size, dtype=znd.default_int)

        # Set dates and durations -- both floats
        for key in self.meta.dates + self.meta.durs:
            self[key] = np.full(pop_size, np.nan, dtype=znd.default_float)

        # Set dates for viral load profile -- floats
        for key in self.meta.ctrl_points:
            self[key] = np.full(pop_size, np.nan, dtype=znd.default_float)

        # Store the dtypes used in a flat dict
        self._dtypes = {key:self[key].dtype for key in self.keys()} # Assign all to float by default
        if strict:
            self.lock() # If strict is true, stop further keys from being set (does not affect attributes)

        # Store flows to be computed during simulation
        self.init_flows()

        # Although we have called init(), we still need to call initialize()
        self.initialized = False


        # Handle all other values, e.g. age
        for key,value in kwargs.items():
            if strict:
                self.set(key, value)
            else:
                self[key] = value

        self._pending_quarantine = defaultdict(list)  # Internal cache to record people that need to be quarantined on each timestep {t:(person, quarantine_end_day)}

        return


    def init_flows(self):
        ''' Initialize flows to be zero '''
        self.flows = {key:0 for key in znd.new_human_flows}
        self.flows_variant = {}
        for key in znd.new_human_flows_by_variant:
            self.flows_variant[key] = np.zeros(self.pars['n_variants'], dtype=znd.default_float)
        return

    def initialize(self, agents_pars=None):
        ''' Perform initializations '''
        self.validate(roster_pars=agents_pars) # First, check that essential-to-match parameters match
        self.set_pars(agents_pars) # Replace the saved parameters with this simulation's
        self.set_prognoses()
        self.initialized = True
        return


    def set_prognoses(self):
        '''
        Set the prognoses for each person based on age during initialization. Need
        to reset the seed because viral loads are drawn stochastically.
        '''

        pars = self.pars # Shorten
        if 'prognoses' not in pars or 'rand_seed' not in pars:
            errormsg = 'This Humans object does not have the required parameters ("prognoses" and "rand_seed"). Create a sim (or parameters), then do e.g. people.set_pars(sim.pars).'
            raise sc.KeyNotFoundError(errormsg)

        def find_cutoff(age_cutoffs, age):
            '''
            Find which age bin each person belongs to -- e.g. with standard
            age bins 0, 10, 20, etc., ages [5, 12, 4, 58] would be mapped to
            indices [0, 1, 0, 5]. Age bins are not guaranteed to be uniform
            width, which is why this can't be done as an array operation.
            '''
            return np.nonzero(age_cutoffs <= age)[0][-1]  # Index of the age bin to use

        znu.set_seed(pars['rand_seed'])

        progs = pars['prognoses']['human']
        inds = np.fromiter((find_cutoff(progs['age_cutoffs'], this_age) for this_age in self.age), dtype=znd.default_int, count=len(self)) # Convert ages to indices
        self.symp_prob[:]   = progs['symp_probs'][inds] # Probability of developing symptoms
        self.severe_prob[:] = progs['severe_probs'][inds]*progs['comorbidities'][inds] # Severe disease probability is modified by comorbidities
        self.death_prob[:]  = progs['death_probs'][inds] # Probability of death
        self.rel_sus[:]     = progs['sus_ORs'][inds]  # Default susceptibilities
        self.rel_trans[:]   = progs['trans_ORs'][inds] * znu.sample(**self.pars['transmission_pars']['human']['beta_dist'], size=len(inds))  # Default transmissibilities, with viral load drawn from a distribution

        return




    def update_states_pre(self, t):
        ''' Perform all state updates at the current timestep '''

        # Initialize
        self.t = t
        self.is_exp = self.true('exposed') # For storing the interim values since used in every subsequent calculation

        # Perform updates
        self.init_flows()
        self.flows['new_infectious']    += self.check_infectious() # For people who are exposed and not infectious, check if they begin being infectious
        self.flows['new_symptomatic']   += self.check_symptomatic()
        self.flows['new_severe']        += self.check_severe()
        self.flows['new_recovered']    += self.check_recovery()
        new_dead, new_known_dead     = self.check_death()
        self.flows['new_dead']        += new_dead
        self.flows['new_known_dead']  += new_known_dead
        return


    def update_states_post(self):
        ''' Perform post-timestep updates '''


        self.flows['new_diagnosed'] += self.check_diagnosed()
        #self.flows['new_quarantined'] += self.check_quar()

        del self.is_exp  # Tidy up

        return



    
    def schedule_behaviour(self):
        ''' Schedules events on the basis of results received today '''


        return

    
    def get_key_to_use(self, policy_dict):
        ''' Helper function used to determine what policy to use when scheduling events '''

        keys = np.array(list(policy_dict.keys()))
        key_to_use = keys[znu.true(keys <= self.t)[-1]]

        return key_to_use


    #%% Methods for updating state



    def check_inds_diagnosed(self, current, date, filter_inds=None):
        '''
        Modified version of check_inds specifically for diagnoses.
        Replacing >= with == ensures that individuals who leave the diagnosed state do not immediately become diagnosed again.
        '''
        if filter_inds is None:
            not_current = znu.false(current)
        else:
            not_current = znu.ifalsei(current, filter_inds)
        has_date = znu.idefinedi(date, not_current)
        inds     = znu.itrue(self.t == date[has_date], has_date)
        return inds
    
    def check_diagnosed(self):
        '''
        Check for new diagnoses. Since most data are reported with diagnoses on
        the date of the test, this function reports counts not for the number of
        people who received a positive test result on a day, but rather, the number
        of people who were tested on that day who are schedule to be diagnosed in
        the future.
        '''

        # Handle people who tested today who will be diagnosed in future (i.e., configure them to have have finished being tested)
        test_pos_inds = self.check_inds(self.diagnosed, self.date_pos_test, filter_inds=None) # Find people who are not diagnosed and have a date of a positive test that is today or earlier
        self.date_pos_test[test_pos_inds] = np.nan # Clear date of having will-be-positive test


        # Handle people who were actually diagnosed today (i.e., set them as diagnosed and remove any of them that were quarantining from quarantine)
        # diag_inds  = self.check_inds(self.diagnosed, self.date_diagnosed, filter_inds=None) # Find who are not diagnosed and have a date of diagnosis that is today or earlier
        diag_inds  = self.check_inds_diagnosed(self.diagnosed, self.date_diagnosed, filter_inds=None) # Find who are not diagnosed and have a date of diagnosis that is today or earlier

        self.diagnosed[diag_inds]   = True # Set these people to be diagnosed

        quarantined = znu.itruei(self.quarantined, diag_inds) # Find individuals who were just diagnosed who are in quarantine
        self.date_end_quarantine[quarantined] = self.t # Set end quarantine date to match when the person left quarantine (and entered isolation)
        self.quarantined[diag_inds] = False # If you are diagnosed, you are isolated, not in quarantine

        # Remove diagnosed individuals if they were diagnosed more than 14 days ago
        diag_expired_inds = znu.true(self.t - self.date_diagnosed > 14)
        self.diag_expired_inds = diag_expired_inds
        self.diagnosed[diag_expired_inds] = False

        return len(test_pos_inds)

    def check_quar(self):
        ''' Update quarantine state '''

        n_quarantined = 0 # Number of people entering quarantine
        for ind,end_day in self._pending_quarantine[self.t]:
            if self.quarantined[ind]: # Update when quarantine should be finished (in case schedule_quarantine is called on someone already in quarantine)
                self.date_end_quarantine[ind] = max(self.date_end_quarantine[ind], end_day) # Extend quarantine if required
            elif not (self.dead[ind] or self.recovered[ind] or self.diagnosed[ind]): # Unclear whether recovered should be included here # elif not (self.dead[ind] or self.diagnosed[ind]):
                self.quarantined[ind] = True
                self.date_quarantined[ind] = self.t
                self.date_end_quarantine[ind] = end_day
                n_quarantined += 1

        # If someone on quarantine has reached the end of their quarantine, release them
        end_inds = self.check_inds(~self.quarantined, self.date_end_quarantine, filter_inds=None) # Note the double-negative here (~)
        self.quarantined[end_inds] = False # Release from quarantine

        # Update the counter for consecutive days in quarantine
        self.cons_days_in_quar[self.quarantined] += 1
        self.cons_days_in_quar[~self.quarantined] = 0

        return n_quarantined


    def check_infectious(self):
        ''' Check if they become infectious '''
        inds = self.check_inds(self.infectious, self.date_infectious, filter_inds=self.is_exp)
        self.infectious[inds] = True
        self.infectious_variant[inds] = self.exposed_variant[inds]
        

        for variant in range(self.pars['n_variants']):
            this_variant_inds = znu.itrue(self.infectious_variant[inds] == variant, inds)
            n_this_variant_inds = len(this_variant_inds)

            self.flows_variant['new_infectious_by_variant'][variant] += n_this_variant_inds
            self.infectious_by_variant[variant, this_variant_inds] = True
        return len(inds)


    def check_symptomatic(self):
        ''' Check for new progressions to symptomatic '''
        inds = self.check_inds(self.symptomatic, self.date_symptomatic, filter_inds=self.is_exp)
        self.symptomatic[inds] = True
        return len(inds)


    def check_severe(self):
        ''' Check for new progressions to severe '''
        inds = self.check_inds(self.severe, self.date_severe, filter_inds=self.is_exp)
        
        self.severe[inds] = True
        return len(inds)



    def check_recovery(self, inds=None, filter_inds='is_exp'):
        '''
        Check for recovery.

        More complex than other functions to allow for recovery to be manually imposed
        for a specified set of indices.
        '''

        # Handle more flexible options for setting indices
        if filter_inds == 'is_exp':
            filter_inds = self.is_exp
        if inds is None:
            inds = self.check_inds(self.recovered, self.date_recovered, filter_inds=filter_inds)

        # Now reset all disease states
        self.exposed[inds]          = False
        self.infectious[inds]       = False
        self.symptomatic[inds]      = False
        self.severe[inds]           = False
        self.recovered[inds]        = True
        self.recovered_variant[inds] = self.exposed_variant[inds]
        self.infectious_variant[inds] = np.nan
        self.exposed_variant[inds]    = np.nan
        self.exposed_by_variant[:, inds] = False
        self.infectious_by_variant[:, inds] = False

        # Handle immunity aspects
        if self.pars['immunity_pars']['human']['use_waning']:

            # Reset additional states
            self.susceptible[inds] = True
            self.diagnosed[inds]   = False # Reset their diagnosis state because they might be reinfected

        return len(inds)


    def check_death(self):
        ''' Check whether or not this agent died on this timestep  '''
        inds = self.check_inds(self.dead, self.date_dead, filter_inds=self.is_exp)  
        
        self.dead[inds]             = True
        diag_inds = inds[self.diagnosed[inds]] # Check whether the person was diagnosed before dying
        self.known_dead[diag_inds]  = True
        self.susceptible[inds]      = False
        self.exposed[inds]          = False
        self.infectious[inds]       = False
        self.symptomatic[inds]      = False
        self.severe[inds]           = False
        #self.known_contact[inds]    = False
        self.quarantined[inds]      = False
        self.recovered[inds]        = False
        self.infectious_variant[inds] = np.nan
        self.exposed_variant[inds]    = np.nan
        self.recovered_variant[inds]  = np.nan

        return len(inds), len(diag_inds)



    #%% Methods to make events occur (infection and diagnosis)




    def infect(self, inds, hosp_max=None, source=None, layer=None, variant=0):
        '''
        Infect people and determine their eventual outcomes.

            * Every infected person can infect other people, regardless of whether they develop symptoms
            * Infected people that develop symptoms are disaggregated into mild vs. severe (=requires hospitalization) vs. critical (=requires ICU)
            * Every asymptomatic, mildly symptomatic, and severely symptomatic person recovers
            * Critical cases either recover or die
            * If the simulation is being run with waning, this method also sets/updates agents' neutralizing antibody levels

        Method also deduplicates input arrays in case one agent is infected many times
        and stores who infected whom in infection_log list.

        Args:
            inds     (array): array of people to infect
            hosp_max (bool): are all the hospital beds occupied
            source   (array): source indices of the people who transmitted this infection (None if an importation or seed infection)
            layer    (str):   contact layer this infection was transmitted on
            variant  (int):   the variant people are being infected by

        Returns:
            count (int): number of people infected
        '''

        if len(inds) == 0:
            return 0

        # Remove duplicates
        inds, unique = np.unique(inds, return_index=True)
        if source is not None:
            source = source[unique]

        # Keep only susceptibles
        keep = self.susceptible[inds] # Unique indices in inds and source that are also susceptible
        inds = inds[keep]
        if source is not None:
            source = source[keep]

        # Deal with variant parameters
        variant_keys = ['rel_symp_prob', 'rel_severe_prob', 'rel_death_prob']
        infect_pars = {k:self.pars[k]['human'] for k in variant_keys}
        variant_label = self.pars['variant_map'][variant]
        if variant:
            for k in variant_keys:
                infect_pars[k] *= self.pars['variant_pars'][variant_label]['human'][k]

        n_infections = len(inds)
        durpars      = self.pars['dur']['human']

        # Retrieve those with a breakthrough infection (defined nabs)
        breakthrough_inds = inds[znu.true(self.peak_nab[inds])]
        if len(breakthrough_inds):
            no_prior_breakthrough = (self.n_breakthroughs[breakthrough_inds] == 0) # We only adjust transmissibility for the first breakthrough
            new_breakthrough_inds = breakthrough_inds[no_prior_breakthrough]
            self.rel_trans[new_breakthrough_inds] *= self.pars['immunity_pars']['human']['trans_redux']

        # Update states, variant info, and flows
        self.susceptible[inds]    = False
        self.recovered[inds]      = False
        self.diagnosed[inds]      = False
        self.exposed[inds]        = True
        self.n_infections[inds]  += 1
        self.n_breakthroughs[breakthrough_inds] += 1
        self.exposed_variant[inds] = variant
        self.exposed_by_variant[variant, inds] = True
        self.flows['new_exposed']   += len(inds)
        self.flows['new_infections'] += len(inds) # Record new infections
        self.flows['new_reinfections'] += len(znu.defined(self.date_recovered[inds])) # Record reinfections
        self.flows_variant['new_infections_by_variant'][variant] += len(inds)

        # Record transmissions
        for i, target in enumerate(inds):
            entry = dict(source=source[i] if source is not None else None, target=target, date=self.t, layer=layer, variant=variant_label)
            self.infection_log.append(entry)

        # Calculate how long before this person can infect other people
        self.dur_exp2inf[inds] = znu.sample(**durpars['exp2inf'], size=n_infections)
        self.date_exposed[inds] = self.t
        self.date_infectious[inds] = self.dur_exp2inf[inds] + self.t

        # Reset all other dates
        for key in ['date_symptomatic', 'date_severe', 'date_diagnosed', 'date_recovered']:
            self[key][inds] = np.nan

        # Use prognosis probabilities to determine what happens to them
        symp_probs = infect_pars['rel_symp_prob']*self.symp_prob[inds]*(1-self.symp_imm[variant, inds]) # Calculate their actual probability of being symptomatic
        is_symp = znu.binomial_arr(symp_probs) # Determine if they develop symptoms
        symp_inds = inds[is_symp]
        asymp_inds = inds[~is_symp] # Asymptomatic
        self.flows_variant['new_symptomatic_by_variant'][variant] += len(symp_inds)

        # CASE 1: Asymptomatic: may infect others, but have no symptoms and do not die
        dur_asym2rec = znu.sample(**durpars['asym2rec'], size=len(asymp_inds))
        self.date_recovered[asymp_inds] = self.date_infectious[asymp_inds] + dur_asym2rec  # Date they recover. Nx1.
        self.dur_disease[asymp_inds] = self.dur_exp2inf[asymp_inds] + dur_asym2rec  # Store how long this person had COVID-19. Nx1. 

        # CASE 2: Symptomatic: can either be mild or severe.
        n_symp_inds = len(symp_inds)
        self.dur_inf2sym[symp_inds] = znu.sample(**durpars['inf2sym'], size=n_symp_inds) # Store how long this person took to develop symptoms
        self.date_symptomatic[symp_inds] = self.date_infectious[symp_inds] + self.dur_inf2sym[symp_inds] # Date they become symptomatic
        sev_probs = infect_pars['rel_severe_prob'] * self.severe_prob[symp_inds]*(1-self.sev_imm[variant, symp_inds]) # Probability of these people being severe
        is_sev = znu.binomial_arr(sev_probs) # See if they're a severe or mild case
        sev_inds = symp_inds[is_sev]
        mild_inds = symp_inds[~is_sev] # Not severe
        self.flows_variant['new_severe_by_variant'][variant] += len(sev_inds)

        # CASE 2.1: Mild symptoms, no hospitalization required and no probability of death
        dur_mild2rec = znu.sample(**durpars['mild2rec'], size=len(mild_inds))
        self.date_recovered[mild_inds] = self.date_symptomatic[mild_inds] + dur_mild2rec  # Date they recover
        self.dur_disease[mild_inds] = self.dur_exp2inf[mild_inds] + self.dur_inf2sym[mild_inds] + dur_mild2rec  # Store how long this person had COVID-19

        # CASE 2.2: Severe cases: hospitalization required, may die
        self.dur_sym2sev[sev_inds] = znu.sample(**durpars['sym2sev'], size=len(sev_inds)) # Store how long this person took to develop severe symptoms
        self.date_severe[sev_inds] = self.date_symptomatic[sev_inds] + self.dur_sym2sev[sev_inds]  # Date symptoms become severe
        death_probs = infect_pars['rel_death_prob'] * self.death_prob[sev_inds] * (self.pars['no_hosp_factor'] if hosp_max else 1.)# Probability they'll die
        is_dead = znu.binomial_arr(death_probs)  # Death outcome
        dead_inds = sev_inds[is_dead]
        alive_inds = sev_inds[~is_dead]
        

        # CASE 2.2.1: Did not die
        dur_sev2rec = znu.sample(**durpars['sev2rec'], size=len(alive_inds))
        self.date_recovered[alive_inds] = self.date_severe[alive_inds] + dur_sev2rec # Date they recover
        self.dur_disease[alive_inds] = self.dur_exp2inf[alive_inds] + self.dur_inf2sym[alive_inds] + self.dur_sym2sev[alive_inds] + self.dur_sev2crit[alive_inds] + dur_sev2rec  # Store how long this person had COVID-19

        # CASE 2.2.2: Did die
        dur_sev2die = znu.sample(**durpars['sev2die'], size=len(dead_inds))
        self.date_dead[dead_inds] = self.date_severe[dead_inds] + dur_sev2die # Date of death
        self.dur_disease[dead_inds] = self.dur_exp2inf[dead_inds] + self.dur_inf2sym[dead_inds] + self.dur_sym2sev[dead_inds] + self.dur_sev2crit[dead_inds] + dur_sev2die   # Store how long this person had COVID-19
        self.date_recovered[dead_inds] = np.nan # If they did die, remove them from recovered

        # HANDLE VIRAL LOAD CONTROL POINTS
        
        # Get P_inf: where viral load crosses 10^6 cp/mL
        self.x_p_inf[inds] = self.dur_exp2inf[inds]
        self.y_p_inf[inds] = 6

        # Get P1: where viral load crosses 10^3 cp/mL;  dummy value for now
        self.x_p1[inds] = np.maximum(self.x_p_inf[inds] - (np.random.gamma(2, 0.35, size=len(inds)) + 0.25), 0)
        self.y_p1[inds] = 3

        # Get P2: where viral load peaks; dummy value for now
        self.x_p2[inds] = self.x_p_inf[inds] + (np.random.gamma(3, 0.26, size=len(inds)) + 0.1)
        self.y_p2[inds] = ((self.y_p_inf[inds] - self.y_p1[inds])*(self.x_p2[inds] - self.x_p1[inds])/(self.x_p_inf[inds] - self.x_p1[inds])) + self.y_p1[inds]

        # Align P1, P_inf, and P2 to current time
        self.x_p1[inds] = self.x_p1[inds] + self.t
        self.x_p_inf[inds] = self.x_p_inf[inds] + self.t
        self.x_p2[inds] = self.x_p2[inds] + self.t

        # Get P3: where viral load drops below 10^6 cp/mL
        time_recovered = np.ones(len(self.date_recovered), dtype=znd.default_float)*self.date_recovered # This is needed to make a copy
        inds_dead = ~np.isnan(self.date_dead)
        time_recovered[inds_dead] = self.date_dead[inds_dead]
        self.x_p3[inds] = np.maximum(time_recovered[inds], self.x_p2[inds])
        self.y_p3[inds] = 6

            # # For testing purposes
            # if self.t < self.pars['x_p1'].shape[1]:
            #     self.pars['x_p1'][:, self.t] = self.x_p1
            #     self.pars['x_p_inf'][:, self.t] = self.x_p_inf
            #     self.pars['x_p2'][:, self.t] = self.x_p2
            #     self.pars['x_p3'][:, self.t] = self.x_p3
            #     self.pars['y_p1'][:, self.t] = self.y_p1
            #     self.pars['y_p_inf'][:, self.t] = self.y_p_inf
            #     self.pars['y_p2'][:, self.t] = self.y_p2
            #     self.pars['y_p3'][:, self.t] = self.y_p3

        # Handle immunity aspects
        # if self.pars['use_waning']:
        #    symp = dict(asymp=asymp_inds, mild=mild_inds, sev=sev_inds)
        #    cvi.update_peak_nab(self, inds, nab_pars=self.pars, symp=symp)

        return n_infections # For incrementing counters


    def test(self, inds, test_sensitivity=1.0, loss_prob=0.0, test_delay=0):
        '''
        Method to test people. Typically not to be called by the user directly;
        see the test_num() and test_prob() interventions.

        Args:
            inds: indices of who to test
            test_sensitivity (float): probability of a true positive
            loss_prob (float): probability of loss to follow-up
            test_delay (int): number of days before test results are ready
        '''

        inds = np.unique(inds)
        self.tested[inds] = True
        self.date_tested[inds] = self.t # Only keep the last time they tested

        is_infectious = znu.itruei(self.infectious, inds)
        pos_test      = znu.n_binomial(test_sensitivity, len(is_infectious))
        is_inf_pos    = is_infectious[pos_test]

        not_diagnosed = is_inf_pos[np.isnan(self.date_diagnosed[is_inf_pos])]
        not_lost      = znu.n_binomial(1.0-loss_prob, len(not_diagnosed))
        final_inds    = not_diagnosed[not_lost]

        # Store the date the person will be diagnosed, as well as the date they took the test which will come back positive
        self.date_diagnosed[final_inds] = self.t + test_delay
        self.date_pos_test[final_inds] = self.t

        return final_inds

