import numpy as np
import sciris as sc
from .. import version as znv
from .. import utils as znu
from .. import defaults as znd
from .Subroster import Subroster

__all__ = ['Flocks']

class FlocksMeta(sc.prettyobj):
    ''' Defines all keys that are used by flocks '''

    def __init__(self):
        
        self.agent = [
            'uid', # int
            'breed', # e.g. breeder, layer, broiler
            'barn', # uid of the barn where the flock is located
            'headcount', # Number of live birds in the flock
            'rel_sus', # Relative susceptibility
            'rel_trans', # Relative Transmissibility
            'infected_headcount',
            'symptomatic_headcount',
            'baseline_symptomatic_rate', 
            'infected_symptomatic_rate', 
            'dead_headcount',
            'baseline_mortality_rate',
            'infected_mortality_rate', 
            'water_consumption', # Water consumption (L/day)
            'baseline_water_rate', # Water consumption rate (L/bird/day)
            'infected_water_rate', # Water consumption rate (L/bird/day)
        ]

        self.states = [
            'susceptible',
            'exposed',
            'infectious',
            'suspected', # Producer suspects flock is infected
            'quarantined', # Quarantined by CFIA agent
        ]

        # NOTE: I'm not sure if we need to record the cfia risk level assessment directly or if it is implied by the other states
        # self.risk_level = [
        #     'negligible', # CFIA agent assesses flock as negligible risk
        #     'tbc_negative', # CFIA agent assesses flock as negative but to be confirmed
        #     'tbc_positive',# CFIA agent assesses flock as positive but to be confirmed
        # ]

        self.variant_states = [
            'exposed_variant',
            'infectious_variant',
            'symptomatic_variant',
        ]

        self.by_variant_states = [
            'exposed_by_variant',
            'infectious_by_variant',
            'symptomatic_by_variant',
        ]

        # Set the dates various events took place: these are floats per agent
        self.state_dates = [f'date_{state}' for state in self.states] # Convert each state into a date

        # Dates for the cfia protocols: these are floats per flock
        self.protocol_dates = [
            'date_result', # Date of lab test result
        ]

        self.dates = self.state_dates + self.protocol_dates

        # Duration of different states: these are floats per flock.
        self.durs = [
            'dur_exp2inf', # Mean time from exposed to infectious
            'dur_inf2symp', # Mean time from infectious to symptomatic
            'dur_inf2peak', # Mean time from infectious to peak infection
            'dur_peak2eq', # Mean time from peak infection to equilibrium infection
            'dur_susp2res', # Mean time from inspected to confirmation of status 
        ]

        # Control points for infection progression
        self.ctrl_points = [
            'x_p1',
            'y_p1',
            'x_p2',
            'y_p2',
            'x_p3',
            'y_p3',
        ]

        self.all_states = self.agent + self.states + self.variant_states + self.by_variant_states + self.dates + self.durs

        # Validate
        self.state_types = ['agent', 'states', 'variant_states', 'by_variant_states']
        for state_type in self.state_types:
            states = getattr(self, state_type)
            n_states        = len(states)
            n_unique_states = len(set(states))
            if n_states != n_unique_states: # pragma: no cover
                errormsg = f'In {state_type}, only {n_unique_states} of {n_states} state names are unique'
                raise ValueError(errormsg)

        return
    

class Flocks(Subroster):
    '''
    A class to perform all the operations on the flock agents -- usually not invoked directly.

    Note that this class handles the mechanics of updating the actual poultry, while
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
        self.meta = FlocksMeta() # Store list of keys and dtypes
        self.contacts = None
        # self.init_contacts() # Initialize the contacts
        self.infection_log = [] # Record of infections - keys for ['source','target','date','layer']
        
        pop_size = pars['pop_size_by_type']['flock']

        # Set person properties -- all floats except for UID
        for key in self.meta.agent:
            if key == 'uid':
                self[key] = np.zeros(pop_size, dtype=znd.default_int)# NOTE: values get passed as kwargs by make_flock
            elif key == 'barn':
                self[key] = np.zeros(pop_size, dtype=znd.default_int)# NOTE: values get passed as kwargs by make_flock
            elif key == 'breed':
                self[key] = np.zeros(pop_size, dtype=znd.default_str) # NOTE: values get passed as kwargs by make_flock
            else:
                self[key] = np.full(pop_size, np.nan, dtype=znd.default_float)

        # Set health states -- only susceptible is true by default -- booleans except exposed by variant which should return the variant that ind is exposed to
        for key in self.meta.states:
            val = (key in ['susceptible']) # Default value is True for susceptible, False otherwise
            self[key] = np.full(pop_size, val, dtype=bool)

        # Set variant states, which store info about which variant an agent is exposed to
        for key in self.meta.variant_states:
            self[key] = np.full(pop_size, np.nan, dtype=znd.default_float)
        for key in self.meta.by_variant_states:
            self[key] = np.full((self.pars['n_variants'], pop_size), False, dtype=bool)

        # Set dates and durations -- both floats
        for key in self.meta.dates + self.meta.durs:
            self[key] = np.full(pop_size, np.nan, dtype=znd.default_float)

        # Set dates for infection profile -- floats
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

        # self._pending_quarantine = defaultdict(list)  # Internal cache to record people that need to be quarantined on each timestep {t:(person, quarantine_end_day)}

        return


    def init_flows(self):
        ''' Initialize flows to be zero '''
        self.flows = {key:0 for key in znd.new_flock_flows}
        self.flows_variant = {}
        for key in znd.new_flock_flows_by_variant:
            self.flows_variant[key] = np.zeros(self.pars['n_variants'], dtype=znd.default_float)
        return

    def initialize(self, sim_pars=None):
        ''' Perform initializations '''
        self.validate(sim_pars=sim_pars) # First, check that essential-to-match parameters match
        self.set_pars(sim_pars) # Replace the saved parameters with this simulation's
        self.set_prognoses()
        self.initialized = True
        return


    def set_prognoses(self):
        '''
        Set the prognoses for each flock based on breed. 
        '''
        pars = self.pars # Shorten
        if 'prognoses' not in pars or 'rand_seed' not in pars:
            errormsg = 'This flock object does not have the required parameters ("prognoses" and "rand_seed").'
            raise sc.KeyNotFoundError(errormsg)


        znu.set_seed(pars['rand_seed'])

        progs = pars['prognoses']['flock']
        breed_to_index = {breed: index for index, breed in enumerate(progs['breeds'])}
        #inds = np.fromiter((breed_to_index[this_breed] for this_breed in self.breed), dtype=znd.default_str)
        inds = np.array([breed_to_index[this_breed] for this_breed in self.breed])
        self.baseline_symptomatic_rate[:] = progs['baseline_symptomatic_rate'][inds]
        self.baseline_mortality_rate[:] = progs['baseline_mortality_rate'][inds]
        self.baseline_water_rate[:] = progs['baseline_water_rate'][inds]
        self.rel_sus[:] = progs['sus_ORs'][inds]
        self.rel_trans[:] = progs['trans_ORs'][inds]

        return




    def update_states_pre(self, t):
        ''' Perform all state updates at the current timestep '''

        # Initialize
        self.t = t


        # Perform updates
        self.init_flows()
        self.flows['new_infectious'] = self.check_infectious() # For flocks that are exposed and not infectious, check if they begin being infectious
        self.update_headcounts() # Update the headcounts and water consumption of the flocks
        return


    def update_states_post(self):
        ''' Perform post-timestep updates '''

        # Update the status of flocks

        return



    #%% Methods for updating state

    def check_inds(self, current, date, filter_inds=None):
        ''' Return indices for which the current state is false and which are assigned a date on or before the current date

        Args:
            current (array): list of boolean values that represent a current state
            date (array): list that contains either a date or a Nan
        '''
        if filter_inds is None:
            not_current = znu.false(current)
        else:
            not_current = znu.ifalsei(current, filter_inds)
        has_date = znu.idefinedi(date, not_current)
        inds     = znu.itrue(self.t >= date[has_date], has_date)
        return inds


    def check_infectious(self):
        ''' Check if they become infectious '''
        inds = self.check_inds(self.infectious, self.date_infectious)
        self.infectious[inds] = True
        self.infectious_variant[inds] = self.exposed_variant[inds]


        for variant in range(self.pars['n_variants']):
            this_variant_inds = znu.itrue(self.infectious_variant[inds] == variant, inds)
            n_this_variant_inds = len(this_variant_inds)
            self.flows_variant['new_infectious_by_variant'][variant] += n_this_variant_inds
            self.infectious_by_variant[variant, this_variant_inds] = True
        return len(inds)
    


    def check_quarantined(self):
        ''' Check for new progressions to quarantined '''
        inds = self.check_inds(self.quarantined, self.date_quarantined)
        self.quarantined[inds] = True
        return len(inds)
    
    
    def update_headcounts(self):
        ''' Update the headcounts and water consumption of the flocks '''
        uninfected_headcount = self.headcount - self.infected_headcount
        dead_infected = self.infected_headcount * self.infected_symptomatic_rate
        dead_uninfected = uninfected_headcount * self.baseline_symptomatic_rate
        self.dead_headcount += dead_infected + dead_uninfected
        self.infected_headcount = self.infected_headcount - dead_infected
        self.headcount -= dead_infected - dead_uninfected
        self.symptomatic_headcount = self.infected_headcount * self.infected_symptomatic_rate + uninfected_headcount * self.baseline_symptomatic_rate
        self.water_consumption = self.infected_headcount * self.infected_water_rate + uninfected_headcount * self.baseline_water_rate
        return

    #%% Methods to make events occur (infection and diagnosis)


    def infect(self, inds, source=None, layer=None, variant=0):
        '''
        Infect agents and determine their eventual outcomes.

            * Every infected agent can infect other agent, regardless of whether they develop symptoms
            * Every asymptomatic, mildly symptomatic, and severely symptomatic person recovers
            * Critical cases either recover or die
            * If the simulation is being run with waning, this method also sets/updates agents' neutralizing antibody levels

        Method also deduplicates input arrays in case one agent is infected many times
        and stores who infected whom in infection_log list.

        Args:
            inds     (array): array of agents to infect
            source   (array): source indices of the agent who transmitted this infection (None if an importation or seed infection)
            layer    (str):   contact layer this infection was transmitted on
            variant  (int):   the variant agents are being infected by

        Returns:
            count (int): number of flocks infected
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
        variant_keys = ['rel_symp_prob', 'rel_death_prob']
        infect_pars = {k:self.pars[k] for k in variant_keys}
        variant_label = self.pars['variant_map'][variant]
        if variant:
            for k in variant_keys:
                infect_pars[k] *= self.pars['variant_pars'][variant_label][k]

        n_infections = len(inds)
        durpars      = self.pars['dur']['flock']

        # Update states, variant info, and flows
        self.susceptible[inds]    = False
        self.exposed[inds]        = True
        self.exposed_variant[inds] = variant
        self.exposed_by_variant[variant, inds] = True
        self.flows['new_infectious']   += len(inds)
        self.flows_variant['new_infectious_by_variant'][variant] += len(inds)

        # Record transmissions
        for i, target in enumerate(inds):
            entry = dict(source=source[i] if source is not None else None, target=target, date=self.t, layer=layer, variant=variant_label)
            self.infection_log.append(entry)

        # Calculate how long before this flock can infect other flocks
        self.dur_exp2inf[inds] = np.maximum(znu.sample(**durpars['exp2inf'], size=n_infections), 0) # Ensure that this is not negative
        self.dur_inf2peak[inds] = np.maximum(znu.sample(**durpars['inf2peak'], size=n_infections), 0) # Ensure that this is not negative
        self.dur_peak2eq[inds] = np.maximum(znu.sample(**durpars['peak2eq'], size=n_infections), 0) # Ensure that this is not negative
        self.date_exposed[inds] = self.t
        self.date_infectious[inds] = self.dur_exp2inf[inds] + self.t

        #Update water rate, symptomatic rate, and mortality rate.
        pars = self.pars # Shorten
        if 'prognoses' not in pars or 'rand_seed' not in pars:
            errormsg = 'This flock object does not have the required parameters ("prognoses" and "rand_seed").'
            raise sc.KeyNotFoundError(errormsg)
        progs = pars['prognoses']['flock']
        breed_to_index = {breed: index for index, breed in enumerate(progs['breeds'])}
        #breed_inds = np.fromiter((breed_to_index[this_breed] for this_breed in self.breed[inds]), dtype=znd.default_str)
        breed_inds = np.array([breed_to_index[this_breed] for this_breed in self.breed[inds]])
        breed, frequency = np.unique(breed_inds, return_counts=True)
        breed_freq = zip(breed, frequency)
        for breed, frequency in breed_freq:
            # NOTE: I'm just guessing at the distribution of these parameters.
            self.infected_symptomatic_rate[inds[breed_inds == breed]] = self.baseline_symptomatic_rate[inds[breed_inds == breed]] + np.maximum(znu.sample('lognormal', progs['mean_symptomatic_rate_increase'][breed], 1.0, size=frequency), 0)*infect_pars['rel_symp_prob']
            self.infected_mortality_rate[inds[breed_inds == breed]] = self.baseline_mortality_rate[inds[breed_inds == breed]] + np.maximum(znu.sample('lognormal', progs['mean_mortality_rate_increase'][breed], 1.0, size=frequency), 0)*infect_pars['rel_death_prob']
            self.infected_water_rate[inds[breed_inds == breed]] = self.baseline_water_rate[inds[breed_inds == breed]] + np.maximum(znu.sample('lognormal', progs['mean_water_rate_increase'][breed], 1.0, size=frequency), 0)



        # HANDLE INFECTION LEVEL CONTROL POINTS

        # Get P1: 
        self.x_p1[inds] = self.date_infectious[inds] # Date of first infections.
        self.y_p1[inds] = 3 # Initial number of infections. TODO: This should be sampled from a distribution, but I don't know what the best distribution is.

        # Get P2: 
        self.x_p2[inds] = self.x_p1[inds] + self.dur_inf2peak[inds] # Date of peak infections.
        self.y_p2[inds] = 0.7*self.headcount[inds] # Peak number of infections. TODO: This should be sampled from a distribution, but I don't know what the best distribution is.

        # Get P3: 

        self.x_p3[inds] = self.x_p2[inds] + self.dur_peak2eq[inds] # Date of equilibrium infections.
        self.y_p3[inds] = 0.1*self.headcount[inds] # Equilibrium number of infections. TODO: This should be sampled from a distribution, but I don't know what the best distribution is.

        return n_infections # For incrementing counters


    def test(self, inds, test_sensitivity=1.0, loss_prob=0.0, test_delay=0):
        '''
        Method to test poultry. Typically not to be called by the user directly;
        see the test_num() and test_prob() interventions.

        Args:
            inds: indices of who to test
            test_sensitivity (float): probability of a true positive
            loss_prob (float): probability of loss to follow-up
            test_delay (int): number of days before test results are ready
        '''
        # TODO: Implement this method
        return 

