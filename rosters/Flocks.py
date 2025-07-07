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
            'exposed_headcount',
            'infectious_headcount',
            'symptomatic_headcount',
            'daily_dead_headcount', # Daily number of dead birds in the flock
            'total_dead_headcount', # Total number of dead birds in the flock
            'baseline_symptomatic_rate', 
            'infected_symptomatic_rate',
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

        self.variant_states = [
            'exposed_variant',
            'infectious_variant',
            'beta_variant', # Variant-specific beta
            'gamma_variant', # Variant-specific gamma
        ]

        self.by_variant_states = [
            'exposed_by_variant',
            'infectious_by_variant',
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
            'dur_inf2out', # Mean time from infectious to outcome (recovered/removed)
            'dur_susp2res', # Mean time from inspected to confirmation of status 
        ]


        self.all_recordable_states = self.agent + self.states + self.variant_states + self.dates + self.durs
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
        self.event_log = [] # Record of events - keys for ['target', 'event_type', 'date']
        self.infection_log = [] # Record of infections - keys for ['source','target','date','layer']
        
        pop_size = pars['pop_size_by_type']['flock']

        # Set flock properties
        for key in self.meta.agent:
            if key in ['uid', 'barn']:
                self[key] = np.zeros(pop_size, dtype=znd.default_int)# NOTE: values get passed as kwargs by make_flock
            elif key == 'breed':
                self[key] = np.zeros(pop_size, dtype=znd.default_str) # NOTE: values get passed as kwargs by make_flock
            else:
                self[key] = np.zeros(pop_size,dtype=znd.default_float)

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

    def initialize(self, agents_pars=None):
        ''' Perform initializations '''
        self.validate(roster_pars=agents_pars) # First, check that essential-to-match parameters match
        self.set_pars(agents_pars) # Replace the saved parameters with this simulation's
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
        self.flows['new_suspected'] = self.check_suspected()
        self.flows['new_quarantined'] = self.check_quarantined() # 
        self.update_water_consumption() # Update the headcounts and water consumption of the flocks
        return


    def update_states_post(self):
        ''' Perform post-timestep updates '''

        # Update the status of flocks

        return

    def update_event_log(self, target_inds, event_type):
        '''
        Add an entry to the event log

        args:
            target_inds: array of indices of flocks that experienced a recordable event
            event (str): The specific event in question
        '''

        if target_inds is None:
            return
        
        for ind in target_inds:
            entry = dict(target = self.uid[ind], event_type = event_type, date = self.t)
            self.event_log.append(entry)

        return

    #%% Methods for updating state

    def check_suspected(self):
        ''' Check for new progressions to suspected '''

        unsuspected_inds = np.where((self.suspected == False) & (self.headcount>0))[0]
        if len(unsuspected_inds) == 0:
            return 0
        actual_symptomatic_rate = self.symptomatic_headcount[unsuspected_inds] / self.headcount[unsuspected_inds]
        actual_mortality_rate = self.daily_dead_headcount[unsuspected_inds] / self.headcount[unsuspected_inds]
        actual_water_rate = self.water_consumption[unsuspected_inds] / self.headcount[unsuspected_inds]

        suspicious_symptomatic_inds = np.where(actual_symptomatic_rate > znd.default_suspicious_symptomatic_rate)[0]
        suspicious_mortality_inds = np.where(actual_mortality_rate > znd.default_suspicious_mortality_rate)[0]
        suspicious_water_inds = np.where(actual_water_rate > znd.default_suspicious_consumption_rate)[0]
        suspicious_inds = np.unique(np.concatenate((suspicious_symptomatic_inds, suspicious_mortality_inds, suspicious_water_inds)))
        if len(suspicious_inds) == 0:
            return 0
        new_suspicious_inds = unsuspected_inds[suspicious_inds]
        self.suspected[new_suspicious_inds] = True
        self.date_suspected[new_suspicious_inds] = self.t
        self.dur_susp2res[new_suspicious_inds] = znu.sample(**self.pars['dur']['flock']['susp2res'], size=len(new_suspicious_inds))
        self.date_result[new_suspicious_inds] = self.date_suspected[new_suspicious_inds] + self.dur_susp2res[new_suspicious_inds]
        self.date_quarantined[new_suspicious_inds] = self.t

        self.update_event_log(new_suspicious_inds, 'suspected')


        return len(suspicious_inds)


    def check_infectious(self):
        ''' Check if they become infectious '''

        inds = np.where((self.infectious_headcount >= 1) & (self.infectious == False))[0] #Question: what is the cutoff for being considered infectious

        self.infectious[inds] = True
        self.infectious_variant[inds] = self.exposed_variant[inds]
        self.date_infectious[inds] = self.t
        self.update_event_log(inds, 'infectious')


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
        self.update_event_log(inds, 'quarantined')
        return len(inds)


    def update_water_consumption(self):
        ''' Update the water consumption of the flocks '''
        infected_headcount = self.exposed_headcount + self.infectious_headcount
        uninfected_headcount = self.headcount - infected_headcount
        self.water_consumption = infected_headcount * self.infected_water_rate + uninfected_headcount * self.baseline_water_rate
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
        infect_pars = {k:self.pars['variant_pars']['wild']['flock'][k] for k in variant_keys}
        variant_label = self.pars['variant_map'][variant]
        if variant:
            for k in variant_keys:
                infect_pars[k] *= self.pars['variant_pars'][variant_label]['flock'][k]

        n_infections = len(inds)
        durpars      = self.pars['dur']['flock']

        # Determine how many birds are initially exposed for each flock
        initial_exposed = znu.n_multinomial(np.ones(6)/6.0, len(inds)) + 1 # All flocks start with between 1 and 6 initial exposures

        # Update states, variant info, and flows
        self.susceptible[inds]    = False
        self.exposed[inds]        = True
        self.exposed_headcount[inds] = initial_exposed
        self.exposed_variant[inds] = variant
        self.exposed_by_variant[variant, inds] = True
        self.flows['new_exposed']   += len(inds)
        self.flows_variant['new_exposed_by_variant'][variant] += len(inds)

        # Record transmissions
        for i, target_ind in enumerate(inds):
            entry = dict(source=source[i] if source is not None else None, target=self.uid[target_ind], date=self.t, layer=layer, variant=variant_label)
            self.infection_log.append(entry)

        # Calculate how long before this flock can infect other flocks
        self.dur_exp2inf[inds] = np.maximum(znu.sample(**durpars['exp2inf'], size=n_infections), 0) # Ensure that this is not negative
        self.dur_inf2out[inds] = np.maximum(znu.sample(**durpars['inf2out'], size=n_infections), 0)
        self.date_exposed[inds] = self.t

        #Update water rate, symptomatic rate, and mortality rate.
        pars = self.pars # Shorten
        if 'prognoses' not in pars or 'rand_seed' not in pars:
            errormsg = 'This flock object does not have the required parameters ("prognoses" and "rand_seed").'
            raise sc.KeyNotFoundError(errormsg)
        progs = pars['prognoses']['flock']
        breed_to_index = {breed: index for index, breed in enumerate(progs['breeds'])}
        breed_inds = np.array([breed_to_index[this_breed] for this_breed in self.breed[inds]])
        breed, frequency = np.unique(breed_inds, return_counts=True)
        breed_freq = zip(breed, frequency)
        for breed, frequency in breed_freq:
            # NOTE: I'm just guessing at the distribution of these parameters.
            self.infected_symptomatic_rate[inds[breed_inds == breed]] = self.baseline_symptomatic_rate[inds[breed_inds == breed]] + np.maximum(znu.sample('lognormal', progs['mean_symptomatic_rate_increase'][breed], 1.0, size=frequency), 0)*infect_pars['rel_symp_prob']
            self.infected_mortality_rate[inds[breed_inds == breed]] = self.baseline_mortality_rate[inds[breed_inds == breed]] + np.maximum(znu.sample('lognormal', progs['mean_mortality_rate_increase'][breed], 1.0, size=frequency), 0)*infect_pars['rel_death_prob']
            self.infected_water_rate[inds[breed_inds == breed]] = self.baseline_water_rate[inds[breed_inds == breed]] + np.maximum(znu.sample('lognormal', progs['mean_water_rate_increase'][breed], 1.0, size=frequency), 0)


        return n_infections # For incrementing counters


    def test(self, inds, sample_size):
        '''
        Method to test poultry. Typically not to be called by the user directly;
        see interventions.

        Args:
            inds: indices of who to test
            sample_size (int): The number of birds to be tested in each flock
        '''
        unsuspected_inds = np.where(self.suspected == False)[0]
        inds_to_test = np.intersect1d(inds, unsuspected_inds)

        already_suspected_inds = np.setdiff1d(inds, inds_to_test)
        self.update_event_log(already_suspected_inds, 'screened_sus')

        prop_infected = (self.exposed_headcount[inds_to_test] + self.infectious_headcount[inds_to_test])/self.headcount[inds_to_test] #exposed + infectious = infected, infected/total=proportion infected
        non_zero_inds = prop_infected.non_zero()[0]
        zero_inds = np.where(prop_infected == 0)[0]
        prob_pos = np.ones(len(non_zero_inds)) - ((np.ones(len(non_zero_inds)) - prop_infected[non_zero_inds])**sample_size)

        positives = (np.random.random(len(non_zero_inds)) < prob_pos[non_zero_inds])
        pos_inds = inds_to_test[non_zero_inds[positives]]
        neg_inds = np.concatenate(inds_to_test[non_zero_inds[positives == False]], inds_to_test[zero_inds])

        if len(pos_inds)>0:
            self.suspected[pos_inds] = True
            self.date_suspected[pos_inds] = self.t
            self.dur_susp2res[pos_inds] = znu.sample(**self.pars['dur']['flock']['susp2res'], size=len(pos_inds))
            self.date_result[pos_inds] = self.date_suspected[pos_inds] + self.dur_susp2res[pos_inds]
            self.date_quarantined[pos_inds] = self.t
            self.update_event_log(pos_inds, 'screened_pos')

        if len(neg_inds)>0:
            self.update_event_log(neg_inds, 'screened_neg')

        return
    
    def story(self, uid, *args):
        '''
        Print out a short history of events in the life of the specified flock.

        Args:
            uid (int/list): the flock or flocks whose story is to be regaled
            args (list): these flocks will tell their stories too

        **Example**::

            sim = cv.Sim()
            sim.run()
            sim.agents.flock.story(12)
            sim.agents.flock.story(795)
        '''

        def label_lkey(lkey):
            ''' Friendly name for common layer keys '''
            if lkey.lower() == 'fb':
                llabel = 'flock-barn contacts'
            if lkey.lower() == 'hf':
                llabel = 'human-flock contacts'
            elif lkey.lower() == 'fw':
                llabel = 'flock-water contacts'
            else:
                llabel = f'"{lkey}"'
            return llabel

        uids = sc.promotetolist(uid)
        uids.extend(args)

        for uid in uids:
            uid_ind = np.where(self.uid == uid)[0]
            breed = self.breed[uid_ind].item()
            barn = self.barn[uid_ind].item()

            intro = f'\nThis is the story of {uid}, a flock of {breed}s housed in barn {barn}.'
            print(f'{intro}')

            events = []

            event_dict = {
                'suspected': 'was suspected of being infected with H5N1',
                'infectious': 'became infectious with H5N1',
                'quarantined': 'was quarantined',
                'cull': 'was culled',
                'cycle_end': 'finished its production cycle',
                'cycle_start':'began its production cycle',
                'negative': 'was confirmed to be negative for H5N1', 
                'screened_pos': 'tested positive during screening',
                'screened_neg': 'tested negative during screening',
                'screened_sus': 'skipped regular screening due to prior suspicion'
            }

            for event in self.event_log:
                if event['target'] == uid:
                    events.append((event['date'], event_dict[event['event_type']]))

            for infection in self.infection_log:
                lkey = infection['layer']
                llabel = label_lkey(lkey)
                if infection['target'] == uid:
                    if lkey:
                        events.append((infection['date'], f'was infected with H5N1 by {infection["source"]} via the {llabel} layer'))
                    else:
                        events.append((infection['date'], 'was infected with H5N1 as a seed infection'))

                if infection['source'] == uid:
                    x = len([a for a in self.infection_log if a['source'] == infection['target']])
                    events.append((infection['date'],f'gave H5N1 to {infection["target"]} via the {llabel} layer ({x} secondary infections)'))

            if len(events):
                for day, event in sorted(events, key=lambda x: x[0]):
                    print(f'On day {day:.0f}, {uid} {event}')
            else:
                print(f'Nothing happened to {uid} during the simulation.')
        return

