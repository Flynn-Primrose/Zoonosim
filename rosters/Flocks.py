import numpy as np
import scipy.stats as stats
import sciris as sc
import scipy.stats as stats
from collections import defaultdict
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
            'barn' # uid of the barn where the flock is located
        ]

        self.states = [
            'susceptible',
            'exposed',
            'infectious',
            'symptomatic',
        ]

        self.variant_states = [
            'exposed_variant',
            'infectious_variant',
            'symptomatic_variant',
        ]

        self.by_variant_states = [
            'exposed_by_variant',
            'infectious_by_variant',
        ]

        # Set the dates various events took place: these are floats per agent
        self.dates = [f'date_{state}' for state in self.states] # Convert each state into a date

        # Duration of different states: these are floats per flock.
        self.durs = [
            'dur_exp2inf',
            'dur_inf2sym',
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
        self.flows = {key:0 for key in znd.new_result_flows}
        self.flows_variant = {}
        for key in znd.new_result_flows_by_variant:
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
        Set the prognoses for each person based on age during initialization. Need
        to reset the seed because viral loads are drawn stochastically.
        '''
        # TODO: Figure out how to model infections in flocks
        pars = self.pars # Shorten
        if 'prognoses' not in pars or 'rand_seed' not in pars:
            errormsg = 'This poultry object does not have the required parameters ("prognoses" and "rand_seed"). Create a sim (or parameters), then do e.g. people.set_pars(sim.pars).'
            raise sc.KeyNotFoundError(errormsg)

        znu.set_seed(pars['rand_seed'])

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

        self.flows['new_critical']      += self.check_critical()
        self.flows['new_recoveries']    += self.check_recovery()
        new_deaths, new_known_deaths     = self.check_death()
        self.flows['new_deaths']        += new_deaths
        self.flows['new_known_deaths']  += new_known_deaths
        return


    def update_states_post(self):
        ''' Perform post-timestep updates '''

        # Update the status of flocks

        del self.is_exp  # Tidy up

        return



    def update_contacts(self):
        ''' Refresh dynamic contacts, e.g. community '''
        # Figure out if anything needs to be done -- e.g. {'h':False, 'c':True}
        for lkey, is_dynam in self.pars['dynam_layer'].items():
            if is_dynam:
                self.contacts[lkey].update(self)

        return self.contacts


    def schedule_behaviour(self, behaviour_pars):
        ''' Schedules events on the basis of results received today '''


        return

    
    def get_key_to_use(self, policy_dict):
        ''' Helper function used to determine what policy to use when scheduling behaviour '''

        keys = np.array(list(policy_dict.keys()))
        key_to_use = keys[znu.true(keys <= self.t)[-1]]

        return key_to_use


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
        
        #  STRATIFICATION
        if self.pars['enable_stratifications']:
            for strat, strat_pars in self.pars['stratification_pars'].items():
                metrics = strat_pars['metrics']

                if 'new_severe' in metrics:
                    bracs = strat_pars['brackets']
                    
                    # update the value for the current day, for each bracket. 
                    for brac in bracs:
                        brac_name = "_".join([ str(brac[0]), str(brac[1])])

                        # Count the number of individuals meeting the criteria. 
                        if strat == 'age':
                            num_inds = np.sum( (self.age[inds] >= brac[0]) & (self.age[inds] < brac[1]) )
                        elif strat == 'income':
                            num_inds = np.sum( (self.income[inds] >= brac[0]) & (self.income[inds] < brac[1]) )
                        else:
                            raise ValueError(f"Stratification {strat} not recognized.")

                        self.stratifications[strat][brac_name]['new_severe'][self.t] += num_inds
        
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
        self.critical[inds]         = False
        self.recovered[inds]        = True
        self.recovered_variant[inds] = self.exposed_variant[inds]
        self.infectious_variant[inds] = np.nan
        self.exposed_variant[inds]    = np.nan
        self.exposed_by_variant[:, inds] = False
        self.infectious_by_variant[:, inds] = False

        # Handle immunity aspects
        if self.pars['use_waning']:

            # Reset additional states
            self.susceptible[inds] = True
            self.diagnosed[inds]   = False # Reset their diagnosis state because they might be reinfected

        # # Handle instances of false positive diagnoses; if someone is not exposed, is diagnosed, and has passed the required time in isolation, un-diagnose them
        # fp_diagnoses = znu.true((~self.exposed) & (self.diagnosed) & (self.t - self.date_diagnosed >= 14))
        # self.diagnosed[fp_diagnoses] = False

        return len(inds)


    def check_death(self):
        ''' Check whether or not this agent died on this timestep  '''
        inds = self.check_inds(self.dead, self.date_dead, filter_inds=self.is_exp) 
        
        self.dead[inds]             = True
        diag_inds = inds[self.diagnosed[inds]] # Check whether the person was diagnosed before dying
        # self.known_dead[diag_inds]  = True
        # self.susceptible[inds]      = False
        # self.exposed[inds]          = False
        # self.infectious[inds]       = False
        # self.symptomatic[inds]      = False
        # self.severe[inds]           = False
        # self.critical[inds]         = False
        # self.known_contact[inds]    = False
        # self.quarantined[inds]      = False
        # self.recovered[inds]        = False
        # self.infectious_variant[inds] = np.nan
        # self.exposed_variant[inds]    = np.nan
        # self.recovered_variant[inds]  = np.nan

        return len(inds), len(diag_inds)



    #%% Methods to make events occur (infection and diagnosis)




    def infect(self, inds, source=None, layer=None, variant=0):
        '''
        Infect agents and determine their eventual outcomes.

            * Every infected person can infect other people, regardless of whether they develop symptoms
            * Infected people that develop symptoms are disaggregated into mild vs. severe (=requires hospitalization) vs. critical (=requires ICU)
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
            count (int): number of poultry agents infected
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
        infect_pars = {k:self.pars[k] for k in variant_keys}
        variant_label = self.pars['variant_map'][variant]
        if variant:
            for k in variant_keys:
                infect_pars[k] *= self.pars['variant_pars'][variant_label][k]

        n_infections = len(inds)
        durpars      = self.pars['dur']

        # Update states, variant info, and flows
        self.susceptible[inds]    = False
        self.recovered[inds]      = False
        self.diagnosed[inds]      = False
        self.exposed[inds]        = True
        self.exposed_variant[inds] = variant
        self.exposed_by_variant[variant, inds] = True
        self.flows['new_infections']   += len(inds)
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

        return # For incrementing counters


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

        inds = np.unique(inds)
        self.tested[inds] = True
        self.date_tested[inds] = self.t # Only keep the last time they tested

        is_infectious = znu.itruei(self.infectious, inds)
        pos_test      = znu.n_binomial(test_sensitivity, len(is_infectious))
        is_inf_pos    = is_infectious[pos_test]

        not_diagnosed = is_inf_pos[np.isnan(self.date_diagnosed[is_inf_pos])]
        not_lost      = znu.n_binomial(1.0-loss_prob, len(not_diagnosed))
        final_inds    = not_diagnosed[not_lost]

        # Store the date the flock will be diagnosed, as well as the date they took the test which will come back positive
        self.date_diagnosed[final_inds] = self.t + test_delay
        self.date_pos_test[final_inds] = self.t

        return final_inds



    #%% Analysis methods

    # def plot(self, *args, **kwargs):
    #     '''
    #     Plot statistics of the population -- age distribution, numbers of contacts,
    #     and overall weight of contacts (number of contacts multiplied by beta per
    #     layer).

    #     Args:
    #         bins      (arr)   : age bins to use (default, 0-100 in one-year bins)
    #         width     (float) : bar width
    #         font_size (float) : size of font
    #         alpha     (float) : transparency of the plots
    #         fig_args  (dict)  : passed to pl.figure()
    #         axis_args (dict)  : passed to pl.subplots_adjust()
    #         plot_args (dict)  : passed to pl.plot()
    #         do_show   (bool)  : whether to show the plot
    #         fig       (fig)   : handle of existing figure to plot into
    #     '''
    #     fig = cvplt.plot_people(people=self, *args, **kwargs)
    #     return fig


    # def story(self, uid, *args):
    #     '''
    #     Print out a short history of events in the life of the specified individual.

    #     Args:
    #         uid (int/list): the person or people whose story is being regaled
    #         args (list): these people will tell their stories too

    #     **Example**::

    #         sim = cv.Sim(pop_type='hybrid', verbose=0)
    #         sim.run()
    #         sim.people.story(12)
    #         sim.people.story(795)
    #     '''

    #     def label_lkey(lkey):
    #         ''' Friendly name for common layer keys '''
    #         if lkey.lower() == 'a':
    #             llabel = 'default contact'
    #         if lkey.lower() == 'h':
    #             llabel = 'household'
    #         elif lkey.lower() == 's':
    #             llabel = 'school'
    #         elif lkey.lower() == 'w':
    #             llabel = 'workplace'
    #         elif lkey.lower() == 'c':
    #             llabel = 'community'
    #         else:
    #             llabel = f'"{lkey}"'
    #         return llabel

    #     uids = sc.promotetolist(uid)
    #     uids.extend(args)

    #     for uid in uids:

    #         p = self[uid]
    #         sex = 'female' if p.sex == 0 else 'male'

    #         intro = f'\nThis is the story of {uid}, a {p.age:.0f} year old {sex}'

    #         if not p.susceptible:
    #             if np.isnan(p.date_symptomatic):
    #                 print(f'{intro}, who had asymptomatic COVID.')
    #             else:
    #                 print(f'{intro}, who had symptomatic COVID.')
    #         else:
    #             print(f'{intro}, who did not contract COVID.')

    #         total_contacts = 0
    #         no_contacts = []
    #         for lkey in p.contacts.keys():
    #             llabel = label_lkey(lkey)
    #             n_contacts = len(p.contacts[lkey])
    #             total_contacts += n_contacts
    #             if n_contacts:
    #                 print(f'{uid} is connected to {n_contacts} people in the {llabel} layer')
    #             else:
    #                 no_contacts.append(llabel)
    #         if len(no_contacts):
    #             nc_string = ', '.join(no_contacts)
    #             print(f'{uid} has no contacts in the {nc_string} layer(s)')
    #         print(f'{uid} has {total_contacts} contacts in total')

    #         events = []

    #         dates = {
    #             'date_critical'       : 'became critically ill and needed ICU care',
    #             'date_dead'           : 'died ☹',
    #             'date_diagnosed'      : 'was diagnosed with COVID',
    #             'date_end_quarantine' : 'ended quarantine',
    #             'date_infectious'     : 'became infectious',
    #             'date_known_contact'  : 'was notified they may have been exposed to COVID',
    #             'date_pos_test'       : 'took a positive test',
    #             'date_quarantined'    : 'entered quarantine',
    #             'date_recovered'      : 'recovered',
    #             'date_severe'         : 'developed severe symptoms and needed hospitalization',
    #             'date_symptomatic'    : 'became symptomatic',
    #             'date_tested'         : 'was tested for COVID',
    #             'date_vaccinated'     : 'was vaccinated against COVID',
    #         }

    #         for attribute, message in dates.items():
    #             date = getattr(p,attribute)
    #             if not np.isnan(date):
    #                 events.append((date, message))

    #         for infection in self.infection_log:
    #             lkey = infection['layer']
    #             llabel = label_lkey(lkey)
    #             if infection['target'] == uid:
    #                 if lkey:
    #                     events.append((infection['date'], f'was infected with COVID by {infection["source"]} via the {llabel} layer'))
    #                 else:
    #                     events.append((infection['date'], 'was infected with COVID as a seed infection'))

    #             if infection['source'] == uid:
    #                 x = len([a for a in self.infection_log if a['source'] == infection['target']])
    #                 events.append((infection['date'],f'gave COVID to {infection["target"]} via the {llabel} layer ({x} secondary infections)'))

    #         if len(events):
    #             for day, event in sorted(events, key=lambda x: x[0]):
    #                 print(f'On day {day:.0f}, {uid} {event}')
    #         else:
    #             print(f'Nothing happened to {uid} during the simulation.')
    #     return