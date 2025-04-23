import numpy as np
import scipy.stats as stats
import sciris as sc
import scipy.stats as stats
from collections import defaultdict
from .. import version as znv
from .. import utils as znu
from .. import defaults as znd
from .Subroster import Subroster

__all__ = ['Barns']

class BarnMeta(sc.prettyobj):
    ''' Defines all keys that are used by barns '''

    def __init__(self):
        
        self.agent = [
            'uid', # int
            'temperature',
            'humidity',
            'flock' # uid of the flock residing here
        ]

        self.states = [
            'contaminated'
        ]

        self.variant_states = [
            'contaminated_variant'
        ]

        self.by_variant_states = [
            'contaminated_by_variant'
        ]

        # Set the dates various events took place: these are floats per agent
        self.dates = [f'date_{state}' for state in self.states] # Convert each state into a date

        # Duration of different states: these are floats per Barn.
        self.durs = [
            'dur_contamination', # Duration of contamination
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
    

class Barns(Subroster):
    '''
    A class to perform all the operations on the agents -- usually not invoked directly.

    This class is usually created automatically by the sim. The only required input
    argument is the population size, but typically the full parameters dictionary
    will get passed instead since it will be needed before the People object is
    initialized. However, ages, contacts, etc. will need to be created separately --
    see ``cv.make_people()`` instead.

    Note that this class handles the mechanics of updating the actual agents, while
    ``cv.BaseRoster`` takes care of housekeeping (saving, loading, exporting, etc.).
    Please see the BaseRoster class for additional methods.

    Args:
        pars (dict): the sim parameters, e.g. sim.pars
        strict (bool): whether or not to only create keys that are already in self.meta.person; otherwise, let any key be set
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
        self.meta = BarnMeta() # Store list of keys and dtypes
  

        self.infection_log = [] # Record of infections - keys for ['source','target','date','layer']

        pop_size = self.pars['pop_size_by_type']['barn']

        # Set person properties -- all floats except for UID
        for key in self.meta.agent:
            if key == 'uid':
                self[key] = np.zeros(pop_size, dtype=znd.default_int) # NOTE: The uid values are passed in kwargs by make_barn()
            elif key == 'flock':
                self[key] = np.zeros(pop_size, dtype=znd.default_int) # NOTE: The flock values are passed in kwargs by make_barn()
            else:
                self[key] = np.full(pop_size, np.nan, dtype=znd.default_float)

        # Set health states -- only susceptible is true by default -- booleans except exposed by variant which should return the variant that ind is exposed to
        for key in self.meta.states:
            val = (key in ['susceptible']) # Default value is True for susceptible and naive, False otherwise
            self[key] = np.full(pop_size, val, dtype=bool)

        # Set variant states, which store info about which variant a person is exposed to
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
        self.flows = {key:0 for key in znd.new_barn_flows}
        self.flows_variant = {}
        for key in znd.new_barn_flows_by_variant:
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

        pars = self.pars # Shorten
        if 'prognoses' not in pars or 'rand_seed' not in pars:
            errormsg = 'This Barns object does not have the required parameters ("prognoses" and "rand_seed").'
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
        self.flows['new_infectious']    += self.check_infectious() # For Barns that are exposed and not infectious, check if they begin being infectious
        return


    def update_states_post(self):
        ''' Perform post-timestep updates '''


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



    def check_cleaned(self):
        ''' Check which barns get cleaned this timestep '''
        # TODO: create this function
        return



    #%% Methods to make events occur (infection and diagnosis)


    def infect(self, inds, source=None, layer=None, variant=0):
        '''
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
        variant_label = self.pars['variant_map'][variant]

        n_infections = len(inds)
        durpars      = self.pars['dur']['barn']


        # Update states, variant info, and flows
        self.susceptible[inds]    = False
        self.contaminated[inds]        = True
        self.contaminated_variant[inds] = variant
        self.contaminated_by_variant[variant, inds] = True
        self.flows['new_barn_contaminated'] += len(inds)
        self.flows_variant['new_barn_contaminated_by_variant'][variant] += len(inds)

        # Record transmissions
        for i, target in enumerate(inds):
            entry = dict(source=source[i] if source is not None else None, target=target, date=self.t, layer=layer, variant=variant_label)
            self.infection_log.append(entry)

        # Calculate how long before this person can infect other people
        self.dur_exp2inf[inds] = znu.sample(**durpars['exp2inf'], size=n_infections)
        self.date_exposed[inds] = self.t
        self.date_infectious[inds] = self.dur_exp2inf[inds] + self.t


        return n_infections # For incrementing counters


    def test(self, inds, test_sensitivity=1.0, loss_prob=0.0, test_delay=0):
        '''
        Method to test agents. Typically not to be called by the user directly;
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