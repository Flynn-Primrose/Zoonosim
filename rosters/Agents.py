'''
Class definition for the Agents object. The Agents object will itself contain a list of all agents as well as all state information that is common to all agent types. 
It will also include a pars object with the parameters needed to generate the population as well as a subroster for each agent type. Each type specific subroster 
will contain a list of all agents of that type as well as all state information that is specific to that data type.

'''

#%% Imports
import numpy as np
import sciris as sc
from .. import version as znv
from .. import utils as znu
from .. import defaults as znd

from .Roster import Roster

#from . import plotting as cvplt
#from . import immunity as cvi
#from . import watches as cvw


__all__ = ['Agents']

class AgentsMeta(sc.prettyobj):
    '''
    For storing all keys that are common across all agent types
    '''

    def __init__(self):

        self.agent = [
            'uid', #Int
            'agent_type', #string? the type of agent, must be one of pars['agent_types']
            'rel_trans', # float - relative transmissibility of the agent
            'rel_sus', # float - relative susceptibility of the agent
        ]

        self.states = [ # all boolean
            'symptomatic', #?
            'susceptible', #?
            'exposed', #?
            'infectious', #?
            'quarantined'

        ]

        self.variant_states =[
            'exposed_variant', # int
            'infectious_variant' # int
        ]

        self.by_variant_states = [
            'exposed_by_variant',
            'infectious_by_variant',
        ]


        self.all_states = self.agent + self.states +self.variant_states + self.by_variant_states

        # Validate
        self.state_types = ['agent', 'states', 'variant_states', 'by_variant_states', 'all_states']
        for state_type in self.state_types:
            states = getattr(self, state_type)
            n_states        = len(states)
            n_unique_states = len(set(states))
            if n_states != n_unique_states: # pragma: no cover
                errormsg = f'In {state_type}, only {n_unique_states} of {n_states} state names are unique'
                raise ValueError(errormsg)

        return

class Agents(Roster):
    '''
    A class to perform all the operations on the agents -- usually not invoked directly.

    Note that this class handles the mechanics of updating the actual agents, while
    ``cv.BaseAgents`` takes care of housekeeping (saving, loading, exporting, etc.).
    Please see the BaseAgents class for additional methods.

    Args:
        pars (dict): the sim parameters, e.g. sim.pars -- alternatively, if a number, interpreted as pop_size
        strict (bool): whether or not to only create keys that are already in self.meta.agents; otherwise, let any key be set
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
        self.meta = AgentsMeta() # Store list of keys and dtypes
        self.contacts = None
        self.human = kwargs['human'] if 'human' in kwargs else None 
        self.flock = kwargs['flock'] if 'flock' in kwargs else None
        self.barn = kwargs['barn'] if 'barn' in kwargs else None
        self.water = kwargs['water'] if 'water' in kwargs else None        
        self.init_contacts() # Initialize the contacts
        self.infection_log = [] # Record of infections - keys for ['source','target','date','layer']

        # Set the population size

        # if self.pars['pop_size'] is None:
        #     pop_size = 0
        #     for type in self.pars['agent_types']:
        #         pop_size += self.pars['pop_size_by_type'][type]
        #     self.pars['pop_size'] = pop_size
        
        # Set agent properties -- all floats except for UID
        for key in self.meta.agent:
            if key == 'uid':
                self[key] = np.arange(self.pars['pop_size'], dtype=znd.default_int)
            elif key == 'agent_type':
                self[key] = np.full(self.pars['pop_size'], 'unknown', dtype=znd.default_str)
            else:
                self[key] = np.full(self.pars['pop_size'], np.nan, dtype=znd.default_float)

        # Set health states -- only susceptible is true by default -- booleans except exposed by variant which should return the variant that ind is exposed to
        for key in self.meta.states:
            val = (key in ['susceptible']) # Default value is True for susceptible, False otherwise
            self[key] = np.full(self.pars['pop_size'], val, dtype=bool)

        # Set variant states, which store info about which variant a person is exposed to
        for key in self.meta.variant_states:
            self[key] = np.full(self.pars['pop_size'], np.nan, dtype=znd.default_float)
        for key in self.meta.by_variant_states:
            self[key] = np.full((self.pars['n_variants'], self.pars['pop_size']), False, dtype=bool)

        # Store the dtypes used in a flat dict
        self._dtypes = {key:self[key].dtype for key in self.keys()} # Assign all to float by default
        if strict:
            self.lock() # If strict is true, stop further keys from being set (does not affect attributes)

        # Store flows to be computed during simulation
        self.init_flows()

        # Although we have called init(), we still need to call initialize()
        self.initialized = False

        # Handle contacts, if supplied (note: they usually are)
        if 'contacts' in kwargs:
            self.add_contacts(kwargs.pop('contacts'))

        # Handle all other values, e.g. age
        for key,value in kwargs.items():
            if key not in ['human', 'flock', 'barn', 'water']: # These are handled separately
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
        self.initialized = True
        return
    

    def update_states_pre(self, t):
        ''' Perform all state updates at the current timestep '''

        # Initialize
        self.t = t
        self.is_exp = self.true('exposed') # For storing the interim values since used in every subsequent calculation

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
    
    def update_states_from_subrosters(self):
        susceptible_human_uids = np.array(self.human['uid'][znu.true(self.human['susceptible'])])
        exposed_human_uids = np.array(self.human['uid'][znu.true(self.human['exposed'])])
        infectious_human_uids = np.array(self.human['uid'][znu.true(self.human['infectious'])])
        symptomatic_human_uids = np.array(self.human['uid'][znu.true(self.human['symptomatic'])])
        quarantined_human_uids = np.array(self.human['uid'][znu.true(self.human['quarantined'])])

        susceptible_flock_uids = np.array(self.flock['uid'][znu.true(self.flock['susceptible'])])
        exposed_flock_uids = np.array(self.flock['uid'][znu.true(self.flock['exposed'])])
        infectious_flock_uids = np.array(self.flock['uid'][znu.true(self.flock['infectious'])])
        symptomatic_flock_uids = np.array(self.flock['uid'][znu.true(self.flock['symptomatic'])])
        quarantined_flock_uids = np.array(self.flock['uid'][znu.true(self.flock['quarantined'])])

        susceptible_barn_uids = np.array(self.barn['uid'][znu.false(self.barn['contaminated'])])
        exposed_barn_uids = np.array(self.barn['uid'][znu.true(self.barn['contaminated'])])
        infectious_barn_uids = np.array(self.barn['uid'][znu.true(self.barn['contaminated'])])

        susceptible_water_uids = np.array(self.water['uid'][znu.false(self.water['contaminated'])])
        exposed_water_uids = np.array(self.water['uid'][znu.true(self.water['contaminated'])])
        infectious_water_uids = np.array(self.water['uid'][znu.true(self.water['contaminated'])])

        susceptible_uids = np.concatenate((susceptible_human_uids, susceptible_flock_uids, susceptible_barn_uids, susceptible_water_uids))
        exposed_uids = np.concatenate((exposed_human_uids, exposed_flock_uids, exposed_barn_uids, exposed_water_uids))
        infectious_uids = np.concatenate((infectious_human_uids, infectious_flock_uids, infectious_barn_uids, infectious_water_uids))
        symptomatic_uids = np.concatenate((symptomatic_human_uids, symptomatic_flock_uids))
        quarantined_uids = np.concatenate((quarantined_human_uids, quarantined_flock_uids))

        self.susceptible = np.isin(self['uid'], susceptible_uids)
        self.exposed = np.isin(self['uid'], exposed_uids)
        self.infectious = np.isin(self['uid'], infectious_uids)
        self.symptomatic = np.isin(self['uid'], symptomatic_uids)
        self.quarantined = np.isin(self['uid'], quarantined_uids)

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
        inds = self.check_inds(self.infectious, self.date_infectious, filter_inds=self.is_exp)
        self.infectious[inds] = True
        self.infectious_variant[inds] = self.exposed_variant[inds]


        for variant in range(self.pars['n_variants']):
            this_variant_inds = znu.itrue(self.infectious_variant[inds] == variant, inds)
            n_this_variant_inds = len(this_variant_inds)
            self.flows_variant['new_infectious_by_variant'][variant] += n_this_variant_inds
            self.infectious_by_variant[variant, this_variant_inds] = True
        return len(inds)




    #%% Methods to make events occur (infection and diagnosis)

    def infect_type(self, agent_type, inds, update = True):
        self[agent_type].infect(inds = inds, layer = 'seed_infections')
        if update: self.update_states_from_subrosters()
        return


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