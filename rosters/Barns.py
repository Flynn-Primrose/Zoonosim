import numpy as np
import sciris as sc
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
            'flock', # uid of the flock residing here
            'repopulations' # Number of times this barn has been repopulated
        ]

        self.states = [
            'uncontaminated', # bool; whether the barn is contaminated
            'contaminated', # bool; whether the barn is contaminated
            'composting', # bool; whether the barn is currently composting
            'cleaning', # bool; whether the barn is currently cleaning
        ]

        self.biosec_states = [
            'green', # bool; standard biosecurity measures in place
            'yellow', # bool; increased biosecurity measures in place (typically due to an outbreak at a different site within 10km)
            'orange', # bool; emergency biosecurity measures in place (typically due to an outbreak at this site)
        ]

        self.variant_states = [
            'contaminated_variant'
        ]

        self.by_variant_states = [
            'contaminated_by_variant'
        ]

        # Set the dates various events took place: these are floats per agent
        self.state_dates = [f'date_{state}' for state in self.states] # Convert each state into a date
            
        self.production_dates = [
            'date_cycle_end', # Date of the end of the production cycle
            'date_repopulate'
        ]

        self.dates = self.state_dates + self.production_dates

        # Duration of different states: these are floats per Barn.
        self.durs = [
            'dur_contamination', # Duration of contamination
            'dur_composting', # Duration of composting process
            'dur_cleaning', # Duration of the cleaning process
        ]

        self.all_recordable_states = self.agent + self.states + self.variant_states + self.dates + self.durs
        self.all_states = self.agent + self.states + self.biosec_states + self.variant_states + self.by_variant_states + self.dates + self.durs

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
  
        self.event_log = []
        self.infection_log = [] # Record of infections - keys for ['source','target','date','layer']

        pop_size = self.pars['pop_size_by_type']['barn']

        # Set person properties -- all floats except for UID
        for key in self.meta.agent:
            if key in ['uid', 'flock', 'repopulations']:
                self[key] = np.zeros(pop_size, dtype=znd.default_int) # NOTE: The uid values are passed in kwargs by make_barn()
            elif key in ['green', 'yellow', 'orange']:
                val = key in ['green']
                self[key] = np.full(pop_size, val, dtype=znd.default_bool) # NOTE: The flock values are passed in kwargs by make_barn()
            else:
                self[key] = np.full(pop_size, np.nan, dtype=znd.default_float)

        # Set states
        for key in self.meta.states:
            val = (key in ['uncontaminated']) # Default value is True for susceptible and naive, False otherwise
            self[key] = np.full(pop_size, val, dtype=bool)

        # Set biosec states
        for key in self.meta.biosec_states:
            val = (key in ['green']) # Default value is True for susceptible and naive, False otherwise
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

    def initialize(self, agents_pars=None):
        ''' Perform initializations '''
        self.validate(roster_pars=agents_pars) # First, check that essential-to-match parameters match
        self.set_pars(agents_pars) # Replace the saved parameters with this simulation's
        self.initialized = True
        return





    def update_states_pre(self, t):
        ''' Perform all state updates at the current timestep '''

        self.t = t # Set the current time

        # Perform updates
        self.init_flows()

        self.flows['new_cleaned'] = self.check_cleaned()
        self.flows['new_uncontaminated'] = self.check_uncontaminated()
        self.flows['new_composted'] = self.check_composted()
        return


    def update_states_post(self):
        ''' Perform post-timestep updates '''


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

    
    def check_uncontaminated(self):
        ''' Check which barns get uncontaminated this timestep '''
        inds = self.check_inds(~self.contaminated, self.date_uncontaminated) # ~self.contaminated is the same as self.uncontaminated
        if len(inds) > 0:
            self.uncontaminated[inds] = True
            self.contaminated[inds] = False
            self.contaminated_variant[inds] = np.nan
            self.date_contaminated[inds] = np.nan
            self.date_uncontaminated[inds] = np.nan

            self.update_event_log(inds, 'uncontaminated')
        return len(inds)
    
    def check_cleaned(self):
        ''' Check which barns get cleaned this timestep '''
        inds = self.check_inds(~self.cleaning, self.date_cleaning)# ~self.cleaning because self.date_cleaning denotes the date when cleaning ends, not when it starts
        if len(inds) > 0:
            self.cleaning[inds] = False
            self.uncontaminated[inds] = True
            self.contaminated[inds] = False
            self.contaminated_variant[inds] = np.nan
            self.contaminated_by_variant[:, inds] = False
            self.date_repopulate[inds] = self.date_cleaning[inds] + 1
            self.date_cleaning[inds] = np.nan
            self.update_event_log(inds, 'cleaning_finished')
        return len(inds)
    
    def check_composted(self):
        ''' Check which barns finish composting this timestep '''
        inds = self.check_inds(~self.composting, self.date_composting)# ~self.composting because self.date_composting denotes the date when composting ends, not when it starts
        if len(inds) > 0:
            self.composting[inds] = False
            self.cleaning[inds] = True
            self.date_cleaning[inds] = self.date_composting[inds] + znu.sample(**self.pars['dur']['barn']['cleaning'], size=len(inds))
            self.date_composting[inds] = np.nan
            self.update_event_log(inds, 'composting_finished')
            self.update_event_log(inds, 'cleaning_started')
        return len(inds)


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

        # Keep only uncontaminated
        keep = self.uncontaminated[inds] # Unique indices in inds and source that are also susceptible
        inds = inds[keep]
        if source is not None:
            source = source[keep]

        # Deal with variant parameters
        variant_keys = ['rel_dur_contamination']
        contamination_pars = {k:self.pars['variant_pars']['wild']['water'][k] for k in variant_keys}
        variant_label = self.pars['variant_map'][variant]
        if variant:
            for k in variant_keys:
                contamination_pars[k] *= self.pars['variant_pars'][variant_label]['barn'][k]

        n_infections = len(inds)
        durpars      = self.pars['dur']['barn']


        # Update states, variant info, and flows
        self.uncontaminated[inds]    = False
        self.contaminated[inds]        = True
        self.contaminated_variant[inds] = variant
        self.contaminated_by_variant[variant, inds] = True
        self.flows['new_contaminated'] += len(inds)
        self.flows_variant['new_contaminated_by_variant'][variant] += len(inds)

        # Record transmissions
        for i, target_ind in enumerate(inds):
            entry = dict(source=source[i] if source is not None else None, target=self.uid[target_ind], date=self.t, layer=layer, variant=variant_label)
            self.infection_log.append(entry)

        # Calculate how long before this person can infect other people
        self.dur_contamination[inds] = znu.sample(**durpars['contamination'], size=n_infections)*contamination_pars['rel_dur_contamination']
        self.date_contaminated[inds] = self.t
        self.date_uncontaminated[inds] = self.dur_contamination[inds] + self.t


        return n_infections # For incrementing counters
    
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
            if lkey.lower() == 'bw':
                llabel = 'barn-water contacts'
            elif lkey.lower() == 'hb':
                llabel = 'human-barn contacts'
            else:
                llabel = f'"{lkey}"'
            return llabel

        uids = sc.promotetolist(uid)
        uids.extend(args)

        for uid in uids:
            uid_ind = np.where(self.uid == uid)[0]
            flock = self.flock[uid_ind]
            repops = self.repopulations[uid_ind]


            intro = f'\nThis is the story of {uid}, a barn housing flock {flock}. It has repopulated {repops} times.'
            print(f'{intro}')

            events = []

            event_dict = {
                'cleaning_started': 'started the cleaning process',
                'cleaning_finished': 'completed the cleaning process',
                'composting_started': 'started the composting process',
                'composting_finished': 'finished the composting process',
                'uncontaminated': 'naturally lost its contamination',
                'cycle_end': 'finished its production cycle',
                'cycle_start':'repopulated and began its production cycle',
            }

            for event in self.event_log:
                if event['target'] == uid:
                    events.append((event['date'], event_dict[event['event_type']]))

            for infection in self.infection_log:
                lkey = infection['layer']
                llabel = label_lkey(lkey)
                if infection['target'] == uid:
                    if lkey:
                        events.append((infection['date'], f'was contaminated with H5N1 by {infection["source"]} via the {llabel} layer'))
                    else:
                        events.append((infection['date'], 'was contaminated with H5N1 as a seed contamination'))

                if infection['source'] == uid:
                    x = len([a for a in self.infection_log if a['source'] == infection['target']])
                    events.append((infection['date'],f'gave H5N1 to {infection["target"]} via the {llabel} layer ({x} secondary infections)'))

            if len(events):
                for day, event in sorted(events, key=lambda x: x[0]):
                    print(f'On day {day:.0f}, {uid} {event}')
            else:
                print(f'Nothing happened to {uid} during the simulation.')
        return
