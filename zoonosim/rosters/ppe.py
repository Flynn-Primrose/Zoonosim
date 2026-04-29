'''
Defines the subroster and metaroster for personal protective equipment 
'''

import numpy as np
import sciris as sc
from collections import defaultdict
from .. import version as znv
from .. import utils as znu
from .. import defaults as znd
from .rosters import Subroster

__all__ = ['PPE', 'PPEMeta']

class PPEMeta(sc.prettyobj):
    ''' Defines all keys that are used by PPE '''

    def __init__(self):
        
        self.agent = [
            'uid', # int
            #'temperature',
            #'humidity',
            'rel_trans',        # Float
            'rel_sus',          # Float
            'human', # uid of the human agent assigned to this equipment
            'cons_days_in_quar',    # Int
        ]

        self.states = [
            'uncontaminated', # bool; whether the barn is contaminated
            'contaminated', # bool; whether the barn is contaminated
            'quarantined', # bool; whether the barn is under quarantine
        ]

        # self.biosec_states = [
        #     'green', # bool; standard biosecurity measures in place
        #     'yellow', # bool; increased biosecurity measures in place (typically due to an outbreak at a different site within 10km)
        #     'orange', # bool; emergency biosecurity measures in place (typically due to an outbreak at this site)
        # ]

        self.variant_states = [
            'contaminated_variant'
        ]

        self.by_variant_states = [
            'contaminated_by_variant'
        ]

        # Set the dates various events took place: these are floats per agent
        self.state_dates = [f'date_{state}' for state in self.states] # Convert each state into a date

        self.dates = self.state_dates
        self.dates.append('date_end_quarantine')

        # Duration of different states: these are floats per Barn.
        self.durs = [
            'dur_contamination', # Duration of contamination
            'dur_quarantine', # Duration of quarantine
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
    
class PPE(Subroster):
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
        self.meta = PPEMeta() # Store list of keys and dtypes
        self.record_all_events = self.pars['record_all_events'] # Whether or not to record all events in the sim. If false, only transmision events are recorded. We set this to true by default since we have so few agents, but it can be set to false to save memory if desired.
        self.event_log = []
        self.infection_log = [] # Record of infections - keys for ['source','target','date','layer']

        pop_size = self.pars['pop_size_by_type']['ppe'] # Get the population size for this subroster

        # Set person properties -- all floats except for UID
        for key in self.meta.agent:
            if key in ['uid', 'human']:
                self[key] = np.zeros(pop_size, dtype=znd.default_int) # NOTE: The uid values are passed in kwargs by make_PPE()
            else:
                self[key] = np.full(pop_size, np.nan, dtype=znd.default_float)

        # Set states
        for key in self.meta.states:
            val = (key in ['uncontaminated']) # Default value is True for susceptible and naive, False otherwise
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

        self._pending_quarantine = defaultdict(list)  # Internal cache to record ppe that needs to be quarantined on each timestep {t:(person, quarantine_end_day)}

        return


    def init_flows(self):
        ''' Initialize flows to be zero '''
        self.flows = {key:0 for key in znd.new_ppe_flows}
        self.flows_variant = {}
        for key in znd.new_ppe_flows_by_variant:
            self.flows_variant[key] = np.zeros(self.pars['n_variants'], dtype=znd.default_float)

        return

    def initialize(self, agents_pars=None):
        ''' Perform initializations '''
        self.validate(roster_pars=agents_pars) # First, check that essential-to-match parameters match
        self.set_pars(agents_pars) # Replace the saved parameters with this simulation's
        self.set_rel_sus() # Set the relative susceptibility of each PPE based on the parameters
        self.set_rel_trans() # Set the relative transmissibility of each PPE based on the parameters
        self.initialized = True
        return

    def set_rel_sus(self):
        ''' Set the relative susceptibility of each PPE based on the parameters '''
        self.rel_sus = np.full(len(self), self.pars['prognoses']['ppe']['sus_ORs'], dtype=znd.default_float)
        return
    
    def set_rel_trans(self):
        ''' Set the relative transmissibility of each PPE based on the parameters '''
        self.rel_trans = np.full(len(self), self.pars['prognoses']['ppe']['trans_ORs'], dtype=znd.default_float)
        return

    def update_states_pre(self, t):
        ''' Perform all state updates at the current timestep '''

        self.t = t # Set the current time

        # Perform updates
        self.init_flows()
        self.flows['new_uncontaminated'] = self.check_uncontaminated()
        return


    def update_states_post(self):
        ''' Perform post-timestep updates '''

        self.flows['new_quarantined'] += self.check_quar()

        return

    def update_event_log(self, target_inds, event_type):
        '''
        Add an entry to the event log

        args:
            target_inds: array of indices of flocks that experienced a recordable event
            event (str): The specific event in question
        '''
        if self.record_all_events == False:
             return

        if target_inds is None:
            return
        
        for ind in target_inds:
            entry = dict(target = self.uid[ind], event_type = event_type, date = self.t)
            self.event_log.append(entry)

        return
    #%% Methods for updating state

    
    def check_uncontaminated(self):
        ''' Check which PPE get uncontaminated this timestep '''
        inds = self.check_inds(~self.contaminated, self.date_uncontaminated) # ~self.contaminated is the same as self.uncontaminated
        if len(inds) > 0:
            self.uncontaminated[inds] = True
            self.contaminated[inds] = False
            self.contaminated_variant[inds] = np.nan
            self.date_contaminated[inds] = np.nan
            self.date_uncontaminated[inds] = np.nan

            self.update_event_log(inds, 'uncontaminated')
        return len(inds)

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
                self.update_event_log(ind, 'quarantined')

        # If someone on quarantine has reached the end of their quarantine, release them
        end_inds = self.check_inds(~self.quarantined, self.date_end_quarantine, filter_inds=None) # Note the double-negative here (~)
        self.quarantined[end_inds] = False # Release from quarantine
        self.update_event_log(end_inds, 'released from quarantine')

        # Update the counter for consecutive days in quarantine
        self.cons_days_in_quar[self.quarantined] += 1
        self.cons_days_in_quar[~self.quarantined] = 0

        return n_quarantined
    



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
        contamination_pars = {k:self.pars['variant_pars']['wild']['ppe'][k] for k in variant_keys}
        variant_label = self.pars['variant_map'][variant]
        if variant:
            for k in variant_keys:
                contamination_pars[k] *= self.pars['variant_pars'][variant_label]['ppe'][k]

        n_infections = len(inds)
        durpars      = self.pars['dur']['ppe']


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
    
    def schedule_quarantine(self, uids, start_date=None, period=None):
        '''
        Schedule a quarantine. Typically not called by the user directly except
        via a custom intervention; see the contact_tracing() intervention instead.

        This function will create a request to quarantine a person on the start_date for
        a period of time. Whether they are on an existing quarantine that gets extended, or
        whether they are no longer eligible for quarantine, will be checked when the start_date
        is reached.

        Args:
            uids (int/list): UIDs of who to quarantine, specified by check_quar()
            start_date (int): day to begin quarantine (defaults to the current day, `sim.t`)
            period (int): quarantine duration (defaults to ``pars['dur']['ppe']['quar']``)
        '''

        

        start_date = self.t if start_date is None else int(start_date)
        period = self.pars['dur']['ppe']['quar'] if period is None else int(period)

        for uid in uids:
            ind = np.where(self.uid == uid)[0]
            self._pending_quarantine[start_date].append((ind, start_date + period))
        return
    
    def story(self, uid, *args):
        '''
        Print out a short history of events in the life of the specified PPE.

        Args:
            uid (int/list): the PPE whose story is to be regaled
            args (list): these PPE will tell their stories too

        **Example**::

            sim = zn.Sim()
            sim.run()
            sim.agents.ppe.story(12)
            sim.agents.ppe.story(795)
        '''

        def label_lkey(lkey):
            ''' Friendly name for common layer keys '''
            if lkey.lower() == 'hp':
                llabel = 'human-ppe contacts'
            elif lkey.lower() == 'pp':
                llabel = 'ppe-ppe contacts'
            elif lkey.lower() == 'pf':
                llabel = 'ppe-flock contacts'
            elif lkey.lower() == 'pb':
                llabel = 'ppe-barn contacts'
            elif lkey.lower() == 'pw':
                llabel = 'ppe-water contacts'
            else:
                llabel = f'"{lkey}"'
            return llabel

        uids = sc.promotetolist(uid)
        uids.extend(args)

        for uid in uids:
            uid_ind = np.where(self.uid == uid)[0]
            flock = self.human[uid_ind]


            intro = f'\nThis is the story of {uid}, a PPE assigned to human {flock}.'
            print(f'{intro}')

            events = []

            event_dict = {
                'contaminated': 'became contaminated',
                'uncontaminated': 'naturally lost its contamination',
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