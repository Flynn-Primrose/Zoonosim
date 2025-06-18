import numpy as np
import sciris as sc
from .. import version as znv
from .. import utils as znu
from .. import defaults as znd
from .Subroster import Subroster

__all__ = ['Water']

class WaterMeta(sc.prettyobj):
    ''' Defines all keys that are used by barns '''

    def __init__(self):
        
        self.agent = [
            'uid', # int
            'temperature'
        ]

        self.states = [
            'uncontaminated',
            'contaminated'
        ]

        self.variant_states = [
            'contaminated_variant',
        ]

        self.by_variant_states = [
            'contaminated_by_variant',
        ]

        # Set the dates various events took place: these are floats per agent
        self.dates = [f'date_{state}' for state in self.states] # Convert each state into a date

        # Duration of different states: these are floats per Water.
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
    

class Water(Subroster):
    '''
    A class to perform all the operations on the agents -- usually not invoked directly.

    '''

    def __init__(self, pars, strict=True, **kwargs):

        # Handle pars and population size
        self.set_pars(pars)
        self.version = znv.__version__ # Store version info

        # Other initialization
        self.t = 0 # Keep current simulation time
        self._lock = False # Prevent further modification of keys
        self.meta = WaterMeta() # Store list of keys and dtypes
        # self.init_contacts() # Initialize the contacts
        self.infection_log = [] # Record of infections - keys for ['source','target','date','layer']

        pop_size = self.pars['pop_size_by_type']['water']
        
        # Set person properties -- all floats except for UID
        for key in self.meta.agent:
            if key == 'uid':
                self[key] = np.zeros(pop_size, dtype=znd.default_int) # NOTE: The actual uids will be passed in kwargs by make_water()
            else:
                self[key] = np.full(pop_size, np.nan, dtype=znd.default_float)


        # Set health states -- only susceptible is true by default -- booleans except exposed by variant which should return the variant that ind is exposed to
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

        #self._pending_quarantine = defaultdict(list)  # Internal cache to record people that need to be quarantined on each timestep {t:(person, quarantine_end_day)}

        return


    def init_flows(self):
        ''' Initialize flows to be zero '''
        self.flows = {key:0 for key in znd.new_water_flows}
        self.flows_variant = {}
        for key in znd.new_water_flows_by_variant:
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

        # Initialize
        self.t = t

        # Perform updates
        self.init_flows()
        self.flows['new_uncontaminated'] = self.check_uncontaminated()
        return


    def update_states_post(self):
        ''' Perform post-timestep updates '''

        return



    #%% Methods for updating state


    def check_uncontaminated(self):
        ''' Check if waterbodies are contaminated and update their states accordingly '''
        inds = self.check_inds(~self.contaminated, self.date_uncontaminated)
        if len(inds) > 0:
            self.uncontaminated[inds] = True
            self.contaminated[inds]    = False

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
        keep = self.uncontaminated[inds]
        inds = inds[keep]
        if source is not None:
            source = source[keep]



        # Deal with variant parameters
        variant_keys = ['rel_dur_contamination']
        contamination_pars = {k:self.pars[k] for k in variant_keys}
        variant_label = self.pars['variant_map'][variant]
        if variant:
            for k in variant_keys:
                contamination_pars[k] *= self.pars['variant_pars']['water'][variant_label][k]

        n_infections = len(inds)
        durpars      = self.pars['dur']['water']


        # Update states, variant info, and flows
        self.uncontaminated[inds]    = False
        self.contaminated[inds]      = True
        self.contaminated_variant[inds] = variant
        self.contaminated_by_variant[variant, inds] = True
        self.flows['new_contaminated']   += len(inds)
        self.flows_variant['new_contaminated_by_variant'][variant] += len(inds)

        # Record transmissions
        for i, target in enumerate(inds):
            entry = dict(source=source[i] if source is not None else None, target=target, date=self.t, layer=layer, variant=variant_label)
            self.infection_log.append(entry)

        # Calculate how long before this person can infect other people
        self.dur_contamination[inds] = znu.sample(**durpars['contamination'], size=n_infections) * contamination_pars['rel_dur_contamination']
        self.date_contaminated[inds] = self.t
        self.date_uncontaminated[inds] = self.dur_contamination[inds] + self.t

        # Reset all other dates


        return n_infections # For incrementing counters
