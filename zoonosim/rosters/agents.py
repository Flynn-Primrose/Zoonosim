'''
Defines the roster of all agents (regardless of type) that are included in the simulation.
'''

import numpy as np
import sciris as sc
from .. import version as znv
from .. import utils as znu
from .. import defaults as znd

from .contacts import Contacts
from .rosters import Roster
from .humans import Humans
from .ppe import PPE
from .flocks import Flocks
from .barns import Barns
from .waters import Water

__all__ = ['Agents']

class AgentsMeta(sc.prettyobj):
    '''
    For storing all keys that are common across all agent types
    '''

    def __init__(self):

        self.agent = [
            'uid', #Int - unique ID for the agent
            'fid', #Int - farm ID, if applicable
            'agent_type', #string? the type of agent, must be one of pars['agent_types']
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

        self.imm_states = [
            'sus_imm',  # Float, by variant
            'symp_imm', # Float, by variant
            'sev_imm',  # Float, by variant
        ]


        self.all_states = self.agent + self.states +self.variant_states + self.by_variant_states #+ self.imm_states

        # Validate
        #self.state_types = ['agent', 'states', 'variant_states', 'by_variant_states', 'imm_states', 'all_states']
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

        human = kwargs.pop('human', None)
        if isinstance(human, Humans):
            self.human = human
        else:
            raise ValueError("human must be an instance of Humans class")
        
        ppe = kwargs.pop('ppe', None)
        if isinstance(ppe, PPE):
            self.ppe = ppe
        else:
            raise ValueError("ppe must be an instance of the PPE class")
            
        flock = kwargs.pop('flock', None)
        if isinstance(flock, Flocks):
            self.flock = flock
        else:
            raise ValueError("flock must be an instance of Flocks class")
        
        barn = kwargs.pop('barn', None)
        if isinstance(barn, Barns):
            self.barn = barn
        else:
            raise ValueError("barn must be an instance of Barns class")
        
        water = kwargs.pop('water', None)
        if isinstance(water, Water):
            self.water = water
        else:
            raise ValueError("water must be an instance of Water class")
       
        self.contacts = None
        self.init_contacts() # Initialize the contacts

        if 'contacts' in kwargs:
            self.add_contacts(kwargs.pop('contacts'))

        self.infection_log = [] # Record of infections - keys for ['source','target','date','layer']
        
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

        # Set immunity and antibody states
        for key in self.meta.imm_states:  # Everyone starts out with no immunity
            self[key] = np.zeros((self.pars['n_variants'], self.pars['pop_size']), dtype=znd.default_float)

        # Store the dtypes used in a flat dict
        self._dtypes = {key:self[key].dtype for key in self.keys()} # Assign all to float by default
        if strict:
            self.lock() # If strict is true, stop further keys from being set (does not affect attributes)

        # Although we have called init(), we still need to call initialize()
        self.initialized = False



        # Handle all other values, e.g. age
        for key,value in kwargs.items():
            if strict:
                self.set(key, value)
            else:
                    self[key] = value

        return

    #%% Methods for updating state

    def initialize(self, sim_pars=None):
        ''' Perform initializations '''
        self.validate(sim_pars=sim_pars) # First, check that essential-to-match parameters match
        self.set_pars(sim_pars) # Replace the saved parameters with this simulation's
        self.human.initialize(self.pars) # Initialize the human subroster
        self.ppe.initialize(self.pars) # Initialize the ppe subroster
        self.flock.initialize(self.pars) # Initialize the flock subroster
        self.barn.initialize(self.pars) # Initialize the barn subroster
        self.water.initialize(self.pars) # Initialize the water subroster
        self.update_states_from_subrosters() # Update the states of the main roster based on the subrosters
        self.initialized = True
        return
    

    def update_states_pre(self, t):
        ''' Perform all state updates at the current timestep '''

        # Initialize
        self.t = t
        self.is_exp = self.true('exposed') # For storing the interim values since used in every subsequent calculation
        self.human.update_states_pre(t)
        self.ppe.update_states_pre(t)
        self.flock.update_states_pre(t)
        self.barn.update_states_pre(t)
        self.water.update_states_pre(t)

        # self.check_result(t)

        # self.check_cycle_end(t)
        
        #self.check_repopulation(t)

        self.update_states_from_subrosters() # Update the states of the main roster

        return


    def update_states_post(self):
        ''' Perform post-timestep updates '''

        # Update the states of the subrosters
        self.human.update_states_post()
        self.ppe.update_states_post()
        self.flock.update_states_post()
        self.barn.update_states_post()
        self.water.update_states_post()
        # Update the states of the main roster
        self.update_states_from_subrosters()
        del self.is_exp  # Tidy up

        return

    def update_contacts(self):
        ''' Refresh dynamic contacts, e.g. community '''
        # Figure out if anything needs to be done -- e.g. {'h':False, 'c':True}
        for lkey, is_dynam in self.pars['dynam_layer'].items():
            if is_dynam:
                self.contacts[lkey].update(self)

        return self.contacts
    
    def update_human_viral_loads(self):
        '''
        update the viral levels of human agents
        This function just calls the relevent function in the human subroster, but is included here for convenience.

        Args:
            t (float): Current time in simulation
        '''

        return self.human.update_viral_loads()
    
    def update_flock_infection_levels(self):
        '''
        Update the infection levels of the flock subroster. 
        This function just calls the relevant function in the flock subroster, but is included here for convenience.
        NOTE: This is not optimized for speed, so it may be slow for large populations.
        '''

        return self.flock.update_infection_levels()
    
    def update_states_from_subrosters(self):
        susceptible_human_uids = np.array(self.human['uid'][znu.true(self.human['susceptible'])])
        exposed_human_uids = np.array(self.human['uid'][znu.true(self.human['exposed'])])
        infectious_human_uids = np.array(self.human['uid'][znu.true(self.human['infectious'])])
        symptomatic_human_uids = np.array(self.human['uid'][znu.true(self.human['symptomatic'])])
        quarantined_human_uids = np.array(self.human['uid'][znu.true(self.human['quarantined'])])

        susceptible_ppe_uids = np.array(self.ppe['uid'][znu.false(self.ppe['contaminated'])])
        exposed_ppe_uids = np.array(self.ppe['uid'][znu.true(self.ppe['contaminated'])])
        infectious_ppe_uids = np.array(self.ppe['uid'][znu.true(self.ppe['contaminated'])])
        quarantined_ppe_uids = np.array(self.ppe['uid'][znu.true(self.ppe['quarantined'])])

        susceptible_flock_uids = np.array(self.flock['uid'][znu.true(self.flock['susceptible'])])
        exposed_flock_uids = np.array(self.flock['uid'][znu.true(self.flock['exposed'])])
        infectious_flock_uids = np.array(self.flock['uid'][znu.true(self.flock['infectious'])])
        quarantined_flock_uids = np.array(self.flock['uid'][znu.true(self.flock['quarantined'])])

        susceptible_barn_uids = np.array(self.barn['uid'][znu.false(self.barn['contaminated'])])
        exposed_barn_uids = np.array(self.barn['uid'][znu.true(self.barn['contaminated'])])
        infectious_barn_uids = np.array(self.barn['uid'][znu.true(self.barn['contaminated'])])

        susceptible_water_uids = np.array(self.water['uid'][znu.false(self.water['contaminated'])])
        exposed_water_uids = np.array(self.water['uid'][znu.true(self.water['contaminated'])])
        infectious_water_uids = np.array(self.water['uid'][znu.true(self.water['contaminated'])])

        susceptible_uids = np.concatenate((susceptible_human_uids, susceptible_ppe_uids, susceptible_flock_uids, susceptible_barn_uids, susceptible_water_uids))
        exposed_uids = np.concatenate((exposed_human_uids, exposed_ppe_uids, exposed_flock_uids, exposed_barn_uids, exposed_water_uids))
        infectious_uids = np.concatenate((infectious_human_uids, infectious_ppe_uids, infectious_flock_uids, infectious_barn_uids, infectious_water_uids))

        symptomatic_uids = symptomatic_human_uids
        quarantined_uids = np.concatenate((quarantined_human_uids, quarantined_ppe_uids, quarantined_flock_uids))

        self.susceptible = np.isin(self['uid'], susceptible_uids)
        self.exposed = np.isin(self['uid'], exposed_uids)
        self.infectious = np.isin(self['uid'], infectious_uids)
        self.symptomatic = np.isin(self['uid'], symptomatic_uids)
        self.quarantined = np.isin(self['uid'], quarantined_uids)
        self.infectious_variant = np.concatenate((self.human.infectious_variant,
                                                  self.ppe.contaminated_variant,
                                                  self.flock.infectious_variant,
                                                  self.barn.contaminated_variant,
                                                  self.water.contaminated_variant))
        return




    #%% Methods for calculating values from subrosters

    def get_human_rel_sus(self):
        return self.human.rel_sus
    
    def get_ppe_rel_sus(self):
        
        return self.ppe.rel_sus
    
    def get_flock_rel_sus(self):
        return self.flock.rel_sus
    
    
    def get_barn_rel_sus(self):

        return self.barn.rel_sus
    
    def get_water_rel_sus(self):
        return self.water.rel_sus
    
    def get_rel_sus(self):
        human_rel_sus = self.get_human_rel_sus()
        ppe_rel_sus = self.get_ppe_rel_sus()
        flock_rel_sus = self.get_flock_rel_sus()
        barn_rel_sus = self.get_barn_rel_sus()
        water_rel_sus = self.get_water_rel_sus()
        return np.concatenate((human_rel_sus, ppe_rel_sus, flock_rel_sus, barn_rel_sus, water_rel_sus))
    
    def get_human_rel_trans(self):
        return self.human.rel_trans
    
    def get_ppe_rel_trans(self):
        return self.ppe.rel_trans 
    
    def get_flock_rel_trans(self):
        return self.flock.rel_trans
    
    def get_barn_rel_trans(self):
        return self.barn.rel_trans
    
    def get_water_rel_trans(self):
        return self.water.rel_trans
    
    def get_rel_trans(self):
        human_rel_trans = self.get_human_rel_trans()
        ppe_rel_trans = self.get_ppe_rel_trans()
        flock_rel_trans = self.get_flock_rel_trans()
        barn_rel_trans = self.get_barn_rel_trans()
        water_rel_trans = self.get_water_rel_trans()
        return np.concatenate((human_rel_trans, ppe_rel_trans, flock_rel_trans, barn_rel_trans, water_rel_trans))


    #%% Methods to make events occur (infection and diagnosis)

    def infect_type(self, agent_type, inds, update = True):
        '''
        Add infections to a specific subroster. Used by init_infections.

        Args:
            agent_type (str): Which subroster to add the infections to
            inds (array): Indices of the agents to be infected. NOTE: this refers to the indices in the subroster not the roster
            update (bool): True means update the roster based on the subrosters 
        '''
        self[agent_type].infect(inds = inds, layer = 'seed_infections')
        if update: self.update_states_from_subrosters()
        return
    
    def infect(self, inds, hosp_max, source, layer, variant):

        human_inds = np.where(np.isin(self.human.uid, self.uid[inds]))
        #human_inds = np.array([i for i, uid in enumerate(self.human.uid) if uid in set(self.uid[inds])]) # Supposedly faster if self.uid[inds] is large

        ppe_inds = np.where(np.isin(self.ppe.uid, self.uid[inds]))

        flock_inds = np.where(np.isin(self.flock.uid, self.uid[inds])) 
        #flock_inds = np.array([i fo i,uid in enumerate(self.flock.uid) if uid in set(self.uid[inds])])

        barn_inds = np.where(np.isin(self.barn.uid, self.uid[inds]))
        #barn_inds = np.array([i fo i,uid in enumerate(self.barn.uid) if uid in set(self.uid[inds])])

        water_inds = np.where(np.isin(self.water.uid, self.uid[inds]))
        #barn_inds = np.array([i fo i,uid in enumerate(self.barn.uid) if uid in set(self.uid[inds])])

        self.human.infect(human_inds, hosp_max, source, layer, variant)
        self.ppe.infect(ppe_inds, source, layer, variant)
        self.flock.infect(flock_inds, source, layer, variant)
        self.barn.infect(barn_inds, source, layer, variant)
        self.water.infect(water_inds, source, layer, variant)

        self.update_states_from_subrosters()
        return
