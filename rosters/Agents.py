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
            'uid', #Int - unique ID for the agent
            'fid', #Int - farm ID, if applicable
            'agent_type', #string? the type of agent, must be one of pars['agent_types']
            #'rel_trans', # float - relative transmissibility of the agent
            #'rel_sus', # float - relative susceptibility of the agent
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

        # Set immunity and antibody states
        for key in self.meta.imm_states:  # Everyone starts out with no immunity
            self[key] = np.zeros((self.pars['n_variants'], self.pars['pop_size']), dtype=znd.default_float)

        # Store the dtypes used in a flat dict
        self._dtypes = {key:self[key].dtype for key in self.keys()} # Assign all to float by default
        if strict:
            self.lock() # If strict is true, stop further keys from being set (does not affect attributes)

        # Store flows to be computed during simulation
        #self.init_flows()

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

    #%% Methods for updating state

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
        self.human.initialize(self.pars) # Initialize the human subroster
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
        self.flock.update_states_pre(t)
        self.barn.update_states_pre(t)
        self.water.update_states_pre(t)

        self.check_result(t)
        self.check_cycle_end(t)
        
        self.check_repopulation(t)

        self.update_states_from_subrosters() # Update the states of the main roster

        return


    def update_states_post(self):
        ''' Perform post-timestep updates '''

        # Update the states of the subrosters
        self.human.update_states_post()
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
    
    def update_human_viral_loads(self, t):
        '''
        update the viral levels of human agents

        Args:
            t (float): Current time in simulation
        '''
        x_p1, y_p1 = self.human.x_p1, self.human.y_p1
        x_p2, y_p2 = self.human.x_p2, self.human.y_p2
        x_p3, y_p3 = self.human.x_p3, self.human.y_p3
        min_vl = znd.default_float(self.pars['transmission_pars']['human']['viral_levels']['min_vl'])
        max_vl = znd.default_float(self.pars['transmission_pars']['human']['viral_levels']['max_vl'])

        self.human.viral_load, human_viral_load = znu.compute_viral_load(t, x_p1, y_p1, x_p2, y_p2, x_p3, y_p3, min_vl, max_vl)
        return human_viral_load
    
    def update_flock_infection_levels(self, t):
        x_p1, y_p1 = self.flock.x_p1, self.flock.y_p1
        x_p2, y_p2 = self.flock.x_p2, self.flock.y_p2
        x_p3, y_p3 = self.flock.x_p3, self.flock.y_p3
        headcount = self.flock.headcount
        self.flock.infected_headcount, flock_infection_levels = znu.compute_infection_level(t, x_p1, y_p1, x_p2, y_p2, x_p3, y_p3, headcount)
        return flock_infection_levels

    def update_states_from_subrosters(self):
        susceptible_human_uids = np.array(self.human['uid'][znu.true(self.human['susceptible'])])
        exposed_human_uids = np.array(self.human['uid'][znu.true(self.human['exposed'])])
        infectious_human_uids = np.array(self.human['uid'][znu.true(self.human['infectious'])])
        symptomatic_human_uids = np.array(self.human['uid'][znu.true(self.human['symptomatic'])])
        quarantined_human_uids = np.array(self.human['uid'][znu.true(self.human['quarantined'])])

        susceptible_flock_uids = np.array(self.flock['uid'][znu.true(self.flock['susceptible'])])
        exposed_flock_uids = np.array(self.flock['uid'][znu.true(self.flock['exposed'])])
        infectious_flock_uids = np.array(self.flock['uid'][znu.true(self.flock['infectious'])])
        #symptomatic_flock_uids = np.array(self.flock['uid'][znu.true(self.flock['symptomatic'])])
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
        #symptomatic_uids = np.concatenate((symptomatic_human_uids, symptomatic_flock_uids))
        symptomatic_uids = symptomatic_human_uids
        quarantined_uids = np.concatenate((quarantined_human_uids, quarantined_flock_uids))

        self.susceptible = np.isin(self['uid'], susceptible_uids)
        self.exposed = np.isin(self['uid'], exposed_uids)
        self.infectious = np.isin(self['uid'], infectious_uids)
        self.symptomatic = np.isin(self['uid'], symptomatic_uids)
        self.quarantined = np.isin(self['uid'], quarantined_uids)
        self.infectious_variant = np.concatenate((self.human.infectious_variant,
                                                  self.flock.infectious_variant,
                                                  self.barn.contaminated_variant,
                                                  self.water.contaminated_variant))
        # self.infectious_by_variant = np.concatenate((self.human.infectious_by_variant,
        #                                              self.flock.infectious_by_variant,
        #                                              self.barn.contaminated_by_variant,
        #                                              self.water.contaminated_by_variant))
        return




    #%% Methods for calculating values from subrosters

    def get_human_rel_sus(self):
        return self.human.rel_sus
    
    def get_flock_rel_sus(self):
        return self.flock.rel_sus
    
    def get_barn_rel_sus(self):
        #TODO: This should depend on the temperature and humidity of the barn
        return np.repeat([1.0], len(self.barn))
    
    def get_water_rel_sus(self):
        #TODO: This should depend on the temperature of the water
        return np.repeat([1.0], len(self.water))
    
    def get_rel_sus(self):
        human_rel_sus = self.get_human_rel_sus()
        flock_rel_sus = self.get_flock_rel_sus()
        barn_rel_sus = self.get_barn_rel_sus()
        water_rel_sus = self.get_water_rel_sus()
        return np.concatenate((human_rel_sus, flock_rel_sus, barn_rel_sus, water_rel_sus))
    
    def get_human_rel_trans(self):
        return self.human.rel_trans
    
    def get_flock_rel_trans(self):
        return self.flock.rel_trans
    
    def get_barn_rel_trans(self):
        #TODO: This should depend on the temperature and humidity of the barn
        return np.repeat([1.0], len(self.barn))
    
    def get_water_rel_trans(self):
        #TODO: This should depend on the temperature of the water
        return np.repeat([1.0], len(self.water))
    
    def get_rel_trans(self):
        human_rel_trans = self.get_human_rel_trans()
        flock_rel_trans = self.get_flock_rel_trans()
        barn_rel_trans = self.get_barn_rel_trans()
        water_rel_trans = self.get_water_rel_trans()
        return np.concatenate((human_rel_trans, flock_rel_trans, barn_rel_trans, water_rel_trans))


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

        #TODO: Refactor sources so they actually work

        human_inds = np.where(np.isin(self.human.uid, self.uid[inds]))
        #human_inds = np.array([i for i, uid in enumerate(self.human.uid) if uid in set(self.uid[inds])]) # Supposedly faster if self.uid[inds] is large

        flock_inds = np.where(np.isin(self.flock.uid, self.uid[inds])) 
        #flock_inds = np.array([i fo i,uid in enumerate(self.flock.uid) if uid in set(self.uid[inds])])

        barn_inds = np.where(np.isin(self.barn.uid, self.uid[inds]))
        #barn_inds = np.array([i fo i,uid in enumerate(self.barn.uid) if uid in set(self.uid[inds])])

        water_inds = np.where(np.isin(self.water.uid, self.uid[inds]))
        #barn_inds = np.array([i fo i,uid in enumerate(self.barn.uid) if uid in set(self.uid[inds])])

        self.human.infect(human_inds, hosp_max, source, layer, variant)
        self.flock.infect(flock_inds, source, layer, variant)
        self.barn.infect(barn_inds, source, layer, variant)
        self.water.infect(water_inds, source, layer, variant)

        self.update_states_from_subrosters()
        return

    #%% Methods that require access to multiple subrosters
    def check_repopulation(self, t):
        '''
        Check for farms that are scheduled to be repopulated and reincarnate the resident flock with proper initial conditions.
        '''
        prod_pars = self.pars['production_cycle'] 
        progs = self.pars['prognoses']['flock']
        barn_inds = np.where(self.barn.date_repopulate <= t)[0]
        if len(barn_inds) > 0:
            self.barn.repopulations[barn_inds]+= 1
            self.barn.date_repopulate[barn_inds] = np.nan
            flock_inds = np.where(np.isin(self.flock.uid, self.barn.flock[barn_inds]))[0]
            breed_to_index = {breed: index for index, breed in enumerate(prod_pars['breeds'])}
            breed_inds = np.array([breed_to_index[this_breed] for this_breed in self.flock.breed[flock_inds]])

            breed, freq = np.unique(breed_inds, return_counts=True)
            breed_dict = dict(zip(breed, freq))
            for breed, freq in breed_dict.items():
                current_breed_inds = np.where(breed_inds == breed)[0]
                self.barn.date_cycle_end[barn_inds[current_breed_inds]] = znu.sample(**prod_pars['cycle_dur'][breed], size = freq)
                self.flock.headcount[flock_inds[current_breed_inds]] = znu.sample(**prod_pars['flock_size'][breed], size = freq)
            
            self.flock.susceptible[flock_inds] = True
            self.flock.suspected[flock_inds] = False
            self.flock.baseline_symptomatic_rate[flock_inds] = progs['baseline_symptomatic_rate'][breed_inds]
            self.flock.baseline_mortality_rate[flock_inds] = progs['baseline_mortality_rate'][breed_inds]
            self.flock.baseline_water_rate[flock_inds] = progs['baseline_water_rate'][breed_inds]
            self.flock.rel_sus[flock_inds] = progs['sus_ORs'][breed_inds]
            self.flock.rel_trans[flock_inds] = progs['trans_ORs'][breed_inds]

        return len(barn_inds)
    

    
    def check_result(self, t):
        '''
        Check for flocks that have recieved their test results and update the states accordingly.
        '''

        flock_inds = np.where(self.flock.date_result <= t)[0]
        barn_inds = np.where(np.isin(self.barn.flock, self.flock.uid[flock_inds]))[0]


        pos_flock_inds = np.where(self.flock.infectious[flock_inds] == True)[0]
        pos_barn_inds = np.where(np.isin(self.barn.flock, self.flock.uid[pos_flock_inds]))[0]


        neg_flock_inds = np.where(self.flock.infectious[flock_inds] == False)[0]
        #neg_barn_inds = np.isin(self.barn.flock, self.flock.uid[neg_flock_inds])

        self.barn.date_cycle_end[pos_barn_inds] = np.nan
        self.barn.composting[pos_barn_inds] = True
        self.barn.date_composting[pos_barn_inds] = t + znu.sample(**self.pars['dur']['barn']['composting'], size = len(pos_barn_inds))

        self.flock.headcount[pos_flock_inds] = 0
        self.flock.infected_headcount[pos_flock_inds] = 0
        self.flock.symptomatic_headcount[pos_flock_inds] = 0
        self.flock.dead_headcount[pos_flock_inds] = 0
        self.flock.water_consumption[pos_flock_inds] = 0

        self.flock.susceptible[pos_flock_inds] = False
        self.flock.exposed[pos_flock_inds] = False
        self.flock.infectious[pos_flock_inds] = False
        self.flock.quarantined[pos_flock_inds] = False
        #self.flock.suspected[pos_flock_inds] = False
        self.flock.date_infectious[pos_flock_inds] = np.nan
        self.flock.date_exposed[pos_flock_inds] = np.nan
        self.flock.date_suspected[pos_flock_inds] = np.nan

        self.flock.suspected[neg_flock_inds] = False
        self.flock.date_suspected[neg_flock_inds] = np.nan
        self.flock.date_result[neg_flock_inds] = np.nan

        return len(barn_inds)
    
    def check_cycle_end(self, t):
        '''
        Check for flocks that are at the end of their production cycle.
        '''

        barn_inds = np.where(self.barn.date_cycle_end <= t)[0]
        flock_inds = np.isin(self.flock.uid, self.barn.flock[barn_inds])

        self.barn.date_cycle_end[barn_inds] = np.nan
        self.barn.cleaning[barn_inds] = True
        self.barn.date_cleaning[barn_inds] = t + znu.sample(**self.pars['dur']['barn']['cleaning'], size = len(barn_inds))

        self.flock.headcount[flock_inds] = 0
        self.flock.infected_headcount[flock_inds] = 0
        self.flock.symptomatic_headcount[flock_inds] = 0
        self.flock.dead_headcount[flock_inds] = 0
        self.flock.water_consumption[flock_inds] = 0

        self.flock.susceptible[flock_inds] = False
        self.flock.exposed[flock_inds] = False
        self.flock.infectious[flock_inds] = False
        self.flock.suspected[flock_inds] = False
        self.flock.quarantined[flock_inds] = False
        self.flock.date_infectious[flock_inds] = np.nan
        self.flock.date_exposed[flock_inds] = np.nan
        self.flock.date_suspected[flock_inds] = np.nan
        

        return len(barn_inds)
   