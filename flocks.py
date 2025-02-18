'''
Defines the Poultry class and functions associated with making Flocks and handling
the transitions between states (e.g., from susceptible to infected).
'''

#%% Imports
import numpy as np

from . import version as cvv
from . import utils as cvu
from . import config as znconfig
from . import base as cvb




__all__ = ['Flocks']

class Flocks(cvb.BaseRoster):
    '''
    A class to perform all the operations on Chicken flock agents -- usually not invoked directly.

    This class is usually created automatically by the sim. The only required input
    argument is the population size, but typically the full parameters dictionary
    will get passed instead since it will be needed before the Flock object is
    initialized. However, ages, contacts, etc. will need to be created separately --
    see ``zn.make_flock()`` instead.

    Note that this class handles the mechanics of updating the actual people, while
    ``zn.BaseRoster`` takes care of housekeeping (saving, loading, exporting, etc.).
    Please see the BaseRoster class for additional methods.

    Args:
        pars (dict): the sim parameters, e.g. sim.pars -- alternatively, if a number, interpreted as pop_size
        strict (bool): whether or not to only create keys that are already in self.meta.flock; otherwise, let any key be set
        kwargs (dict): the actual data, e.g. from a popdict, being specified

    **Examples**::
        flk1 = zn.Flock(2000)

        sim = zn.Sim()
        flk2 = zn.Flock(sim.pars)
    '''

    def __init__(self, pars, strict=True, **kwargs):

        # Handle pars and population size
        self.set_pars(pars)
        self.version = cvv.__version__ # Store version info

        # Other initialization
        self.t = 0 # Keep current simulation time
        self._lock = False # Prevent further modification of keys
        self.meta = znconfig.FlocksMeta() # Store list of keys and dtypes
        self.contacts = None
        self.init_contacts() # Initialize the contacts
        self.infection_log = [] # Record of infections - keys for ['source','target','date','layer']
        
        # Set flock properties -- all floats except for UID
        for key in self.meta.flock:
            if key == 'uid':
                self[key] = np.arange(self.pars['pop_size'], dtype=znconfig.default_int)
            else:
                self[key] = np.full(self.pars['pop_size'], np.nan, dtype=znconfig.default_float)

        # Set health states -- only susceptible is true by default -- booleans except exposed by variant which should return the variant that ind is exposed to
        for key in self.meta.states:
            val = (key in ['susceptible']) # Default value is True for susceptible, False otherwise
            self[key] = np.full(self.pars['pop_size'], val, dtype=bool)

        # Set variant states, which store info about which variant a person is exposed to
        for key in self.meta.variant_states:
            self[key] = np.full(self.pars['pop_size'], np.nan, dtype=znconfig.default_float)
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
            if strict:
                self.set(key, value)
            else:
                self[key] = value

        return


    def init_flows(self):
        ''' Initialize flows to be zero '''
        self.flows = {key:0 for key in znconfig.new_result_flows}
        self.flows_variant = {}
        for key in znconfig.new_result_flows_by_variant:
            self.flows_variant[key] = np.zeros(self.pars['n_variants'], dtype=znconfig.default_float)
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

        # Perform updates
        self.init_flows()
        self.flows['new_infectious']    += self.check_infectious() # For flocks are exposed and not infectious, check if they begin being infectious
        return


    def update_states_post(self):
        ''' Perform post-timestep updates '''

        return



    def update_contacts(self):
        ''' Refresh dynamic contacts, e.g. community '''
        # Figure out if anything needs to be done -- e.g. {'h':False, 'c':True}
        for lkey, is_dynam in self.pars['dynam_layer'].items():
            if is_dynam:
                self.contacts[lkey].update(self)

        return self.contacts



    #%% Methods for updating state

    def check_inds(self, current, date, filter_inds=None):
        ''' Return indices for which the current state is false and which are assigned a date on or before the current date

        Args:
            current (array): list of boolean values that represent a current state
            date (array): list that contains either a date or a Nan
        '''
        if filter_inds is None:
            not_current = cvu.false(current)
        else:
            not_current = cvu.ifalsei(current, filter_inds)
        has_date = cvu.idefinedi(date, not_current)
        inds     = cvu.itrue(self.t >= date[has_date], has_date)
        return inds


    def check_infectious(self):
        ''' Check if they become infectious '''
        inds = self.check_inds(self.infectious, self.date_infectious, filter_inds=self.is_exp)
        self.infectious[inds] = True
        self.infectious_variant[inds] = self.exposed_variant[inds]

        for variant in range(self.pars['n_variants']):
            this_variant_inds = cvu.itrue(self.infectious_variant[inds] == variant, inds)
            n_this_variant_inds = len(this_variant_inds)
            self.flows_variant['new_infectious_by_variant'][variant] += n_this_variant_inds
            self.infectious_by_variant[variant, this_variant_inds] = True
        return len(inds)
    
    def check_symptomatic(self):
        ''' Check if they become symptomatic '''
        inds = self.check_inds(self.symptomatic, self.date_symptomatic, filter_inds=self.is_infectious)
        self.symptomatic[inds] = True
        for variant in range(self.pars['n_variants']):
            this_variant_inds = cvu.itrue(self.infectious_variant[inds] == variant, inds)
            n_this_variant_inds = len(this_variant_inds)
            self.flows_variant['new_symptomatic_by_variant'][variant] += n_this_variant_inds
            self.symptomatic_by_variant[variant, this_variant_inds] = True
        return len(inds)
    
    def check_culled(self):
        ''' Check if they become culled '''
        inds = self.check_inds(self.culled, self.date_culled, filter_inds=self.is_symptomatic)
        self.culled[inds] = True
        self.infectious[inds] = False
        self.symptomatic[inds] = False
        for variant in range(self.pars['n_variants']):
            this_variant_inds = cvu.itrue(self.infectious_variant[inds] == variant, inds)
            n_this_variant_inds = len(this_variant_inds)
            self.flows_variant['new_culled_by_variant'][variant] += n_this_variant_inds
            self.culled_by_variant[variant, this_variant_inds] = True
        return len(inds)

    #%% Methods to make events occur (infection)

    def infect(self, inds, source, layer=None, variant = 0):
        '''
        Infect the flock
        '''
        # TODO: figure out the best way to handle this