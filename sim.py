'''
Defines the Sim class, ZoonoSim's core class.
'''

#%% Imports
from pickle import NONE
import numpy as np
import pandas as pd
import sciris as sc
from . import utils as znu
from . import misc as znm
from . import base as znb
from . import population as znpop
from . import parameters as znpar
from . import immunity as znimm
from . import people as znppl
from . import pathogens as pat
from. import pathogen_interactions as p_int

from .config import options as zno
from .config import DataTypes as zndt

# Almost everything in this file is contained in the Sim class
__all__ = ['Sim', 'diff_sims', 'demo', 'AlreadyRunError']


class Sim(znb.BaseSim):
    '''
    The Sim class handles the running of the simulation: the creation of the
    population and the dynamics of the epidemic. This class handles the mechanics
    of the actual simulation, while BaseSim takes care of housekeeping (saving,
    loading, exporting, etc.). Please see the BaseSim class for additional methods.

    Args:
        pars     (dict):   parameters to modify from their default values
        datafile (str/df): filename of (Excel, CSV) data file to load, or a pandas dataframe of the data
        datacols (list):   list of column names of the data to load
        label    (str):    the name of the simulation (useful to distinguish in batch runs)
        simfile  (str):    the filename for this simulation, if it's saved
        popfile  (str):    if supplied, load the population from this file (saved People object from infection module, not behaviour module)
        people   (varies): if supplied, use these pre-generated people (either People object, or BehaviourModel object (synthpop population output)) instead of loading or generating new ones
        version  (str):    if supplied, use default parameters from this version of Zoonosim instead of the latest
        kwargs   (dict):   additional parameters; passed to ``cv.make_pars()``

    **Examples**::

        sim = zn.Sim()
        sim = zn.Sim(pop_size=10e3, datafile='my_data.xlsx', label='Sim with data')
    '''

    def __init__(self, pars=None, datafile=None, label=None, simfile=None,
                 popfile=None, people=None, version=None, **kwargs):

        # Set attributes
        self.label         = label    # The label/name of the simulation
        self.created       = None     # The datetime the sim was created
        self.simfile       = simfile  # The filename of the sim
        self.datafile      = datafile # The name of the data file
        self.popfile       = popfile  # The population file
        self.data          = None     # The actual data
        self.popdict       = people   # The population dictionary
        self.people        = None     # TODO: how to we modify this to include multiple species?
        self.t             = None     # The current time in the simulation (during execution); outside of sim.step(), its value corresponds to next timestep to be computed
        self.results       = {}     # For storing results
        self.summary       = None     # For storing a summary of the results
        self.initialized   = False    # Whether or not initialization is complete
        self.complete      = False    # Whether a simulation has completed running
        self.results_ready = False    # Whether or not results are ready
        self._default_ver  = version  # Default version of parameters used
        self._orig_pars    = None     # Store original parameters to optionally restore at the end of the simulation 
        self.initialized_pathogens = False

        # Make default parameters (using values from parameters.py)
        default_pars = znpar.make_pars(version=version) # Start with default pars
        super().__init__(default_pars) # Initialize and set the parameters as attributes
             
             
        # Now update everything
        self.set_metadata(simfile) # Set the simulation date and filename
        self.update_pars(pars, **kwargs) # Update the parameters, if provided
        self.load_data(datafile) # Load the data, if provided
        

        #TODO: this has to change to support multiple species. 
        if not isinstance(self.pars['pathogens'], list):
            self.pars['pathogens'] = [self.pars['pathogens']]
        self.pathogens = self.pars['pathogens']
        self.n_pathogens  = len(self.pathogens)



    def load_data(self, datafile=None, verbose=None, **kwargs):
        ''' Load the data to calibrate against, if provided '''
        if verbose is None:
            verbose = self['verbose']
        self.datafile = datafile # Store this
        if datafile is not None: # If a data file is provided, load it
            self.data = znm.load_data(datafile=datafile, verbose=verbose, start_day=self['start_day'], **kwargs)

        return


    def initialize(self, reset=False, init_infections=True, **kwargs):
        '''
        Perform all initializations on the sim.

        This includes validating the parameters, setting the random number seed,
        creating the results structure, initializing the people, validating the
        layer parameters (which requires the people), and initializing the interventions.

        Note: to create a population to save for later use, use ``init_infections=False``.
        This will then create a fresh People object which other sims can finish
        initializing.

        Args:
            reset (bool): whether or not to reset people even if they already exist
            init_infections (bool): whether to initialize infections (default true so sim is ready, but can't reuse people then)
            kwargs (dict): passed to init_people
        '''
        self.t = 0  # The current time index
        self.validate_pars() # Ensure parameters have valid values
        self.set_seed() # Reset the random seed before the population is created
        self.init_pathogens() #Initialize pathogen information
        self.init_pathogen_interactions() #Validate pathogen-pathogen matrices
        self.init_immunity() # initialize information about immunity (if use_waning=True)
        self.init_results() # After initializing the variant, create the results structure
        self.init_people(reset=reset, init_infections=init_infections, **kwargs) # TODO: this has to change to support multiple species.
        self.init_infections()   
        self.validate_layer_pars() # Once the population is initialized, validate the layer parameters again
        self.set_seed() # Reset the random seed again so the random number stream is consistent
         

        self.initialized   = True
        self.complete      = False
        self.results_ready = False 
        return self

    def layer_keys(self):
        '''
        Attempt to retrieve the current layer keys, in the following order: from
        the people object (for an initialized sim), from the popdict (for one in
        the process of being initialized), from the beta_layer parameter (for an
        uninitialized sim), or by assuming a default (if none of the above are
        available).
        '''
        try:
            keys = list(self['beta_layer'].keys()) # Get keys from beta_layer since the "most required" layer parameter
        except: # pragma: no cover
            keys = []
        return keys


    def reset_layer_pars(self, layer_keys=None, force=False):
        '''
        Reset the parameters to match the population.

        Args:
            layer_keys (list): override the default layer keys (use stored keys by default)
            force (bool): reset the parameters even if they already exist
        '''
        if layer_keys is None:
            if self.people is not None: # If people exist
                layer_keys = self.people.contacts.keys()
            elif self.popdict is not None:
                layer_keys = self.popdict['layer_keys']
        znpar.reset_layer_pars(self.pars, layer_keys=layer_keys, force=force)
        return


    def validate_layer_pars(self):
        '''
        Handle layer parameters, since they need to be validated after the population
        creation, rather than before.
        '''

        # First, try to figure out what the layer keys should be and perform basic type checking
        layer_keys = self.layer_keys()
        layer_pars = znpar.layer_pars # The names of the parameters that are specified by layer
        for lp in layer_pars:
            val = self[lp]
            if sc.isnumber(val): # It's a scalar instead of a dict, assume it's all contacts
                self[lp] = {k:val for k in layer_keys}

        # Handle key mismatches
        for lp in layer_pars:
            lp_keys = set(self.pars[lp].keys())
            if not lp_keys == set(layer_keys):
                errormsg = 'At least one layer parameter is inconsistent with the layer keys; all parameters must have the same keys:'
                errormsg += f'\nsim.layer_keys() = {layer_keys}'
                for lp2 in layer_pars: # Fail on first error, but re-loop to list all of them
                    errormsg += f'\n{lp2} = ' + ', '.join(self.pars[lp2].keys())
                raise sc.KeyNotFoundError(errormsg)

        # Handle mismatches with the population
        if self.people is not None:
            pop_keys = set(self.people.contacts.keys())
            if pop_keys != set(layer_keys): # pragma: no cover
                if not len(pop_keys):
                    errormsg = f'Your population does not have any layer keys, but your simulation does {layer_keys}. If you called cv.People() directly, you probably need cv.make_people() instead.'
                    raise sc.KeyNotFoundError(errormsg)
                else:
                    errormsg = f'Please update your parameter keys {layer_keys} to match population keys {pop_keys}. You may find sim.reset_layer_pars() helpful.'
                    raise sc.KeyNotFoundError(errormsg)

        return


    def validate_pars(self, validate_layers=True):
        '''
        Some parameters can take multiple types; this makes them consistent.

        Args:
            validate_layers (bool): whether to validate layer parameters as well via validate_layer_pars() -- usually yes, except during initialization
        '''

        # Handle population size
        pop_size   = self.pars.get('pop_size')
        scaled_pop = self.pars.get('scaled_pop')
        pop_scale  = self.pars.get('pop_scale')
        if scaled_pop is not None: # If scaled_pop is supplied, try to use it
            if pop_scale in [None, 1.0]: # Normal case, recalculate population scale
                self['pop_scale'] = scaled_pop/pop_size
            else: # Special case, recalculate number of agents
                self['pop_size'] = int(scaled_pop/pop_scale)

        # Handle types
        for key in ['pop_size']:
            try:
                self[key] = int(self[key])
            except Exception as E:
                errormsg = f'Could not convert {key}={self[key]} of {type(self[key])} to integer'
                raise ValueError(errormsg) from E

        # Handle start day
        start_day = self['start_day'] # Shorten
        if start_day in [None, 0]: # Use default start day
            start_day = '2020-03-01'
        self['start_day'] = sc.date(start_day)

        # Handle end day and n_days
        end_day = self['end_day']
        n_days = self['n_days']
        if end_day:
            self['end_day'] = sc.date(end_day)
            n_days = sc.daydiff(self['start_day'], self['end_day'])
            if n_days <= 0:
                errormsg = f"Number of days must be >0, but you supplied start={str(self['start_day'])} and end={str(self['end_day'])}, which gives n_days={n_days}"
                raise ValueError(errormsg)
            else:
                self['n_days'] = int(n_days)
        else:
            if n_days:
                self['n_days'] = int(n_days)
                self['end_day'] = self.date(n_days) # Convert from the number of days to the end day
            else:
                errormsg = f'You must supply one of n_days and end_day, not "{n_days}" and "{end_day}"'
                raise ValueError(errormsg)

        # Handle population data
        popdata_choices = ['random', 'behaviour_module']
        choice = self['pop_type']
        if choice and choice not in popdata_choices: # pragma: no cover
            choicestr = ', '.join(popdata_choices)
            errormsg = f'Population type "{choice}" not available; choices are: {choicestr}'
            raise ValueError(errormsg)



        # Optionally handle layer parameters
        if validate_layers:
            self.validate_layer_pars()

        # Handle versioning
        if self._legacy_trans is None:
            default_ver = self._default_ver if self._default_ver else self.version
            self._legacy_trans = sc.compareversions(default_ver, '<3.1.1') # Handle regression

        # Handle verbose
        if self['verbose'] == 'brief':
            self['verbose'] = -1
        if not sc.isnumber(self['verbose']): # pragma: no cover
            errormsg = f'Verbose argument should be either "brief", -1, or a float, not {type(self["verbose"])} "{self["verbose"]}"'
            raise ValueError(errormsg)

        return

    def init_pathogens(self):
         
        self.pars['n_pathogens'] = len(self.pars['pathogens'])
        self['n_pathogens'] = self.pars['n_pathogens']

        assert self.pars['n_pathogens'] >=1, f"No pathogens are given to the simulation!"
          

        #Initialize the variants 
        for i in range(len(self.pathogens)): 
            self.pathogens[i].initialize(self)
             
        self.initialized_pathogens = True
         

    def init_pathogen_interactions(self):
 
        #Matrices of pathogen-pathogen interaction
        for key in ['Mtrans', 'Miimm', 'Mcimm', 'Mdur', 'Msev']:
            if self.pars[key] is None:
                if key in ['Mcimm']:
                    self[key] = np.matrix(np.full((self.n_pathogens, self.n_pathogens), 0.001, dtype = float)) #we assume a cross-immunity of 0 between pathogens (if no matrix is provided)
                else:
                    self[key] = np.matrix(np.ones((self.n_pathogens, self.n_pathogens),dtype = float))
            else:
                #validate  
                assert len(self.pars[key]) == self.n_pathogens, f'Matrix {key} is not the right size (should be {self.n_pathogens}X{self.n_pathogens})'
                for i in range(len(self.pars[key])):
                    assert len(self.pars[key][i]) == self.n_pathogens, f'Matrix {key} is not the right size (should be {self.n_pathogens}X{self.n_pathogens})'
                    for j in range(len(self.pars[key][i])):
                        if self.pars[key][i][j] == 0:
                            self.pars[key][i][j] = 0.001
                
                self[key] = np.matrix(self.pars[key])
                 
        



    def init_results(self):
        '''
        Create the main results structure.
        We differentiate between flows, stocks, and cumulative results
        The prefix "new" is used for flow variables, i.e. counting new events (infections/deaths/recoveries) on each timestep
        The prefix "n" is used for stock variables, i.e. counting the total number in any given state (sus/inf/rec/etc) on any particular timestep
        The prefix "cum" is used for cumulative variables, i.e. counting the total number that have ever been in a given state at some point in the sim
        Note that, by definition, n_dead is the same as cum_deaths and n_recovered is the same as cum_recoveries, so we only define the cumulative versions
        '''

        def init_res(*args, **kwargs):
            ''' Initialize a single result object '''
            output = znb.Result(*args, **kwargs, npts=self.npts)
            return output

        for i in range(len(self.pathogens)):
            self.results[i] = {}
            # Flows and cumulative flows
            for key,label in znd.result_flows.items():
                self.results[i][f'cum_{key}'] = init_res(f'Cumulative {label}')  # Cumulative variables -- e.g. "Cumulative infections"

            for key,label in znd.result_flows.items(): # Repeat to keep all the cumulative keys together
                self.results[i][f'new_{key}'] = init_res(f'Number of new {label}') # Flow variables -- e.g. "Number of new infections"

            # Stock variables
            for key,label in znd.result_stocks.items():
                self.results[i][f'n_{key}'] = init_res(label)

            # Other variables
            self.results[i]['n_imports']           = init_res('Number of imported infections', scale=True) 
            self.results[i]['n_alive']             = init_res('Number alive', scale=True)
            self.results[i]['n_naive']             = init_res('Number never infected', scale=True)
            self.results[i]['n_preinfectious']     = init_res('Number preinfectious', scale=True)
            self.results[i]['n_removed']           = init_res('Number removed', scale=True)
            self.results[i]['prevalence']          = init_res('Prevalence', scale=False)
            self.results[i]['incidence']           = init_res('Incidence', scale=False)
            self.results[i]['r_eff']               = init_res('Effective reproduction number', scale=False)
            self.results[i]['doubling_time']       = init_res('Doubling time', scale=False)
            self.results[i]['test_yield']          = init_res('Testing yield', scale=False)
            self.results[i]['rel_test_yield']      = init_res('Relative testing yield', scale=False)
            self.results[i]['frac_vaccinated']     = init_res('Proportion vaccinated', scale=False)
            self.results[i]['pop_imm']            = init_res('Population immunity levels', scale=False)
            self.results[i]['pop_nabs']            = init_res('Population nab levels', scale=False)
            self.results[i]['pop_protection']      = init_res('Population immunity protection', scale=False)
            self.results[i]['pop_symp_protection'] = init_res('Population symptomatic protection', scale=False)

            #IgG Level Results: The result for each day is kept in an array len = num_days +1 which contains an array of 
            #people's IgG levels in the population at each day 
            self.results[i]['IgG_level'] = np.full(((self.pars['n_days'] +1), self.pars['pop_size']), 0, dtype = int)
        
        
            # TODO: Need to check this with Andrew. Discussion. 
            self.results[i]['new_diagnoses_custom']      = init_res('Number of new diagnoses with custom testing module')
            self.results[i]['cum_diagnoses_custom']      = init_res('Cumulative diagnoses with custom testing module')

            # Handle variants 
            self.results[i]['variant'] = {}
            self.results[i]['variant']['prevalence_by_variant'] = init_res('Prevalence by variant', scale=False, n_variants=self.pathogens[i].n_variants)
            self.results[i]['variant']['incidence_by_variant']  = init_res('Incidence by variant', scale=False, n_variants=self.pathogens[i].n_variants)
            for key,label in znd.result_flows_by_variant.items():
                self.results[i]['variant'][f'cum_{key}'] = init_res(f'Cumulative {label}', n_variants=self.pathogens[i].n_variants)  # Cumulative variables -- e.g. "Cumulative infections"
            for key,label in znd.result_flows_by_variant.items():
                self.results[i]['variant'][f'new_{key}'] = init_res(f'Number of new {label}', n_variants=self.pathogens[i].n_variants) # Flow variables -- e.g. "Number of new infections"
            for key,label in znd.result_stocks_by_variant.items(): 
                self.results[i]['variant'][f'n_{key}'] = init_res(label, n_variants=self.pathogens[i].n_variants)



        # Populate the rest of the results
        if self['rescale']:
            scale = 1
        else:
            scale = self['pop_scale']
        self.rescale_vec   = scale*np.ones(self.npts) # Not included in the results, but used to scale them
        self.results['date'] = self.datevec
        self.results['t']    = self.tvec
        self.results_ready   = False

        for key,label in znd.result_flows.items():
            self.results[f'cum_{key}'] = init_res(f'Cumulative {label}')  # Cumulative variables -- e.g. "Cumulative infections"
                
        for key,label in znd.result_flows.items(): # Repeat to keep all the cumulative keys together
            self.results[f'new_{key}'] = init_res(f'Number of new {label}') # Flow variables -- e.g. "Number of new infections"

        self.results['co-infections'] = init_res(f'New co-infections') 

        return

    def init_people(self, popdict=None, init_infections=False, reset=False, verbose=None, **kwargs):
        '''
        Create the people.

        Use ``init_infections=False`` for creating a fresh People object for use
        in future simulations

        Args:
            popdict         (any):  pre-generated people of various formats
            init_infections (bool): whether to initialize infections (default false when called directly)
            reset           (bool): whether to regenerate the people even if they already exist
            verbose         (int):  detail to print
            kwargs          (dict): passed to cv.make_people()
        '''

        # Handle inputs
        if verbose is None:
            verbose = self['verbose']
        if popdict is not None:
            self.popdict = popdict
        if verbose > 0:
            resetstr= ''
            if self.people:
                resetstr = ' (resetting people)' if reset else ' (warning: not resetting sim.people)'
            print(f'Initializing sim{resetstr} with {self["pop_size"]:0n} people for {self["n_days"]} days')

        if not self.initialized_pathogens:
            self.init_pathogens()

        if self.popfile and self.popdict is None: # If we passed zn.popfile, load it into self.popdict
            self.popdict = znm.load(self.popfile)

        # Actually make the people 
        if(self.popdict is None and self.popfile is None):
            self.people = znpop.make_people(self, popdict = None, reset=reset, verbose=verbose, **kwargs)  #Random Population
            
        elif (isinstance(self.popdict, znppl.People)): #input is Popfile or People (infection.People)
            self.people = znpop.make_people(self, popdict = self.popdict, reset=reset, verbose=verbose, **kwargs)  

        else: #input is BehaviourModel
                self.people = znpop.make_people(self, popdict = self.popdict.popdict, workplaces = self.popdict.workplaces, n_workplaces = self.popdict.n_workplaces,reset=reset, verbose=verbose, **kwargs)


        self.people.initialize(sim_pars=self.pars, sim = self) # Fully initialize the people
        self.reset_layer_pars(force=False) # Ensure that layer keys match the loaded population
             
         

        return self
 

    def init_immunity(self, create=False):
        ''' Initialize immunity matrices and precompute nab waning for each variant '''
        if self['use_waning']:
            znimm.init_immunity(self, create=create)
        return

    def init_infections(self, force=False, verbose=None):
        if verbose is None:
            verbose = self['verbose']

        # If anyone is non-naive, don't re-initialize
        if self.people.count_not('naive') == 0 or force: # Everyone is naive

            # Handle anyone who isn't susceptible.
            if self['frac_susceptible'] < 1:
                inds = znu.choose(self['pop_size'], np.round((1-self['frac_susceptible'])*self['pop_size']))
                self.people.make_nonnaive(inds=inds) 

            for i in range(len(self.pathogens)):
                # Create the seed infections
                if self.pathogens[i].pop_infected != 0:
                    inds = znu.choose(self['pop_size'], self.pathogens[i].pop_infected)
                    self.people.infect(inds=inds, layer='seed_infection', pathogen_index = i) # Not counted by results since flows are re-initialized during the step

        elif verbose:
            print(f'People already initialized with {self.people.count_not("naive")} people non-naive and {self.people.count1d("exposed")} exposed (with any pathogen); not reinitializing')

        return

    def rescale(self):
        ''' Dynamically rescale the population -- used during step() '''
        if self['rescale']:
            pop_scale = self['pop_scale']
            current_scale = self.rescale_vec[self.t]
            if current_scale < pop_scale: # We have room to rescale
                not_naive_inds = self.people.false('naive') # Find everyone not naive
                n_not_naive = len(not_naive_inds) # Number of people who are not naive
                n_people = self['pop_size'] # Number of people overall
                current_ratio = n_not_naive/n_people # Current proportion not naive
                threshold = self['rescale_threshold'] # Threshold to trigger rescaling
                if current_ratio > threshold: # Check if we've reached point when we want to rescale
                    max_ratio = pop_scale/current_scale # We don't want to exceed the total population size
                    proposed_ratio = max(current_ratio/threshold, self['rescale_factor']) # The proposed ratio to rescale: the rescale factor, unless we've exceeded it
                    scaling_ratio = min(proposed_ratio, max_ratio) # We don't want to scale by more than the maximum ratio
                    self.rescale_vec[self.t:] *= scaling_ratio # Update the rescaling factor from here on
                    n = int(round(n_not_naive*(1.0-1.0/scaling_ratio))) # For example, rescaling by 2 gives n = 0.5*not_naive_inds
                    choices = znu.choose(max_n=n_not_naive, n=n) # Choose who to make naive again
                    new_naive_inds = not_naive_inds[choices] # Convert these back into indices for people
                    self.people.make_naive(new_naive_inds) # Make people naive again
        return


    def step(self):  
        '''
        Step the simulation forward in time. Usually, the user would use sim.run()
        rather than calling sim.step() directly.
        '''
       # print(self.results[0]['n_infectious'])
        #print(self.results[1]['n_infectious'])

        # Set the time and if we have reached the end of the simulation, then do nothing
        if self.complete:
            raise AlreadyRunError('Simulation already complete (call sim.initialize() to re-run)')

        t = self.t 
         
        # If it's the first timestep, infect people
        #if t == 0:
            #self.init_infections(verbose=False)
            
        # Perform initial operations
        #self.rescale() # Check if we need to rescale
        people = self.people # Shorten this for later use
         
        contacts = people.update_contacts() # Compute new contacts. For dynamic contacts. 
        hosp_max = people.count1d('severe')   > self['n_beds_hosp'] if self['n_beds_hosp'] is not None else False # Check for acute bed constraint
        icu_max  = people.count1d('critical') > self['n_beds_icu']  if self['n_beds_icu']  is not None else False # Check for ICU bed constraint
         
         
        people.update_states_pre(t=t) # Update the state of everyone and count the flows. This isn't infecting people nor updating their SEIR's. The date of infection seems to be pre-assigned. 
        #For every pathogen, import pathogens, and import variants
        for current_pathogen in range(len(self.pathogens)):
            # Randomly infect some people (imported infections)
            if self.pathogens[current_pathogen].n_imports:  
                n_imports = 0
                if isinstance(self.pathogens[current_pathogen].n_imports, list):#If we provide an array, then its a timeline of infections, and no random sampling
                    if (self.t < len(self.pathogens[current_pathogen].n_imports)):
                        n_imports = self.pathogens[current_pathogen].n_imports[self.t]
                else: 
                    n_imports = znu.poisson(self.pathogens[current_pathogen].n_imports/self.rescale_vec[self.t])  
                if n_imports>0:
                    importation_inds = znu.choose(max_n=self['pop_size'], n=n_imports)
                    people.infect(inds=importation_inds, hosp_max=hosp_max, icu_max=icu_max, layer='importation', pathogen_index = current_pathogen)
             
            # Add variants
            for variant in self.pathogens[current_pathogen].variants:
                if isinstance(variant, pat.Pathogen.Variant):
                    variant.apply(self) 
    
        
        # Compute viral loads for every pathogen
        viral_load = {}
        for current_pathogen in range(len(self.pathogens)): 
            x_p1, y_p1 = people.x_p1[current_pathogen], people.y_p1[current_pathogen]
            x_p2, y_p2 = people.x_p2[current_pathogen], people.y_p2[current_pathogen]
            x_p3, y_p3 = people.x_p3[current_pathogen], people.y_p3[current_pathogen]
            min_vl = zndt.default_float(self.pathogens[current_pathogen].viral_levels['min_vl'])
            max_vl = zndt.default_float(self.pathogens[current_pathogen].viral_levels['max_vl'])
            people.viral_load[current_pathogen], viral_load[current_pathogen] = znu.compute_viral_load(t, x_p1, y_p1, x_p2, y_p2, x_p3, y_p3, min_vl, max_vl)
        
        
        # Apply interventions
            
        
        for current_pathogen in range(len(self.pathogens)): 
            # Implement state changes relating to quarantine and isolation/diagnosis; compute associated statistics
            people.update_states_post(pathogen = current_pathogen)

            # Shorten useful parameters
            nv = self.pathogens[current_pathogen].n_variants # variants of CURRENT PATHOGEN
            sus = people.p_susceptible[current_pathogen]
            symp = people.p_symptomatic[current_pathogen]
            diag = people.p_diagnosed[current_pathogen]
            prel_trans = people.rel_trans[current_pathogen]
            prel_sus   = people.rel_sus[current_pathogen]
            iso  = people.isolated

            # Iterate through n_variants to calculate infections. The meat of the simulation. 
            for variant in range(nv):

                # Check immunity 
                znimm.check_immunity(people, variant, current_pathogen)

                # Deal with variant parameters
                rel_beta = 1
                asymp_factor = self.pathogens[current_pathogen].asymp_factor
                if variant: 
                    if variant != 0: #if not wild type, multiply by variant rel_beta
                        rel_beta *= self.pathogens[current_pathogen].variants[variant-1].p['rel_beta'] 

                beta = zndt.default_float(self.pathogens[current_pathogen].beta * rel_beta)

                for lkey, layer in contacts.items():
                    p1 = layer['p1']
                    p2 = layer['p2']
                    betas = layer['beta']

                    # Compute relative transmission and susceptibility
                    inf_variant = people.p_infectious[current_pathogen] * (people.p_infectious_variant[current_pathogen] == variant) # TODO: move out of loop?
                    sus_imm = people.sus_imm[current_pathogen,variant,:]
                    iso_factor  = zndt.default_float(self['iso_factor'][lkey]) # Effect of isolating. 
                    quar_factor = zndt.default_float(self['quar_factor'][lkey]) # Ex: 0.2. Probably the effect on beta of quarantining. 
                    beta_layer  = zndt.default_float(self['beta_layer'][lkey]) # A scalar; beta for the layer. Ex: 1.0. 
                    rel_trans, rel_sus = znu.compute_trans_sus(prel_trans, prel_sus, inf_variant, sus, beta_layer, viral_load[current_pathogen], symp, diag, iso, asymp_factor, iso_factor, quar_factor, sus_imm)
                     
                    rel_trans = p_int.mod_rel_trans(current_pathogen, people.p_exposed, self.n_pathogens, rel_trans, self['Mtrans'])
                    rel_sus = p_int.mod_rel_sus(current_pathogen, rel_sus, people.p_exposed, self['Miimm'], self['Mcimm'], people.sus_imm, self.n_pathogens, self.pars['pop_size'])

                    # Calculate actual transmission
                    pairs = [[p1,p2]] if not self._legacy_trans else [[p1,p2], [p2,p1]] # Support slower legacy method of calculation, but by default skip this loop
                    for p1,p2 in pairs:
                        source_inds, target_inds = znu.compute_infections(beta, p1, p2, betas, rel_trans, rel_sus, legacy=self._legacy_trans)  # Calculate transmission! 
                        people.infect(inds=target_inds, hosp_max=hosp_max, icu_max=icu_max, source=source_inds, layer=lkey, variant=variant, pathogen_index = current_pathogen)  # Actually infect people
                                  
        ##### CALCULATE STATISTICS #####
        for current_pathogen in range(len(self.pathogens)):

            nv = self.pathogens[current_pathogen].n_variants
            # Update counts for this time step: stocks.
            for key in znd.result_stocks.keys(): 
                #TODO remove this filtering
                if key not in ['known_dead', 'quarantined', 'vaccinated', 'isolated']: #SELECTING WHICH STATES ARE CURRENTLY IMPLEMENTED FOR PER PATHOGEN TRACKING 
                    self.results[current_pathogen][f'n_{key}'][t] = np.count_nonzero(self.people[f'p_{key}'][current_pathogen][self.stratification_indices])
                else:
                    self.results[current_pathogen][f'n_{key}'][t] = np.count_nonzero(self.people[key][self.stratification_indices])
                   
            for key in znd.result_stocks_by_variant.keys():
                for variant in range(nv): 
                    self.results[current_pathogen]['variant'][f'n_{key}'][variant,t] =np.count_nonzero(self.people[f'p_{key}'][current_pathogen,variant,:][self.stratification_indices]) 
        
            # Update counts for this time step: flows
            for key,count in people.flows[current_pathogen].items(): 
                self.results[current_pathogen][key][t] += count
            for key,count in people.flows_variant[current_pathogen].items():
                for variant in range(nv):
                    self.results[current_pathogen]['variant'][key][variant][t] += count[variant]
 
        for key in znd.new_result_flows:     
            self.results[key][t] += people.flows[key] 


        for current_pathogen in range(len(self.pathogens)): 
            # Update nab, immunity and IgG for this time step
            #if self['use_waning']: 
            if self.pathogens[current_pathogen].use_nab_framework: 
                has_nabs = znu.true(people.peak_nab[current_pathogen])
                if len(has_nabs):
                    znimm.update_nab(people, inds=has_nabs, pathogen = current_pathogen)
                znimm.update_IgG(people, current_pathogen)
            else:
                has_imm = np.where(people['imm_peak'][current_pathogen] > 0)
                if len(has_imm): 
                    znimm.update_imm(people, has_imm[0], current_pathogen, self.people['imm_min'],  self.people['imm_peak'],  self.people['decay_rate'], self.people['growth_rate'])
        

        for current_pathogen in range(len(self.pathogens)): 
                                                                                                                   
            self.results[current_pathogen]['pop_protection'][t]      = np.nanmean(people.sus_imm[current_pathogen])  
            self.results[current_pathogen]['pop_symp_protection'][t] = np.nanmean(people.symp_imm[current_pathogen]) 
            
            if self.pathogens[current_pathogen].use_nab_framework: #For IgG levels in population 
                #TESTING
                #print("Size of self.results[current_pathogen]['IgG_level']: ", len(self.results[current_pathogen]['IgG_level']))
                #print("Value of t: ", t)
                #print("Shape of people['IgG_level']: ", people['IgG_level'].shape)
                #print("people['IgG_level']: ", people['IgG_level'])
                # Verify the content of people['IgG_level']
                
                self.results[current_pathogen]['IgG_level'][t]       = np.array(people['IgG_level'])



        # if this is a day we want to test on then apply 
        # day we want to test on is either: 
        # a) a day that is in the sampling interval for longitudinal studies 
        # b) a day that is in the collection period for cross-sectional studies

        # for program in active population surveillance 
        # self.sampler.apply() --> this will be the same people for longitudinal 
        # tester.apply() --> will be the same test 
        # send results.append to store 


        # Apply analyzers (note this comes with Covasim) -- same syntax as interventions
        for i, analyzer in enumerate(self['analyzers']):
            analyzer(self)

        # Tidy up
        self.t += 1
        if self.t == self.npts:
            self.complete = True
        return
            

    def run(self, do_plot=False, until=None, restore_pars=True, reset_seed=True, verbose=None):
        '''
        Run the simulation.

        Args:
            do_plot (bool): whether to plot
            until (int/str): day or date to run until
            restore_pars (bool): whether to make a copy of the parameters before the run and restore it after, so runs are repeatable
            reset_seed (bool): whether to reset the random number stream immediately before run
            verbose (float): level of detail to print, e.g. -1 = one-line output, 0 = no output, 0.1 = print every 10th day, 1 = print every day

        Returns:
            A pointer to the sim object (with results modified in-place)
        '''

        # Initialization steps -- start the timer, initialize the sim and the seed, and check that the sim hasn't been run
        T = sc.timer() 
        if not self.initialized:
            self.initialize()
            self._orig_pars = sc.dcp(self.pars) # Create a copy of the parameters, to restore after the run, in case they are dynamically modified

        if verbose is None:
            verbose = self['verbose']

        if reset_seed:
            # Reset the RNG. If the simulation is newly created, then the RNG will be reset by sim.initialize() so the use case
            # for resetting the seed here is if the simulation has been partially run, and changing the seed is required
            self.set_seed()

        # Check for AlreadyRun errors
        errormsg = None
        until = self.npts if until is None else self.day(until)
        if until > self.npts:
            errormsg = f'Requested to run until t={until} but the simulation end is t={self.npts}'
        if self.t >= until: # NB. At the start, self.t is None so this check must occur after initialization
            errormsg = f'Simulation is currently at t={self.t}, requested to run until t={until} which has already been reached'
        if self.complete:
            errormsg = 'Simulation is already complete (call sim.initialize() to re-run)'
        if self.people.t not in [self.t, self.t-1]: # Depending on how the sim stopped, either of these states are possible
            errormsg = f'The simulation has been run independently from the people (t={self.t}, people.t={self.people.t}): if this is intentional, manually set sim.people.t = sim.t. Remember to save the people object before running the sim.'
        if errormsg:
            raise AlreadyRunError(errormsg)
         
        # Main simulation loop
        while self.t < until:

            # Check if we were asked to stop
            elapsed = T.toc(output=True)
            if self['timelimit'] and elapsed > self['timelimit']:
                sc.printv(f"Time limit ({self['timelimit']} s) exceeded; call sim.finalize() to compute results if desired", 1, verbose)
                return
            elif self['stopping_func'] and self['stopping_func'](self):
                sc.printv("Stopping function terminated the simulation; call sim.finalize() to compute results if desired", 1, verbose)
                return

            # Print progress
            if verbose:
                simlabel = f'"{self.label}": ' if self.label else ''
                string = f'  Running {simlabel}{self.datevec[self.t]} ({self.t:2.0f}/{self.pars["n_days"]}) ({elapsed:0.2f} s) '
                if verbose >= 2:
                    sc.heading(string)
                elif verbose>0:
                    if not (self.t % int(1.0/verbose)):
                        sc.progressbar(self.t+1, self.npts, label=string, length=20, newline=True)
                        
            # Do the heavy lifting -- actually run the model! 
            self.step()   
            self.validate_people_states()

        # If simulation reached the end, finalize the results
        if self.complete:
            self.finalize(verbose=verbose, restore_pars=restore_pars)
            sc.printv(f'Run finished after {elapsed:0.2f} s.\n', 1, verbose)
        
 
        return self
 
    
    def finalize(self, verbose=None, restore_pars=True):
        ''' Compute final results '''

        if self.results_ready:
            # Because the results are rescaled in-place, finalizing the sim cannot be run more than once or
            # otherwise the scale factor will be applied multiple times
            raise AlreadyRunError('Simulation has already been finalized')
         
        # Scale the results
        for p in range(len(self.pathogens)):
            for reskey in self.result_keys():
                if self.results[p][reskey].scale: # Scale the result dynamically
                    self.results[p][reskey].values *= self.rescale_vec
            for reskey in self.result_keys('variant'):
                if self.results[p]['variant'][reskey].scale: # Scale the result dynamically
                    self.results[p]['variant'][reskey].values = np.einsum('ij,j->ij', self.results[p]['variant'][reskey].values, self.rescale_vec)


            # Calculate cumulative results
            for key in znd.result_flows.keys():

                self.results[p][f'cum_{key}'][:] = np.cumsum(self.results[p][f'new_{key}'][:], axis=0)
            for key in znd.result_flows_by_variant.keys():
                for variant in range(self.pathogens[p].n_variants):
                    self.results[p]['variant'][f'cum_{key}'][variant, :] = np.cumsum(self.results[p]['variant'][f'new_{key}'][variant, :], axis=0)
            if self.pars['enable_multiregion']:
                for res in [self.results[p]['cum_infections'], self.results[p]['variant']['cum_infections_by_variant']]: # Include initially infected people
                    for val in self.pathogens[p].pop_infected.values():
                        res.values += val*self.rescale_vec[0]
            else:
                for res in [self.results[p]['cum_infections'], self.results[p]['variant']['cum_infections_by_variant']]: # Include initially infected people
                    res.values += self.results[p]['n_exposed'][0]*self.rescale_vec[0]  

        for key in znd.result_flows.keys(): 
            self.results[f'cum_{key}'][:] = np.cumsum(self.results[f'new_{key}'][:], axis=0)
       
        for i in range(len(self.results['cum_infections'])):
            for p in range(len(self.pathogens)):
                self.results['cum_infections'][i] += self.results[p]['n_exposed'][0] 

        # Finalize interventions and analyzers
        self.finalize_interventions()
        self.finalize_analyzers()

        # Final settings
        self.results_ready = True # Set this first so self.summary() knows to print the results
        self.t -= 1 # During the run, this keeps track of the next step; restore this be the final day of the sim

        # Perform calculations on results
        self.compute_results(verbose=verbose) # Calculate the rest of the results

        for p in range(len(self.pathogens)):
            self.results[p] = sc.objdict(self.results[p]) # Convert results to a odicts/objdict to allow e.g. sim.results.diagnoses

        if restore_pars and self._orig_pars:
            preserved = ['analyzers', 'interventions']
            orig_pars_keys = list(self._orig_pars.keys()) # Get a list of keys so we can iterate over them
            for key in orig_pars_keys:
                if key not in preserved:
                    self.pars[key] = self._orig_pars.pop(key) # Restore everything except for the analyzers and interventions

        # Optionally print summary output
        if verbose: # Verbose is any non-zero value
            if verbose>0: # Verbose is any positive number
                self.summarize() # Print medium-length summary of the sim
            else:
                self.brief() # Print brief summary of the sim

        return


    def compute_results(self, verbose=None):
        ''' Perform final calculations on the results '''
        self.compute_states()
        self.compute_yield()
        self.compute_doubling()

        for p in range(len(self.pathogens)):
            self.compute_r_eff(pathogen = p)

        self.compute_summary()
        return


    def compute_states(self):
        '''
        Compute prevalence, incidence, and other states. Prevalence is the current
        number of infected people divided by the number of people who are alive.
        Incidence is the number of new infections per day divided by the susceptible
        population. Also calculates the number of people alive, the number preinfectious,
        the number removed, and recalculates susceptibles to handle scaling.
        '''
        for p in range(len(self.pathogens)):
            res = self.results[p]
            count_recov = 1-self['use_waning'] # If waning is on, don't count recovered people as removed
            self.results[p]['n_alive'][:]         = len(self.stratification_indices) - res['cum_deaths'][:] # Number of people still alive
            self.results[p]['n_naive'][:]         = len(self.stratification_indices) - res['cum_deaths'][:] - res['n_recovered'][:] - res['n_exposed'][:] # Number of people naive
            self.results[p]['n_susceptible'][:]   = res['n_alive'][:] - res['n_exposed'][:] - count_recov*res['cum_recoveries'][:] # Recalculate the number of susceptible people, not agents
            self.results[p]['n_preinfectious'][:] = res['n_exposed'][:] - res['n_infectious'][:] # Calculate the number not yet infectious: exposed minus infectious
            self.results[p]['n_removed'][:]       = count_recov*res['cum_recoveries'][:] + res['cum_deaths'][:] # Calculate the number removed: recovered + dead
            self.results[p]['prevalence'][:]      = res['n_exposed'][:]/res['n_alive'][:] # Calculate the prevalence
            self.results[p]['incidence'][:]       = res['new_infections'][:]/res['n_susceptible'][:] # Calculate the incidence
            self.results[p]['frac_vaccinated'][:] = res['n_vaccinated'][:]/res['n_alive'][:] # Calculate the fraction vaccinated

            self.results[p]['variant']['incidence_by_variant'][:] = np.einsum('ji,i->ji',res['variant']['new_infections_by_variant'][:], 1/res['n_susceptible'][:]) # Calculate the incidence
            self.results[p]['variant']['prevalence_by_variant'][:] = np.einsum('ji,i->ji',res['variant']['new_infections_by_variant'][:], 1/res['n_alive'][:])  # Calculate the prevalence

        return


    def compute_yield(self):
        '''
        Compute test yield -- number of positive tests divided by the total number
        of tests, also called test positivity rate. Relative yield is with respect
        to prevalence: i.e., how the yield compares to what the yield would be from
        choosing a person at random from the population.
        '''
        # Absolute yield
        for p in range(len(self.pathogens)):
            res = self.results[p]
            inds = znu.true(res['new_tests'][:]) # Pull out non-zero numbers of tests
            self.results[p]['test_yield'][inds] = res['new_diagnoses'][inds]/res['new_tests'][inds] # Calculate the yield

            # Relative yield
            inds = znu.true(res['n_infectious'][:]) # To avoid divide by zero if no one is infectious
            denom = res['n_infectious'][inds] / (res['n_alive'][inds] - res['cum_diagnoses'][inds]) # Alive + undiagnosed people might test; infectious people will test positive
            self.results[p]['rel_test_yield'][inds] = self.results[p]['test_yield'][inds]/denom # Calculate the relative yield
        return


    def compute_doubling(self, window=3, max_doubling_time=30):
        '''
        Calculate doubling time using exponential approximation -- a more detailed
        approach is in utils.py. Compares infections at time t to infections at time
        t-window, and uses that to compute the doubling time. For example, if there are
        100 cumulative infections on day 12 and 200 infections on day 19, doubling
        time is 7 days.

        Args:
            window (float): the size of the window used (larger values are more accurate but less precise)
            max_doubling_time (float): doubling time could be infinite, so this places a bound on it

        Returns:
            doubling_time (array): the doubling time results array
        '''
        for p in range(len(self.pathogens)):
            cum_infections = self.results[p]['cum_infections'].values
            infections_now = cum_infections[window:]
            infections_prev = cum_infections[:-window]
            use = (infections_prev > 0) & (infections_now > infections_prev)
            doubling_time = window * np.log(2) / np.log(infections_now[use] / infections_prev[use])
            self.results[p]['doubling_time'][:] = np.nan
            self.results[p]['doubling_time'][window:][use] = np.minimum(doubling_time, max_doubling_time)
        return self.results[p]['doubling_time'].values


    def compute_r_eff(self, method='daily', smoothing=2, window=7, pathogen = 0):
        '''
        Effective reproduction number based on number of people each person infected.

        Args:
            method (str): 'daily' uses daily infections, 'infectious' counts from the date infectious, 'outcome' counts from the date recovered/dead
            smoothing (int): the number of steps to smooth over for the 'daily' method
            window (int): the size of the window used for 'infectious' and 'outcome' calculations (larger values are more accurate but less precise)

        Returns:
            r_eff (array): the r_eff results array
        '''

        # Initialize arrays to hold sources and targets infected each day
        sources = np.zeros(self.npts)
        targets = np.zeros(self.npts)
        window = int(window)

        # Default method -- calculate the daily infections
        if method == 'daily':

            # Find the dates that everyone became infectious and recovered, and hence calculate infectious duration
            recov_inds   = self.people.defined('date_p_recovered',0)
            dead_inds    = self.people.defined('date_p_dead',0)
            date_recov   = self.people.date_p_recovered[pathogen, recov_inds]
            date_dead    = self.people.date_dead[dead_inds]
            date_outcome = np.concatenate((date_recov, date_dead))
            inds         = np.concatenate((recov_inds, dead_inds))
            date_inf     = self.people.date_p_infectious[pathogen, inds]
            if len(date_outcome):
                mean_inf     = date_outcome.mean() - date_inf.mean()
            else:
                warnmsg ='There were no infections during the simulation'
                znm.warn(warnmsg)
                mean_inf = 0 # Doesn't matter since r_eff is 0

            # Calculate R_eff as the mean infectious duration times the number of new infectious divided by the number of infectious people on a given day
            new_infections = self.results[pathogen]['new_infections'].values - self.results[pathogen]['n_imports'].values
            n_infectious = self.results[pathogen]['n_infectious'].values
            raw_values = mean_inf*np.divide(new_infections, n_infectious, out=np.zeros(self.npts), where=n_infectious>0)

            # Handle smoothing, including with too-short arrays
            len_raw = len(raw_values) # Calculate the number of raw values
            dur_pars = self.pathogens[pathogen].dur[0] if isinstance(self.pathogens[pathogen].dur, list) else self.pathogens[pathogen].dur # Note: does not take variants into account
            if len_raw >= 3: # Can't smooth arrays shorter than this since the default smoothing kernel has length 3
                initial_period = dur_pars['exp2inf']['par1'] + dur_pars['asym2rec']['par1'] # Approximate the duration of the seed infections for averaging
                initial_period = int(min(len_raw, initial_period)) # Ensure we don't have too many points
                for ind in range(initial_period): # Loop over each of the initial inds
                    raw_values[ind] = raw_values[ind:initial_period].mean() # Replace these values with their average
                values = sc.smooth(raw_values, smoothing)
                values[:smoothing] = raw_values[:smoothing] # To avoid numerical effects, replace the beginning and end with the original
                values[-smoothing:] = raw_values[-smoothing:]
            else:
                values = raw_values

        # Alternate (traditional) method -- count from the date of infection or outcome
        elif method in ['infectious', 'outcome']:

            # Store a mapping from each source to their date
            source_dates = {}

            for t in self.tvec:

                # Sources are easy -- count up the arrays for all the people who became infections on that day
                if method == 'infectious':
                    inds = znu.true(t == self.people.date_p_infectious[pathogen]) # Find people who became infectious on this timestep
                elif method == 'outcome':
                    recov_inds = znu.true(t == self.people.date_p_recovered[pathogen]) # Find people who recovered on this timestep
                    dead_inds  = znu.true(t == self.people.date_dead)  # Find people who died on this timestep
                    inds       = np.concatenate((recov_inds, dead_inds))
                sources[t] = len(inds)

                # Create the mapping from sources to dates
                for ind in inds:
                    source_dates[ind] = t

            # Targets are hard -- loop over the transmission tree
            for transdict in self.people.infection_log:
                source = transdict['source']
                if source is not None and source in source_dates: # Skip seed infections and people with e.g. recovery after the end of the sim
                    source_date = source_dates[source]
                    targets[source_date] += 1

                # for ind in inds:
                #     targets[t] += len(self.people.transtree.targets[ind])

            # Populate the array -- to avoid divide-by-zero, skip indices that are 0
            r_eff = np.divide(targets, sources, out=np.full(self.npts, np.nan), where=sources > 0)

            # Use stored weights calculate the moving average over the window of timesteps, n
            num = np.nancumsum(r_eff * sources)
            num[window:] = num[window:] - num[:-window]
            den = np.cumsum(sources)
            den[window:] = den[window:] - den[:-window]
            values = np.divide(num, den, out=np.full(self.npts, np.nan), where=den > 0)

        # Method not recognized
        else: # pragma: no cover
            errormsg = f'Method must be "daily", "infectious", or "outcome", not "{method}"'
            raise ValueError(errormsg)

        # Set the values and return
        self.results[pathogen]['r_eff'].values[:] = values

        return self.results[pathogen]['r_eff'].values


    def compute_gen_time(self, pathogen = 0):
        '''
        Calculate the generation time (or serial interval). There are two
        ways to do this calculation. The 'true' interval (exposure time to
        exposure time) or 'clinical' (symptom onset to symptom onset).

        Returns:
            gen_time (dict): the generation time results
        '''

        intervals1 = np.zeros(len(self.people))
        intervals2 = np.zeros(len(self.people))
        pos1 = 0
        pos2 = 0
        date_exposed = self.people.date_p_exposed[pathogen]
        date_symptomatic = self.people.date_p_symptomatic[pathogen]

        for infection in self.people.infection_log:
            if infection['source'] is not None:
                source_ind = infection['source']
                target_ind = infection['target']
                intervals1[pos1] = date_exposed[target_ind] - date_exposed[source_ind]
                pos1 += 1
                if np.isfinite(date_symptomatic[source_ind]) and np.isfinite(date_symptomatic[target_ind]):
                    intervals2[pos2] = date_symptomatic[target_ind] - date_symptomatic[source_ind]
                    pos2 += 1

        self.results[pathogen]['gen_time'] = {
                'true':         np.mean(intervals1[:pos1]),
                'true_std':     np.std(intervals1[:pos1]),
                'clinical':     np.mean(intervals2[:pos2]),
                'clinical_std': np.std(intervals2[:pos2])}
        return self.results[pathogen]['gen_time']


    def compute_summary(self, full=None, t=None, update=True, output=False, require_run=False):
        '''
        Compute the summary dict and string for the sim. Used internally; see
        sim.summarize() for the user version.

        Args:
            full (bool): whether or not to print all results (by default, only cumulative)
            t (int/str): day or date to compute summary for (by default, the last point)
            update (bool): whether to update the stored sim.summary
            output (bool): whether to return the summary
            require_run (bool): whether to raise an exception if simulations have not been run yet
        '''
        if t is None:
            t = self.day(self.t)

        # Compute the summary
        if require_run and not self.results_ready:
            errormsg = 'Simulation not yet run'
            raise RuntimeError(errormsg)

        summary = {}
        
        for p in range(len(self.pathogens)):
            summary[p] = sc.objdict()
            for key in self.result_keys():
                summary[p][key] = self.results[p][key][t]

        for key in ['cum_infections', 'cum_reinfections', 'cum_infectious','cum_symptomatic', 'cum_severe', 'cum_isolated', 'cum_critical','cum_recoveries', 'cum_deaths']:
            summary[key] = self.results[key][t]
        # Update the stored state
        if update:
            self.summary = summary

        # Optionally return
        if output:
            return summary
        else:
            return


    def summarize(self, full=False, t=None, sep=None, output=False):
        '''
        Print a medium-length summary of the simulation, drawing from the last time
        point in the simulation by default. Called by default at the end of a sim run.
        See also sim.disp() (detailed output) and sim.brief() (short output).

        Args:
            full   (bool):    whether or not to print all results (by default, only cumulative)
            t      (int/str): day or date to compute summary for (by default, the last point)
            sep    (str):     thousands separator (default ',')
            output (bool):    whether to return the summary instead of printing it

        **Examples**::

            sim = cv.Sim(label='Example sim', verbose=0) # Set to run silently
            sim.run() # Run the sim
            sim.summarize() # Print medium-length summary of the sim
            sim.summarize(t=24, full=True) # Print a "slice" of all sim results on day 24
        '''

        # Compute the summary
        summary = self.compute_summary(full=full, t=t, update=False, output=True)

        # Construct the output string
        if sep is None: sep = zno.sep # Default separator
        labelstr = f' "{self.label}"' if self.label else ''
        string = f'\nSimulation{labelstr} summary:\n'
        
        for p in range(len(self.pathogens)):
            string += f'   {"":5}Summary of pathogen {self.pathogens[p].label}:\n'
            for key in self.result_keys():
                if key in ['cum_quarantined', 'cum_isolated']:
                    continue
                if full or key.startswith('cum_'): 
                        val = np.round(summary[p][key])
                        string += f'   {val:15,.0f} {self.results[p][key].name.lower()}\n'.replace(',', sep) # Use replace since it's more flexible
            string += '\n'

         
        coinfections = self.results['co-infections'] 
        quarantines = self.results['cum_quarantined'][self.t]
        iso = self.results['cum_isolated'][self.t]
        string += f'   {sum(coinfections):15,.0f} cumulative co-infections\n'
        string += f'   {iso:15,.0f} cumulative isolations started\n'
        string += f'   {quarantines:15,.0f} cumulative quarantines started\n'
        #string += f'   {sum(coinfections_deaths):15,.0f} cumulative deaths when co-infected'
        string += '\n'
        # Print or return string
        if not output:
            print(string)
        else:
            return string


    def disp(self, output=False):
        '''
        Display a verbose description of a sim. See also sim.summarize() (medium
        length output) and sim.brief() (short output).

        Args:
            output (bool): if true, return a string instead of printing output

        **Example**::

            sim = cv.Sim(label='Example sim', verbose=0) # Set to run silently
            sim.run() # Run the sim
            sim.disp() # Displays detailed output
        '''
        string = self._disp()
        if not output:
            print(string)
        else:
            return string


    def brief(self, output=False):
        '''
        Print a one-line description of a sim. See also sim.disp() (detailed output)
        and sim.summarize() (medium length output). The symbol "⚙" is used to show
        infections, and "☠" is used to show deaths.

        Args:
            output (bool): if true, return a string instead of printing output

        **Example**::

            sim = cv.Sim(label='Example sim', verbose=0) # Set to run silently
            sim.run() # Run the sim
            sim.brief() # Prints one-line output
        '''
        string = self._brief()
        if not output:
            print(string)
        else:
            return string

     

def diff_sims(sim1, sim2, skip_key_diffs=False, skip=None, output=False, die=False):
    '''
    Compute the difference of the summaries of two simulations, and print any
    values which differ.

    Args:
        sim1 (sim/dict): either a simulation object or the sim.summary dictionary
        sim2 (sim/dict): ditto
        skip_key_diffs (bool): whether to skip keys that don't match between sims
        skip (list): a list of values to skip
        output (bool): whether to return the output as a string (otherwise print)
        die (bool): whether to raise an exception if the sims don't match
        require_run (bool): require that the simulations have been run

    **Example**::

        s1 = cv.Sim(beta=0.01)
        s2 = cv.Sim(beta=0.02)
        s1.run()
        s2.run()
        cv.diff_sims(s1, s2)
    '''

    if isinstance(sim1, Sim):
        sim1 = sim1.compute_summary(update=False, output=True, require_run=True)
    if isinstance(sim2, Sim):
        sim2 = sim2.compute_summary(update=False, output=True, require_run=True)
    for sim in [sim1, sim2]:
        if not isinstance(sim, dict): # pragma: no cover
            errormsg = f'Cannot compare object of type {type(sim)}, must be a sim or a sim.summary dict'
            raise TypeError(errormsg)

    # Compare keys
    keymatchmsg = ''
    sim1_keys = set(sim1.keys())
    sim2_keys = set(sim2.keys())
    if sim1_keys != sim2_keys and not skip_key_diffs: # pragma: no cover
        keymatchmsg = "Keys don't match!\n"
        missing = list(sim1_keys - sim2_keys)
        extra   = list(sim2_keys - sim1_keys)
        if missing:
            keymatchmsg += f'  Missing sim1 keys: {missing}\ns'
        if extra:
            keymatchmsg += f'  Extra sim2 keys: {extra}\n'

    # Compare values
    valmatchmsg = ''
    mismatches = {}
    skip = sc.tolist(skip)
    for key in sim2.keys(): # To ensure order
        if key in sim1_keys and key not in skip: # If a key is missing, don't count it as a mismatch
            sim1_val = sim1[key] if key in sim1 else 'not present'
            sim2_val = sim2[key] if key in sim2 else 'not present'
            both_nan = sc.isnumber(sim1_val, isnan=True) and sc.isnumber(sim2_val, isnan=True)
            if sim1_val != sim2_val and not both_nan:
                mismatches[key] = {'sim1': sim1_val, 'sim2': sim2_val}

    if len(mismatches):
        valmatchmsg = '\nThe following values differ between the two simulations:\n'
        df = pd.DataFrame.from_dict(mismatches).transpose()
        diff   = []
        ratio  = []
        change = []
        small_change = 1e-3 # Define a small change, e.g. a rounding error
        for mdict in mismatches.values():
            old = mdict['sim1']
            new = mdict['sim2']
            numeric = sc.isnumber(sim1_val) and sc.isnumber(sim2_val)
            if numeric and old>0:
                this_diff  = new - old
                this_ratio = new/old
                abs_ratio  = max(this_ratio, 1.0/this_ratio)

                # Set the character to use
                if abs_ratio<small_change:
                    change_char = '≈'
                elif new > old:
                    change_char = '↑'
                elif new < old:
                    change_char = '↓'
                else:
                    errormsg = f'Could not determine relationship between sim1={old} and sim2={new}'
                    raise ValueError(errormsg)

                # Set how many repeats it should have
                repeats = 1
                if abs_ratio >= 1.1:
                    repeats = 2
                if abs_ratio >= 2:
                    repeats = 3
                if abs_ratio >= 10:
                    repeats = 4

                this_change = change_char*repeats
            else: # pragma: no cover
                this_diff   = np.nan
                this_ratio  = np.nan
                this_change = 'N/A'

            diff.append(this_diff)
            ratio.append(this_ratio)
            change.append(this_change)

        df['diff'] = diff
        df['ratio'] = ratio
        for col in ['sim1', 'sim2', 'diff', 'ratio']:
            df[col] = df[col].round(decimals=3)
        df['change'] = change
        valmatchmsg += str(df)

    # Raise an error if mismatches were found
    mismatchmsg = keymatchmsg + valmatchmsg
    if mismatchmsg: # pragma: no cover
        if die:
            raise ValueError(mismatchmsg)
        elif output:
            return mismatchmsg
        else:
            print(mismatchmsg)
    else:
        if not output:
            print('Sims match')
    return


class AlreadyRunError(RuntimeError):
    '''
    This error is raised if a simulation is run in such a way that no timesteps
    will be taken. This error is a distinct type so that it can be safely caught
    and ignored if required, but it is anticipated that most of the time, calling
    sim.run() and not taking any timesteps, would be an inadvertent error.
    '''
    pass