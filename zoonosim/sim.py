'''
Defines the Sim class, zoonosims core class.
'''



import numpy as np
import sciris as sc

from . import options as zno

from . import utils as znu
from . import misc as znm
from . import base as znb
from . import defaults as znd
from . import interventions as zni
from . import analysis as zna
from . import parameters as znpar
from . import population as znpop
from . import immunity as znimm
from . import plotting as znplt
from . import testing as znt

__all__ = ['Sim']

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
        popfile  (str):    if supplied, load the population from this file
        agents   (varies): if supplied, use these pre-generated agents (as an agents object) instead of loading or generating new ones
        version  (str):    if supplied, use default parameters from this version of Covasim instead of the latest
        kwargs   (dict):   additional parameters; passed to ``cv.make_pars()``

    **Examples**::

        sim = zn.Sim()
        sim = zn.Sim(n_farms=100, datafile='my_data.xlsx', label='Sim with data')
    '''

    def __init__(self, pars=None, datafile=None, label=None, simfile=None,
                 popfile=None, agents=None, version=None, **kwargs):

        """

        """
        # Set attributes
        self.label         = label    # The label/name of the simulation
        self.created       = None     # The datetime the sim was created
        self.simfile       = simfile  # The filename of the sim
        self.datafile      = datafile # The name of the data file
        self.popfile       = popfile  # The population file
        self.data          = None     # The actual data
        self.popdict       = agents   # The population dictionary
        self.agents        = None     # Initialize these here so methods that check their length can see they're empty
        self.t             = None     # The current time in the simulation (during execution); outside of sim.step(), its value corresponds to next timestep to be computed
        self.results       = {}       # For storing results
        self.summary       = None     # For storing a summary of the results
        self.initialized   = False    # Whether or not initialization is complete
        self.complete      = False    # Whether a simulation has completed running
        self.results_ready = False    # Whether or not results are ready
        self._default_ver  = version  # Default version of parameters used
        self._legacy_trans = None     # Whether to use the legacy transmission calculation method (slower; for reproducing earlier results)
        self._orig_pars    = None     # Store original parameters to optionally restore at the end of the simulation

        
        # Make default parameters (using values from parameters.py)
        default_pars = znpar.make_pars(version=version) # Start with default pars
        super().__init__(default_pars) # Initialize and set the parameters as attributes

        # Now update everything
        self.set_metadata(simfile) # Set the simulation date and filename
        self.update_pars(pars, **kwargs) # Update the parameters, if provided
        self.load_data(datafile) # Load the data, if provided

                # Process test object pars
        if self.pars['enable_testobjs']: 
            self.process_testobj_pars()
            print("Testing is enabled!")


        return


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
        creating the results structure, initializing the agents, validating the
        layer parameters (which requires the agents), and initializing the interventions.

        Note: to create a population to save for later use, use ``init_infections=False``.
        This will then create a fresh Agents object which other sims can finish
        initializing.

        Args:
            reset (bool): whether or not to reset people even if they already exist
            init_infections (bool): whether to initialize infections (default true so sim is ready, but can't reuse agents then)
            kwargs (dict): passed to init_agents
        '''
        self.t = 0  # The current time index

        self.validate_pars() # Ensure parameters have valid values

        self.set_seed() # Reset the random seed before the population is created

        self.init_variants() # Initialize the variants
 
        self.init_immunity() # initialize information about immunity (if use_waning=True)

        self.init_results() # After initializing the variant, create the results structure

        self.init_agents(reset=reset, init_infections=init_infections, **kwargs)

        self.init_interventions()  # Initialize the interventions...

        self.init_testobjs() # ...and the testing objects...

        self.init_analyzers() # ...and the analyzers...

        

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
            if self.agents is not None: # If agents exist
                layer_keys = self.agents.contacts.keys()
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
        if self.agents is not None:
            pop_keys = set(self.agents.contacts.keys())
            if pop_keys != set(layer_keys): # pragma: no cover
                if not len(pop_keys):
                    errormsg = f'Your population does not have any layer keys, but your simulation does {layer_keys}. If you called zn.Agents() directly, you probably need zn.make_agents() instead.'
                    raise sc.KeyNotFoundError(errormsg)
                else:
                    errormsg = f'Please update your parameter keys {layer_keys} to match population keys {pop_keys}. You may find sim.reset_layer_pars() helpful.'
                    raise sc.KeyNotFoundError(errormsg)

        return


    def validate_pars(self, validate_layers=False):
        '''
        Some parameters can take multiple types; this makes them consistent.

        Args:
            validate_layers (bool): whether to validate layer parameters as well via validate_layer_pars() -- usually yes, except during initialization
        '''

        # Handle population size
        # Can't handle population size here as population hasn't been created yet.
        # recorded_pop_size = self['pop_size']
        # actual_pop_size = 0
        # for key in self.pars['agent_types']:
        #     if key in self.pars['pop_size_by_type']:
        #         actual_pop_size += self.pars['pop_size_by_type'][key]
        #     else:
        #         errormsg = f'You have specified agent_types {self.pars["agent_types"]} but not provided a pop_size_by_type entry for "{key}"'
        #         raise KeyError(errormsg)
        # if recorded_pop_size != actual_pop_size:
        #     if recorded_pop_size != 0:
        #         warnmsg = f'pop_size ({recorded_pop_size}) does not match sum of pop_size_by_type ({actual_pop_size}); using the latter'
        #         znm.warn(warnmsg)
        #     self['pop_size'] = actual_pop_size

        # Handle types


        # Handle start day
        start_day = self['start_day'] # Shorten
        if start_day in [None, 0]: # Use default start day
            start_day = '2025-01-01'
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


        # Handle interventions, analyzers, and variants
        for key in ['interventions', 'analyzers', 'variants']: # Ensure all of them are lists
            self[key] = sc.dcp(sc.tolist(self[key], keepnone=False)) # All of these have initialize functions that run into issues if they're reused
        for i,interv in enumerate(self['interventions']):
            if isinstance(interv, dict): # It's a dictionary representation of an intervention
                self['interventions'][i] = zni.InterventionDict(**interv)
        self['variant_map'] = {int(k):v for k,v in self['variant_map'].items()} # Ensure keys are ints, not strings of ints if loaded from JSON

        # Optionally handle layer parameters
        if validate_layers:
            self.validate_layer_pars()

        # Handle versioning


        # Handle verbose
        if self['verbose'] == 'brief':
            self['verbose'] = -1
        if not sc.isnumber(self['verbose']): # pragma: no cover
            errormsg = f'Verbose argument should be either "brief", -1, or a float, not {type(self["verbose"])} "{self["verbose"]}"'
            raise ValueError(errormsg)

        return

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

        human_dcols = znd.get_default_colors('human') # Get default human colors
        flock_dcols = znd.get_default_colors('flock') # Get default flock colors
        breed_dcols = znd.get_default_colors('breed') # Get default breed colors
        barn_dcols  = znd.get_default_colors('barn')  # Get default barn colors
        water_dcols = znd.get_default_colors('water') # Get default water colors
        misc_dcols = znd.get_default_colors('misc') # Get default misc colors
        sw_dcols = znd.get_default_colors('smartwatch')

        # Flows and cumulative flows

        for key,label in znd.human_flows.items():
            self.results[f'cum_human_{key}'] = init_res(f'Cumulative {label}', color=human_dcols[key])  # Cumulative variables -- e.g. "Cumulative infections"

        for key,label in znd.human_flows.items():
            self.results[f'new_human_{key}'] = init_res(f'Number of new {label}', color=human_dcols[key]) # Flow variables -- e.g. "Number of new infections"

        for key,label in znd.flock_flows.items():
            self.results[f'cum_flock_{key}'] = init_res(f'Cumulative {label}', color=flock_dcols[key])

        for key,label in znd.flock_flows.items():
            self.results[f'new_flock_{key}'] = init_res(f'Number of new {label}', color=flock_dcols[key])

        for breed in znd.default_flock_breeds:
            for key,label in znd.flock_flows.items():
                self.results[f'cum_flock_{breed}_{key}'] = init_res(f'Cumulative {label} ({breed})', color=breed_dcols[breed])

            for key,label in znd.flock_flows.items():
                self.results[f'new_flock_{breed}_{key}'] = init_res(f'Number of new {label} ({breed})', color=breed_dcols[breed])

        for key,label in znd.barn_flows.items():
            self.results[f'cum_barn_{key}'] = init_res(f'Cumulative {label}', color=barn_dcols[key])

        for key,label in znd.barn_flows.items():
            self.results[f'new_barn_{key}'] = init_res(f'Number of new {label}', color=barn_dcols[key])

        for key,label in znd.water_flows.items():
            self.results[f'cum_water_{key}'] = init_res(f'Cumulative {label}', color=water_dcols[key])

        for key,label in znd.water_flows.items():
            self.results[f'new_water_{key}'] = init_res(f'Number of new {label}', color=water_dcols[key])

        # Stock variables

        for key,label in znd.human_stocks.items():
            self.results[f'n_human_{key}'] = init_res(label, color=human_dcols[key])

        for key,label in znd.flock_stocks.items():
            self.results[f'n_flock_{key}'] = init_res(label, color=flock_dcols[key])

        for key,label in znd.barn_stocks.items():
            self.results[f'n_barn_{key}'] = init_res(label, color=barn_dcols[key])

        for key,label in znd.water_stocks.items():
            self.results[f'n_water_{key}'] = init_res(label, color=water_dcols[key])


        # Other variables
        self.results['n_human_imports'] = init_res('Number of imported human infections', scale = True)
        self.results['n_flock_imports'] = init_res('Number of imported flock infections', scale = True)
        self.results['n_barn_imports'] = init_res('Number of imported barn infections', scale = True)
        self.results['n_water_imports'] = init_res('Number of imported water infections', scale = True)

        # Nab levels and population immunity
        self.results['pop_nabs']            = init_res('nab levels in human population', scale=False, color=human_dcols.pop_nabs)
        self.results['pop_protection']      = init_res('Population immunity protection', scale=False, color=human_dcols.pop_protection)
        self.results['pop_symp_protection'] = init_res('Population symptomatic protection', scale=False, color=human_dcols.pop_symp_protection)

        # Testing
        self.results['new_tests']                     = init_res('Number of new tests administered', color=misc_dcols.misc1) # This is specific for tests UNRELATED to the custom testing module
        self.results['cum_tests']                     = init_res('Cumulative tests administered', color=misc_dcols.misc2) # This is specific for tests UNRELATED to the custom testing module
        if self.pars['enable_testobjs']:
            self.results['new_PCR_tests']                 = init_res('Number of new PCR tests administered', color=misc_dcols.misc1)
            self.results['cum_PCR_tests']                 = init_res('Cumulative PCR tests administered', color=misc_dcols.misc2)
            self.results['new_RAT_tests']                 = init_res('Number of new RAT tests administered', color=misc_dcols.misc3)
            self.results['cum_RAT_tests']                 = init_res('Cumulative RAT tests administered', color=misc_dcols.misc4)
            self.results['new_diagnoses_custom']      = init_res('Number of new diagnoses with custom testing module', color=misc_dcols.misc5)
            self.results['cum_diagnoses_custom']      = init_res('Cumulative diagnoses with custom testing module', color=misc_dcols.misc6)

        # Smartwatches
        if self.pars['enable_smartwatches']:
            for key,label in znd.smartwatch_flows.items():
                color = sw_dcols[key]
                self.results[f'new_{key}'] = init_res(f'Number of new {label}', color=color)
                self.results[f'cum_{key}'] = init_res(f'Cumulative {label}', color=color)

        # Handle variants

        nv = self['n_variants']
        self.results['variant'] = {}
        #self.results['variant']['prevalence_by_variant'] = init_res('Prevalence by variant', scale=False, n_variants=nv)
        #self.results['variant']['incidence_by_variant']  = init_res('Incidence by variant', scale=False, n_variants=nv)
        for key,label in znd.human_flows_by_variant.items():
            self.results['variant'][f'cum_human_{key}'] = init_res(f'Cumulative {label}', color=human_dcols[key], n_variants=nv)  # Cumulative variables -- e.g. "Cumulative infections"

        for key,label in znd.flock_flows_by_variant.items():
            self.results['variant'][f'cum_flock_{key}'] = init_res(f'Cumulative {label}', color=flock_dcols[key], n_variants=nv )

        for key,label in znd.barn_flows_by_variant.items():
            self.results['variant'][f'cum_barn_{key}'] = init_res(f'Cumulative {label}', color=barn_dcols[key], n_variants=nv )

        for key,label in znd.water_flows_by_variant.items():
            self.results['variant'][f'cum_water_{key}'] = init_res(f'Cumulative {label}', color=water_dcols[key], n_variants=nv )



        for key,label in znd.human_flows_by_variant.items():
            self.results['variant'][f'new_human_{key}'] = init_res(f'Number of new {label}', color=human_dcols[key], n_variants=nv) # Flow variables -- e.g. "Number of new infections"

        for key,label in znd.flock_flows_by_variant.items():
            self.results['variant'][f'new_flock_{key}'] = init_res(f'Number of new {label}', color=flock_dcols[key], n_variants=nv)

        for key,label in znd.barn_flows_by_variant.items():
            self.results['variant'][f'new_barn_{key}'] = init_res(f'Number of new {label}', color=barn_dcols[key], n_variants=nv)

        for key,label in znd.water_flows_by_variant.items():
            self.results['variant'][f'new_water_{key}'] = init_res(f'Number of new {label}', color=water_dcols[key], n_variants=nv)


        for key,label in znd.human_stocks_by_variant.items():
            self.results['variant'][f'n_human_{key}'] = init_res(label, color=human_dcols[key], n_variants=nv)

        for key,label in znd.flock_stocks_by_variant.items():
            self.results['variant'][f'n_flock_{key}'] = init_res(label, color=flock_dcols[key], n_variants=nv)

        for key,label in znd.barn_stocks_by_variant.items():
            self.results['variant'][f'n_barn_{key}'] = init_res(label, color=barn_dcols[key], n_variants=nv)

        for key,label in znd.water_stocks_by_variant.items():
            self.results['variant'][f'n_water_{key}'] = init_res(label, color=water_dcols[key], n_variants=nv)

        # Populate the rest of the results

        self.results['date'] = self.datevec
        self.results['t']    = self.tvec
        self.results_ready   = False

        return


    def load_population(self, popfile=None, init_agents=True, **kwargs):
        '''
        Load the population dictionary from file -- typically done automatically
        as part of ``sim.initialize()``.

        Supports loading either saved population dictionaries (popdicts, file ending
        .pop by convention), or ready-to-go Agents objects (file ending .agnts by
        convention). Either object an also be supplied directly. Once a population
        file is loaded, it is removed from the Sim object.

        Args:
            popfile (str or obj): if a string, name of the file; otherwise, the popdict or Agents object to load
            init_people (bool): whether to immediately convert the loaded popdict into an initialized Agents object
            kwargs (dict): passed to ``sim.init_people()``
        '''
        # Set the file path if not is provided
        if popfile is None and self.popfile is not None:
            popfile = self.popfile

        # Load the population into the popdict
        self.popdict = znm.load(popfile) # NOTE: I'm honestly not sure if this will work given the changes to how populations are handled.
        if self['verbose']:
            print(f'Loading population from {popfile}')

        if init_agents:
            self.init_agents(**kwargs)

        return


    def init_agents(self, popdict=None, init_infections=False, reset=False, verbose=None, **kwargs):
        '''
        Create the agents.

        Use ``init_infections=False`` for creating a fresh Agents object for use
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
            if self.agents:
                resetstr = ' (resetting agents)' if reset else ' (warning: not resetting sim.agents)'
            print(f'Initializing sim{resetstr} with {self["n_farms"]:0n} farms for {self["n_days"]} days')

        if self.popfile and self.popdict is None: # If there's a popdict, we initialize it via cvpop.make_people()
            self.load_population(init_agents=False)
        
        # Actually make the people
        self.agents = znpop.make_agents(self, reset=reset, verbose=verbose, **kwargs)
        self.agents.initialize(sim_pars=self.pars) # Fully initialize the people
        self.reset_layer_pars(force=False) # Ensure that layer keys match the loaded population

        if self.pars['enable_smartwatches']:
            self.agents.human.init_watches(self.pars['smartwatch_pars'])


        if init_infections:
            self.init_infections(verbose=verbose)

        return self

    def init_interventions(self, popdict=None, init_infections=False, reset=False, verbose=None, **kwargs):
        ''' Initialize and validate the interventions '''

        # Initialization
        if self._orig_pars and 'interventions' in self._orig_pars:
            self['interventions'] = self._orig_pars.pop('interventions') # Restore

        for i,intervention in enumerate(self['interventions']):
            if isinstance(intervention, zni.Intervention):
                intervention.initialize(self)

        # Validation
        # trace_ind = np.nan # Index of the tracing intervention(s)
        # test_ind = np.nan # Index of the tracing intervention(s)
        # for i,intervention in enumerate(self['interventions']):
        #     if isinstance(intervention, (cvi.contact_tracing)):
        #         trace_ind = np.fmin(trace_ind, i) # Find the earliest-scheduled tracing intervention
        #     elif isinstance(intervention, (cvi.test_num, cvi.test_prob)):
        #         test_ind = np.fmax(test_ind, i) # Find the latest-scheduled testing intervention

        # if not np.isnan(trace_ind): # pragma: no cover
        #     warnmsg = ''
        #     if np.isnan(test_ind):
        #         warnmsg = 'Note: you have defined a contact tracing intervention but no testing intervention was found. Unless this is intentional, please define at least one testing intervention.'
        #     elif trace_ind < test_ind:
        #         warnmsg = f'Note: contact tracing (index {trace_ind:.0f}) is scheduled before testing ({test_ind:.0f}); this creates a 1-day delay. Unless this is intentional, please reorder the interentions.'
        #     if warnmsg:
        #         znm.warn(warnmsg)

        return

    def finalize_interventions(self):
        for intervention in self['interventions']:
            if isinstance(intervention, zni.Intervention):
                intervention.finalize(self)

    def init_testobjs(self):
        if self._orig_pars and 'testing' in self._orig_pars: 
            self['testing'] = self._orig_pars.pop('testing')
        
        for testobj in self['testing']:
            if isinstance(testobj, znt.TestObj):
                testobj.init_common_trackers(self)
                testobj.initialize(self)

    def finalize_testobjs(self): 
        for testobj in self['testing']: 
            if isinstance(testobj, znt.TestObj):
                testobj.finalize(self)


    def init_analyzers(self):
        ''' Initialize the analyzers '''
        if self._orig_pars and 'analyzers' in self._orig_pars:
            self['analyzers'] = self._orig_pars.pop('analyzers') # Restore

        for analyzer in self['analyzers']:
            if isinstance(analyzer, zna.Analyzer):
                analyzer.initialize(self)
        return


    def finalize_analyzers(self):
        for analyzer in self['analyzers']:
            if isinstance(analyzer, zna.Analyzer):
                analyzer.finalize(self)


    def init_variants(self):
        ''' Initialize the variants '''
        if self._orig_pars and 'variants' in self._orig_pars:
            self['variants'] = self._orig_pars.pop('variants') # Restore

        for i,variant in enumerate(self['variants']):
            if isinstance(variant, znimm.variant):
                if not variant.initialized:
                    variant.initialize(self)
            else: # pragma: no cover
                errormsg = f'Variant {i} ({variant}) is not a zn.variant object; please create using zn.variant()'
                raise TypeError(errormsg)

        len_pars = len(self['variant_pars'])
        len_map = len(self['variant_map'])
        assert len_pars == len_map, f"variant_pars and variant_map must be the same length, but they're not: {len_pars} ≠ {len_map}"
        self['n_variants'] = len_pars # Each variant has an entry in variant_pars

        return


    def init_immunity(self, create=False):
        ''' Initialize immunity matrices and precompute nab waning for each variant '''
        znimm.init_immunity(self, create=create)
        return



    def init_infections(self, force=False, verbose=None):
        if verbose is None:
            verbose = self['verbose']
        
        human_pop_size = self.pars['pop_size_by_type']['human']
        flock_pop_size = self.pars['pop_size_by_type']['flock']
        barn_pop_size = self.pars['pop_size_by_type']['barn']
        water_pop_size = self.pars['pop_size_by_type']['water']


        requested_human_exposures = self.pars['initial_conditions']['human']
        requested_flock_exposures = self.pars['initial_conditions']['flock']
        requested_barn_contaminations = self.pars['initial_conditions']['barn']
        requested_water_contaminations = self.pars['initial_conditions']['water']


        if (self.agents.count('exposed') == 0 and self.agents.count('infectious')==0) or force:
            if requested_human_exposures > 0:
                if human_pop_size >= requested_human_exposures:
                    inds = znu.choose(human_pop_size, requested_human_exposures)
                    self.agents.infect_type('human', inds, update=False)
                else:
                    errormsg = (f'requested number of exposed humans ({requested_human_exposures}) '
                    f'is greater than the human population size ({human_pop_size})')
                    ValueError(errormsg)
            if requested_flock_exposures > 0:
                if flock_pop_size >= requested_flock_exposures:
                    inds = znu.choose(flock_pop_size, requested_flock_exposures)
                    self.agents.infect_type('flock', inds, update=False)
                else:
                    errormsg = (f'requested number of exposed flocks ({requested_flock_exposures}) '
                                f'is greater than the flock population ({flock_pop_size})')
                    ValueError(errormsg)
            if requested_barn_contaminations > 0:
                if barn_pop_size >= requested_barn_contaminations:
                    inds = znu.choose(barn_pop_size, requested_barn_contaminations)
                    self.agents.infect_type('barn', inds, update = False)
                else:
                    errormsg = (f'requested number of contaminated barns ({requested_barn_contaminations}) '
                                f'is greater than the barn population ({barn_pop_size})')
                    ValueError(errormsg)
            if requested_water_contaminations > 0:
                if water_pop_size >= requested_water_contaminations:
                    inds = znu.choose(water_pop_size, requested_water_contaminations)
                    self.agents.infect_type('water', inds, update = False)
                else:
                    errormsg = (f'requested number of contaminated waterbodies ({requested_water_contaminations}) '
                                f'is greater than the waterbody population ({water_pop_size})')
                    ValueError(errormsg)
            
            self.agents.update_states_from_subrosters()
        elif verbose:
            n_exposed = self.agents.count('exposed')
            n_infectious = self.agents.count('infectious')
            print(f'agents are already initialized with {n_exposed} exposed agents and {n_infectious} infectious agents')

        return


    def step(self):
        '''
        Step the simulation forward in time. Usually, the user would use sim.run()
        rather than calling sim.step() directly.
        '''

        # Set the time and if we have reached the end of the simulation, then do nothing
        if self.complete:
            raise AlreadyRunError('Simulation already complete (call sim.initialize() to re-run)')

        t = self.t

        # If it's the first timestep, infect people
        if t == 0:
            self.init_infections(verbose=False)

        # Perform initial operations
        #self.rescale() # Check if we need to rescale
        agents = self.agents # Shorten this for later use

        agents.update_states_pre(t=t) # Update the state of everyone and count the flows. This isn't infecting people nor updating their SEIR's. The date of infection seems to be pre-assigned. 
        contacts = agents.update_contacts() # Compute new contacts. For dynamic contacts. 
        hosp_max = agents.type_count('human', 'severe')   > self['n_beds_hosp'] if self['n_beds_hosp'] is not None else False # Check for acute bed constraint.
        
        # Randomly infect some people (imported infections)
        if self['n_imports']['human']>0:
            n_human_imports = znu.poisson(self['n_imports']['human']) # imported human cases
            if n_human_imports>0:
                human_inds = znu.choose(max_n=len(self.agents.human), n=n_human_imports)
                self.agents.infect_type('human', human_inds)
                self.results['n_human_imports'][t] += n_human_imports
        if self['n_imports']['flock']>0:
            n_flock_imports = znu.poisson(self['n_imports']['flock']) # imported flock cases
            if n_flock_imports>0:
                flock_inds = znu.choose(max_n=len(self.agents.flock), n=n_flock_imports)
                self.agents.infect_type('flock', flock_inds)
                self.results['n_flock_imports'][t] += n_flock_imports
        if self['n_imports']['barn']>0:
            n_barn_imports = znu.poisson(self['n_imports']['barn']) # imported barn contaminations
            if n_barn_imports>0:
                barn_inds = znu.choose(max_n=len(self.agents.barn), n=n_barn_imports)
                self.agents.infect_type('barn', barn_inds)
                self.results['n_barn_imports'][t] += n_barn_imports
        if self['n_imports']['water']>0:
            n_water_imports = znu.poisson(self['n_imports']['water']) # imported water contaminations
            if n_water_imports>0:
                water_inds = znu.choose(max_n = len(self.agents.water), n=n_water_imports)
                self.agents.infect_type('water', water_inds)
                self.results['n_water_imports'][t] += n_water_imports

        # Add variants
        for variant in self['variants']:
            if isinstance(variant, znimm.variant):
                variant.apply(self)

        # Sent smartwatch alerts and update alert histories
        if self.pars['enable_smartwatches']:
            self.agents.human.watches.send_alerts_baseline(self)
            self.agents.human.watches.send_alerts_infected(self)
            self.agents.human.watches.update_alert_histories(self)

        # Compute viral loads in humans
        human_viral_load = self.agents.update_human_viral_loads(t=t)


        # Compute infection levels in flocks
        flock_infection_levels = self.agents.update_flock_infection_levels()

        
        # Compute modifiers for barns
        # In principle the relative transmission and susceptibility of barns should vary based on temperature and humidity
        # For now, I'm just setting them to 1
        barn_modifiers = np.repeat(1.0, len(agents.barn)) # TODO: Implement barn modifiers based on temperature and humidity

        # Compute modifiers for water
        # In principle the relative transmission and susceptibility of water should vary based on temperature. 
        # For now, I'm just setting them to 1
        water_modifiers = np.repeat(1.0, len(agents.water)) # TODO: Implement water modifiers based on temperature

        # Set modifiers for all agent types
        # NOTE: Currently barn and water have no modifiers, I'm setting them to 1.0 for now.
        misc_modifiers = np.concatenate((human_viral_load, flock_infection_levels, barn_modifiers, water_modifiers)) 


                # Toggle testing. Background ILI and symptoms are only simulated for testing, so they are enclosed here too.
        if self.pars['enable_testobjs']:

            # Determine background ILI infections. It is justifiable to call this before testing, because ILI should exist even on the first day of testing.
            znt.infect_ILI(self)

            # Apply testing and record diagnoses.
            for i, testobj in enumerate(self['testing']):

                testobj(self)
                
                # Interface with quarantine behaviour.
                self.agents.human.date_pos_test[testobj.date_pos_test == t] = t  # Update date_pos_test with people who received a test today that will return positive
                self.agents.human.date_diagnosed[testobj.date_positive == t] = t  # Update date_diagnosed with people who received at least one positive test today
                
                self.results['new_diagnoses_custom'][t] += sum(testobj.date_positive == t)
                self.results['cum_diagnoses_custom'][t] += sum(self.results['new_diagnoses_custom'][:t])

        # Apply interventions
        for i,intervention in enumerate(self['interventions']):
            intervention(self) # If it's a function, call it directly


        # Implement state changes relating to quarantine and isolation/diagnosis; compute associated statistics
        agents.update_states_post()

        # Shorten useful parameters
        nv = self['n_variants'] # Shorten number of variants
        sus = agents.susceptible
        symp = agents.symptomatic
        quar = agents.quarantined

        prel_trans = agents.get_rel_trans()
        prel_sus   = agents.get_rel_sus()





        # Iterate through n_variants to calculate infections. The meat of the simulation. 
        for variant in range(nv):

            # Check immunity
            
            znimm.check_immunity(agents.human, variant) # NOTE: Only affects human agents

            # The 'wild' variant defines the baseline for all these values
            asymp_factor = np.repeat([self['variant_pars']['wild']['human']['rel_asymp_fact'], 1.0], 
                                  [len(agents.human), (len(agents.flock) + len(agents.barn) + len(agents.water))])
            human_rel_beta = self['variant_pars']['wild']['human']['rel_beta']
            flock_rel_beta = self['variant_pars']['wild']['flock']['rel_beta']
            barn_rel_beta  = self['variant_pars']['wild']['barn']['rel_beta']
            water_rel_beta = self['variant_pars']['wild']['human']['rel_beta']

            if variant:
                variant_label = self.pars['variant_map'][variant]
                asymp_factor *= np.repeat([self['variant_pars'][variant_label]['human']['rel_asymp_fact'], 1.0], 
                                  [len(agents.human), (len(agents.flock) + len(agents.barn) + len(agents.water))])
                human_rel_beta *= self['variant_pars'][variant_label]['human']['rel_beta']
                flock_rel_beta *= self['variant_pars'][variant_label]['flock']['rel_beta']
                barn_rel_beta  *= self['variant_pars'][variant_label]['barn']['rel_beta']
                water_rel_beta *= self['variant_pars'][variant_label]['human']['rel_beta']



            beta = np.repeat([znd.default_float(self['beta']['human']*human_rel_beta),
                              znd.default_float(self['beta']['flock']*flock_rel_beta), 
                              znd.default_float(self['beta']['barn']*barn_rel_beta), 
                              znd.default_float(self['beta']['water']*water_rel_beta)], 
                              [len(agents.human), len(agents.flock), len(agents.barn), len(agents.water)])
            

            for lkey, layer in contacts.items():
                p1 = layer['p1']
                p2 = layer['p2']
                betas = layer['beta']

                # Compute relative transmission and susceptibility
                inf_variant = agents.infectious * (agents.infectious_variant == variant) # TODO: move out of loop?
                sus_imm = agents.sus_imm[variant,:]
                quar_factor = znd.default_float(self['quar_factor'][lkey]) # Ex: 0.2. Probably the effect on beta of quarantining. 
                beta_layer  = znd.default_float(self['beta_layer'][lkey]) # A scalar; beta for the layer. Ex: 1.0. 

                rel_trans, rel_sus = znu.compute_trans_sus(prel_trans, prel_sus, inf_variant, sus, beta_layer, misc_modifiers, symp, quar, asymp_factor, quar_factor, sus_imm)
                

                # Calculate actual transmission
                pairs = [[p1,p2]] 
                for p1,p2 in pairs:
                    source_inds, target_inds = znu.compute_infections(beta, p1, p2, betas, rel_trans, rel_sus)  # Calculate transmission!
                    agents.infect(inds=target_inds, hosp_max=hosp_max, source=source_inds, layer=lkey, variant=variant)  # Actually infect people

        ##### CALCULATE STATISTICS #####
        
        # Update counts for this time step: Human Stocks
        for key in znd.human_stocks:
            self.results[f'n_human_{key}'][t] = agents.human.count(key)
        
        # Update counts for this time step: Flock stocks
        for key in znd.flock_stocks:
            self.results[f'n_flock_{key}'][t] = agents.flock.count(key)

        # Update counts for this time step: Barn stocks
        for key in znd.barn_stocks:
            self.results[f'n_barn_{key}'][t] = agents.barn.count(key)

        # Update counts for this time step: Water Stocks
        for key in znd.water_stocks:
            self.results[f'n_water_{key}'][t] = agents.water.count(key)

        # Update stocks_by_variant: Human
        for key in znd.human_stocks_by_variant:
            for variant in range(nv):
                self.results['variant'][f'n_human_{key}'][variant, t] = agents.human.count_by_variant(key, variant)

        # Update stocks_by_variant: Flock
        for key in znd.flock_stocks_by_variant:
            for variant in range(nv):
                self.results['variant'][f'n_flock_{key}'][variant, t] = agents.flock.count_by_variant(key, variant)

        # Update stocks_by_variant: Barn
        for key in znd.barn_stocks_by_variant:
            for variant in range(nv):
                self.results['variant'][f'n_barn_{key}'][variant, t] = agents.barn.count_by_variant(key, variant)
        
        # Update stocks_by_variant: Water
        for key in znd.water_stocks_by_variant:
            for variant in range(nv):
                self.results['variant'][f'n_water_{key}'][variant, t] = agents.water.count_by_variant(key, variant)
        
        # Update counts for this time step: Human flows
        for key in znd.human_flows:
            self.results[f'new_human_{key}'][t] = agents.human.flows[f'new_{key}']
        
        if self.pars['enable_smartwatches']:
            for key in znd.smartwatch_flows:
                self.results[f'new_{key}'][t] = agents.human.flows[f'new_{key}']

        # Update counts for this time step: Flock flows
        for key in znd.flock_flows:
            self.results[f'new_flock_{key}'][t] = agents.flock.flows[f'new_{key}']
            for breed in znd.flock_breeds:
                self.results[f'new_flock_{key}_{breed}'][t] = agents.flock.flows_breed[f'new_{key}'][breed]

        # Update counts for this time step: Barn flows
        for key in znd.barn_flows:
            self.results[f'new_barn_{key}'][t] = agents.barn.flows[f'new_{key}']

        # Update counts for this time step: Water flows
        for key in znd.water_flows:
            self.results[f'new_water_{key}'][t] = agents.water.flows[f'new_{key}']

        # Update counts for this time step: Human flows by variant
        for key in znd.human_flows_by_variant:
            for variant in range(nv):
                self.results['variant'][f'new_human_{key}'][variant, t] = agents.human.flows_variant[f'new_{key}'][variant]
        
        # Update counts for this time step: Flock flows by variant
        for key in znd.flock_flows_by_variant:
            for variant in range(nv):
                self.results['variant'][f'new_flock_{key}'][variant, t] = agents.flock.flows_variant[f'new_{key}'][variant]

        
        # Update counts for this time step: Barn flows by variant
        for key in znd.barn_flows_by_variant:
            for variant in range(nv):
                self.results['variant'][f'new_barn_{key}'][variant, t] = agents.barn.flows_variant[f'new_{key}'][variant]

        # Update counts for this time step: Water flows by variant
        for key in znd.water_flows_by_variant:
            for variant in range(nv):
                self.results['variant'][f'new_water_{key}'][variant, t] = agents.water.flows_variant[f'new_{key}'][variant]

        # # Update nab and immunity for this time step
        
        has_nabs = znu.true(agents.human.peak_nab)
        if len(has_nabs):
            znimm.update_nab(agents.human, inds=has_nabs)

        inds_alive = znu.false(agents.human.dead)
        self.results['pop_nabs'][t]            = np.sum(agents.human.nab[inds_alive[znu.true(agents.human.nab[inds_alive])]])/len(inds_alive)
        self.results['pop_protection'][t]      = np.nanmean(agents.human.sus_imm)
        self.results['pop_symp_protection'][t] = np.nanmean(agents.human.symp_imm)

        ##### DONE CALCULATING STATISTICS #####

        # Apply analyzers (note this comes with Covasim) -- same syntax as interventions
        for i, analyzer in enumerate(self['analyzers']):
            analyzer(self)

        # Tidy up
        self.t += 1
        if self.t == self.npts:
            self.complete = True
        return


    
    def process_testobj_pars(self):
        '''
        Build test objects using supplied test object parameter dictionary. The dictionary must adhere to documented format. 
        Then, assign the test objects to self.pars['testing'].
        Alternatively, test objects can be built outside and externally supplied to sim.pars['testing'] parameter. 
        '''
        test_params = self.pars['testobjs']
        options = ['PCR_disc', 'RAT_disc'] # TESTING
        testobjs = []
        if not(test_params is None):
            for testname, pars in test_params.items(): 
                if not(testname in options): 
                    raise RuntimeError(f'Test parameter dictionary specifies a non-existent test object. Options are {options}')
                
                # These are the only possible test objects
                if testname == 'RAT_disc': 
                    testobjs.append(znt.RAT_disc(**pars))  # Unpack dictionary parameters as keyword arguments. Keys need to match argument name. 
                elif testname == 'PCR_disc': 
                    testobjs.append(znt.PCR_disc(**pars))
            self.pars['testing'] = testobjs
        else: 
            if len(self.pars['testing']) == 0: 
                print("Warning: Test objects enabled but no test object parameters, or already-built test objects, were supplied")

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
        if self.agents and self.agents.t not in [self.t, self.t-1]: # Depending on how the sim stopped, either of these states are possible 
            errormsg = f'The simulation has been run independently from the agents (t={self.t}, agents.t={self.agents.t}): if this is intentional, manually set sim.agents.t = sim.t. Remember to save the agents object before running the sim.'
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
        
        # Calculate cumulative results for smartwatches
        if self.pars['enable_smartwatches']:
            for key in znd.smartwatch_flows.keys():
                self.results[f'cum_{key}'][:] = np.cumsum(self.results[f'new_{key}'][:], axis=0)

        # Calculate cumulative results: Human
        for key in znd.human_flows:
            self.results[f'cum_human_{key}'][:] = np.cumsum(self.results[f'new_human_{key}'][:], axis=0)
        for key in znd.human_flows_by_variant:
            for variant in range(self['n_variants']):
                self.results['variant'][f'cum_human_{key}'][variant, :] = np.cumsum(self.results['variant'][f'new_human_{key}'][variant, :], axis=0)
        
        # Calculate cumulative results: Flock
        for key in znd.flock_flows:
            self.results[f'cum_flock_{key}'][:] = np.cumsum(self.results[f'new_flock_{key}'][:], axis=0)
            for breed in znd.flock_breeds:
                self.results[f'cum_flock_{key}_{breed}'][:] = np.cumsum(self.results[f'new_flock_{key}_{breed}'][:], axis=0)
        for key in znd.flock_flows_by_variant:
            for variant in range(self['n_variants']):
                self.results['variant'][f'cum_flock_{key}'][variant, :] = np.cumsum(self.results['variant'][f'new_flock_{key}'][variant, :], axis=0)
        
        # Calculate cumulative results: Barn
        for key in znd.barn_flows:
            self.results[f'cum_barn_{key}'][:] = np.cumsum(self.results[f'new_barn_{key}'][:], axis=0)
        for key in znd.barn_flows_by_variant:
            for variant in range(self['n_variants']):
                self.results['variant'][f'cum_barn_{key}'][variant, :] = np.cumsum(self.results['variant'][f'new_barn_{key}'][variant, :], axis=0)
        
        # Calculate cumulative results: Water
        for key in znd.water_flows:
            self.results[f'cum_water_{key}'][:] = np.cumsum(self.results[f'new_water_{key}'][:], axis=0)
        for key in znd.water_flows_by_variant:
            for variant in range(self['n_variants']):
                self.results['variant'][f'cum_water_{key}'][variant, :] = np.cumsum(self.results['variant'][f'new_water_{key}'][variant, :], axis=0)


        # Finalize interventions and analyzers
        self.finalize_interventions()
        #self.finalize_testobjs() # TODO: Ritchie toggle. 
        self.finalize_analyzers()

        # Final settings
        self.results_ready = True # Set this first so self.summary() knows to print the results
        self.t -= 1 # During the run, this keeps track of the next step; restore this be the final day of the sim

        # Perform calculations on results
        self.compute_results(verbose=verbose) # Calculate the rest of the results
        self.results = sc.objdict(self.results) # Convert results to a odicts/objdict to allow e.g. sim.results.diagnoses

        if restore_pars and self._orig_pars:
            preserved = ['analyzers', 'interventions', 'testing']
            orig_pars_keys = list(self._orig_pars.keys()) # Get a list of keys so we can iterate over them
            for key in orig_pars_keys:
                if key not in preserved:
                    self.pars[key] = self._orig_pars.pop(key) # Restore everything except for the analyzers, interventions, and testing

        # Optionally print summary output
        if verbose: # Verbose is any non-zero value
            if verbose>0: # Verbose is any positive number
                self.summarize() # Print medium-length summary of the sim
            else:
                self.brief() # Print brief summary of the sim

        return


    def compute_results(self, verbose=None):
         ''' Perform final calculations on the results '''
        # self.compute_states()
        # self.compute_yield()
        # self.compute_doubling()
        # self.compute_r_eff()
        # self.compute_summary()
         return


    # def compute_states(self):
    #     '''
    #     Compute prevalence, incidence, and other states. Prevalence is the current
    #     number of infected people divided by the number of people who are alive.
    #     Incidence is the number of new infections per day divided by the susceptible
    #     population. Also calculates the number of people alive, the number preinfectious,
    #     the number removed, and recalculates susceptibles to handle scaling.
    #     '''
    #     res = self.results
    #     count_recov = 1-self['use_waning'] # If waning is on, don't count recovered people as removed
    #     self.results['n_alive'][:]         = self.scaled_pop_size - res['cum_deaths'][:] # Number of people still alive
    #     self.results['n_naive'][:]         = self.scaled_pop_size - res['cum_deaths'][:] - res['n_recovered'][:] - res['n_exposed'][:] # Number of people naive
    #     self.results['n_susceptible'][:]   = res['n_alive'][:] - res['n_exposed'][:] - count_recov*res['cum_recoveries'][:] # Recalculate the number of susceptible people, not agents
    #     self.results['n_preinfectious'][:] = res['n_exposed'][:] - res['n_infectious'][:] # Calculate the number not yet infectious: exposed minus infectious
    #     self.results['n_removed'][:]       = count_recov*res['cum_recoveries'][:] + res['cum_deaths'][:] # Calculate the number removed: recovered + dead
    #     self.results['prevalence'][:]      = res['n_exposed'][:]/res['n_alive'][:] # Calculate the prevalence
    #     self.results['incidence'][:]       = res['new_infections'][:]/res['n_susceptible'][:] # Calculate the incidence
    #     self.results['frac_vaccinated'][:] = res['n_vaccinated'][:]/res['n_alive'][:] # Calculate the fraction vaccinated

    #     self.results['variant']['incidence_by_variant'][:] = np.einsum('ji,i->ji',res['variant']['new_infections_by_variant'][:], 1/res['n_susceptible'][:]) # Calculate the incidence
    #     self.results['variant']['prevalence_by_variant'][:] = np.einsum('ji,i->ji',res['variant']['new_infections_by_variant'][:], 1/res['n_alive'][:])  # Calculate the prevalence

    #     return


    # def compute_yield(self):
    #     '''
    #     Compute test yield -- number of positive tests divided by the total number
    #     of tests, also called test positivity rate. Relative yield is with respect
    #     to prevalence: i.e., how the yield compares to what the yield would be from
    #     choosing a person at random from the population.
    #     '''
    #     # Absolute yield
    #     res = self.results
    #     inds = cvu.true(res['new_tests'][:]) # Pull out non-zero numbers of tests
    #     self.results['test_yield'][inds] = res['new_diagnoses'][inds]/res['new_tests'][inds] # Calculate the yield

    #     # Relative yield
    #     inds = cvu.true(res['n_infectious'][:]) # To avoid divide by zero if no one is infectious
    #     denom = res['n_infectious'][inds] / (res['n_alive'][inds] - res['cum_diagnoses'][inds]) # Alive + undiagnosed people might test; infectious people will test positive
    #     self.results['rel_test_yield'][inds] = self.results['test_yield'][inds]/denom # Calculate the relative yield
    #     return


    # def compute_doubling(self, window=3, max_doubling_time=30):
    #     '''
    #     Calculate doubling time using exponential approximation -- a more detailed
    #     approach is in utils.py. Compares infections at time t to infections at time
    #     t-window, and uses that to compute the doubling time. For example, if there are
    #     100 cumulative infections on day 12 and 200 infections on day 19, doubling
    #     time is 7 days.

    #     Args:
    #         window (float): the size of the window used (larger values are more accurate but less precise)
    #         max_doubling_time (float): doubling time could be infinite, so this places a bound on it

    #     Returns:
    #         doubling_time (array): the doubling time results array
    #     '''

    #     cum_infections = self.results['cum_infections'].values
    #     infections_now = cum_infections[window:]
    #     infections_prev = cum_infections[:-window]
    #     use = (infections_prev > 0) & (infections_now > infections_prev)
    #     doubling_time = window * np.log(2) / np.log(infections_now[use] / infections_prev[use])
    #     self.results['doubling_time'][:] = np.nan
    #     self.results['doubling_time'][window:][use] = np.minimum(doubling_time, max_doubling_time)
    #     return self.results['doubling_time'].values

    def plot(self, *args, **kwargs):
        '''
        Plot the results of a single simulation.

        Args:
            to_plot      (dict): Dict of results to plot; see get_default_plots() for structure
            do_save      (bool): Whether or not to save the figure
            fig_path     (str):  Path to save the figure
            fig_args     (dict): Dictionary of kwargs to be passed to ``pl.figure()``
            plot_args    (dict): Dictionary of kwargs to be passed to ``pl.plot()``
            scatter_args (dict): Dictionary of kwargs to be passed to ``pl.scatter()``
            axis_args    (dict): Dictionary of kwargs to be passed to ``pl.subplots_adjust()``
            legend_args  (dict): Dictionary of kwargs to be passed to ``pl.legend()``; if show_legend=False, do not show
            date_args    (dict): Control how the x-axis (dates) are shown (see below for explanation)
            show_args    (dict): Control which "extras" get shown: uncertainty bounds, data, interventions, ticks, the legend; additionally, "outer" will show the axes only on the outer plots
            style_args   (dict): Dictionary of kwargs to be passed to Matplotlib; options are dpi, font, fontsize, plus any valid key in ``pl.rcParams``
            n_cols       (int):  Number of columns of subpanels to use for subplot
            font_size    (int):  Size of the font
            font_family  (str):  Font face
            grid         (bool): Whether or not to plot gridlines
            commaticks   (bool): Plot y-axis with commas rather than scientific notation
            setylim      (bool): Reset the y limit to start at 0
            log_scale    (bool): Whether or not to plot the y-axis with a log scale; if a list, panels to show as log
            do_show      (bool): Whether or not to show the figure
            colors       (dict): Custom color for each result, must be a dictionary with one entry per result key in to_plot
            sep_figs     (bool): Whether to show separate figures for different results instead of subplots
            fig          (fig):  Handle of existing figure to plot into
            ax           (axes): Axes instance to plot into
            kwargs       (dict): Parsed among figure, plot, scatter, date, and other settings (will raise an error if not recognized)

        The optional dictionary "date_args" allows several settings for controlling
        how the x-axis of plots are shown, if this axis is dates. These options are:

            - ``as_dates``:   whether to format them as dates (else, format them as days since the start)
            - ``dateformat``: string format for the date (if not provided, choose based on timeframe)
            - ``rotation``:   whether to rotate labels
            - ``start``:      the first day to plot
            - ``end``:        the last day to plot
            - ``outer``:      only show the date labels on the outer (bottom) plots

        The ``show_args`` dictionary allows several other formatting options, such as:

            - ``tight``:    use tight layout for the figure (default false)
            - ``maximize``: try to make the figure full screen (default false)
            - ``outer``:    only show outermost (bottom) date labels (default false)

        Date, show, and other arguments can also be passed directly, e.g. ``sim.plot(tight=True)``.

        For additional style options, see ``cv.options.with_style()``, which is the
        final refuge of arguments that are not picked up by any of the other parsers,
        e.g. ``sim.plot(**{'ytick.direction':'in'})``.

        Returns:
            fig: Figure handle

        **Examples**::

            sim = cv.Sim().run()
            sim.plot() # Default plotting
            sim.plot('overview') # Show overview
            sim.plot('overview', maximize=True, outer=True, rotation=15) # Make some modifications to make plots easier to see
            sim.plot(style='seaborn-whitegrid') # Use a built-in Matplotlib style
            sim.plot(style='simple', font='Rosario', dpi=200) # Use the other house style with several customizations

        | New in version 2.1.0: argument passing, date_args, and mpl_args
        | New in version 3.1.2: updated date arguments; mpl_args renamed style_args
        '''
        fig = znplt.plot_sim(sim=self, *args, **kwargs)
        return fig


    def plot_result(self, key, *args, **kwargs):
        '''
        Simple method to plot a single result. Useful for results that aren't
        standard outputs. See sim.plot() for explanation of other arguments.

        Args:
            key (str): the key of the result to plot

        Returns:
            fig: Figure handle

        **Example**::

            sim = cv.Sim().run()
            sim.plot_result('r_eff')
        '''
        fig = znplt.plot_result(sim=self, key=key, *args, **kwargs)
        return fig

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

        summary = sc.objdict()
        for key in self.result_keys():
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

            sim = zn.Sim(label='Example sim', verbose=0) # Set to run silently
            sim.run() # Run the sim
            sim.summarize() # Print medium-length summary of the sim
            sim.summarize(t=24, full=True) # Print a "slice" of all sim results on day 24
        '''
        # Compute the summary
        summary = self.compute_summary(full=full, t=t, update=False, output=True)

        # Construct the output string
        if sep is None: sep = zno.options.sep # Default separator
        labelstr = f' "{self.label}"' if self.label else ''
        string = f'Simulation{labelstr} summary:\n'
        for key in self.result_keys():
            if full or key.startswith('cum_'):
                val = np.round(summary[key])
                string += f'   {val:10,.0f} {self.results[key].name.lower()}\n'.replace(',', sep) # Use replace since it's more flexible

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

            sim = zn.Sim(label='Example sim', verbose=0) # Set to run silently
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

            sim = zn.Sim(label='Example sim', verbose=0) # Set to run silently
            sim.run() # Run the sim
            sim.brief() # Prints one-line output
        '''
        string = self._brief()
        if not output:
            print(string)
        else:
            return string
        
class AlreadyRunError(RuntimeError):
    '''
    This error is raised if a simulation is run in such a way that no timesteps
    will be taken. This error is a distinct type so that it can be safely caught
    and ignored if required, but it is anticipated that most of the time, calling
    sim.run() and not taking any timesteps, would be an inadvertent error.
    '''
    pass