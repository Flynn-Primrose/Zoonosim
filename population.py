'''
Defines functions for creating agents
'''

import numpy as np # Needed for a few things not provided by pl
import sciris as sc

from . import misc as znm
from . import rosters as znr
from . import defaults as znd
from . import utils as znu

__all__ = ['make_agents']

def make_agents(sim, popdict=None, reset = False, **kwargs):
    '''
    Create agents for the simulation.

    Args:
        sim             (Sim)  : the simulation object; population parameters are taken from the sim object
        popdict         (any)  : either None, pop.popdict, popfile, or People object
        reset           (bool) : whether to force population creation even if self.popdict/self.people exists
        kwargs          (dict) : passed to make_randpop() or make_synthpop()

    Returns:
        agents (Agents): an instance of the Agents class containing the agents and their contacts
    '''

        # If a people object or popdict is supplied, use it
    if sim.agents and not reset:
        sim.agents.initialize(sim_pars=sim.pars)
        return sim.agents # If it's already there, just return
    elif sim.popdict and popdict is None:
        popdict = sim.popdict # Use stored one
        sim.popdict = None # Once loaded, remove
    else:
        popdict = make_popdict(sim, **kwargs) # Create a new one

    validate_popdict(popdict, sim.pars) # Check the popdict is valid

    agents = znr.Agents(sim, uid=popdict['uid'], agent_type=popdict['agent_type'], contacts=popdict['contacts']) # Create the agents object

    return agents

def validate_popdict(popdict, pars, verbose=True):
    '''
    Check that the popdict is the correct type, has the correct keys, and has
    the correct length
    '''

    # Check it's the right type
    try:
        popdict.keys() # Although not used directly, this is used in the error message below, and is a good proxy for a dict-like object
    except Exception as E:
        errormsg = f'The popdict should be a dictionary or zn.Agents object, but instead is {type(popdict)}'
        raise TypeError(errormsg) from E

    # Check keys and lengths
    required_keys = ['uid', 'agent_type']
    popdict_keys = popdict.keys()
    for key in required_keys:

        if key not in popdict_keys:
            errormsg = f'Could not find required key "{key}" in popdict; available keys are: {sc.strjoin(popdict.keys())}'
            sc.KeyNotFoundError(errormsg)

        isnan = np.isnan(popdict[key]).sum()
        if isnan:
            errormsg = f'Population not fully created: {isnan:,} NaNs found in {key}. This can be caused by calling zn.Agents() instead of zn.make_agents().'
            raise ValueError(errormsg)

    if ('contacts' not in popdict_keys) and (not hasattr(popdict, 'contacts')) and verbose:
        warnmsg = 'No contacts found. Please remember to add contacts before running the simulation.'
        znm.warn(warnmsg)

    return

def make_popdict(sim, **kwargs):
    '''
    Create a population dictionary for the simulation.

    Args:
        sim             (Sim)  : the simulation object; population parameters are taken from the sim object
        kwargs          (dict) : passed to make_contacts()

    Returns:
        popdict (dict): a dictionary containing the population data
    '''

    n_farms = sim.pars['n_farms']
    n_water = round(n_farms * pop_pars['avg_water_per_farm']) # Number of water sources

    # Set pop_pars, these are required parameters for creating the population. For the most part they
    # define the means of the distributions used to create the population. The default values are set in defaults.py
    # and can be overridden by passing them in as kwargs.
    # The default values are currently just dummy values and should be changed to reflect the actual population parameters.
    if 'pop_pars' in kwargs:
        pop_pars = kwargs['pop_pars']
    else:
        pop_pars = znd.default_pop_pars

    farm_id = np.arange(n_farms) # Farm IDs
    n_barns = znu.n_poisson(pop_pars['avg_barns_per_farm'], n_farms) # Number of barns per farm
    n_occupied_barns = np.zeros(n_farms) # Number of occupied barns per farm
    n_workers = np.zeros(n_farms) # Number of workers per farm
    
    water_id = np.arange(n_water) # Water IDs

    for farm in range(n_farms):
        n_occupied_barns[farm] = sum(znu.n_binomial(pop_pars['avg_barn_occupancy'], n_barns[farm])) # Number of occupied barns per farm
        n_workers[farm] = sum(znu.n_binomial(pop_pars['avg_workers_per_barn'], n_barns[farm])) # Number of workers per barn


    water_id = znu.choose_r(n_water, n_farms) # Water IDs





    return popdict
