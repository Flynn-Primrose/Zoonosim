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
        popdict         (any)  : either None, pop.popdict, popfile, or Agents object
        reset           (bool) : whether to force population creation even if self.popdict/self.agents exists
        kwargs          (dict) : passed to make_randpop() or make_synthpop()

    Returns:
        agents (Agents): an instance of the Agents class containing the agents and their contacts
    '''

        # If a agents object or popdict is supplied, use it
    if sim.agents and not reset:
        sim.agents.initialize(sim_pars=sim.pars)
        return sim.agents # If it's already there, just return
    elif sim.popdict and popdict is None:
        popdict = sim.popdict # Use stored one
        sim.popdict = None # Once loaded, remove
    else:
        popdict = make_popdict(sim, **kwargs) # Create a new one

    validate_popdict(popdict, sim.pars) # Check the popdict is valid

    humans = make_humans(sim.pars, popdict['human_uids'])
    flocks = make_flocks(sim.pars, popdict['flock_uids'], popdict['flock2barn'])
    barns = make_barns(sim.pars, popdict['barn_uids'], popdict['flock2barn'], popdict['barn2water'])
    water = make_water(sim.pars, popdict['water_uids'], popdict['barn2water'])

    contacts = make_contacts(popdict['contactdict'])

    agents = znr.Agents(sim.pars, 
                        uid = popdict['uid'], 
                        agent_type = popdict['agent_type'], 
                        humans = humans,
                        flocks = flocks,
                        barns = barns,
                        water = water,
                        contacts = contacts)

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
    required_keys = ['uid',
                     'agent_type',
                     'human_uids',
                     'barn_uids',
                     'flock_uids',
                     'water_uids',
                     'contactdict',
                     'barn2water',
                     'flock2barn',
                     ]
    popdict_keys = popdict.keys()
    for key in required_keys:

        if key not in popdict_keys:
            errormsg = f'Could not find required key "{key}" in popdict; available keys are: {sc.strjoin(popdict.keys())}'
            sc.KeyNotFoundError(errormsg)

        isnan = np.isnan(popdict[key]).sum()
        if isnan:
            errormsg = f'Population not fully created: {isnan:,} NaNs found in {key}. This can be caused by calling zn.Agents() instead of zn.make_agents().'
            raise ValueError(errormsg)

    # if ('contacts' not in popdict_keys) and (not hasattr(popdict, 'contacts')) and verbose:
    #     warnmsg = 'No contacts found. Please remember to add contacts before running the simulation.'
    #     znm.warn(warnmsg)

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

    popdict = {}
    contactdict = {}

    # Set pop_pars, these are required parameters for creating the population. For the most part they
    # define the means of the distributions used to create the population. The default values are set in defaults.py
    # and can be overridden by passing them in as kwargs.
    # The default values are currently just dummy values and should be changed to reflect the actual population parameters.
    if 'pop_pars' in kwargs:
        pop_pars = kwargs['pop_pars']
    else:
        pop_pars = znd.default_pop_pars


    # Farms
    n_farms = sim.pars['n_farms']



    # Create water sources
    n_water = round(n_farms * pop_pars['avg_water_per_farm']) # Number of water sources

    # Create barns
    n_barns_by_farm = znu.n_poisson(pop_pars['avg_barns_per_farm'], n_farms) # Number of barns per farm
    n_barns = sum(n_barns_by_farm)
    n_occupied_barns_by_farm = np.zeros(n_farms) # Number of occupied barns per farm

    # Create workers
    n_humans_by_farm = np.zeros(n_farms) # Number of workers per farm
    for farm in range(n_farms):
        n_occupied_barns_by_farm[farm] = sum(znu.n_binomial(pop_pars['avg_barn_occupancy'], n_barns_by_farm[farm])) # Number of occupied barns per farm
        n_humans_by_farm[farm] = sum(znu.n_binomial(pop_pars['avg_humans_per_barn'], n_barns_by_farm[farm])) # Number of workers per barn
    n_humans = sum(n_humans_by_farm) # Total number of workers in the simulation

    # Create flocks
    n_flocks_by_farm = n_occupied_barns_by_farm # 
    n_flocks = sum(n_flocks_by_farm)


    

    n_agents = n_humans + n_barns + n_flocks + n_water

    popdict['uid'] = np.arange(n_agents)
    popdict['agent_type'] = np.concatenate(np.repeat('human', n_humans), 
                                           np.repeat('flock', n_flocks), 
                                           np.repeat('barn', n_barns), 
                                           np.repeat('water', n_water))
    popdict['human_uids'] = popdict['uid'][popdict['agent_type'] == 'human']
    popdict['barn_uids'] = popdict['uid'][popdict['agent_type'] == 'barn']
    popdict['flock_uids'] = popdict['uid'][popdict['agent_type'] == 'flock']
    popdict['water_uids'] = popdict['uid'][popdict['agent_type'] == 'water']

    human_index = 1
    barn_index = 1
    flock_index = 1
    water_index = znu.choose_r(n_water, n_farms) # Randomly assign water sources to farms
    flock2barn = {}
    barn2water = {}

    for farm in range(n_farms):
        contactdict[farm] = {
            'humans':popdict['human_uids'][human_index:human_index + n_humans_by_farm[farm]],
            'barns':popdict['barn_uids'][barn_index:barn_index + n_barns_by_farm[farm]],
            'flocks':popdict['flock_uids'][flock_index:flock_index + n_flocks_by_farm[farm]],
            'water':popdict['water_uids'][water_index[farm]],
        }
        human_index += n_humans_by_farm[farm]
        barn_index += n_barns_by_farm[farm]
        flock_index += n_flocks_by_farm[farm]

        occupied_barns = znu.choose(n_barns_by_farm[farm], n_occupied_barns_by_farm[farm])
        contactdict[farm][flock2barn] = dict(zip(contactdict[farm]['flocks'], contactdict[farm]['barns'][occupied_barns]))# Map flocks to barns for this farm

        contactdict[farm][barn2water] = {barn :contactdict[farm]['water'][0] for barn in contactdict[farm]['barns']} # Map barns to water sources for this farm


        flock2barn.update(contactdict[farm][flock2barn]) # Map flocks to barns for all farms
        barn2water.update(contactdict[farm][barn2water]) # Map barns to water sources for all farms
    
    popdict['contactdict'] = contactdict # Add the contact dictionary to the population dictionary
    popdict['barn2water'] = barn2water # Add the barn to water mapping to the population dictionary
    popdict['flock2barn'] = flock2barn # Add the flock to barn mapping to the population dictionary
    return popdict

def make_humans():

    return
def make_flocks():

    return

def make_barns():

    return

def make_water():

    return

def make_contacts():

    return
