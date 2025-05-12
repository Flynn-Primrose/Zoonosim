'''
Defines functions for creating agents
'''

import numpy as np # Needed for a few things not provided by pl
import sciris as sc

from . import base as znb
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

    new_pars = {
        'pop_size': len(popdict['uid']),
        'pop_size_by_type': {
            'human': len(popdict['human_uids']),
            'barn': len(popdict['barn_uids']),
            'flock': len(popdict['flock_uids']),
            'water': len(popdict['water_uids']),
        }
    }

    sim.pars.update(new_pars) # Update the simulation parameters with the new population size
    sim.popdict = popdict # Store the population dictionary in the simulation object

    human = make_humans(sim.pars, popdict['human_uids'])
    flock = make_flocks(sim.pars, popdict['flock_uids'], popdict['flock2barn'])
    barn = make_barns(sim.pars, popdict['barn_uids'], popdict['barn2flock'])
    water = make_water(sim.pars, popdict['water_uids'])

    contacts = make_contacts(popdict['contactdict'])

    agents = znr.Agents(sim.pars, 
                        uid = popdict['uid'], 
                        agent_type = popdict['agent_type'], 
                        human = human,
                        flock = flock,
                        barn = barn,
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

        #isnan = np.isnan(popdict[key]).sum()
        #if isnan:
            #errormsg = f'Population not fully created: {isnan:,} NaNs found in {key}. This can be caused by calling zn.Agents() instead of zn.make_agents().'
            #raise ValueError(errormsg)

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
    n_barns_by_farm = np.maximum(np.array(znu.n_poisson(pop_pars['avg_barns_per_farm'], n_farms), dtype=znd.default_int), np.repeat(1, n_farms)) # Number of barns per farm
    n_barns = sum(n_barns_by_farm)
    #n_occupied_barns_by_farm = np.zeros(n_farms) # Number of occupied barns per farm

    # Create workers
    n_humans_by_farm = np.zeros(n_farms, dtype=znd.default_int) # Number of workers per farm
    for farm in range(n_farms):
        #n_occupied_barns_by_farm[farm] = sum(znu.n_binomial(pop_pars['avg_barn_occupancy'], n_barns_by_farm[farm])) # Number of occupied barns per farm
        n_humans_by_farm[farm] = int(sum(znu.n_binomial(pop_pars['avg_humans_per_barn'], n_barns_by_farm[farm]))) # Number of workers per barn
    n_humans = sum(n_humans_by_farm) # Total number of workers in the simulation

    # Create flocks
    n_flocks_by_farm = n_barns_by_farm # 
    n_flocks = sum(n_flocks_by_farm)


    

    n_agents = n_humans + n_barns + n_flocks + n_water

    popdict['uid'] = np.arange(n_agents, dtype=znd.default_int) # Create a list of unique IDs for each agent
    popdict['agent_type'] = np.repeat(['human', 'barn', 'flock', 'water'], [n_humans, n_barns, n_flocks, n_water]) # Create a list of agent types
    popdict['human_uids'] = popdict['uid'][popdict['agent_type'] == 'human']
    popdict['barn_uids'] = popdict['uid'][popdict['agent_type'] == 'barn']
    popdict['flock_uids'] = popdict['uid'][popdict['agent_type'] == 'flock']
    popdict['water_uids'] = popdict['uid'][popdict['agent_type'] == 'water']

    human_index = 0
    barn_index = 0
    flock_index = 0
    water_index = znu.choose_r(n_water, n_farms) # Randomly assign water sources to farms
    flock2barn = {} # Dictionary to hold the mapping of flocks to barns
    barn2water = {} # Dictionary to hold the mapping of barns to water sources


    for farm in range(n_farms):
        contactdict[farm] = {
            'humans':popdict['human_uids'][human_index:(human_index + n_humans_by_farm[farm])],
            'barns':popdict['barn_uids'][barn_index:(barn_index + n_barns_by_farm[farm])],
            'flocks':popdict['flock_uids'][flock_index:(flock_index + n_flocks_by_farm[farm])],
            'water':popdict['water_uids'][water_index[farm]],
        }
        human_index += n_humans_by_farm[farm]
        barn_index += n_barns_by_farm[farm]
        flock_index += n_flocks_by_farm[farm]

        #occupied_barns = znu.choose(n_barns_by_farm[farm], n_occupied_barns_by_farm[farm])
        contactdict[farm]['flock2barn'] = dict(zip(contactdict[farm]['flocks'], contactdict[farm]['barns']))# Map flocks to barns for this farm

        contactdict[farm]['barn2water'] = {barn :contactdict[farm]['water'] for barn in contactdict[farm]['barns']} # Map barns to water sources for this farm


        flock2barn.update(contactdict[farm]['flock2barn']) # Map flocks to barns for all farms
        barn2water.update(contactdict[farm]['barn2water']) # Map barns to water sources for all farms
    
    popdict['contactdict'] = contactdict # Add the contact dictionary to the population dictionary
    popdict['barn2water'] = barn2water # Add the barn to water mapping to the population dictionary
    popdict['flock2barn'] = flock2barn # Add the flock to barn mapping to the population dictionary
    popdict['barn2flock'] = {v: k for k, v in flock2barn.items()}
    return popdict

def make_humans(sim_pars, uid):
    sex = znu.n_binomial(0.5, len(uid))
    age  = np.maximum(18, znu.n_poisson(40, len(uid))) # NOTE: Dummy values (assume average worker age of 40)

    humans = znr.Humans(sim_pars, strict = False, uid = uid, age = age, sex = sex)

    return humans

def make_flocks(sim_pars, uid, flock2barn):
    breed_index = znu.n_multinomial(znd.default_flock_breed_freqs, len(uid))
    breed = np.empty(len(uid), dtype = object)
    barn = np.empty(len(uid), dtype=znd.default_int)
    for index in range(len(uid)):
        breed[index] = znd.default_flock_breeds[breed_index[index]] # Get the breed for this flock
        barn[index] = flock2barn[uid[index]]

    flocks = znr.Flocks(sim_pars, strict = False, uid=uid, breed = breed, barn = barn)
    return flocks

def make_barns(sim_pars, uid, barn2flock):
    temperature = znu.n_poisson(22.5, len(uid)) # NOTE: Dummy values
    humidity = znu.n_poisson(45, len(uid)) # NOTE: Dummy values
    flock = np.empty(len(uid), dtype=znd.default_int)
    for index in range(len(uid)): # NOTE: There is probably a better way of doing this
        flock[index] = barn2flock.get(uid[index])
        
    barns = znr.Barns(sim_pars, strict = False, uid=uid, temperature = temperature, humidity = humidity, flock = flock)
    return barns

def make_water(sim_pars, uid):

    temperature = znu.n_poisson(22.5, len(uid)) # NOTE: Dummy values
    water = znr.Water(sim_pars, strict = False, uid = uid, temperature = temperature)
    
    return water

def make_contacts(contactdict):
    
    fb_layer = make_fb_contacts(contactdict) # Flock-barn contacts
    bw_layer = make_bw_contacts(contactdict) # Barn-water contacts
    fw_layer = make_fw_contacts(contactdict) # Flock-water contacts
    hb_layer = make_hb_contacts(contactdict) # Human-barn contacts
    hf_layer = make_hf_contacts(contactdict) # Human-flock contacts
    hh_layer = make_hh_contacts(contactdict) # Human-human contacts

    return znb.Contacts(fb = fb_layer,
                        bw = bw_layer,
                        fw = fw_layer,
                        hb = hb_layer,
                        hf = hf_layer,
                        hh = hh_layer,
                        )

def make_fb_contacts(contactdict):
    '''
    Create flock-barn contacts for the simulation.

    Args:
        sim_pars        (Sim)  : the simulation object; population parameters are taken from the sim object
        contactdict     (dict) : dictionary containing the contacts between agents

    Returns:
        fb_layer: a Layer object containing the flock-barn contacts
    '''

    fb_p1 = []

    fb_p2 = []

    for farm, farm_contacts in contactdict.items():
        for flock in farm_contacts['flocks']:
            fb_p1.append(flock)
            fb_p2.append(farm_contacts['flock2barn'][flock]) # Get the barn for this flock

    beta = np.repeat(0.25, len(fb_p1)) # NOTE: Dummy values

    fb_layer = znb.Layer(p1 = fb_p1,
                         p2 = fb_p2,
                         beta = beta,
                         label = 'flock-barn contacts'
                        )

    return fb_layer

def make_bw_contacts(contactdict):
    '''
    Create barn-water contacts for the simulation.

    Args:
        sim_pars        (Sim)  : the simulation object; population parameters are taken from the sim object
        contactdict     (dict) : dictionary containing the contacts between agents

    Returns:
        bw_layer: a Layer object containing the barn-water contacts
    '''

    bw_p1 = []

    bw_p2 = []

    for farm, farm_contacts in contactdict.items():
        for barn in farm_contacts['barns']:
            bw_p1.append(barn)
            bw_p2.append(farm_contacts['barn2water'][barn]) # Get the water source for this barn

    beta = np.repeat(0.25, len(bw_p1)) # NOTE: Dummy values

    bw_layer = znb.Layer(p1 = bw_p1,
                         p2 = bw_p2,
                         beta = beta,
                         label = 'barn-water contacts'
                        )

    return bw_layer

def make_fw_contacts(contactdict):
    '''
    Create flock-water contacts for the simulation.

    Args:
        sim_pars        (Sim)  : the simulation object; population parameters are taken from the sim object
        contactdict     (dict) : dictionary containing the contacts between agents

    Returns:
        fw_layer: a Layer object containing the flock-water contacts
    '''

    fw_p1 = []

    fw_p2 = []

    for farm, farm_contacts in contactdict.items():
        for flock in farm_contacts['flocks']:
            fw_p1.append(flock) # Get the barn for this flock
            fw_p2.append(farm_contacts['barn2water'][farm_contacts['flock2barn'][flock]]) # Get the water source for this flock

    beta = np.repeat(0.25, len(fw_p1)) # NOTE: Dummy values

    fw_layer = znb.Layer(p1 = fw_p1,
                         p2 = fw_p2,
                         beta = beta,
                         label = 'flock-water contacts'
                        )

    return fw_layer

def make_hb_contacts(contactdict):
    '''
    Create human-barn contacts for the simulation.

    Args:
        sim_pars        (Sim)  : the simulation object; population parameters are taken from the sim object
        contactdict     (dict) : dictionary containing the contacts between agents

    Returns:
        hb_layer: a Layer object containing the human-barn contacts
    '''

    hb_p1 = []

    hb_p2 = []

    for farm, farm_contacts in contactdict.items():
        for human in farm_contacts['humans']:
            for barn in farm_contacts['barns']:
                hb_p1.append(human)
                hb_p2.append(barn) # Get the barn for this human
                # NOTE: This assumes that humans have contact with all barns on the farm

    beta = np.repeat(0.25, len(hb_p1)) # NOTE: Dummy values

    hb_layer = znb.Layer(p1 = hb_p1,
                         p2 = hb_p2,
                         beta = beta,
                         label = 'human-barn contacts'
                        )

    return hb_layer

def make_hf_contacts(contactdict):
    '''
    Create human-flock contacts for the simulation.

    Args:
        sim_pars        (Sim)  : the simulation object; population parameters are taken from the sim object
        contactdict     (dict) : dictionary containing the contacts between agents

    Returns:
        hf_layer: a Layer object containing the human-flock contacts
    '''

    hf_p1 = []

    hf_p2 = []

    for farm, farm_contacts in contactdict.items():
        for human in farm_contacts['humans']:
            for flock in farm_contacts['flocks']:
                hf_p1.append(human)
                hf_p2.append(flock) # Get the barn for this flock
                # NOTE: This assumes that humans have contact with all flocks on the farm

    beta = np.repeat(1, len(hf_p1)) # NOTE: Dummy values

    hf_layer = znb.Layer(p1 = hf_p1,
                         p2 = hf_p2,
                         beta = beta,
                         label = 'human-flock contacts'
                        )

    return hf_layer

def make_hh_contacts(contactdict):
    '''
    Create human-human contacts for the simulation.

    Args:
        sim_pars        (Sim)  : the simulation object; population parameters are taken from the sim object
        contactdict     (dict) : dictionary containing the contacts between agents

    Returns:
        hh_layer: a Layer object containing the human-human contacts
    '''

    hh_p1 = []

    hh_p2 = []

    for farm, farm_contacts in contactdict.items():
        for human1 in farm_contacts['humans']:
            for human2 in farm_contacts['humans']:
                if human1 != human2:
                    hh_p1.append(human1)
                    hh_p2.append(human2) # Get the barn for this human
                    # NOTE: This assumes that humans have contact with all humans on the farm

    beta = np.repeat(0.25, len(hh_p1)) # NOTE: Dummy values

    hh_layer = znb.Layer(p1 = hh_p1,
                         p2 = hh_p2,
                         beta = beta,
                         label = 'human-human contacts'
                        )

    return hh_layer