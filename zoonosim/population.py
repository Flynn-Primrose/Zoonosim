'''
Defines functions for creating agents
'''

import numpy as np # Needed for a few things not provided by pl
import sciris as sc

from . import defaults as znd
from . import utils as znu

from .rosters.agents import Agents
from .rosters.humans import Humans
from .rosters.ppe import PPE
from .rosters.barns import Barns
from .rosters.flocks import Flocks
from .rosters.waters import Water
from .rosters.contacts import Contacts, Layer


__all__ = ['make_agents']

def make_agents(sim, popdict=None, reset = False, **kwargs):
    '''
    Create agents for the simulation.

    Args:
        sim             (Sim)  : the simulation object; population parameters are taken from the sim object
        popdict         (any)  : either None, pop.popdict, popfile, or Agents object
        reset           (bool) : whether to force population creation even if self.popdict/self.agents exists
        kwargs          (dict) : 'skip_layers' is passed as an argument to make_contacts. All other keywords are passed to make_popdict() which looks for pop_pars and poultry_pars.

    Returns:
        agents (Agents): an instance of the Agents class containing the agents and their contacts
    '''

    skip_layers = kwargs.pop('skip_layers', None) # Retrieve any layers to be skipped, these will be passed to make_contacts().

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
            'ppe': len(popdict['ppe_uids']),
            'flock': len(popdict['flock_uids']),
            'barn': len(popdict['barn_uids']),
            'water': len(popdict['water_uids']),
        }
    }

    sim.pars.update(new_pars) # Update the simulation parameters with the new population size
    sim.popdict = popdict # Store the population dictionary in the simulation object

    ppe = make_ppe(sim.pars, popdict['ppe_uids'], popdict['ppe2human'])
    human = make_humans(sim.pars, popdict['human_uids'], popdict['transient_uids'], popdict['human2ppe'], ppe.schedule_quarantine) # We pass the ppe quarantine function to the human roster so that when a human is quarantined, their assigned ppe can also be quarantined
    flock = make_flocks(sim.pars, popdict['flock_uids'], popdict['flock2barn'], popdict['flock_breed_index'], barn.schedule_cycle_end, barn.schedule_composting) # We pass the barn schedule_cycle_end function to the flock roster so that when a flock is repopulated, the assigned barn can be provided with a new cycle end date
    barn = make_barns(sim.pars, popdict['barn_uids'], popdict['barn2type'], popdict['barn2resident'], popdict['barn2breed'], flock.repopulate, flock.end_cycle) # We pass the flock.repopulate function to the barn roster so that when a barn is scheduled for repopulation, the assigned flock can also be repopulated
    water = make_water(sim.pars, popdict['water_uids'])
    contacts = make_contacts(sim.pars, popdict, skip_layers)

    agents = Agents(sim.pars, 
                    uid = popdict['uid'],
                    fid = popdict['fid'], 
                    agent_type = popdict['agent_type'], 
                    human = human,
                    ppe = ppe,
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
                     'fid',
                     'agent_type',
                     'human_uids',
                     'barn_uids',
                     'flock_uids',
                     'water_uids',
                     'farmdict',
                     'barn2water',
                     'flock2barn',
                     'flock_breed_index',
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




    # Set pop_pars, these are required parameters for creating the population. For the most part they
    # define the means of the distributions used to create the population.
    # and can be overridden by passing them in as kwargs.
    if 'pop_pars' in kwargs:
        pop_pars = kwargs['pop_pars']
    else:
        pop_pars = sim.pars['pop_pars']

    # Set poultry_pars, these are parameters used for creating the poultry flocks.
    # Specifically we will need the frequency of each breed type
    # Can be overridden by passing them as kwargs.
    if 'poultry_pars' in kwargs:
        poultry_pars = kwargs['poultry_pars']
    else:
        poultry_pars = sim.pars['poultry_pars']

    # Farms
    n_farms = sim.pars['n_farms']
    farm_ids = np.arange(n_farms, dtype=znd.default_int) # Create a list of unique IDs for each farm


    # Create water sources
    n_water = round(n_farms * pop_pars['avg_water_per_farm']) # Number of water sources

    # Create barns
    n_barns_by_farm = np.maximum(np.array(znu.n_poisson(pop_pars['avg_barns_per_farm'], n_farms), dtype=znd.default_int), np.repeat(1, n_farms)) # Number of barns per farm
    n_barns = sum(n_barns_by_farm)
    

    # Create workers
    n_humans_by_farm = np.zeros(n_farms, dtype=znd.default_int) # Number of workers per farm
    n_ppe_by_farm = np.zeros(n_farms, dtype = znd.default_int) # Number of PPE per farm
    for farm in range(n_farms):
        #n_occupied_barns_by_farm[farm] = sum(znu.n_binomial(pop_pars['avg_barn_occupancy'], n_barns_by_farm[farm])) # Number of occupied barns per farm
        n_humans_by_farm[farm] = int(sum(znu.n_binomial(pop_pars['avg_humans_per_barn'], n_barns_by_farm[farm]))) # Number of workers per barn
        n_ppe_by_farm[farm] = n_humans_by_farm[farm] # Every farm will have the same number of ppe as humans
    n_humans = sum(n_humans_by_farm) + pop_pars['number_of_transients'] # Total number of workers in the simulation
    n_ppe = n_humans # We assume one ppe per human

    # Create flocks
    n_flocks = sum(n_barns_by_farm)


    n_agents = n_humans + n_ppe + n_barns + n_flocks + n_water


    popdict = {}
    popdict['uid'] = np.arange(n_agents, dtype=znd.default_int) # Create a list of unique IDs for each agent
    popdict['fid'] = np.zeros(n_agents, dtype=znd.default_int) # Create a list of farm IDs for each agent
    popdict['agent_type'] = np.repeat(['human', 'ppe', 'barn', 'flock', 'water'], [n_humans, n_humans, n_barns, n_flocks, n_water]) # Create a list of agent types

    popdict['human_uids'] = popdict['uid'][popdict['agent_type'] == 'human']
    popdict['ppe_uids'] = popdict['uid'][popdict['agent_type'] == 'ppe']
    popdict['barn_uids'] = popdict['uid'][popdict['agent_type'] == 'barn']
    popdict['flock_uids'] = popdict['uid'][popdict['agent_type'] == 'flock']
    popdict['water_uids'] = popdict['uid'][popdict['agent_type'] == 'water']

    popdict['visits_per_day'] = pop_pars['visits_per_day']

    human_index = 0
    ppe_index = 0
    barn_index = 0
    flock_index = 0

    water_index = znu.choose_r(n_water, n_farms) # Randomly assign water sources to farms

    flock_breed_index = znu.n_multinomial(poultry_pars['breed_freqs'], len(popdict['flock_uids']))# Assign each flock a breed
    flock2breed = dict(zip(popdict['flock_uids'], flock_breed_index))


    flock2barn = {} # Dictionary to hold the mapping of flocks to barns
    human2ppe = {} # Dictionary to hold the mapping of humans to ppe
    barn2water = {} # Dictionary to hold the mapping of barns to water sources

    farmdict = {}
    for farm in range(n_farms):
        flock_value = popdict['flock_uids'][flock_index:(flock_index + n_barns_by_farm[farm])]
        farmdict[farm] = {
            'humans':popdict['human_uids'][human_index:(human_index + n_humans_by_farm[farm])],
            'ppe':popdict['ppe_uids'][ppe_index:(ppe_index + n_ppe_by_farm[farm])],
            'barns':popdict['barn_uids'][barn_index:(barn_index + n_barns_by_farm[farm])],
            'water':popdict['water_uids'][water_index[farm]],
            'flocks':flock_value
        }
        present_uids = np.concatenate((farmdict[farm]['humans'], farmdict[farm]['ppe'], farmdict[farm]['barns'], farmdict[farm]['flocks'])) # Combine all uids for this farm excluding water sources as they are shared among multiple farms
        present_inds = np.isin(popdict['uid'], present_uids) # Get the indices of the agents in this farm
        popdict['fid'][present_inds] =  farm_ids[farm] # Assign the farm ID to the agents in this farm
        human_index += n_humans_by_farm[farm]
        ppe_index += n_ppe_by_farm[farm]
        barn_index += n_barns_by_farm[farm]
        flock_index += n_barns_by_farm[farm]

        if farmdict[farm]['flocks'].size > 0: # Only create the flock to barn mapping if there are flocks on this farm
            farmdict[farm]['flock2barn'] = dict(zip(farmdict[farm]['flocks'], farmdict[farm]['barns']))# Map flocks to barns for this farm
            farmdict[farm]['flock2breed'] = {flock : flock2breed[flock] for flock in farmdict[farm]['flocks']} # Map flocks to breeds for this farm
            flock2barn.update(farmdict[farm]['flock2barn']) # Map flocks to barns for all farms

        farmdict[farm]['human2ppe'] = dict(zip(farmdict[farm]['humans'], farmdict[farm]['ppe']))# Map people to their ppe

        farmdict[farm]['barn2water'] = {barn :farmdict[farm]['water'] for barn in farmdict[farm]['barns']} # Map barns to water sources for this farm




        human2ppe.update(farmdict[farm]['human2ppe']) # map humans to ppe for all farms
        barn2water.update(farmdict[farm]['barn2water']) # Map barns to water sources for all farms

    transient_uids = popdict['human_uids'][human_index:]
    transient_ppe = popdict['ppe_uids'][ppe_index:]
    transient_human2ppe = dict(zip(transient_uids, transient_ppe))

    human2ppe.update(transient_human2ppe)

    popdict['transient_uids'] = transient_uids

    popdict['farmdict'] = farmdict # Add the contact dictionary to the population dictionary
    popdict['flock_breed_index'] = flock_breed_index

    popdict['barn2water'] = barn2water # Add the barn to water mapping to the population dictionary

    popdict['flock2barn'] = flock2barn # Add the flock to barn mapping to the population dictionary
    barn2flock = {v: k for k, v in flock2barn.items()}
    popdict['barn2flock'] = barn2flock


    popdict['human2ppe'] = human2ppe
    ppe2human = {v: k for k, v in human2ppe.items()}
    popdict['ppe2human'] = ppe2human
    popdict['barn2breed'] = {k: flock2breed[v] for k, v in barn2flock.items()}
    popdict['barn2type'] = {k: 'flock' for k in barn2flock.keys()}
    return popdict

def make_humans(sim_pars, uid, transient_uid, human2ppe, schedule_ppe_quarantine=None):
    permenent_uids = np.setdiff1d(uid, transient_uid) # Get the permanent human uids by taking the set difference of all human uids and transient human uids
    sex = znu.n_binomial(0.5, len(uid))
    age  = np.maximum(18, znu.n_poisson(40, len(uid))) # NOTE: Dummy values (assume average worker age of 40)
    ppe = np.empty(len(uid), dtype=znd.default_int)
    for index in range(len(uid)):
        ppe[index] = human2ppe[uid[index]]
    if sim_pars['enable_smartwatches']: # If smartwatches are enabled we randomly assign them to 25% of the human population
        if sim_pars['smartwatch_pars']['who'] == 'all':
            n_true = int(len(uid)*sim_pars['smartwatch_pars']['participation_rate']) # Number of humans with smartwatches
            n_false = len(uid) - n_true
            has_watch = np.array([True]*n_true + [False]*n_false)
            np.random.shuffle(has_watch)
        elif sim_pars['smartwatch_pars']['who'] == 'permanent':
            n_true = int(len(permenent_uids)*sim_pars['smartwatch_pars']['participation_rate']) # Number of permanent humans with smartwatches
            n_false = len(permenent_uids) - n_true
            has_watch_perm = np.array([True]*n_true + [False]*n_false)
            np.random.shuffle(has_watch_perm)
            has_watch = np.zeros(len(uid), dtype=bool)
            has_watch[np.isin(uid, permenent_uids)] = has_watch_perm
        elif sim_pars['smartwatch_pars']['who'] == 'transient':
            n_true = int(len(transient_uid)*sim_pars['smartwatch_pars']['participation_rate']) # Number of transient humans with smartwatches
            n_false = len(transient_uid) - n_true
            has_watch_transient = np.array([True]*n_true + [False]*n_false)
            np.random.shuffle(has_watch_transient)
            has_watch = np.zeros(len(uid), dtype=bool)
            has_watch[np.isin(uid, transient_uid)] = has_watch_transient
        humans = Humans(sim_pars, schedule_ppe_quarantine=schedule_ppe_quarantine, strict = False, uid = uid, age = age, sex = sex, ppe = ppe, has_watch = has_watch)
    else: humans = Humans(sim_pars, schedule_ppe_quarantine=schedule_ppe_quarantine, strict = False, uid = uid, age = age, sex = sex, ppe = ppe)

    return humans

def make_ppe(sim_pars, uid, ppe2human):
    human = np.empty(len(uid), dtype=znd.default_int)
    for index in range(len(uid)):
        human[index] = ppe2human[uid[index]]
    ppe = PPE(sim_pars, strict = False, uid=uid, human=human)
    return ppe

def make_flocks(sim_pars, uid, flock2barn, breed_index, schedule_cycle_end=None, schedule_composting=None):
    prod_pars = sim_pars['poultry_pars']

    breed = np.empty(len(uid), dtype = object)
    headcount = np.empty(len(uid), dtype=znd.default_float)
    barn = np.empty(len(uid), dtype=znd.default_int)
    for index in range(len(uid)):
        breed[index] = prod_pars['breeds'][breed_index[index]] # Get the breed for this flock
        #breed[index] = znd.default_flock_breeds[breed_index[index]] # Get the breed for this flock
        barn[index] = flock2barn[uid[index]]
    # breed = sim_pars['flock_breeds'][breed_index]
    # barn = flock2barn[uid]

    this_breed, freq = np.unique(breed_index, return_counts=True)
    breed_dict = dict(zip(this_breed, freq))
    for this_breed, freq in breed_dict.items():
        headcount[breed_index == this_breed] = znu.sample(**prod_pars['flock_size'][this_breed], size = freq)
    
    flocks = Flocks(sim_pars, schedule_cycle_end=schedule_cycle_end, schedule_composting=schedule_composting, strict = False, uid=uid, breed = breed, barn = barn, headcount=headcount)
    return flocks

def make_barns(sim_pars, uid, barn2type, barn2resident, barn2breed, repopulate=None, end_cycle=None):
    temperature = znu.n_poisson(22.5, len(uid)) # NOTE: Dummy values
    humidity = znu.n_poisson(45, len(uid)) # NOTE: Dummy values
    resident_uid = np.empty(len(uid), dtype=znd.default_int)
    breed_index = np.empty(len(uid), dtype=znd.default_int)
    type_index = np.empty(len(uid), dtype=znd.default_str)
    date_cycle_end = np.empty(len(uid), dtype=znd.default_float)
    for index in range(len(uid)):
        type_index[index] = barn2type[uid[index]]
        resident_uid[index] = barn2resident[uid[index]]
        breed_index[index] = barn2breed[uid[index]]


    poultry_barn_uids = [uid for uid in uid if barn2type[uid] == 'flock']
    poultrybarn2breed = {uid: barn2breed[uid] for uid in poultry_barn_uids}
    breed, freq = np.unique(list(poultrybarn2breed.values()), return_counts=True)
    breed_dict = dict(zip(breed, freq))
    for breed, freq in breed_dict.items():
        date_cycle_end[(type_index == 'flock') & (breed_index == breed)] = znu.sample(**sim_pars['poultry_pars']['cycle_dur'][breed], size = freq)

    barns = Barns(sim_pars, repopulate=repopulate, end_cycle=end_cycle, strict = False, uid=uid, temperature = temperature, humidity = humidity, resident_uid = resident_uid, date_cycle_end = date_cycle_end)
    return barns

def make_water(sim_pars, uid):

    temperature = znu.n_poisson(22.5, len(uid)) # NOTE: Dummy values
    water = Water(sim_pars, strict = False, uid = uid, temperature = temperature)
    
    return water

def make_contacts(sim_pars, popdict, skip_layers=None, layers_to_make=None):
    
    '''
    Create contacts for the simulation.

    Args:
        popdict     (dict) : dictionary containing details regarding the population and its contacts. Expected to be generated by 'make_popdict'
                            or loaded from a saved file
        skip_layers   (list) : list of layer names to skip when creating contacts
        layers_to_make (list) : list of layer names to create contacts for
    '''

    if skip_layers is None:
        skip_layers = []

    if layers_to_make is None:
        all_layers = sim_pars['beta_layer'].keys()
    else:
        all_layers = layers_to_make


    data = {}
    for layer in all_layers:
        if layer not in skip_layers:
            if layer in sim_pars['beta_layer'].keys(): # Only make the layer if it's specified in the beta_layer dictionary, otherwise skip it
                data[layer] = make_layer_contacts(layer, popdict, sim_pars['beta_layer'][layer])
            else:
                warnmsg = f'Layer "{layer}" specified in layers_to_make but not found in sim_pars["beta_layer"]; skipping this layer.'
                znu.warn(warnmsg)
    return Contacts(data=data)

def make_layer_contacts(layer, popdict, beta):
    '''
    Create contacts for a specific layer.

    Args:
        layer (str): the name of the layer to create contacts for
        popdict (dict): a dictionary containing the details needed to make the contacts. Expected to be generated by 'make_popdict'
        beta (float): the weight of this layer

    Returns:
        layer (Layer): A Layer object containing the contacts for this layer
    '''
    human_human_layer = None
    human_flock_layer = None
    human_barn_layer = None
    human_water_layer = None

    match layer:
        case 'human_ppe':
            return human_ppe_contacts(popdict, beta)
        case 'human_human':
            if human_human_layer is None: # If the human-human layer has not already been created, we need to create it in order to create the human-human layer
                human_human_layer = human_human_contacts(popdict, beta) # Needed to create the ppe-ppe layer, as we assume that ppe-ppe contacts are determined by the human-human contacts (i.e. if two humans have contact, their ppe also have contact)
            return human_human_layer
        case 'human_flock':
            if human_flock_layer is None: # If the human-flock layer has not already been created, we need to create it in order to create the human-flock layer
                human_flock_layer = human_flock_contacts(popdict, beta) # Needed to create the ppe-flock layer, as we assume that ppe-flock contacts are determined by the human-flock contacts (i.e. if a human has contact with a flock, their ppe also has contact with that flock)
            return human_flock_layer
        case 'human_barn':
            if human_barn_layer is None: # If the human-barn layer has not already been created, we need to create it in order to create the human-barn layer
                human_barn_layer = human_barn_contacts(popdict, beta) # Needed to create the ppe-barn layer, as we assume that ppe-barn contacts are determined by the human-barn contacts (i.e. if a human has contact with a barn, their ppe also has contact with that barn)
            return human_barn_layer
        case 'human_water':
            if human_water_layer is None: # If the human-water layer has not already been created, we need to create it in order to create the human-water layer
                human_water_layer = human_water_contacts(popdict, beta) # Needed to create the ppe-water layer, as we assume that ppe-water contacts are determined by the human-water contacts (i.e. if a human has contact with water, their ppe also has contact with that water)
            return human_water_layer
        case 'ppe_ppe':
            if human_human_layer is None: # If the human-human layer has not already been created, we need to create it in order to create the ppe-ppe layer
                human_human_layer = human_human_contacts(popdict, beta)
            return ppe_ppe_contacts(human_human_layer, popdict, beta)
        case 'ppe_flock':
            if human_flock_layer is None: # If the human-flock layer has not already been created, we need to create it in order to create the ppe-flock layer
                human_flock_layer = human_flock_contacts(popdict, beta)
            return ppe_flock_contacts(human_flock_layer, popdict, beta)
        case 'ppe_barn':
            if human_barn_layer is None: # If the human-barn layer has not already been created, we need to create it in order to create the ppe-barn layer
                human_barn_layer = human_barn_contacts(popdict, beta)
            return ppe_barn_contacts(human_barn_layer, popdict, beta)
        case 'ppe_water':
            if human_water_layer is None: # If the human-water layer has not already been created, we need to create it in order to create the ppe-water layer
                human_water_layer = human_water_contacts(popdict, beta)
            return ppe_water_contacts(human_water_layer, popdict, beta)
        case 'flock_barn':
            return flock_barn_contacts(popdict, beta)
        case 'flock_water':
            return flock_water_contacts(popdict, beta)
        case 'barn_water':
            return barn_water_contacts(popdict, beta)
        case 'transient':
            return transient_contacts(popdict, beta)

def human_ppe_contacts(popdict, beta):
    '''
    Create contacts between all humans (including transients) and their assigned PPE

    Args:
        popdict (dict): A dictionary containing the details needed to make the contacts. Expected to be generated by 'make_popdict'
        beta (float): the weight of this layer
    Returns:
        hp_layer (Layer): A layer object containing all human-ppe contacts
    '''
    hp_1 = popdict['human_uids']
    hp_2 = np.vectorize(popdict['human2ppe'].get)(popdict['human_uids']) # This is a more efficient way to get the ppe for each human than looping through them
    hp_beta = np.repeat(beta, len(hp_1))


    return Layer(p1 = hp_1, p2 = hp_2, beta = hp_beta, label = 'human-PPE contacts')

def human_human_contacts(popdict, beta):
    '''
    Create human-human contacts for the simulation.

    Args:
        farmdict     (dict) : dictionary containing the contacts between agents

    Returns:
        hh_layer: a Layer object containing the human-human contacts
    '''

    hh_p1 = []

    hh_p2 = []

    farmdict = popdict['farmdict']

    for farm, farm_contacts in farmdict.items():
        for human1 in farm_contacts['humans']:
            for human2 in farm_contacts['humans']:
                if human1 != human2:
                    hh_p1.append(human1)
                    hh_p2.append(human2) # Get the barn for this human
                    # NOTE: This assumes that humans have contact with all humans on the farm

    hh_beta = np.repeat(beta, len(hh_p1)) 
    return Layer(p1 = hh_p1, p2 = hh_p2, beta = hh_beta, label = 'human-human contacts')

def human_flock_contacts(popdict, beta):
    '''
    Create human-flock contacts for the simulation.

    Args:
        popdict     (dict) :
        beta (float) :

    Returns:
        hf_layer: a Layer object containing the human-flock contacts
    '''

    hf_p1 = []

    hf_p2 = []

    farmdict = popdict['farmdict']

    for farm, farm_contacts in farmdict.items():
        for human in farm_contacts['humans']:
            for flock in farm_contacts['flocks']:
                hf_p1.append(human)
                hf_p2.append(flock) # Get the barn for this flock
                # NOTE: This assumes that humans have contact with all flocks on the farm

    hf_beta = np.repeat(beta, len(hf_p1)) 

    return Layer(p1 = hf_p1, p2 = hf_p2, beta = hf_beta, label = 'human-flock contacts')


def human_barn_contacts(popdict, beta):
    '''
    Create human-barn contacts for the simulation.

    Args:
        popdict     (dict) : 
        beta    (float) :

    Returns:
        hb_layer: a Layer object containing the human-barn contacts
    '''

    hb_p1 = []

    hb_p2 = []

    farmdict = popdict['farmdict']

    for farm, farm_contacts in farmdict.items():
        for human in farm_contacts['humans']:
            for barn in farm_contacts['barns']:
                hb_p1.append(human)
                hb_p2.append(barn) # Get the barn for this human
                # NOTE: This assumes that humans have contact with all barns on the farm

    hb_beta = np.repeat(beta, len(hb_p1)) 

    return Layer(p1 = hb_p1, p2 = hb_p2, beta = hb_beta, label = 'human-barn contacts')

def human_water_contacts(popdict, beta):
    '''
    Create human-water contacts

    Args:
        popdict (dict):
        beta    (float):
    returns:
        hw_layer    (Layer):
    '''

    hw_1 = []
    hw_2 = []
    farmdict = popdict['farmdict']

    for farm, farm_contacts in farmdict.items():
        for human in farm_contacts['humans']:
            hw_1.append(human)
            hw_2.append(farm_contacts['water']) # Get the water source for this human
            # NOTE: This assumes all humans on a farm come into contact with all water sources on that farm
    
    hw_beta = np.repeat(beta, len(hw_1))
 
    return Layer(p1 = hw_1, p2 = hw_2, beta = hw_beta, label = 'human-water contacts')

def ppe_ppe_contacts(human_human_layer, popdict, beta):
    '''
    Create PPE-PPE contacts

    Args:
        hh_layer    (Layer): The Layer object with human-human contacts
        popdict (dict):
        beta    (float): 
    Returns:
        pp_layer    (Layer):
    '''
    pp_1 = np.vectorize(popdict['human2ppe'].get)(human_human_layer['p1']) # This is a more efficient way to get the ppe for each human than looping through them
    pp_2 = np.vectorize(popdict['human2ppe'].get)(human_human_layer['p2']) # This is a more efficient way to get the ppe for each human than looping through them
    pp_beta = np.repeat(beta, len(pp_1))

    return Layer(p1 = pp_1, p2 = pp_2, beta = pp_beta, label = 'PPE-PPE contacts')

def ppe_flock_contacts(human_flock_layer, popdict, beta):
    '''
    Create PPE-flock contacts

    Args:
        hf_layer    (Layer): The Layer object with human-flock contacts
        popdict (dict):
        beta    (float): 
    Returns:
        pf_layer    (Layer):
    '''

    pf_1 = np.vectorize(popdict['human2ppe'].get)(human_flock_layer['p1']) # This is a more efficient way to get the ppe for each human than looping through them
    pf_2 = human_flock_layer['p2']
    pf_beta = np.repeat(beta, len(pf_1))
    return Layer(p1 = pf_1, p2 = pf_2, beta = pf_beta, label = 'PPE-Flock contacts')


def ppe_barn_contacts(human_barn_layer, popdict, beta):
    '''
    Create PPE-barn contacts

    Args:
        hb_layer    (Layer): The Layer object with human-barn contacts
        popdict (dict):
        beta    (float): 
    Returns:
        pb_layer    (Layer):
    '''
    # pb_1 = human_barn_layer['p1'][popdict['human2ppe']]
    pb_1 = np.vectorize(popdict['human2ppe'].get)(human_barn_layer['p1']) # This is a more efficient way to get the ppe for each human than looping through them
    pb_2 = human_barn_layer['p2']
    pb_beta = np.repeat(beta, len(pb_1))
    pb_layer = Layer(p1 = pb_1,
                     p2 = pb_2,
                     beta = pb_beta,
                     label = 'PPE-Barn contacts')
    return pb_layer

def ppe_water_contacts(human_water_layer, popdict, beta):
    '''
    Create PPE-water contacts

    Args:
        hw_layer    (Layer): The Layer object with human-water contacts
        popdict (dict):
        beta    (float): 
    Returns:
        pb_layer    (Layer):
    '''
    # pw_1 = human_water_layer['p1'][popdict['human2ppe']]
    pw_1 = np.vectorize(popdict['human2ppe'].get)(human_water_layer['p1']) # This is a more efficient way to get the ppe for each human than looping through them
    pw_2 = human_water_layer['p2']
    pw_beta = np.repeat(beta, len(pw_1))
    pw_layer = Layer(p1 = pw_1,
                     p2 = pw_2,
                     beta = pw_beta,
                     label = 'PPE-Water contacts')
    return pw_layer

def flock_barn_contacts(popdict, beta):
    '''
    Create flock-barn contacts for the simulation.

    Args:
        popdict     (dict) : 
        beta    (float) :

    Returns:
        fb_layer: a Layer object containing the flock-barn contacts
    '''

    fb_p1 = []

    fb_p2 = []

    farmdict = popdict['farmdict']

    for farm, farm_contacts in farmdict.items():
        for flock in farm_contacts['flocks']:
            fb_p1.append(flock)
            fb_p2.append(farm_contacts['flock2barn'][flock]) # Get the barn for this flock

    fb_beta = np.repeat(beta, len(fb_p1))

    return Layer(p1 = fb_p1, p2 = fb_p2, beta = fb_beta, label = 'flock-barn contacts')


def flock_water_contacts(popdict, beta):
    '''
    Create flock-water contacts for the simulation.

    Args:
        popdict     (dict) : 
        beta    (float) :

    Returns:
        fw_layer: a Layer object containing the flock-water contacts
    '''

    fw_p1 = []

    fw_p2 = []

    farmdict = popdict['farmdict']

    for farm, farm_contacts in farmdict.items():
        for flock in farm_contacts['flocks']:
            fw_p1.append(flock) # Get the barn for this flock
            fw_p2.append(farm_contacts['barn2water'][farm_contacts['flock2barn'][flock]]) # Get the water source for this flock

    fw_beta = np.repeat(beta, len(fw_p1))

    return Layer(p1 = fw_p1, p2 = fw_p2, beta = fw_beta, label = 'flock-water contacts')

def barn_water_contacts(popdict, beta):
    '''
    Create barn-water contacts for the simulation.

    Args:
        popdict     (dict) : 
        beta    (float) :

    Returns:
        bw_layer: a Layer object containing the barn-water contacts
    '''

    bw_p1 = []

    bw_p2 = []

    farmdict = popdict['farmdict']

    for farm, farm_contacts in farmdict.items():
        for barn in farm_contacts['barns']:
            bw_p1.append(barn)
            bw_p2.append(farm_contacts['barn2water'][barn]) # Get the water source for this barn

    bw_beta = np.repeat(beta, len(bw_p1))

    return Layer(p1 = bw_p1, p2 = bw_p2, beta = bw_beta, label = 'barn-water contacts')


class transient_Layer(Layer):
    '''
    A subclass of Layer to represent transient contacts. The key difference is the override of update() that allows the contacts to change at each timestep.
    '''
    def __init__(self, p1, p2, beta, label, popdict):
        super().__init__(p1=p1, p2=p2, beta=beta, label=label)
        self.popdict = popdict

    def update(self):
        '''
        Update the contacts for this layer.
        '''

        p1 = np.array([], dtype=znd.default_int)
        p2 = np.array([], dtype=znd.default_int)

        for transient in self.popdict['transient_humans']:
            farms = np.random.choice(list(self.popdict['farmdict'].keys()), self.popdict['visits_per_day']) # Randomly assign 'visits_per_day' farms to this transient
            for farm in farms:
                humans = self.popdict['farmdict'][farm]['humans']
                flocks = self.popdict['farmdict'][farm]['flocks']
                barns = self.popdict['farmdict'][farm]['barns']
                water = self.popdict['farmdict'][farm]['water']
                transient_p1 = np.repeat(transient, len(humans) + len(flocks) + len(barns) + 1) # The +1 is for the water source.
                transient_p2 = np.concatenate((humans, flocks, barns, [water]))
                ppe_p1 = [self.popdict['human2ppe'][transient] for transient in transient_p1]
                ppe_p2 = np.concatenate(([self.popdict['human2ppe'][human] for human in humans], flocks, barns, [water]))
                p1 = np.concatenate((p1, transient_p1), dtype=transient_p1.dtype)
                p2 = np.concatenate((p2, transient_p2), dtype=transient_p2.dtype)
                p1 = np.concatenate((p1, ppe_p1), dtype=ppe_p1.dtype)
                p2 = np.concatenate((p2, ppe_p2), dtype=ppe_p2.dtype)
        self.p1 = p1
        self.p2 = p2
        self.beta = np.repeat(self.beta, len(self.p1)) # Assuming the same beta for all transient contacts, we just repeat the first value to match the new length of p1 and p2
        return

def transient_contacts(popdict, beta):
    '''
    Create transient contacts for the simulation. These are contacts that change at each timestep, representing the fact that transient workers may visit different farms on different days.

    Args:
        popdict     (dict) : dictionary containing details regarding the population and its contacts. Expected to be generated by 'make_popdict'
        beta    (float) : the weight of this layer

    Returns:
        transient_layer: a transient_Layer object containing the transient contacts
    '''

    p1 = np.array([], dtype=znd.default_int)
    p2 = np.array([], dtype=znd.default_int)

    for transient in popdict['transient_uids']:
        farms = np.random.choice(list(popdict['farmdict'].keys()), popdict['visits_per_day']) # Randomly assign 'visits_per_day' farms to this transient
        for farm in farms:
            humans = popdict['farmdict'][farm]['humans']
            flocks = popdict['farmdict'][farm]['flocks']
            barns = popdict['farmdict'][farm]['barns']
            water = popdict['farmdict'][farm]['water']
            transient_p1 = np.repeat(transient, len(humans) + len(flocks) + len(barns) + 1) # The +1 is for the water source.
            transient_p2 = np.concatenate((humans, flocks, barns, [water]))
            ppe_p1 = np.array([popdict['human2ppe'][transient] for transient in transient_p1])
            ppe_p2 = np.concatenate(([popdict['human2ppe'][human] for human in humans], flocks, barns, [water]))
            p1 = np.concatenate((p1, transient_p1), dtype=transient_p1.dtype)
            p2 = np.concatenate((p2, transient_p2), dtype=transient_p2.dtype)
            p1 = np.concatenate((p1, ppe_p1), dtype=ppe_p1.dtype)
            p2 = np.concatenate((p2, ppe_p2), dtype=ppe_p2.dtype)
    transient_beta = np.repeat(beta, len(p1)) # Assuming the same beta for all transient contacts, we just repeat the first value to match the new length of p1 and p2
    transient_layer = transient_Layer(p1, p2, transient_beta, label = 'transient contacts', popdict = popdict)

    return transient_layer