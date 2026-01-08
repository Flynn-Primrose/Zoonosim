'''
Utilities for working with nested parameter dictionaries.
'''
import numpy as np
from .. import defaults as znd

def compare_pars(dict_new, dict_orig):
    """
    Compare two dictionaries recursively.
    Returns:
        common  -> keys (and subkeys) present in both dict_new and dict_orig
        missing -> keys (and subkeys) present in dict_new but not in dict_orig
    """
    common = {}
    missing = {}
    for key, val in dict_new.items():
        if key in dict_orig:
            if isinstance(val, dict) and isinstance(dict_orig[key], dict):
                # Recurse into subdicts
                sub_common, sub_missing = compare_pars(val, dict_orig[key])
                if sub_common:
                    common[key] = sub_common
                if sub_missing:
                    missing[key] = sub_missing
            else:
                # Key exists in both, keep dict_new's value
                common[key] = val
        else:
            # Key missing from dict_orig
            missing[key] = val
    return common, missing

def is_empty(d):
    ''' Check if a dictionary (possibly nested) is empty '''
    if not isinstance(d, dict):
        return False
    return all(is_empty(v) for v in d.values()) if d else True

def pars_sampler(trial, calib_pars, par_samplers, past_keys=None):
    '''
    A helper function to sample parameters from Optuna trials

    Args:
        trial        (optuna.trial.Trial): the Optuna trial object
        calib_pars   (dict): the calibration parameters dictionary
        par_samplers (dict): mapping from parameters to the Optuna sampler to use for choosing new points for each.
        past_keys    (list): used for recursion to keep track of nested keys
    '''
    sampled_pars = {}
    if past_keys is None:
        past_keys = []
    for key, value in calib_pars.items():
        if key == 'prognoses':
            sampled_pars[key] = {}
            # errormsg = 'Calibration of prognoses parameters is not supported'
            # raise ValueError(errormsg)
            for agent_type, agent_prog in value.items():
                    matches = [(k, v) for k, v in agent_prog.items() if isinstance(v, np.ndarray)]
                    if len(matches) == 0:
                        raise KeyError("No numpy array found in prognosis calibration dictionary")
                    if len(matches) > 1:
                        raise ValueError("More than one numpy array found in prognosis calibration dictionary")
                    strat_key, strat_array = matches[0] # this key value pair should contain the stratification information for the remaining prognosis parameters
                    prog_pars = {strat_key: strat_array}
                    for subkey, subval in agent_prog.items():
                        if subkey == strat_key:
                            continue
                        par_vector = np.array(strat_array.shape[0], dtype=znd.default_float)
                        for i in range(strat_array.shape[0]):
                            par_name = f'{key}.{agent_type}.{subkey}.{strat_array[i]}'
                            best, low, high = (arr[i] for arr in subval)
                            if par_name in par_samplers:
                                if not isinstance(par_samplers[par_name], str):
                                    errormsg = f'Parameter sampler for "{par_name}" must be a string representing an Optuna Trial method'
                                    raise ValueError(errormsg)
                                try:
                                    sampler_fn = getattr(trial, par_samplers[par_name])
                                except Exception as E:
                                    errormsg = 'The requested sampler function is not found: ensure it is a valid attribute of an Optuna Trial object'
                                    raise AttributeError(errormsg) from E
                            else:
                                sampler_fn = trial.suggest_float
                            par_vector[i] = sampler_fn(".".join(past_keys + [par_name]), low, high)
                        prog_pars[subkey] = par_vector
                    sampled_pars[key][agent_type] = prog_pars
            continue
        if isinstance(value, list) and len(value) == 3:
            best, low, high = value
            if key in par_samplers: # If a custom sampler is used, get it now
                if not isinstance(par_samplers[key], str):
                    errormsg = f'Parameter sampler for "{".".join(past_keys + [key])}" must be a string representing an Optuna Trial method'
                    raise ValueError(errormsg)
                try:
                    sampler_fn = getattr(trial, par_samplers[key])
                except Exception as E:
                    errormsg = 'The requested sampler function is not found: ensure it is a valid attribute of an Optuna Trial object'
                    raise AttributeError(errormsg) from E
            else:
                sampler_fn = trial.suggest_float
            sampled_pars[key] = sampler_fn(".".join(past_keys + [key]), low, high) # Sample from values within this range
        elif isinstance(value, list):
            errormsg = f'Parameter "{".".join(past_keys + [key])}" must be a list of [best, low, high]'
            raise ValueError(errormsg)
        elif isinstance(value, dict):
            sampled_pars[key] = pars_sampler(trial, value, par_samplers.get(key, {}), past_keys=past_keys + [key]) # Recurse into sub-dictionaries
        else:
            errormsg = f'Parameter "{".".join(past_keys + [key ])}" must be a list of [best, low, high] or a dictionary for nested parameters'
            raise ValueError(errormsg)
    return sampled_pars
            
def pars_parser(calib_pars):
    '''
    A helper function to parse calibration parameters into initial, and bounds.
    '''
    initial_pars = {}
    par_bounds = {}
    for key, value in calib_pars.items():
        if isinstance(value, list) and len(value) == 3:
            best, low, high = value
            initial_pars[key] = best
            par_bounds[key] = [low, high]
        elif isinstance(value, dict):
            sub_best, sub_bounds = pars_parser(value) # Recurse into sub-dictionaries
            initial_pars[key] = sub_best
            par_bounds[key] = sub_bounds
        elif isinstance(value, np.ndarray):
            initial_pars[key] = value
            par_bounds[key] = value
        else:
            errormsg = f'Parameter "{key}" must be a list of [best, low, high] or a dictionary for nested parameters'
            raise ValueError(errormsg)
    return initial_pars, par_bounds

def unflatten_dict(flat_dict):
    '''
    A helper function to unflatten a dictionary with dot-separated keys into a nested dictionary.

    Args:
        flat_dict (dict): the flattened dictionary

    Returns:
        A nested dictionary
    '''
    nested_dict = {}
    for flat_key, value in flat_dict.items():
        keys = flat_key.split('.')
        d = nested_dict
        for key in keys[:-1]:
            if key not in d:
                d[key] = {}
            d = d[key]
        d[keys[-1]] = value
    return nested_dict

def unflatten_progs(init_pars, flat_progs):
    '''
    A helper function to unflatten a dictionary of prognosis parameters with dot-separated keys into a nested dictionary.

    Args:
        flat_progs (dict): the flattened prognosis parameters dictionary

    Returns:
        A nested prognosis parameters dictionary
    '''
    # First we need to discover what agent types are present in the flat_progs and retrieve the correct stratification array for each one.
    agent_types = []
    for flat_key in flat_progs.keys():
        keys = flat_key.split('.')
        if len(keys) != 4 or keys[0] != 'prognoses':
            errormsg = f'Invalid prognosis parameter key format: "{flat_key}". Expected format: "prognoses.<agent_type>.<par_name>.<strat_value>"'
            raise ValueError(errormsg)
        if keys[1] not in agent_types:
            agent_types.append(keys[1])
    strat_vectors = {}
    for agent_type in agent_types:
        agent_strat_keys = []
        for flat_key in flat_progs.keys():
            keys = flat_key.split('.')
            if keys[0] == 'prognoses' and keys[1] == agent_type:
                strat_key = keys[3]
                if strat_key not in agent_strat_keys:
                    agent_strat_keys.append(strat_key)
        for par_name, par_array in init_pars.get(agent_type, {}).items():
            if all(strat in par_array for strat in agent_strat_keys):
                strat_vectors[agent_type] = par_array
                break
        if agent_type not in strat_vectors:
            errormsg = f'Could not find stratification array for agent type "{agent_type}" in init_pars'
            raise ValueError(errormsg)
    # AI generated below here
    prognoses = {}
    for flat_key, value in flat_progs.items():
        keys = flat_key.split('.')
        if keys[0] != 'prognoses':
            continue
        agent_type = keys[1]
        par_name = keys[2]
        strat_value = keys[3]
        if agent_type not in prognoses:
            prognoses[agent_type] = {}
        if par_name not in prognoses[agent_type]:
            prognoses[agent_type][par_name] = np.empty_like(strat_vectors[agent_type], dtype=znd.default_float)
        # Find index of strat_value in the stratification array
        try:
            index = np.where(strat_vectors[agent_type] == strat_value)[0][0]
        except IndexError:
            errormsg = f'Stratification value "{strat_value}" not found in prognosis for agent type "{agent_type}"'
            raise ValueError(errormsg)
        prognoses[agent_type][par_name][index] = value
    return prognoses

def unflatten_pars(init_pars, flat_pars):
    '''
    A helper function to unflatten a dictionary of parameters with dot-separated keys into a nested dictionary.

    Args:
        flat_pars (dict): the flattened parameters dictionary

    Returns:
        A nested parameters dictionary
    '''
    final_pars = {}
    flat_prognoses = {k: v for k, v in flat_pars.items() if k.startswith('prognoses.')}
    flat_other = {k: v for k, v in flat_pars.items() if not k.startswith('prognoses.')}
    if flat_prognoses:
        final_pars['prognoses'] = unflatten_progs(init_pars.get('prognoses', {}), flat_prognoses)
    final_pars.update(unflatten_dict(flat_other))
    return final_pars