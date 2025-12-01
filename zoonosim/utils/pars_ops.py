'''
Utilities for working with nested parameter dictionaries.
'''


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
            errormsg = 'Calibration of prognoses parameters is not supported'
            raise ValueError(errormsg)
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

