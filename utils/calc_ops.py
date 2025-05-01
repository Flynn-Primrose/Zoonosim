import numba as nb # For faster computations
import numpy as np # For numerics
from .. import defaults as znd

__all__ = ['compute_viral_load', 'compute_infection_level', 'compute_trans_sus', 'compute_infections', 'find_contacts']

# Set dtypes -- note, these cannot be changed after import since Numba functions are precompiled
nbbool  = znd.nbbool
nbint   = znd.nbint
nbfloat = znd.nbfloat

# Specify whether to allow parallel Numba calculation -- 10% faster for safe and 20% faster for random, but the random number stream becomes nondeterministic for the latter

safe_parallel = znd.numba_parallel in znd.safe_opts + znd.full_opts
rand_parallel = znd.numba_parallel in znd.full_opts
if znd.numba_parallel not in [0, 1, 2, '0', '1', '2', 'none', 'safe', 'full']:
    errormsg = f'Numba parallel must be "none", "safe", or "full", not "{znd.numba_parallel}"'
    raise ValueError(errormsg)
cache = znd.numba_cache # Turning this off can help switching parallelization options

@nb.njit(             (nbint,   nbfloat[:], nbfloat[:], nbfloat[:], nbfloat[:], nbfloat[:], nbfloat[:], nbfloat,    nbfloat), cache=cache, parallel=safe_parallel)
def compute_viral_load(t,       x_p1,       y_p1,       x_p2,       y_p2,       x_p3,       y_p3,       min_vl,     max_vl): # pragma: no cover
    '''
    Calculate viral load for infectious humans?

    Args:
        t: (int) timestep
        x_p1: (float[]) date of onset of infectiousness
        y_p1: (float[]) viral load (3 cp/mL)
        x_p2: (float[]) date of peak viral load
        y_p2: (float[]) peak viral load
        x_p3: (float[]) date of hitting viral load of 6 cp/mL
        y_p3: (float[]) viral load (6 cp/mL)

    Returns:
        viral_load (float): viral load
    '''

    # Set up arrays
    N = len(x_p2)
    vl = np.zeros(N, dtype=znd.default_float)
    vl_rescaled = np.zeros(N, dtype=znd.default_float)

    # Calculate viral load for those for whom it is rising
    rising_vl = t < x_p2
    vl[rising_vl] = y_p1[rising_vl] + (y_p2[rising_vl] - y_p1[rising_vl])*(t - x_p1[rising_vl])/(x_p2[rising_vl] - x_p1[rising_vl])

    # Calculate viral load for those for whom it is falling
    falling_vl = t >= x_p2
    vl[falling_vl] = y_p2[falling_vl] + (y_p3[falling_vl] - y_p2[falling_vl])*(t - x_p2[falling_vl])/(x_p3[falling_vl] - x_p2[falling_vl])

    # Rescale viral load for Covasim computation
    # NOTE: I think I need to talk to Ritchie about this -- I don't understand how this works
    infected = ~np.isnan(x_p2)
    infected_vl = vl[infected]
    infected_vl = min_vl + (max_vl - min_vl)*(infected_vl - 6)/(11 - 6) 
    infected_vl[infected_vl < min_vl] = 0
    vl_rescaled[infected] = infected_vl

    # Clip viral load when it falls below 10^0 cp/mL to reflect LoD
    vl[vl <= 0] = 0

    return vl, vl_rescaled

@nb.njit(                   (nbint,   nbfloat[:], nbfloat[:], nbfloat[:], nbfloat[:], nbfloat[:], nbfloat[:], nbfloat[:]), cache=cache, parallel=safe_parallel)
def compute_infection_level(t,       x_p1,       y_p1,       x_p2,       y_p2,       x_p3,       y_p3,       headcount): # pragma: no cover
    '''
    Calculate infection levels for infected flocks?

        Args:
        t: (int) timestep
        x_p1: (float[]) date of first infection
        y_p1: (float[]) initial number of infected birds
        x_p2: (float[]) date of peak infection
        y_p2: (float[]) number of infected birds at peak infection
        x_p3: (float[]) date of equilibrium infection
        y_p3: (float[]) number of birds infected at equilibrium
        headcount: (int[]) total number of birds in the flock

    Returns:
        infected_headcount (int): the number of infected birds in the flock
        infected_headcount_rescaled (int): the number of infected birds in the flock, rescaled to a fraction of total headcount
    '''
        # Set up arrays
    N = len(x_p2)
    il = np.zeros(N, dtype=znd.default_float)
    il_rescaled = np.zeros(N, dtype=znd.default_float)

    # Calculate viral load for those for whom it is rising
    rising_il = t < x_p2
    il[rising_il] = y_p1[rising_il] + (y_p2[rising_il] - y_p1[rising_il])*(t - x_p1[rising_il])/(x_p2[rising_il] - x_p1[rising_il])

    # Calculate viral load for those for whom it is falling
    falling_il = t >= x_p2
    il[falling_il] = y_p2[falling_il] + (y_p3[falling_il] - y_p2[falling_il])*(t - x_p2[falling_il])/(x_p3[falling_il] - x_p2[falling_il])

    il_rescaled = il / headcount # Rescale to fraction of total headcount
    il_rescaled[il_rescaled < 0] = 0 # Clip to 0
    il_rescaled[il_rescaled > 1] = 1 # Clip to 1

    return il, il_rescaled

# jit means you let Numba's combiler optimize this function. 
@nb.njit(            (nbfloat[:], nbfloat[:], nbbool[:], nbbool[:], nbfloat,    nbfloat[:], nbbool[:], nbbool[:], nbfloat,  nbfloat,     nbfloat[:]), cache=cache, parallel=safe_parallel)
def compute_trans_sus(rel_trans,  rel_sus,    inf,       sus,       beta_layer, viral_load, symp,      quar,      asymp_factor,  quar_factor, immunity_factors): # pragma: no cover
    ''' Calculate relative transmissibility and susceptibility '''
    f_asymp   =  symp + ~symp * asymp_factor # Asymptomatic factor, changes e.g. [0,1] with a factor of 0.8 to [0.8,1.0]
    f_quar    = ~quar +  quar * quar_factor # Quarantine, changes e.g. [0,1] with a factor of 0.5 to [1,0.5] 
    rel_trans = rel_trans * inf * f_quar * f_asymp * beta_layer * viral_load # Recalculate transmissibility
    rel_sus   = rel_sus * sus * f_quar * (1-immunity_factors) # Recalculate susceptibility 
    return rel_trans, rel_sus


@nb.njit(             (nbfloat,  nbint[:],  nbint[:], nbfloat[:],   nbfloat[:], nbfloat[:]), cache=cache, parallel=rand_parallel)
def compute_infections(beta,     p1,        p2,       layer_betas,  rel_trans,  rel_sus): # pragma: no cover
    '''
    Compute who infects whom

    The heaviest step of the model -- figure out who gets infected on this timestep.
    Cannot be easily parallelized since random numbers are used. Loops over contacts
    in both directions (i.e., targets become sources).

    Args:
        beta: overall transmissibility
        p1: person 1
        p2: person 2
        layer_betas: per-contact transmissibilities
        rel_trans: the source's relative transmissibility
        rel_sus: the target's relative susceptibility
    '''
    slist = np.empty(0, dtype=nbint)
    tlist = np.empty(0, dtype=nbint)
    pairs = [[p1,p2], [p2,p1]]
    for sources,targets in pairs:
        source_trans     = rel_trans[sources] # Pull out the transmissibility of the sources (0 for non-infectious people). Cx1 array. 
        inf_inds         = source_trans.nonzero()[0] # Infectious indices -- remove noninfectious people. Smaller array of the non-zero inds.
        # betas: I x 1 array, I is # contacts where P1 is infectious. 
        betas            = beta * layer_betas[inf_inds] * source_trans[inf_inds] * rel_sus[targets[inf_inds]] # Calculate the raw transmission probabilities
        # There will be zeros intoduced, for example if someone isn't susciptible due to an intervention. 
        nonzero_inds     = betas.nonzero()[0] # Find nonzero entries
        nonzero_inf_inds = inf_inds[nonzero_inds] # Map onto original indices
        nonzero_betas    = betas[nonzero_inds] # Remove zero entries from beta
        nonzero_sources  = sources[nonzero_inf_inds] # Remove zero entries from the sources
        nonzero_targets  = targets[nonzero_inf_inds] # Remove zero entries from the targets
        transmissions    = (np.random.random(len(nonzero_betas)) < nonzero_betas).nonzero()[0] # Compute the actual infections!
        source_inds      = nonzero_sources[transmissions]
        target_inds      = nonzero_targets[transmissions] # Filter the targets on the actual infections
        slist = np.concatenate((slist, source_inds), axis=0)
        tlist = np.concatenate((tlist, target_inds), axis=0)
    return slist, tlist


@nb.njit((nbint[:], nbint[:], nb.int64[:]), cache=cache)
def find_contacts(p1, p2, inds): # pragma: no cover
    """
    Numba for Layer.find_contacts()

    A set is returned here rather than a sorted array so that custom tracing interventions can efficiently
    add extra people. For a version with sorting by default, see Layer.find_contacts(). Indices must be
    an int64 array since this is what's returned by true() etc. functions by default.
    """
    pairing_partners = set()
    inds = set(inds)
    for i in range(len(p1)):
        if p1[i] in inds:
            pairing_partners.add(p2[i])
        if p2[i] in inds:
            pairing_partners.add(p1[i])
    return pairing_partners