import numba as nb # For faster computations
import numpy as np # For numerics
from .. import defaults as znd
from ..settings import options

__all__ = ['compute_viral_load', 'compute_trans_sus', 'compute_infections', 'find_contacts']

# Set dtypes -- note, these cannot be changed after import since Numba functions are precompiled
nbbool  = znd.nbbool
nbint   = znd.nbint
nbfloat = znd.nbfloat

# Specify whether to allow parallel Numba calculation -- 10% faster for safe and 20% faster for random, but the random number stream becomes nondeterministic for the latter

safe_parallel = options.numba_parallel in znd.safe_opts + znd.full_opts
rand_parallel = options.numba_parallel in znd.full_opts
if options.numba_parallel not in [0, 1, 2, '0', '1', '2', 'none', 'safe', 'full']:
    errormsg = f'Numba parallel must be "none", "safe", or "full", not "{options.numba_parallel}"'
    raise ValueError(errormsg)
cache = options.numba_cache # Turning this off can help switching parallelization options

@nb.njit(             (nbint,   nbfloat[:],     nbfloat[:],  nbfloat[:],    nbfloat,                    nbfloat,    nbfloat,    nbfloat), cache=cache, parallel=safe_parallel)
def compute_viral_load(t,       t_detectable,   t_peak,      t_recovered,   minimum_detectable_load,    peak_load,  min_scl,     max_scl): # pragma: no cover
    '''
    Calculate viral load for infectious humans?

    Args:
        t: (int) timestep
        t_detectable: (float) time when viral load becomes detectable
        t_peak: (float) time when viral load reaches its peak
        t_recovered: (float) time when viral load drops to recovered level
        minimum_detectable_load: (float) minimum detectable viral load
        peak_load: (float) peak viral load
        min_scl: (float) relative transmissibility scaling factor for minimum viral load
        max_scl: (float) relative transmissibility scaling factor for peak viral load

    Returns:
        viral_load (float): viral load
    '''

    # Set up arrays
    vl = np.zeros(len(t_detectable), dtype=znd.default_float)
    vl_rescaled = np.zeros(len(t_detectable), dtype=znd.default_float)

    # Calculate viral load for those for whom it is rising
    rising_vl = t < t_peak
    vl[rising_vl] = minimum_detectable_load + (peak_load - minimum_detectable_load)*(t - t_detectable[rising_vl])/(t_peak[rising_vl] - t_detectable[rising_vl])

    # Calculate viral load for those for whom it is falling
    falling_vl = t >= t_peak
    vl[falling_vl] = minimum_detectable_load + (peak_load - minimum_detectable_load)*(t - t_peak[falling_vl])/(t_recovered[falling_vl] - t_peak[falling_vl])

    # Rescale viral load to represent relative transmissibility
    vl_rescaled = min_scl + (max_scl - min_scl)*(vl - minimum_detectable_load)/(peak_load - minimum_detectable_load)
    vl_rescaled[vl_rescaled < min_scl] = min_scl #Not sure if this is necessary, but it seems like a good idea to prevent negative values.
    vl_rescaled[vl_rescaled > max_scl] = max_scl #Not sure if this is necessary, but it seems like a good idea to prevent values above max_scl.

    # Clip viral load when it falls below 10^0 cp/mL to reflect LoD
    vl[vl <= 0] = 0

    return vl, vl_rescaled


# jit means you let Numba's combiler optimize this function. 
@nb.njit(            (nbfloat[:], nbfloat[:], nbbool[:], nbbool[:], nbfloat,    nbfloat[:], nbbool[:], nbbool[:], nbfloat[:],  nbfloat,     nbfloat[:]), cache=cache, parallel=safe_parallel)
def compute_trans_sus(rel_trans,  rel_sus,    inf,       sus,       beta_layer, viral_load, symp,      quar,      asymp_factor,  quar_factor, immunity_factors): # pragma: no cover
    ''' Calculate relative transmissibility and susceptibility '''
    f_asymp   =  symp + ~symp * asymp_factor # Asymptomatic factor, changes e.g. [0,1] with a factor of 0.8 to [0.8,1.0]
    f_quar    = ~quar +  quar * quar_factor # Quarantine, changes e.g. [0,1] with a factor of 0.5 to [1,0.5] 
    rel_trans = rel_trans * inf * f_quar * f_asymp * beta_layer * viral_load # Recalculate transmissibility
    rel_sus   = rel_sus * sus * f_quar * (1-immunity_factors) # Recalculate susceptibility 
    return rel_trans, rel_sus


@nb.njit(             (nbfloat[:],  nbint[:],  nbint[:], nbfloat[:],   nbfloat[:], nbfloat[:]), cache=cache, parallel=rand_parallel)
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
        betas            = beta[targets[inf_inds]] * layer_betas[inf_inds] * source_trans[inf_inds] * rel_sus[targets[inf_inds]] # Calculate the raw transmission probabilities
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