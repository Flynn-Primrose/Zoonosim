import numba as nb # For faster computations
import numpy as np # For numerics
import random # Used only for resetting the seed
import sciris as sc # For additional utilities
import pandas as pd
import itertools
import bisect
import scipy
from scipy import stats as st
import warnings
from .. import defaults as znd
from ..settings import options
from . import stats_ops as znso # For additional statistical operations

__all__ = ['sample', 'get_pdf', 'set_seed', 'fast_choice',
           'sample_single_dict', 'sample_single_arr', 'resample_age',
           'sample_from_range', 'check_dist', 'check_normal', 'check_poisson',
           'check_truncated_poisson', 'statistic_test']

# Set dtypes -- note, these cannot be changed after import since Numba functions are precompiled
nbbool  = nb.bool_
nbint   = znd.nbint
nbfloat = znd.nbfloat

#safe_opts = [1, '1', 'safe'] # TODO: Move this to config
#full_opts = [2, '2', 'full'] # TODO: Move this to config
safe_parallel = options.numba_parallel in znd.safe_opts + znd.full_opts
rand_parallel = options.numba_parallel in znd.full_opts
if options.numba_parallel not in [0, 1, 2, '0', '1', '2', 'none', 'safe', 'full']:
    errormsg = f'Numba parallel must be "none", "safe", or "full", not "{options.numba_parallel}"'
    raise ValueError(errormsg)
cache = options.numba_cache # Turning this off can help switching parallelization options

def sample(dist=None, par1=None, par2=None, size=None, **kwargs):
    '''
    Draw a sample from the distribution specified by the input. The available
    distributions are:

    - 'uniform'       : uniform distribution from low=par1 to high=par2; mean is equal to (par1+par2)/2
    - 'normal'        : normal distribution with mean=par1 and std=par2
    - 'lognormal'     : lognormal distribution with mean=par1 and std=par2 (parameters are for the lognormal distribution, *not* the underlying normal distribution)
    - 'normal_pos'    : right-sided normal distribution (i.e. only positive values), with mean=par1 and std=par2 *of the underlying normal distribution*
    - 'normal_int'    : normal distribution with mean=par1 and std=par2, returns only integer values
    - 'lognormal_int' : lognormal distribution with mean=par1 and std=par2, returns only integer values
    - 'poisson'       : Poisson distribution with rate=par1 (par2 is not used); mean and variance are equal to par1
    - 'neg_binomial'  : negative binomial distribution with mean=par1 and k=par2; converges to Poisson with k=∞

    Args:
        dist (str):   the distribution to sample from
        par1 (float): the "main" distribution parameter (e.g. mean)
        par2 (float): the "secondary" distribution parameter (e.g. std)
        size (int):   the number of samples (default=1)
        kwargs (dict): passed to individual sampling functions

    Returns:
        A length N array of samples

    **Examples**::

        cv.sample() # returns Unif(0,1)
        cv.sample(dist='normal', par1=3, par2=0.5) # returns Normal(μ=3, σ=0.5)
        cv.sample(dist='lognormal_int', par1=5, par2=3) # returns a lognormally distributed set of values with mean 5 and std 3

    Notes:
        Lognormal distributions are parameterized with reference to the underlying normal distribution (see:
        https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.random.lognormal.html), but this
        function assumes the user wants to specify the mean and std of the lognormal distribution.

        Negative binomial distributions are parameterized with reference to the mean and dispersion parameter k
        (see: https://en.wikipedia.org/wiki/Negative_binomial_distribution). The r parameter of the underlying
        distribution is then calculated from the desired mean and k. For a small mean (~1), a dispersion parameter
        of ∞ corresponds to the variance and standard deviation being equal to the mean (i.e., Poisson). For a
        large mean (e.g. >100), a dispersion parameter of 1 corresponds to the standard deviation being equal to
        the mean.
    '''

    # Some of these have aliases, but these are the "official" names
    choices = [
        'uniform',
        'normal',
        'normal_pos',
        'normal_int',
        'lognormal',
        'lognormal_int',
        'poisson',
        'neg_binomial',
    ]

    # Ensure it's an integer
    if size is not None:
        size = int(size)

    # Compute distribution parameters and draw samples
    # NB, if adding a new distribution, also add to choices above
    if   dist in ['unif', 'uniform']: samples = np.random.uniform(low=par1, high=par2, size=size, **kwargs)
    elif dist in ['norm', 'normal']:  samples = np.random.normal(loc=par1, scale=par2, size=size, **kwargs)
    elif dist == 'normal_pos':        samples = np.abs(np.random.normal(loc=par1, scale=par2, size=size, **kwargs))
    elif dist == 'normal_int':        samples = np.round(np.abs(np.random.normal(loc=par1, scale=par2, size=size, **kwargs)))
    elif dist == 'poisson':           samples = znso.n_poisson(rate=par1, n=size, **kwargs) # Use Numba version below for speed
    elif dist == 'neg_binomial':      samples = znso.n_neg_binomial(rate=par1, dispersion=par2, n=size, **kwargs) # Use custom version below
    elif dist in ['lognorm', 'lognormal', 'lognorm_int', 'lognormal_int']:
        if par1>0:
            mean  = np.log(par1**2 / np.sqrt(par2**2 + par1**2)) # Computes the mean of the underlying normal distribution
            sigma = np.sqrt(np.log(par2**2/par1**2 + 1)) # Computes sigma for the underlying normal distribution
            samples = np.random.lognormal(mean=mean, sigma=sigma, size=size, **kwargs)
        else:
            samples = np.zeros(size)
        if '_int' in dist:
            samples = np.round(samples)
    else:
        errormsg = f'The selected distribution "{dist}" is not implemented; choices are: {sc.newlinejoin(choices)}'
        raise NotImplementedError(errormsg)

    return samples



def get_pdf(dist=None, par1=None, par2=None):
    '''
    Return a probability density function for the specified distribution. This
    is used for example by test_num to retrieve the distribution of times from
    symptom-to-swab for testing. For example, for Washington State, these values
    are dist='lognormal', par1=10, par2=170.
    '''
    import scipy.stats as sps # Import here since slow

    choices = [
        'none',
        'uniform',
        'lognormal',
    ]

    if dist in ['None', 'none', None]:
        return None
    elif dist == 'uniform':
        pdf = sps.uniform(loc=par1, scale=par2)
    elif dist == 'lognormal':
        mean  = np.log(par1**2 / np.sqrt(par2 + par1**2)) # Computes the mean of the underlying normal distribution
        sigma = np.sqrt(np.log(par2/par1**2 + 1)) # Computes sigma for the underlying normal distribution
        pdf   = sps.lognorm(sigma, loc=-0.5, scale=np.exp(mean))
    else:
        choicestr = '\n'.join(choices)
        errormsg = f'The selected distribution "{dist}" is not implemented; choices are: {choicestr}'
        raise NotImplementedError(errormsg)

    return pdf


def set_seed(seed=None):
    '''
    Reset the random seed -- complicated because of Numba, which requires special
    syntax to reset the seed. This function also resets Python's built-in random
    number generated.

    Args:
        seed (int): the random seed
    '''

    @nb.njit((nbint,), cache=cache)
    def set_seed_numba(seed):
        return np.random.seed(seed)

    def set_seed_regular(seed):
        return np.random.seed(seed)

    # Dies if a float is given
    if seed is not None:
        seed = int(seed)

    set_seed_regular(seed) # If None, reinitializes it
    if seed is None: # Numba can't accept a None seed, so use our just-reinitialized Numpy stream to generate one
        seed = np.random.randint(1e9)
    set_seed_numba(seed)
    random.seed(seed) # Finally, reset Python's built-in random number generator, just in case (used by SynthPops)

    return

def fast_choice(weights):
    """
    Choose an option -- quickly -- from the provided weights. Weights do not need
    to be normalized.

    Reimplementation of random.choices(), removing everything inessential.

    Example:
        fast_choice([0.1,0.2,0.3,0.2,0.1]) # might return 2
    """
    cum_weights = list(itertools.accumulate(weights))
    return bisect.bisect(cum_weights, random.random()*(cum_weights[-1]), 0, len(cum_weights)-1)


# @nb.njit(cache=True)
def sample_single_dict(distr_keys, distr_vals):
    """
    Sample from a distribution.

    Args:
        distr (dict or np.ndarray): distribution

    Returns:
        A single sampled value from a distribution.

    """
    return distr_keys[fast_choice(distr_vals)]


def sample_single_arr(distr):
    """
    Sample from a distribution.

    Args:
        distr (dict or np.ndarray): distribution

    Returns:
        A single sampled value from a distribution.
    """
    return fast_choice(distr)


def resample_age(age_dist_vals, age):
    """
    Resample age from single year age distribution.

    Args:
        single_year_age_distr (arr) : age distribution, ordered by age
        age (int)                   : age as an integer
    Returns:
        Resampled age as an integer.
    """
    if age == 0:
        age_min = 0
        age_max = 1
    elif age == 1:
        age_min = 0
        age_max = 2
    elif age >= 2 and age <= 98:
        age_min = age - 2
        age_max = age + 2
    elif age == 99:
        age_min = 97
        age_max = 99
    else:
        age_min = 98
        age_max = 100

    age_distr = age_dist_vals[age_min:age_max + 1]  # create an array of the values, not yet normalized
    age_range = np.arange(age_min, age_max + 1)
    return age_range[fast_choice(age_distr)]


__all__ += ['norm_dic', 'norm_age_group']


def norm_dic(dic):
    """
    Normalize the dictionary ``dic``.

    Args:
        dic (dict): A dictionary with numerical values.

    Returns:
        A normalized dictionary.
    """
    total = float(sum(dic.values()))
    if total == 0.0:
        return dic
    return {k: v / total for k, v in dic.items()}

def norm_age_group(age_dic, age_min, age_max):
    """
    Create a normalized dictionary for the range ``age_min`` to ``age_max``, inclusive.

    Args:
        age_dic (dict) : A dictionary with numerical values.
        age_min (int)  : The minimum value of the range for the dictionary.
        age_max (int)  : The maximum value of the range for the dictionary.

    Returns:
        A normalized dictionary for keys in the range ``age_min`` to ``age_max``, inclusive.
    """
    dic = {a: age_dic[a] for a in range(age_min, age_max + 1)}
    return norm_dic(dic)

def sample_from_range(distr, min_val, max_val):
    """
    Sample from a distribution from min_val to max_val, inclusive.

    Args:
        distr (dict)  : distribution with integer keys
        min_val (int) : minimum of the range to sample from
        max_val (int) : maximum of the range to sample from
    Returns:
        A sampled number from the range min_val to max_val in the distribution distr.
    """
    new_distr = norm_age_group(distr, min_val, max_val)
    distr_keys = np.array(list(new_distr.keys()), dtype=np.int64)
    distr_vals = np.array(list(new_distr.values()), dtype=np.float64)
    return sample_single_dict(distr_keys, distr_vals)


def check_dist(actual, expected, std=None, dist='norm', check='dist', label=None, alpha=0.05, size=10000, verbose=True, die=False, stats=False):
    """
    Check whether counts match the expected distribution. The distribution can be
    any listed in scipy.stats. The parameters for the distribution should be supplied
    via the "expected" argument. The standard deviation for a normal distribution is
    a special case; it can be supplied separately or calculated from the (actual) data.

    Args:
        actual (int, float, or array) : the observed value, or distribution of values
        expected (int, float, tuple)  : the expected value; or, a tuple of arguments
        std (float)                   : for normal distributions, the standard deviation of the expected value (taken from data if not supplied)
        dist (str)                    : the type of distribution to use
        check (str)                   : what to check: 'dist' = entire distribution (default), 'mean' (equivalent to supplying np.mean(actual)), or 'median'
        label (str)                   : the name of the variable being tested
        alpha (float)                 : the significance level at which to reject the null hypothesis
        size (int)                    : the size of the sample from the expected distribution to compare with if distribution is discrete

        verbose (bool)                : print a warning if the null hypothesis is rejected
        die (bool)                    : raise an exception if the null hypothesis is rejected
        stats (bool)                  : whether to return statistics

    Returns:
        If stats is True, returns statistics: whether null hypothesis is
        rejected, pvalue, number of samples, expected quintiles, observed
        quintiles, and the observed quantile.

    **Examples**::

        sp.check_dist(actual=[3,4,4,2,3], expected=3, dist='poisson')
        sp.check_dist(actual=[0.14, -3.37,  0.59, -0.07], expected=0, std=1.0, dist='norm')
        sp.check_dist(actual=5.5, expected=(1, 5), dist='lognorm')
    """
    # Handle inputs
    label = f' "{label}"' if label else ''
    is_dist = sc.isiterable(actual)

    # Set distribution
    if dist.lower() in ['norm', 'normal', 'gaussian']:
        if std is None:
            if is_dist:
                std = np.std(actual)  # Get standard deviation from the data
            else: # pragma: no cover
                std = 1.0
        args = (expected, std)
        scipydist = getattr(scipy.stats, 'norm')
        truedist = scipy.stats.norm(expected, std)
    else:
        try:
            if sc.isnumber(expected):
                args = (expected, )
            else:
                args = tuple(expected)
            scipydist = getattr(scipy.stats, dist)
            truedist = scipydist(*args)
        except Exception as E:
            errormsg = f'Distribution "{dist}" not supported with the expected values supplied; valid distributions are those in scipy.stats'
            raise NotImplementedError(errormsg) from E

    # Calculate stats
    if is_dist and check == 'dist':
        quantile = truedist.cdf(np.median(actual))

        # only if distribution is continuous
        if isinstance(scipydist, scipy.stats.rv_continuous):
            teststat, pvalue = scipy.stats.kstest(rvs=actual, cdf=dist, args=args)  # Use the K-S test to see if came from the same distribution

        # ks test against large sample from the theoretical distribution
        elif isinstance(scipydist, scipy.stats.rv_discrete):
            expected_r = truedist.rvs(size=size)
            teststat, pvalue = scipy.stats.ks_2samp(actual, expected_r)

        else: # pragma: no cover
            errormsg = 'Distribution is neither continuous or discrete and so not supported at this time.'
            raise NotImplementedError(errormsg)
        null = pvalue > alpha

    else:
        if check == 'mean':
            value = np.mean(actual)
        elif check == 'median':
            value = np.median(actual)
        else:
            value = actual
        quantile = truedist.cdf(value)  # If it's a single value, see where it lands on the Poisson CDF
        pvalue = 1.0-2*abs(quantile-0.5)  # E.g., 0.975 maps on to p=0.05
        minquant = alpha/2  # e.g., 0.025 for alpha=0.05
        maxquant = 1-alpha/2  # e.g., 0.975 for alpha=0.05
        minval = truedist.ppf(minquant)
        maxval = truedist.ppf(maxquant)
        quant_check = (minquant <= quantile <= maxquant)  # True if above minimum and below maximum
        val_check = (minval <= value <= maxval)  # Check values
        null = quant_check or val_check  # Consider it to pass if either passes

    # Additional stats
    n_samples = len(actual) if is_dist else 1
    eps = 1.0/n_samples if n_samples > 4 else 1e-2  # For small number of samples, use default limits
    quintiles = [eps, 0.25, 0.5, 0.75, 1-eps]
    obvs_quin = np.quantile(actual, quintiles) if is_dist else actual
    expect_quin = truedist.ppf(quintiles)

    # If null hypothesis is rejected, print a warning or error
    if not null:
        msg = f''''
Variable{label} with n={n_samples} samples is out of range using the distribution:
    {dist}({args}) →
    p={pvalue} < α={alpha}
Expected quintiles are: {expect_quin}
Observed quintiles are: {obvs_quin}
Observed median is in quantile: {quantile}'''
        if die:
            raise ValueError(msg)
        elif verbose:
            warnings.warn(msg)

    # If null hypothesis is not rejected, under verbose, print a confirmation
    if null and verbose:
        print(f'Check passed. Null hypothesis with expected distribution: {dist}{args} not rejected.')
        if is_dist and check == 'dist':
            print(f'Test statistic: {teststat}, pvalue: {pvalue}')

    if not stats:
        return null
    else:
        s = sc.objdict()
        s.null = null
        s.pvalue = pvalue
        s.n_samples = n_samples
        s.expected_quintiles = expect_quin
        s.observed_quintiles = obvs_quin
        s.observed_quantile = quantile
        return s


def check_normal(*args, **kwargs):
    ''' Alias to check_dist(dist='normal') '''
    dist = kwargs.pop('dist', 'norm')
    return check_dist(*args, **kwargs, dist=dist)


def check_poisson(*args, **kwargs):
    ''' Alias to check_dist(dist='poisson') '''
    dist = kwargs.pop('dist', 'poisson')
    return check_dist(*args, **kwargs, dist=dist)


def check_truncated_poisson(testdata, mu, lowerbound=None, upperbound=None, skipcheck=False, **kwargs):
    """
    test if data fits in truncated poisson distribution between upperbound and lowerbound using kstest
    Args:
        testdata (array) : data to be tested
        mu (float) : expected mean for the poisson distribution
        lowerbound (float) : lowerbound for truncation
        upperbound (float) : upperbound for truncation

    Returns:
        (bool) return True if statistic check passed, else return False
    """
    sample_size = len(testdata)
    # need to exclude any value below or equal to lowerbound and any value above or equal to upperbound, so we first find the quantile location for
    # lowerbound and upperbound then only generate poisson cdf values in between these 2 locations
    minquantile = st.poisson.cdf(lowerbound, mu=mu) if lowerbound else 0
    maxquantile = st.poisson.cdf(upperbound, mu=mu) if upperbound else 1
    # create uniformly distributed number between minquantile and maxquantile (in the cdf quantile space)
    q = np.random.uniform(low=minquantile, high=maxquantile, size=sample_size)
    # use percent point function to get inverse of cdf
    expected_data = st.poisson.ppf(q=q, mu=mu)
    result = True
    if not skipcheck:
        try:
            statistic_test(expected_data, testdata, test=st.kstest, verbose=True, die=True)
            result = True
        except ValueError as e:
            if 'reject the hypothesis' in str(e):
                result = False
            else:
                raise Exception(e)

    #plot comparison
    bins_count = int(min(10, max(expected_data)-min(expected_data)+1))
    print(f"bins:{bins_count}")

    expected, bins = np.histogram(expected_data, bins=bins_count)
    actual = np.histogram(testdata, bins=bins)[0]
    kwargs["generated"] = actual
    #merge 11 bins to 10 for bar plot align at center
    merged_bins = [round((bins[idx] + bins[idx+1])/2,2) for idx, val in enumerate(bins) if idx < len(bins)-1]
    kwargs["xvalue"] = merged_bins
    return result


def statistic_test(expected, actual, test=st.chisquare, verbose=True, die=False, **kwargs):
    """
    Perform statistical checks for expected and actual data based on the null
    hypothesis that expected and actual distributions are identical. Throw
    assertion if the expected and actual data differ significantly based on the
    test selected.
    See https://docs.scipy.org/doc/scipy/reference/stats.html#statistical-tests.

    Args:
        expected (array)    : the expected value; or, a tuple of arguments
        actual (array)      : the observed value, or distribution of values
        test (scipy.stats)  : scipy statistical tests functions, for example scipy.stats.chisquare
        verbose (bool)      : print a warning if the null hypothesis is rejected
        die (bool)          : raise an exception if the null hypothesis is rejected
        **kwargs (dict)     : optional arguments for statistical tests

    Returns:
        None.
    """
    # data = {'expected': expected, 'actual': actual}
    df_expected = pd.DataFrame(expected, columns=['expected'])
    df_actual = pd.DataFrame(actual, columns=['actual'])
    details = pd.merge(df_expected.describe(), df_actual.describe(), left_index=True, right_index=True, suffixes=('expected', 'actual'))

    print(f"use {str(test.__name__)} to check actual distribution")
    s, p = test(expected, actual, **kwargs)
    print(f"statistics: {s} pvalue:{p}")
    if p < 0.05: # pragma: no cover
        msg = f"Under the null hypothesis the expected/actual distributions are identical." \
              f"If statistics is small or the p-value is high (greater than the significance level 5%)," \
              f" then we cannot reject the hypothesis. But we got p = {p} and s = {s}."
        if die: # pragma: no cover
            raise ValueError(msg)
        elif verbose: # pragma: no cover
            warnings.warn(f"data: \n{details}")
            warnings.warn(msg)