import numpy as np
import copy

from itertools import combinations
from functools import reduce
from .background import *

'''
This file implements functions for checking test-allocation.
'''

def allocate_tests(sim, testobj, seekers_dict=None, mode='proportional'): 

    if mode == 'proportional': 
        if isinstance(testobj.capacity, np.ndarray) or isinstance(testobj.capacity, list):
            capacity = int(testobj.capacity[sim.t] * sim.pars['pop_size_by_type']['human'])
        else: 
            capacity = int(testobj.capacity * sim.pars['pop_size_by_type']['human'])
        return proportional_allocation(capacity, seekers_dict)

    elif mode == 'basic': 
        if isinstance(testobj.capacity, np.ndarray) or isinstance(testobj.capacity, list):
            capacity = int(testobj.capacity[sim.t] * sim.pars['pop_size_by_type']['human'])
        else: 
            capacity = int(testobj.capacity * sim.pars['pop_size_by_type']['human'])
        return basic_allocation(capacity, seekers_dict)

    else: 
        raise RuntimeError("Invalid allocation policy supplied. Options are: 'proportional', 'basic'")

def basic_allocation(test_capacity, seekers_dict): 
    '''
    Each day, get people who are seeking a PCR test. Randomly determine who gets a test. There is no temporal dimension, people do not stay in any sort of queue if 
    they do not get a test today.  

    Args: 
        test_capacity    (int)       : Number of tests available to allocate
        seekers_dict     (dict)      : Dictionary of people from each criterion who is seeking a test

    Returns: 
        testers          (np.array)  : Array of tester UIDs, where max(length) = test capacity
    '''
    # Unpack test seekers
    test_seekers = set()  # Need uniqueness since some test seekers may satisfy multiple criteria, and we do not want to count them twice
    for value in seekers_dict.values():
        test_seekers.update(value)

    # Determine who actually gets a test
    if len(test_seekers) <= test_capacity:  # All test-seekers can get tested (very unlikely)
        testers = np.fromiter(test_seekers, dtype=int)
    else: 
        testers = np.random.choice(list(test_seekers), test_capacity, replace=False)  # Insufficient capacity, randomly choose a lucky few test-seekers

    return testers

def proportional_allocation(test_capacity, seekers_dict): 
    '''
    Compute the power set of criteria.

    For each power set element:
    (1) Get N, the number of people who precisely meet these criteria.
    (2) Let p = N / S, where S is population size.
    (3) S is the MAXIMUM proportion of tests which will be allocated to those associated with the power set element.

    Note that the power set contains the empty set. This corresponds to people who did not satisfy any criteria, 
    which is untrue of any test seeker.

    Args: 
        test_capacity   (int) : Number of tests to distribute on this day 
        seekers_dict    (dict) : Keys are criteria, values are test-seekers who meet the criteria
    '''
    valid_criteria = set([criterion for criterion in seekers_dict.keys()])

    criteria_powset    = compute_powset(valid_criteria)
    seeker_crit_subset = compute_members(valid_criteria, criteria_powset, seekers_dict)  # get test seekers and the powerset elements they correspond to
    test_proportions   = compute_proportions(seeker_crit_subset)  # compute proportions and allocate

    testers = set()
    for sname in test_proportions.keys(): 
        subset_seekers = np.array(list(seeker_crit_subset[sname]))  # Test seekers who meet the criteria described by subset
        subset_n_tests = int(test_proportions[sname] * test_capacity)

        if subset_n_tests < len(subset_seekers):
            subset_testers = set(np.random.choice(subset_seekers, subset_n_tests, replace=False).tolist())
        else: 
            subset_testers = subset_seekers

        testers.update(subset_testers)

    # Redistribute extras by randomly sampling subsets according to proportions
    extra_tests = test_capacity - len(testers)

    while extra_tests > 0: 

        total_not_tested = sum([len(seeker_crit_subset[sname].difference(testers)) for sname in test_proportions.keys()])

        if total_not_tested == 0: 
            break

        probs = [test_proportions[sname] for sname in test_proportions.keys()]
        extra_sname = np.random.choice(list(test_proportions.keys()), 1, p=probs)[0]
        not_tested = seeker_crit_subset[extra_sname].difference(testers)

        if len(not_tested) > 0: 
            testers.add(np.random.choice(list(not_tested)))
            extra_tests -= 1

    return np.fromiter(testers, dtype=int)

def compute_members(valid_criteria, criteria_powset, criteria_dict): 
    '''
    Given a powerset of criteria, figure out which people precisely meet each subset of criteria (no more or less). 

    Args: 
        valid_criteria      (set)  : Set of valid criteria that was used to compute the power set
        criteria_powset     (dict) : Keys are subset names (e.g. 'subset1'), values are the actual criteria subsets
        criteria_dict       (dict) : Keys are criteria, values are sets of people meeting those criteria

    Returns: 
        meets_crit_subsets  (dict) : Keys are subset names, values are people who precisely meet the criteria described in the subset
    '''
    # Get the number of people who precisely meet the criteria specified by each power set element.
    meets_crit_subset = dict()
    for subset_key in criteria_powset.keys(): 
        crit_subset_copy = copy.deepcopy(criteria_powset[subset_key])  # Criteria we want seekers to meet. A copied object avoids bugs with referencing when we do pop(). 
        complement = valid_criteria.difference(crit_subset_copy)  # Criteria we don't want seekers to meet. crit_subset_copy should be a subset of valid_criteria.

        members = criteria_dict[crit_subset_copy.pop()]  # Choose any criterion to start with, pop() can only happen after line above 

        # Get people who meet at least the prescribed subset criteria 
        if len(crit_subset_copy) > 0:  # May be zero if it only contained one element initially, in that case then we already have all members 
            for criterion in crit_subset_copy: 
                members = members.intersection(criteria_dict[criterion]) 

        # Remove anyone who meets other criteria
        if len(complement) > 0:  # May be zero when the crit_subset contains all criteria, in that case then there is nothing to remove
            for criterion in complement: 
                members = members - members.intersection(criteria_dict[criterion])

        meets_crit_subset[subset_key] = members  # Meets all and only the subset criteria

    return meets_crit_subset

    
def compute_proportions(crit_subset_dict): 
    '''
    Determine which proportion of tests should be delegated to each subset element of the criteria powerset

    Args: 
        eligible_crit_subset  (dict) : Keys are powerset subset names, values are corresponding test-eligible people who meet the criteria subset
    '''
    proportions = dict()
    n_eligible = 0
    for subset_ppl in crit_subset_dict.values(): 
        n_eligible += len(subset_ppl)

    for sname in crit_subset_dict.keys(): 
        if n_eligible != 0: 
            proportions[sname] = len(crit_subset_dict[sname]) / n_eligible
        else: 
            proportions[sname] = 0  # If no one is eligible, then no one should be given a test

    return proportions
    

def compute_powset(valid_criteria): 
    ''' 
    Given a set of criteria, compute and return the powerset.
    
    Args: 
        valid_criteria  (set)   : Set of valid criteria to compute the power set from. 

    Returns: 
        criteria_powset (dict)  : Keys identify subsets, values are all subsets of valid_criteria
    '''
    criteria_powset = dict()
    for i in range(1, len(valid_criteria) + 1):  # We skip the empty set, since this corresponds to someone satisfying no criteria
        subsets_size_i = list(combinations(valid_criteria, i))
        for j in range(len(subsets_size_i)): 
            criteria_powset['subset' + '_' + str(i) + '_' + str(j)] = set(subsets_size_i[j])

    return criteria_powset