import numpy as np 

def LOD_test(testobj, sim, inds): 
    '''
    Administer tests to a list of users.
    Tests fail with some probability p_fail, and return the opposite truth value.
    Assume people have Avian Influenza only if they are exposed (i.e. any time between infection to recovery).
    Conduct analysis if desired. Record whether they were symptomatic at the time of testing.

    Args:
        testobj        (testobj): Test object
        sim                (sim): Simulation object
        inds          (np.array): Human agent indices to test

    Returns: 
        test_results  (np.array): Boolean array of test results
    '''
    LOD = testobj.LOD
    test_results = sim.agents.human.viral_load[inds] >= LOD
    failed = np.random.binomial(1, testobj.p_fail, len(test_results))  # draw Bernoulli trials to determine which test results failed
    test_results[failed.nonzero()[0]] = ~test_results[failed.nonzero()[0]]  # for the failed results, return opposite truth value

    return test_results

def basic_test(testobj, sim, inds): 
    '''
    Administer tests to a list of users. Assume people have COVID-19 only if they are exposed (i.e. any time between infection to recovery). 
    Conduct analysis if desired. Record whether they were symptomatic at the time of testing.

    Args: 
        sim                (sim): Simulation object
        inds          (np.array): Human agent indices to test
        se               (float): Test sensitivity 
        sp               (float): Test specificity

    Returns: 
        test_results  (np.array): Boolean array of test results
    '''
    se = testobj.sensitivity
    sp = testobj.specificity

    test_results = np.full(len(inds), False)
    test_results[sim.agents.human.exposed[inds]] = np.random.uniform(0, 1, sum(sim.agents.human.exposed[inds])) < se
    test_results[~sim.agents.human.exposed[inds]] = np.random.uniform(0, 1, sum(~sim.agents.human.exposed[inds])) > sp

    return test_results