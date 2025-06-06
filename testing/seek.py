import numpy as np

'''
This file implements functions for checking test-seeking.
'''

def get_seekers(candidates, prob): 
        ''' 
        Takes as input a set of potential test-seekers. Determines if they actually will seek a test, probabilistically. 

        Args: 
            candidates  (set): Set of people who satisfy a test-seeking criterion

        Returns: 
            seekers     (set): Set of people who will seek a test
        '''
        if not(prob is None): 
            test_seeker_uids = np.fromiter(candidates, dtype=int)[np.random.binomial(1, prob, len(candidates)).nonzero()[0]] 
            return set(test_seeker_uids.tolist())

        return set()

def get_all_seekers(eligible_dict, seekprobs):
        '''
        Look at who satisfies the test criteria, and determine which of these people seek a test

        Args: 
            seekprobs        (dict)          : Dictionary containing test seeking probability associated with each criterion
            eligible_dict    (dict)          : Dictionary containing sets of people who meet each criterion

        Returns: 
            test_seekers  (dict)    : Dictionary of test seekers for each criterion
        '''
        test_seekers = dict()

        for criterion in eligible_dict.keys():
            test_seekers[criterion] = get_seekers(eligible_dict[criterion], seekprobs[criterion])

        return test_seekers