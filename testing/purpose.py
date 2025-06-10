'''
Defines classes related to the purpose of testing (i.e. Diagnostic, Screening, and Surveillance)
'''

from . import eligible as eg
from . import seek as sk

import numpy as np

class Diagnostic:
    # Used to conduct diagnostic testing
    # Takes as input:
    # - Type of test to be used (PCR or RAT)    (mandatory)
    #     - Test-specific details: SE/SP, 
    # - Eligibility criteria                    (optional)
    # - Test seeking probabilities              (optional)

    def __init__(self, criteria, seekprobs): 
        self.system = 'Diagnostic'
        self.criteria = criteria
        self.seekprobs = seekprobs

    def initialize(self): 
        return

    def apply(self, testobj, sim): 
        # Apply test criteria
        if not(testobj.crit_kwargs is None):
            all_eligible = eg.check_criteria(self, testobj, sim, **testobj.crit_kwargs)
        else:
            all_eligible = eg.check_criteria(self, testobj, sim)

        # Apply test seeking probabilities
        all_seekers = sk.get_all_seekers(all_eligible, self.seekprobs)

        return all_eligible, all_seekers


class Screening: 

    def __init__(self, criteria, seekprobs): 
        self.system = 'Screening'
        self.criteria = criteria
        self.seekprobs = seekprobs

    def initialize(self, agents): 
        if not(agents is None):
            fids = np.unique(agents.fid)
            fid2uid = dict()
            for fid in fids: 
                fid2uid[fid] = np.where((agents.fid == fid) & (agents.type == 'human'))[0]
            self.fid2uid = fid2uid
            self.n_fids = len(self.fid2uid.keys())
        return

    def apply(self, testobj, sim): 
        # Apply test criteria
        if not(testobj.crit_kwargs is None):
            all_eligible = eg.check_criteria(self, testobj, sim, **testobj.crit_kwargs)
        else:
            all_eligible = eg.check_criteria(self, testobj, sim)

        # Apply test seeking probabilities
        all_seekers = sk.get_all_seekers(all_eligible, self.seekprobs)

        return all_eligible, all_seekers


class Surveillance: 

    def __init__(self): 
        self.system = 'Surveillance'
        return
    
    def initialize(self): 
        return 

    def apply(self, testobj, sim): 
        return