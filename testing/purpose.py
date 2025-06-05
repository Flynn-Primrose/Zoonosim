'''
Defines classes related to the purpose of testing (i.e. Diagnostic, Screening, and Surveillance)
'''




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

    def initialize(self, pop_type, people): 
        if not(people is None) and (pop_type == 'behaviour_module'):
            self.workplaces = people.workplaces
            self.n_workplaces = people.n_workplaces
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