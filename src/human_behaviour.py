from .. import utils as znu

class HumanCategories:
    """
    OBJECTIVE: Categories human agents, for example if they have a sw alert. 
    SKELE: Basic categories: 
        - Tested positive, i.e. union of RAT and PCR positive testers. 
        - Have Avian Influenza Specific Symptoms
        - Have SW alert
    TASKS: 
        - Refactor and populate with additional categories. 
        - complete testers() function
    """
    def __init__(self, sim, criteria=None):
        if criteria is None:
            self.criteria = ['pos_test', 
                             'ai_sx', 
                             'sw_alert']
        else:
            # May want to only specify a few relevant criteria. 
            self.criteria = criteria

        self.cat_dict = self.check_criteria(sim)

    def check_criteria(self, sim): 
        '''
        Given the criterion names, call respective criterion functions to identify who meets each criterion. 

        Args: 
            sim     (simulation)    : Initialized simulation object
            surv    (surveillance)  : The PCR surveillance system being used 
        '''
        meets_criteria = dict()  # Will be a dictionary of sets

        for name in self.criteria:  # Not sure if there's a cleaner way to do this.
            if name == 'pos_test': 
                meets_criteria['pos_test'] = self.testers(sim, mode='pos')
            if name == 'ai_sx': 
                meets_criteria['ai_sx'] = self.sx(sim, mode='specific')
            if name == 'sw_alert': 
                meets_criteria['sw_alert'] = self.sw_alerts(sim, mode='pos')
        
        return meets_criteria

    def testers(self, sim, mode='pos'):
        '''
        Get people who have tested, + we can apply different filters. 

        Args: 
            sim          (sim) : Simulation object 
            mode         (str) : Choose return option 

        Returns: 
            sw_alerted  (set) : Set of people who have tested, with the following filters: 
                                    If mode == 'pos', returns set of all people who received a positive result on that day. 
                                    If mode == 'neg', returns set of all people who received a negative result on that day.
                                    If mode == 'RAT_pos' returns set of all people who received a RAT positive.
                                    If mode == 'RAT_neg' returns set of all people who received a RAT negative.
                                    If mode == 'PCR_pos' returns set of all people who received a PCR positive.
                                    If mode == 'PCR_neg' returns set of all people who received a PCR negative.
        '''
        return_options = ['pos'] # all others are to be implemented.

        if not(mode in return_options): 
            raise RuntimeError('Return choice is invalid, the options are: "all", "specific, "nonspecific"')

        if mode == 'pos':
            return sim.agents.human.date_diagnosed == sim.t

    def sw_alerts(self, sim, mode='pos'):
        '''
        Get people who have positive smartwatch alerts. 

        Args: 
            sim          (sim) : Simulation object 
            mode         (str) : Choose return option 

        Returns: 
            sw_alerted  (set) : Set of people who have smartwatches and ...
                                    If mode == 'pos', returns set of people who have positive alerts on the day. 
                                    If mode == 'neg', returns set of watch users without alerts on the day. 
        '''
        return_options = ['pos'] # neg is to be implemented. 

        if not(mode in return_options): 
            raise RuntimeError('Return choice is invalid, the options are: "all", "specific, "nonspecific"')

        alerted = sim.agents.human.true('sw_alarmed') # Should this be 'alerted'?
        if mode == 'pos':
            return alerted

    def sx(self, sim, mode='all'): 
        '''
        Get people who are currently displaying symptoms. These people either have COVID specific or nonspecific symptoms.

        Args: 
            sim          (sim) : Simulation object 
            mode         (str) : Choose return option 

        Returns: 
            symptomatic  (set) : Set of people who are symptomatic with COVID or ILI. 
                                    If mode == 'all', returns set of people who have any symptoms
                                    If mode == 'specific', returns set of people who have Avian Influenza specific symptoms
                                    If mode == 'nonspecific', returns set of people who have non COVID specific symptoms
        '''
        return_options = ['all', 'specific', 'nonspecific']

        if not(mode in return_options): 
            raise RuntimeError('Return choice is invalid, the options are: "all", "specific, "nonspecific"')

        sx_AI = set(znu.true(sim.agents.human.symptomatic))
        
        if sim.pars['enable_testObj']: # 
            sx_ILI = set(znu.true(sim.agents.human.symptomatic_ILI))
        else: sx_ILI = set()
        all_sx = sx_AI.union(sx_ILI)
        cs_sx, ncs_sx = self.stratify_sx(sim, all_sx) #stratify_sx does not exist

        if mode == 'all': 
            return all_sx

        elif mode == 'specific': 
            return cs_sx

        elif mode == 'nonspecific': 
            return ncs_sx
        

def send_to_quarantine(sim, inds, quar_period):
    """
    Full name: State transition - quarantine
    Args: 
        length of quarantine
    """

    sim.agents.human.schedule_quarantine(inds, start_date=sim.t, period=quar_period)

#### Code some mock-up functions to simplify development. ####
def random_tests(sim):
    # Randomly select 5% of people to seek tests; assume instant results. 
    test_inds = znu.choose(sim.n, 0.05)

    # Filter them by who's exposed. 100% SP. 
    test_inds = test_inds[sim.agents.human.exposed[test_inds]]

    # Set their test date to today. 100% SE. Assume 0 day delay. 
    sim.agents.human.date_pos_test[test_inds] = sim.t
    sim.agents.human.date_diagnosed[test_inds] = sim.t

class BehaviourUpdater:
    def __init__(self, **kwargs):
        self.make_pars()
        self.update_pars(**kwargs)
    
    def make_pars(self):
        self.quarantine_probs = {
            'ai_sx': 0.3,
            'pos_test': 0.8,
            'sw_alert': 0.05
        }
        self.quar_period = 14 # This should be sim.pars['quar_period], or at least equal to it.

    def update_behaviour(self, sim):
        """
        Simple model: each data point gives people a certain probability of quarantine. 
        """
        ppl_cats = HumanCategories(sim, criteria=self.quarantine_probs.keys())
        all_quar_inds = set()

        for cat in ppl_cats.cat_dict:
            # Draw quarantines from the category. Coin flip. 
            quar_draws = znu.n_binomial(self.quarantine_probs[cat], len(ppl_cats.cat_dict[cat]))

            # Get the indicies corresponding to quar_draws being true.
            quar_inds = ppl_cats.cat_dict[cat][quar_draws.nonzero()[0]]
            all_quar_inds = all_quar_inds.update(quar_inds)

        # Schedule quarantines.
        send_to_quarantine(sim, all_quar_inds, self.l_Q)

    def update_pars(self, **kwargs):
        """
        Args: 
            **kwargs: keyword arguments to update parameters.
        Body:
            Update fields with parameters. New fields cannot be initialized in this way. 
        """
        for k, v in kwargs.items():
            # Check if k is already a field. 
            if hasattr(self, k):
                setattr(self, k, v)
            else:
                raise ValueError('Invalid parameter: %s' % k)