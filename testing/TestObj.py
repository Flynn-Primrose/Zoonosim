import pickle as pkl
import numpy as np
import pandas as pd 
from scipy import stats
from . import purpose as tp
from . import allocate as ac
from . import mechanism as tm

__all__ = ['TestObj', 'PCR_disc', 'RAT_disc']

class TestObj: 
    def __init__(self): 
        self.system = self.__class__.__name__
        self.initialized = False
        self.finalized = False

        return
    
    def __call__(self, *args, **kwargs):
        if not self.initialized:  # pragma: no cover
            errormsg = f'TestObj (label={self.system}, {type(self)}) has not been initialized'
            raise RuntimeError(errormsg)
        return self.apply(*args, **kwargs)
    
    def init_common_trackers(self, sim):
        '''
        Initialize common trackers between test objects (PCR and RAT).
        '''
        # Initialize count trackers
        self.pos_history             = np.zeros(sim.pars['n_days'] + 1)  # The simulation always goes on for one extra day, hence +1
        self.neg_history             = np.zeros(sim.pars['n_days'] + 1)
        self.true_pos                = np.zeros(sim.pars['n_days'] + 1)
        self.false_pos               = np.zeros(sim.pars['n_days'] + 1)
        self.true_neg                = np.zeros(sim.pars['n_days'] + 1)
        self.false_neg               = np.zeros(sim.pars['n_days'] + 1)
        self.test_count_split        = {k : np.zeros((sim.pars['n_days'] + 1)) for k in ['I_C', 'I_nC', 'nI_C', 'nI_nC']}

        # Distributions for time from exposed to receipt of true positive result
        self.rel_date_TP             = np.zeros(sim.pars['n_days'] + 1)  # time since exposure (true positives only)
        self.rel_date_TP_I_C         = np.zeros(sim.pars['n_days'] + 1)
        self.rel_date_TP_I_nC        = np.zeros(sim.pars['n_days'] + 1)
        self.rel_date_TP_nI_C        = np.zeros(sim.pars['n_days'] + 1)
        self.rel_date_TP_nI_nC       = np.zeros(sim.pars['n_days'] + 1)

        # Distributions for time from exposed to receipt of false negative result
        self.rel_date_FN             = np.zeros(sim.pars['n_days'] + 1)  # time since exposure (false negatives only)
        self.rel_date_FN_I_C         = np.zeros(sim.pars['n_days'] + 1)
        self.rel_date_FN_I_nC        = np.zeros(sim.pars['n_days'] + 1)
        self.rel_date_FN_nI_C        = np.zeros(sim.pars['n_days'] + 1)
        self.rel_date_FN_nI_nC       = np.zeros(sim.pars['n_days'] + 1)

        # Distributions for time from exposed to administration of test
        self.rel_date_test           = np.zeros(sim.pars['n_days'] + 1)  # time since exposure (true positive and false negative only)
        self.rel_date_test_I_C       = np.zeros(sim.pars['n_days'] + 1)
        self.rel_date_test_I_nC      = np.zeros(sim.pars['n_days'] + 1)
        self.rel_date_test_nI_C      = np.zeros(sim.pars['n_days'] + 1)
        self.rel_date_test_nI_nC     = np.zeros(sim.pars['n_days'] + 1)

        # Distributions for eligible waiting time (time from eligible to administration of test)
        self.rel_date_egwait         = np.zeros(sim.pars['n_days'] + 1)
        self.rel_date_egwait_I_C     = np.zeros(sim.pars['n_days'] + 1)
        self.rel_date_egwait_I_nC    = np.zeros(sim.pars['n_days'] + 1)
        self.rel_date_egwait_nI_C    = np.zeros(sim.pars['n_days'] + 1)
        self.rel_date_egwait_nI_nC   = np.zeros(sim.pars['n_days'] + 1)

        # Distributions for seek waiting time (time from eligible to administration of test)
        self.rel_date_skwait         = np.zeros(sim.pars['n_days'] + 1)
        self.rel_date_skwait_I_C     = np.zeros(sim.pars['n_days'] + 1)
        self.rel_date_skwait_I_nC    = np.zeros(sim.pars['n_days'] + 1)
        self.rel_date_skwait_nI_C    = np.zeros(sim.pars['n_days'] + 1)
        self.rel_date_skwait_nI_nC   = np.zeros(sim.pars['n_days'] + 1)

        # Distributions for consecutive days test-eligible before end of eligibility
        self.dist_egcons             = np.zeros(sim.pars['n_days'] + 1)
        self.dist_egcons_I_C         = np.zeros(sim.pars['n_days'] + 1)
        self.dist_egcons_I_nC        = np.zeros(sim.pars['n_days'] + 1)
        self.dist_egcons_nI_C        = np.zeros(sim.pars['n_days'] + 1)
        self.dist_egcons_nI_nC       = np.zeros(sim.pars['n_days'] + 1)

        # Distributions for consecutive days test-seeking before end of eligibility
        self.dist_skcons             = np.zeros(sim.pars['n_days'] + 1)
        self.dist_skcons_I_C         = np.zeros(sim.pars['n_days'] + 1)
        self.dist_skcons_I_nC        = np.zeros(sim.pars['n_days'] + 1)
        self.dist_skcons_nI_C        = np.zeros(sim.pars['n_days'] + 1)
        self.dist_skcons_nI_nC       = np.zeros(sim.pars['n_days'] + 1)

        # Initialize date trackers -- important for integrating with Covasim quarantine behaviour
        self.date_positive           = np.full(sim.pars['pop_size_by_type']['human'], -1)
        self.date_negative           = np.full(sim.pars['pop_size_by_type']['human'], -1)
        self.date_pos_test           = np.full(sim.pars['pop_size_by_type']['human'], -1)

        # Initialize overflow trackers
        self.pos_overflow            = 0
        self.neg_overflow            = 0

        # Initialize other trackers
        self.retest_tracker        = np.zeros(sim.pars['pop_size_by_type']['human'])     # Track whether an individual is allowed to retest
        self.cons_days_eg          = np.zeros(sim.pars['pop_size_by_type']['human'])     # Track the number of consecutive days individual is eligible for a test
        self.cons_days_sk          = np.zeros(sim.pars['pop_size_by_type']['human'])     # Track the number of consecutive days individual is seeking a test
        self.eventual_test         = np.zeros(sim.pars['n_days'] + 1)   # Track the number of individuals infected on a given day who eventually receive a test DURING the associated infection

    def update_pars(self, pars_dict): 
        for param_name in pars_dict.keys():
            if hasattr(self, param_name): 
                setattr(self, pars_dict[param_name])
            else: 
                raise RuntimeError(f'Tried to assign {param_name} to testobj, but testobj does not have that attribute')
            
    def save_pars(self, save_path):
        result = {k:v for k, v in self.__dict__.items() if not(k.startswith('__'))}
        with open(save_path, 'wb') as f:
            pkl.dump(result, f)

    def finalize(self):
        '''
        Finalize the TestObj object
        '''
        if self.finalized: # pragma: no cover
            raise RuntimeError('TestObj already finalized')  # Raise an error because finalizing multiple times has a high probability of producing incorrect results e.g. applying rescale factors twice
        self.finalized = True
        return


class PCR_disc(TestObj): 

    def __init__(self, 
                 di_criteria=None,
                 di_seekprobs=None,
                 sc_criteria=None,
                 sc_seekprobs=None,
                 crit_kwargs=None,
                 sensitivity=0.85,
                 specificity=1.00,
                 LOD=3,
                 p_fail=0.01,
                 capacity=0.015,
                 pos_retest_delay=90,  
                 neg_retest_delay=14,  
                 mean_latency=3,
                 static_latency=None,
                 RAT_ind=None,
                 timeline=None,
                 test_mode='basic',
                 name='PCR'):
        '''
        PCR-based 'di'agnostic testing and 'sc'reening.

        Args: 
            system             (str)    : Test system name
            di_criteria        (list)   : List of eligibility criteria used for diagnostic testing
            di_seekprobs       (dict)   : Dictionary of test seeking probabilities (value) associated with each criteria (key)
            sc_criteria        (list)   : List of eligibiltiy criteria used for screening 
            sc_seekprobs       (dict)   : Dictionary of test seeking probabilities (value) associated with each criteria (key) 
            sensitivity        (float)  : Test sensitivity
            specificity        (float)  : Test specificity 
            capacity           (n/s)    : Test capacity as a proportion of total population (float or array)
            pos_retest_delay   (int)    : Number of days before an individual can retest, after receiving a positive test
            neg_retest_delay   (int)    : Number of days before an individual can retest, after receiving a negative test
            latency_shape      (int)    : Shape of gamma distribution for test reporting latency 
            latency_scale      (int)    : Scale of gamma distribution for test reporting latency 
            static_latency     (int)    : Specificy a static test reporting latency (optional)
            RAT_ind            (int)    : Index of RAT testing system, which also inherits from TestObj. Usually RAT_disc, can also be RAT_surv. 
            hybrid             (bool)   : Whether the testing system is a PCR-RAT hybrid 
            timeline           (dict)   : Parameter updates at different times 
        '''
        super().__init__()
        self.di_criteria = di_criteria
        self.di_seekprobs = di_seekprobs
        self.sc_criteria = sc_criteria
        self.sc_seekprobs = sc_seekprobs
        self.crit_kwargs = crit_kwargs
        self.sensitivity = sensitivity
        self.specificity = specificity
        self.LOD = LOD
        self.p_fail = p_fail
        self.capacity = capacity 
        self.pos_retest_delay = pos_retest_delay 
        self.neg_retest_delay = neg_retest_delay 
        self.mean_latency = mean_latency
        self.latency_shape = 2  # fixed
        self.static_latency = static_latency
        self.RAT_ind = RAT_ind  # For ease of use, we could remove this and hide it from user, force RAT_ind=0 always. 
        self.RAT = None
        self.hybrid = False
        self.timeline = timeline
        self.test_mode = test_mode
        self.name = name

    def initialize(self, sim): 
        self.latency_scale = self.mean_latency / self.latency_shape

        self.diagnostic = tp.Diagnostic(self.di_criteria, self.di_seekprobs)
        self.screening = tp.Screening(self.sc_criteria, self.sc_seekprobs)

        # Initialize base test objects
        self.diagnostic.initialize()
        self.screening.initialize(sim.agents)

        # Initialize two-way interface with another testing system - assumed to be RAT currently, can it be PCR? 
        if not(self.RAT_ind is None): 
            self.RAT = sim.pars['testing'][self.RAT_ind] 
            self.RAT.PCR = self

            # Activate
            self.RAT.hybrid = True
            self.hybrid = True

        self.initialized = True
        

    def apply(self, sim): 
        '''
        Every day: 
            1) Get people eligible for diagnostic and screening
            2) Allocate tests to test-seekers
            3) Conduct testing
        '''
        if not(self.timeline is None): 
            update_parameters(self, sim.t)

            # Refresh diagnostic and screening objects with any updates to criteria/probabilities
            self.diagnostic.criteria = self.di_criteria 
            self.diagnostic.seekprobs = self.di_seekprobs
            self.screening.criteria = self.sc_criteria
            self.screening.seekprobs = self.sc_seekprobs

        # --------------------------
        # Update retest restrictions
        # --------------------------
        self.retest_tracker -= 1

        # -----------------------------------------------
        # Determine test-eligible and test-seeking agents
        # -----------------------------------------------
        di_eg, di_sk = self.diagnostic.apply(self, sim)
        sc_eg, sc_sk = self.screening.apply(self, sim)
        
        # Merge dictionaries
        all_eg = {**di_eg, **sc_eg}
        all_sk = {**di_sk, **sc_sk}

        # --------------------------
        # Allocate and conduct tests
        # --------------------------
        test_uids = ac.allocate_tests(sim, self, all_sk, mode='proportional')
        if self.test_mode == 'basic':
            test_results = tm.basic_test(self, sim, test_uids)
        elif self.test_mode == 'LOD':
            test_results = tm.LOD_test(self, sim, test_uids)
        else: 
            raise RuntimeError('Invalid test mode supplied. Options are "basic" and "LOD"')
        
        #TODO: Why are there still so many people receiving several diagnoses on the same day?
        # --------------------------
        # Update retest restrictions
        # --------------------------
        self.retest_tracker[test_uids[test_results]] = self.pos_retest_delay
        self.retest_tracker[test_uids[~test_results]] = self.neg_retest_delay

        # ------------------
        # Record information
        # ------------------
        # TODO: Optimize record_pipeline
        record_pipeline(self, sim, all_eg, all_sk, test_uids)
        record_test_results(self, sim, test_uids, test_results)  # Latencies get recorded here
        record_consumed(self, sim, len(test_uids))

        # -----------------------------------
        # Deal with consecutive days eligible
        # -----------------------------------
        eg = list(set().union(*all_eg.values()))
        eg_ind = np.zeros(sim.pars['pop_size_by_type']['human'], dtype=bool)
        eg_ind[eg] = True  # corresponds to which agents are eligible

        # Track distribution of days test-eligible before end of eligibility
        eg_end = np.logical_and.reduce((~eg_ind, self.cons_days_eg>=1))  # eligibility has ended if not eligible today + was eligible yesterday
        eg_end_uids = eg_end.nonzero()[0]
        record_cons_days_eg_dist(self, sim, eg_end_uids)  # record BEFORE tidying up

        # ----------------------------------
        # Deal with consecutive days seeking
        # ----------------------------------
        sk = list(set().union(*all_sk.values()))
        sk_ind = np.zeros(sim.pars['pop_size_by_type']['human'], dtype=bool)
        sk_ind[sk] = True  # corresponds to which agents are seeking

        # Track distribution of days test-seeking before end of eligibility
        sk_end = np.logical_and.reduce((~sk_ind, self.cons_days_sk>=1, eg_end))  # test-seeking has ended if not seeking today + was seeking yesterday, then filter those whose eligibility ended
        sk_end_uids = sk_end.nonzero()[0]
        record_cons_days_sk_dist(self, sim, sk_end_uids)

        # -------
        # Tidy up
        # -------
        self.cons_days_eg[eg_ind] += 1
        self.cons_days_eg[~eg_ind] = 0  # reset everyone else
        self.cons_days_eg[test_uids] = 0

        self.cons_days_sk[sk_ind] += 1
        self.cons_days_sk[~sk_ind] = 0  # reset everyone else
        self.cons_days_sk[test_uids] = 0  # we want every period of test-seeking to only be associated with receipt of one test

        sim.results['new_PCR_tests'][sim.t] += len(test_uids)
        sim.results['cum_PCR_tests'][sim.t] += sum(sim.results['new_PCR_tests'][:sim.t])

        return

    def finalize(self, sim):
        '''
        Do any finalizing work here (e.g. resetting to save memory)
        '''
        super().finalize()
        return


class RAT_disc(TestObj): 

    def __init__(self, 
                di_criteria=None, 
                di_seekprobs=None,
                sc_criteria=None, 
                sc_seekprobs=None,
                crit_kwargs=None,
                sensitivity=0.70,
                specificity=0.99,
                LOD=6,  # Log scale
                p_fail=0.10,
                capacity=0.05,
                pos_retest_delay=14,
                neg_retest_delay=0,  
                mean_latency=0,
                static_latency=None,
                timeline=None,
                test_mode='basic',
                name='RAT'):
        '''
        RAT-based 'di'agnostic testing and 'sc'reening.
        '''
        super().__init__()
        self.di_criteria = di_criteria
        self.di_seekprobs = di_seekprobs
        self.sc_criteria = sc_criteria
        self.sc_seekprobs = sc_seekprobs
        self.crit_kwargs = crit_kwargs
        self.sensitivity = sensitivity
        self.specificity = specificity
        self.LOD = LOD
        self.p_fail = p_fail
        self.capacity = capacity
        self.pos_retest_delay = pos_retest_delay
        self.neg_retest_delay = neg_retest_delay
        self.mean_latency = mean_latency
        self.latency_shape = 2
        self.static_latency = static_latency
        self.PCR = None
        self.RAT = None
        self.hybrid = False
        self.timeline = timeline
        self.test_mode = test_mode
        self.name = name

    def initialize(self, sim): 
        self.latency_scale = self.mean_latency / self.latency_shape

        self.diagnostic = tp.Diagnostic(self.di_criteria, self.di_seekprobs)
        self.screening = tp.Screening(self.sc_criteria, self.sc_seekprobs)

        # Initialize base test objects
        self.diagnostic.initialize()
        self.screening.initialize(sim.agents)

        # Initialize trackers specific to RATs
        self.cons_days_neg_rat = np.zeros(sim.pars['pop_size_by_type']['human'])  # Track the number of consecutive days individual has received negative RAT

        self.initialized = True

    def apply(self, sim): 
        '''
        Every day: 
            1) Get people eligible for diagnostic and screening
            2) Allocate tests to test-seekers
            3) Conduct testing
        '''
        if not(self.timeline is None): 
            update_parameters(self, sim.t)

            # Refresh diagnostic and screening objects with any updates to criteria/probabilities
            self.diagnostic.criteria = self.di_criteria
            self.diagnostic.seekprobs = self.di_seekprobs
            self.screening.criteria = self.sc_criteria
            self.screening.seekprobs = self.sc_seekprobs

        # --------------------------
        # Update retest restrictions 
        # --------------------------
        self.retest_tracker -= 1

        # -----------------------------------------------
        # Determine test-eligible and test-seeking agents
        # -----------------------------------------------
        di_eg, di_sk = self.diagnostic.apply(self, sim)
        sc_eg, sc_sk = self.screening.apply(self, sim)
        
        # Merge dictionaries
        all_eg = {**di_eg, **sc_eg}
        all_sk = {**di_sk, **sc_sk}
        
        # --------------------------
        # Allocate and conduct tests
        # --------------------------
        test_uids = ac.allocate_tests(sim, self, all_sk, mode='proportional')
        if self.test_mode == 'basic':
            test_results = tm.basic_test(self, sim, test_uids)
        elif self.test_mode == 'LOD':
            test_results = tm.LOD_test(self, sim, test_uids)
        else: 
            raise RuntimeError('Invalid test mode supplied. Options are "basic" and "LOD"')

        #TODO: Why are there still so many people receiving several diagnoses on the same day?
        # --------------------------
        # Update retest restrictions
        # --------------------------
        self.retest_tracker[test_uids[test_results]] = self.pos_retest_delay
        self.retest_tracker[test_uids[~test_results]] = self.neg_retest_delay

        # ------------------
        # Record information
        # ------------------
        record_pipeline(self, sim, all_eg, all_sk, test_uids)
        record_test_results(self, sim, test_uids, test_results)
        record_consumed(self, sim, len(test_uids))

        # -------------------------------
        # Track consecutive days negative
        # -------------------------------
        # Assume agents instantly know test result
        tested_negative = np.zeros(sim.pars['pop_size_by_type']['human'], dtype=bool) 
        tested_negative[test_uids[~test_results]] = True
        self.cons_days_neg_rat[tested_negative] += 1
        self.cons_days_neg_rat[~tested_negative] = 0  # reset everyone else

        # -----------------------------------
        # Deal with consecutive days eligible
        # -----------------------------------
        eg = list(set().union(*all_eg.values()))
        eg_ind = np.zeros(sim.pars['pop_size_by_type']['human'], dtype=bool)
        eg_ind[eg] = True  # corresponds to which agents are eligible

        # Track distribution of days test-eligible before end of eligibility
        eg_end = np.logical_and.reduce((~eg_ind, self.cons_days_eg>=1))  # eligibility has ended if not eligible today + was eligible yesterday
        eg_end_uids = eg_end.nonzero()[0]
        record_cons_days_eg_dist(self, sim, eg_end_uids)  # record BEFORE tidying up

        # ----------------------------------
        # Deal with consecutive days seeking
        # ----------------------------------
        sk = list(set().union(*all_sk.values()))
        sk_ind = np.zeros(sim.pars['pop_size_by_type']['human'], dtype=bool)
        sk_ind[sk] = True  # corresponds to which agents are seeking

        # Track distribution of days test-seeking before end of eligibility
        sk_end = np.logical_and.reduce((~sk_ind, self.cons_days_sk>=1, eg_end))  # test-seeking has ended if not seeking today + was seeking yesterday, then filter those whose eligibility ended
        sk_end_uids = sk_end.nonzero()[0]
        record_cons_days_sk_dist(self, sim, sk_end_uids)
        
        # -------
        # Tidy up
        # -------
        # Must only tidy up after record functions, since for example the waiting time from initial eligibility to test is allowed to be 0
        self.cons_days_eg[eg_ind] += 1
        self.cons_days_eg[~eg_ind] = 0  # reset everyone else
        self.cons_days_eg[test_uids] = 0

        self.cons_days_sk[sk_ind] += 1
        self.cons_days_sk[~sk_ind] = 0  # reset everyone else
        self.cons_days_sk[test_uids] = 0  # we want every period of test-seeking to only be associated with receipt of one test

        sim.results['new_RAT_tests'][sim.t] += len(test_uids)
        sim.results['cum_RAT_tests'][sim.t] += sum(sim.results['new_RAT_tests'][:sim.t])

    def finalize(self, sim):
        '''
        Do any finalizing work here (e.g. resetting to save memory)
        '''
        super().finalize()
        return


# ---------
# Utilities
# ---------
def find_seed_infections(sim):
    seed_infections = set()
    for entry in sim.people.infection_log:
        if entry['layer'] == 'seed_infection':
            seed_infections.add(entry['target'])
    return np.array(list(seed_infections))

def update_dist(testobj, times, attr):
    '''
    Given an array of times, find the unique values and their associated counts.
    Update the associated distribution.
    '''
    unique, counts = np.unique(times, return_counts=True)
    if not(len(unique) == 0):
        getattr(testobj, attr)[unique.astype(int)] += counts

def time_since_exposure(sim, testobj, end, uids, attr):
    '''
    Given agents who test (uids), retrieve the time elapsed between exposure and some end condition (e.g., administration of test).
    Update the associated "rel_date" distribution.
    '''
    rel_date = end - sim.agents.human.date_exposed[uids]
    update_dist(testobj, rel_date, attr)

def time_waiting(testobj, cons_days, uids, attr):
    '''
    Given agents who test (uids), retrieve the number of consecutive days (cons_days) they have been eligible, or have sought a test.
    Update the associated waiting time distribution.
    '''
    wait_time = cons_days[uids]
    update_dist(testobj, wait_time, attr)

def record_test_results(testobj, sim, test_uids, test_results): 
    cur = sim.t
    
    # Record the positive results with a result latency term (swab delay is not included in this term, as we are not delaying test administration)
    if not(testobj.static_latency is None): 
        test_latencies = np.full(len(test_results), testobj.static_latency)
        unique_latencies = [testobj.static_latency]
    else:
        test_latencies = np.trunc(stats.gamma.rvs(a=testobj.latency_shape, scale=testobj.latency_scale, size=len(test_results)))  # latencies are always sampled as zero if scale is zero
        unique_latencies = np.unique(test_latencies).astype(np.intc)

    # Record date of tests that will eventually return positive (i.e., the date before latency is introduced)
    testobj.date_pos_test[test_uids[test_results]] = cur

    # Stratify population into ILI and COVID, ILI only, COVID only, healthy
    I_C_uids   = test_uids[np.logical_and( sim.agents.human.symptomatic_ILI[test_uids],  sim.agents.human.exposed[test_uids])]
    I_nC_uids  = test_uids[np.logical_and( sim.agents.human.symptomatic_ILI[test_uids], ~sim.agents.human.exposed[test_uids])]
    nI_C_uids  = test_uids[np.logical_and(~sim.agents.human.symptomatic_ILI[test_uids],  sim.agents.human.exposed[test_uids])]
    nI_nC_uids = test_uids[np.logical_and(~sim.agents.human.symptomatic_ILI[test_uids], ~sim.agents.human.exposed[test_uids])]

    # Track stratification of test count by those having COVID/ILI at time of administration
    testobj.test_count_split['I_C'][cur]   += len(I_C_uids)
    testobj.test_count_split['I_nC'][cur]  += len(I_nC_uids)
    testobj.test_count_split['nI_C'][cur]  += len(nI_C_uids)
    testobj.test_count_split['nI_nC'][cur] += len(nI_nC_uids)

    # Record number of people infected on a given day who eventually got tested
    exposure_dates = sim.agents.human.date_exposed[test_uids[sim.agents.human.exposed[test_uids]]]  # get exposure dates for those those that tested and are exposed
    exposure_dates = exposure_dates[~np.isnan(exposure_dates)]  # remove NaNs (TODO: why are there NaNs? If they are exposed, they should have an exposure date)
    update_dist(testobj, exposure_dates, 'eventual_test')
    
    # unique_dates, date_counts = np.unique(exposure_dates, return_counts=True)  # replaced with 'update_dist' function
    # if len(unique_dates > 0):
    #     testobj.eventual_test[unique_dates.astype(int)] += date_counts

    # Update waiting time distributions
    time_waiting(testobj, testobj.cons_days_eg, test_uids,  'rel_date_egwait')
    time_waiting(testobj, testobj.cons_days_eg, I_C_uids,   'rel_date_egwait_I_C')
    time_waiting(testobj, testobj.cons_days_eg, I_nC_uids,  'rel_date_egwait_I_nC')
    time_waiting(testobj, testobj.cons_days_eg, nI_C_uids,  'rel_date_egwait_nI_C')
    time_waiting(testobj, testobj.cons_days_eg, nI_nC_uids, 'rel_date_egwait_nI_nC')

    time_waiting(testobj, testobj.cons_days_sk, test_uids,  'rel_date_skwait')
    time_waiting(testobj, testobj.cons_days_sk, I_C_uids,   'rel_date_skwait_I_C')
    time_waiting(testobj, testobj.cons_days_sk, I_nC_uids,  'rel_date_skwait_I_nC')
    time_waiting(testobj, testobj.cons_days_sk, nI_C_uids,  'rel_date_skwait_nI_C')
    time_waiting(testobj, testobj.cons_days_sk, nI_nC_uids, 'rel_date_skwait_nI_nC')

    # Record dates that involve a latency - these dates are when test results become known to the individual
    for l in unique_latencies:  # This should be a short loop
        if cur + l <= sim.pars['n_days']:

            # Track case counts
            testobj.pos_history[cur + l] += np.sum(test_results[test_latencies == l])
            testobj.neg_history[cur + l] += np.sum(~test_results[test_latencies == l])  # Evaluates ~test_results before boolean indexing

            # Track dates
            testobj.date_positive[test_uids[np.logical_and(test_results, test_latencies == l)]] = cur + l  # Pick people who returned positive result and was assigned latency l 
            testobj.date_negative[test_uids[np.logical_and(~test_results, test_latencies == l)]] = cur + l  # Pick people who returned negative result and was assigned latency l

            # Track boolean states
            TP_uids       = test_uids[np.logical_and.reduce((test_results,  test_latencies==l, sim.agents.human.viral_load[test_uids]>=testobj.LOD))]
            FP_uids       = test_uids[np.logical_and.reduce((test_results,  test_latencies==l, sim.agents.human.viral_load[test_uids]<testobj.LOD))]
            TN_uids       = test_uids[np.logical_and.reduce((~test_results, test_latencies==l, sim.agents.human.viral_load[test_uids]<testobj.LOD))]
            FN_uids       = test_uids[np.logical_and.reduce((~test_results, test_latencies==l, sim.agents.human.viral_load[test_uids]>=testobj.LOD))]

            TP_uids_I_C   = TP_uids[np.logical_and( sim.agents.human.symptomatic_ILI[TP_uids],  sim.agents.human.exposed[TP_uids])]  # has ILI and has COVID
            TP_uids_I_nC  = TP_uids[np.logical_and( sim.agents.human.symptomatic_ILI[TP_uids], ~sim.agents.human.exposed[TP_uids])]  # has ILI and no COVID
            TP_uids_nI_C  = TP_uids[np.logical_and(~sim.agents.human.symptomatic_ILI[TP_uids],  sim.agents.human.exposed[TP_uids])]  # no ILI and has COVID
            TP_uids_nI_nC = TP_uids[np.logical_and(~sim.agents.human.symptomatic_ILI[TP_uids], ~sim.agents.human.exposed[TP_uids])]  # no ILI and no COVID

            FN_uids_I_C   = FN_uids[np.logical_and( sim.agents.human.symptomatic_ILI[FN_uids],  sim.agents.human.exposed[FN_uids])]  
            FN_uids_I_nC  = FN_uids[np.logical_and( sim.agents.human.symptomatic_ILI[FN_uids], ~sim.agents.human.exposed[FN_uids])]  
            FN_uids_nI_C  = FN_uids[np.logical_and(~sim.agents.human.symptomatic_ILI[FN_uids],  sim.agents.human.exposed[FN_uids])]  
            FN_uids_nI_nC = FN_uids[np.logical_and(~sim.agents.human.symptomatic_ILI[FN_uids], ~sim.agents.human.exposed[FN_uids])] 

            testobj.true_pos[cur + l]  += len(TP_uids)
            testobj.false_pos[cur + l] += len(FP_uids)
            testobj.true_neg[cur + l]  += len(TN_uids)
            testobj.false_neg[cur + l] += len(FN_uids)

            # Track time since exposure to ADMINISTRATION of test, whether true positive or false negative
            # NOTE: Since we are only ever supplying TP and FN uids, we never get NaN values for date_exposed
            time_since_exposure(sim, testobj, end=cur,   uids=np.union1d(TP_uids,        FN_uids),        attr='rel_date_test')
            time_since_exposure(sim, testobj, end=cur,   uids=np.union1d(TP_uids_I_C,    FN_uids_I_C),    attr='rel_date_test_I_C')
            time_since_exposure(sim, testobj, end=cur,   uids=np.union1d(TP_uids_I_nC,   FN_uids_I_nC),   attr='rel_date_test_I_nC')
            time_since_exposure(sim, testobj, end=cur,   uids=np.union1d(TP_uids_nI_C,   FN_uids_nI_C),   attr='rel_date_test_nI_C')
            time_since_exposure(sim, testobj, end=cur,   uids=np.union1d(TP_uids_nI_nC,  FN_uids_nI_nC),  attr='rel_date_test_nI_nC')

            # Track time since exposure to RECEPTION of true positive result
            time_since_exposure(sim, testobj, end=cur+l, uids=TP_uids,                                    attr='rel_date_TP')
            time_since_exposure(sim, testobj, end=cur+l, uids=TP_uids_I_C,                                attr='rel_date_TP_I_C')
            time_since_exposure(sim, testobj, end=cur+l, uids=TP_uids_I_nC,                               attr='rel_date_TP_I_nC')
            time_since_exposure(sim, testobj, end=cur+l, uids=TP_uids_nI_C,                               attr='rel_date_TP_nI_C')
            time_since_exposure(sim, testobj, end=cur+l, uids=TP_uids_nI_nC,                              attr='rel_date_TP_nI_nC')

            # Track time since exposure to RECEPTION of false negative result
            time_since_exposure(sim, testobj, end=cur+l, uids=FN_uids,                                    attr='rel_date_FN')
            time_since_exposure(sim, testobj, end=cur+l, uids=FN_uids_I_C,                                attr='rel_date_FN_I_C')
            time_since_exposure(sim, testobj, end=cur+l, uids=FN_uids_I_nC,                               attr='rel_date_FN_I_nC')
            time_since_exposure(sim, testobj, end=cur+l, uids=FN_uids_nI_C,                               attr='rel_date_FN_nI_C')
            time_since_exposure(sim, testobj, end=cur+l, uids=FN_uids_nI_nC,                              attr='rel_date_FN_nI_nC')

        else:  # Pool tests that are recorded after the outbreak ends
            testobj.pos_overflow += np.sum(test_results[test_latencies == l])
            testobj.neg_overflow += np.sum(~test_results[test_latencies == l])
    return 

def record_pipeline(testobj, sim, all_eligible, all_seekers, test_uids):
    criteria = all_eligible.keys()
    record_eligible_groups(testobj, sim, criteria, all_eligible)
    record_seek_groups(testobj, sim, criteria, all_seekers)
    record_allocate_groups(testobj, sim, criteria, all_eligible, test_uids)
    return 

def record_eligible_groups(testobj, sim, criteria, all_eligible): 
    '''
    Get daily number of test-eligible agents, stratified by:
        (1) Eligibility criteria
        (2) Disease status
    '''
    if sim.t == 0: 
        setattr(testobj, 'crit_groups', dict())
        for criterion in criteria: 
            testobj.crit_groups[criterion] = []

        setattr(testobj, 'crit_disease_groups', dict())
        for status in ['I_C', 'I_nC', 'nI_C', 'nI_nC']:
            testobj.crit_disease_groups[status] = []

    # Record daily test-eligible agents stratified by criteria
    for criterion in criteria:
        testobj.crit_groups[criterion].append(len(all_eligible[criterion])) 
    
    # Record daily test-eligible agents stratified by disease status
    eg_uids           = list(set().union(*all_eligible.values()))
    eligible          = np.zeros(sim.pars['pop_size_by_type']['human'], dtype=bool)
    eligible[eg_uids] = True

    testobj.crit_disease_groups['I_C'].append(np.sum(np.logical_and.reduce((eligible,    sim.agents.human.symptomatic_ILI,  sim.agents.human.exposed))))
    testobj.crit_disease_groups['I_nC'].append(np.sum(np.logical_and.reduce((eligible,   sim.agents.human.symptomatic_ILI, ~sim.agents.human.exposed))))
    testobj.crit_disease_groups['nI_C'].append(np.sum(np.logical_and.reduce((eligible,  ~sim.agents.human.symptomatic_ILI,  sim.agents.human.exposed))))
    testobj.crit_disease_groups['nI_nC'].append(np.sum(np.logical_and.reduce((eligible, ~sim.agents.human.symptomatic_ILI, ~sim.agents.human.exposed))))

def record_seek_groups(testobj, sim, criteria, all_seekers): 
    '''
    Get daily number of test-seeking agents, stratified by:
        (1) Eligibility criteria
        (2) Disease status
    '''
    if sim.t == 0: 
        setattr(testobj, 'seek_groups', dict())
        for criterion in criteria: 
            testobj.seek_groups[criterion] = []

        setattr(testobj, 'seek_disease_groups', dict())
        for status in ['I_C', 'I_nC', 'nI_C', 'nI_nC']:
            testobj.seek_disease_groups[status] = []

    # Record daily test-seeking agents stratified by criteria
    for criterion in criteria:
        testobj.seek_groups[criterion].append(len(all_seekers[criterion]))

    # Record daily test-seeking agents stratified by disease status
    sk_uids          = list(set().union(*[all_seekers[criterion] for criterion in criteria]))
    seekers          = np.zeros(sim.pars['pop_size_by_type']['human'], dtype=bool)
    seekers[sk_uids] = True

    testobj.seek_disease_groups['I_C'].append(np.sum(np.logical_and.reduce((seekers,    sim.agents.human.symptomatic_ILI,  sim.agents.human.exposed))))
    testobj.seek_disease_groups['I_nC'].append(np.sum(np.logical_and.reduce((seekers,   sim.agents.human.symptomatic_ILI, ~sim.agents.human.exposed))))
    testobj.seek_disease_groups['nI_C'].append(np.sum(np.logical_and.reduce((seekers,  ~sim.agents.human.symptomatic_ILI,  sim.agents.human.exposed))))
    testobj.seek_disease_groups['nI_nC'].append(np.sum(np.logical_and.reduce((seekers, ~sim.agents.human.symptomatic_ILI, ~sim.agents.human.exposed))))

def record_allocate_groups(testobj, sim, criteria, all_eligible, testers): 
    '''
    Among those that received a PCR test, identify which criteria each person satisfied - this is their "reason" for testing. 
    '''
    if sim.t == 0: 
        setattr(testobj, 'alloc_groups', dict())
        for criterion in criteria: 
            testobj.alloc_groups[criterion] = []
    
    testers = set(testers)

    for criterion in criteria:
        testobj.alloc_groups[criterion].append(len(testers.intersection(all_eligible[criterion])))

    return

def record_consumed(testobj, sim, tests): 
    if sim.t == 0: 
        setattr(testobj, 'tests_consumed', [])

    testobj.tests_consumed.append(tests)

def record_cons_days_eg_dist(testobj, sim, uids):
    '''
    Tracks the distribution of consecutive days test-eligible before eligibility ends.
    Stratifies the distribution by disease status.
    '''
    # Stratify identified agents by disease status
    I_C_uids   = uids[np.logical_and( sim.agents.human.symptomatic_ILI[uids],  sim.agents.human.exposed[uids])]
    I_nC_uids  = uids[np.logical_and( sim.agents.human.symptomatic_ILI[uids], ~sim.agents.human.exposed[uids])]
    nI_C_uids  = uids[np.logical_and(~sim.agents.human.symptomatic_ILI[uids],  sim.agents.human.exposed[uids])]
    nI_nC_uids = uids[np.logical_and(~sim.agents.human.symptomatic_ILI[uids], ~sim.agents.human.exposed[uids])]

    # if len(nI_C_uids) > 0:
    #     sample_uid = nI_C_uids[0]
    #     print(f"Day {sim.t}: Is agent {sample_uid} symptomatic? {sim.people.symptomatic[sample_uid]}")
    #     print(f"Day {sim.t}: Is agent {sample_uid} exposed? {sim.people.exposed[sample_uid]}")
    #     print(f"Day {sim.t}: Agent {sample_uid} got ILI on {sim.people.date_symptomatic_ILI[sample_uid]}")
    #     print(f"Day {sim.t}: Agent {sample_uid} got COVID on {sim.people.date_exposed[sample_uid]}")
    #     print(f"Day {sim.t}: Agent {sample_uid} recovered from ILI on {sim.people.date_recovered_ILI[sample_uid]}")
    #     print(f"Day {sim.t}: Agent {sample_uid} recovered from COVID on {sim.people.date_recovered[sample_uid]}")
    #     print(f"Day {sim.t}: Agent {sample_uid} started to have ILI symptoms on {sim.people.date_symptomatic_ILI[sample_uid]}")
    #     print(f"Day {sim.t}: Agent {sample_uid} started to have COVID symptoms on {sim.people.date_symptomatic[sample_uid]}")
    #     print()

    # Retrieve how many days they have been eligible
    cons_days       = testobj.cons_days_eg[uids]
    cons_days_I_C   = testobj.cons_days_eg[I_C_uids]
    cons_days_I_nC  = testobj.cons_days_eg[I_nC_uids]
    cons_days_nI_C  = testobj.cons_days_eg[nI_C_uids]
    cons_days_nI_nC = testobj.cons_days_eg[nI_nC_uids]

    # Update the associated distributions
    update_dist(testobj, cons_days,       'dist_egcons')
    update_dist(testobj, cons_days_I_C,   'dist_egcons_I_C')
    update_dist(testobj, cons_days_I_nC,  'dist_egcons_I_nC')
    update_dist(testobj, cons_days_nI_C,  'dist_egcons_nI_C')
    update_dist(testobj, cons_days_nI_nC, 'dist_egcons_nI_nC')

def record_cons_days_sk_dist(testobj, sim, uids):
    '''
    Tracks the distribution of consecutive days test-seeking before eligibility ends.
    Stratifies the distribution by disease status.
    '''
    # Stratify identified agents by disease status
    I_C_uids   = uids[np.logical_and( sim.agents.human.symptomatic_ILI[uids],  sim.agents.human.exposed[uids])]
    I_nC_uids  = uids[np.logical_and( sim.agents.human.symptomatic_ILI[uids], ~sim.agents.human.exposed[uids])]
    nI_C_uids  = uids[np.logical_and(~sim.agents.human.symptomatic_ILI[uids],  sim.agents.human.exposed[uids])]
    nI_nC_uids = uids[np.logical_and(~sim.agents.human.symptomatic_ILI[uids], ~sim.agents.human.exposed[uids])]

    # Retrieve how many days they have been eligible
    cons_days       = testobj.cons_days_sk[uids]
    cons_days_I_C   = testobj.cons_days_sk[I_C_uids]
    cons_days_I_nC  = testobj.cons_days_sk[I_nC_uids]
    cons_days_nI_C  = testobj.cons_days_sk[nI_C_uids]
    cons_days_nI_nC = testobj.cons_days_sk[nI_nC_uids]

    # Update the associated distributions
    update_dist(testobj, cons_days,       'dist_skcons')
    update_dist(testobj, cons_days_I_C,   'dist_skcons_I_C')
    update_dist(testobj, cons_days_I_nC,  'dist_skcons_I_nC')
    update_dist(testobj, cons_days_nI_C,  'dist_skcons_nI_C')
    update_dist(testobj, cons_days_nI_nC, 'dist_skcons_nI_nC')

def update_parameters(testobj, t):
    '''
    Check if parameter updates have been scheduled for the current day, and update accordingly. 
    '''
    # Throw error if first element of keys is a string. (only numbers work)
    if len(testobj.timeline.keys()) > 0:
        if isinstance(list(testobj.timeline.keys())[0], str):
            raise TypeError('Timeline keys must be integers.')

    if t in testobj.timeline.keys(): 
        for param, param_val in testobj.timeline[t].items():
            if not(param in testobj.__dict__.keys()):
                raise RuntimeError(f'The {param} parameter of {testobj.__class__} does not exist, update attempt failed.')
            else: 
                setattr(testobj, param, param_val)
        print(f'Parameters for {testobj.__class__} have been updated!')