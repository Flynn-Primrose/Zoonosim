import numpy as np
from .. import utils as znu

'''
This file implements functions for checking test-eligibility.
'''

def check_criteria(base_testobj, testobj, sim, **crit_kwargs): 
        '''
        Given the criterion names, call respective criterion functions to identify who meets each criterion. Apply restrictions at the end.
        '''
        all_eligible = dict()  # Will be a dictionary of sets
        RAT = testobj.RAT
        criteria = base_testobj.criteria

        for name in criteria: 

            # Diagnostic
            # if name == 'cont_conf': 
            #     all_eligible['cont_conf'] = contact_conf(sim)
            # if name == 'cont_conf_quar':
            #     all_eligible['cont_conf_quar'] = contact_conf_quar(sim)
            # if name == 'cont_vuln': 
            #     all_eligible['cont_vuln'] = contact_vuln(sim)
            # if name == 'cont_sx': 
            #     all_eligible['cont_sx'] = contact_sx(sim)
            if name == 'sx': 
                all_eligible['sx'] = sx(sim)
            if name == 'sx_quar':
                all_eligible['sx_quar'] = sx_quar(sim)
            if name == 'sev_crit':
                all_eligible['sev_crit'] = sev_crit(sim)
            if name == 'alerted':
                all_eligible['alerted'] = alerted(sim)
            if name == 'other':
                if 'other' in crit_kwargs.keys(): 
                    all_eligible['other'] = other(sim, **crit_kwargs['other'])
                else:
                    all_eligible['other'] = other(sim)
            if name == 'pos_RAT': 
                all_eligible['pos_RAT'] = RAT_testers(sim, RAT, mode='pos_RAT')
            if name == 'neg_RAT': 
                all_eligible['neg_RAT'] = RAT_testers(sim, testobj, mode='neg_RAT')

            # Screening
            if name == 'work': 
                if 'work' in crit_kwargs.keys():
                    all_eligible['work'] = work(sim, base_testobj, **crit_kwargs['work'])
                else:
                    all_eligible['work'] = work(sim, base_testobj)

        # Apply retest restrictions to those seeking diagnostic testing and finally return
        return apply_restrictions(testobj, base_testobj, all_eligible)

# def contact_conf(sim):
#     '''
#     Select people who contacted a confirmed case in the previous day
#     '''
#     contacted_inds = cvu.true(sim.people.date_known_contact == sim.t - 1)
#     contacts = set(contacted_inds.tolist())
#     return contacts

# def contact_conf_quar(sim, d=7):
#     '''
#     Select people who have contacted confirmed case in last 7 days and they are not quarantined)
#     '''
#     prev_contacted = np.logical_and(sim.t - sim.people.date_known_contact <= d, sim.t - sim.people.date_known_contact >= 0)  # logical_and not strictly necessary
#     not_quarantined = sim.people.quarantined == False
#     return set(cvu.true(np.logical_and(prev_contacted, not_quarantined)).tolist())

# def contact_vuln(sim):  
#     '''
#     Get people who are in contact with vulnerable people in the household. Vulnerability defined in terms of odds ratio of susceptibility. 
#     '''
#     # Get people who are vulnerable (and not dead) 
#     vuln_uid = cvu.true(np.logical_and(sim.people.rel_sus > 1.0, ~sim.people.dead))  # See parameters.py for assignment of susceptibility odds ratios by age bracket

#     # Get their household contacts. It may be more realistic to only get some subset of the household contacts.
#     return sim.people.contacts['h'].find_contacts(vuln_uid, as_array=False)  # Sped up by Numba

# def contact_sx(sim, d=3): 
#     '''
#     Get people who were in contact with a symptomatic person, in the last d days. 
    
#     Args: 
#         sim            (sim) : Simulation object 
#         d              (int) : Contacts of symptomatic people, in the last d days, satisfy this criterion
#         mode           (str) : Choose return option

#     Returns: 
#         contacts       (set) : Set of contacts. 
#     '''
#     all_sx = cvu.true(np.logical_or(sim.people.symptomatic, sim.people.symptomatic_ILI))  # Anyone who is currently symptomatic with COVID or ILI

#     all_contacts = set() 
#     for layer in sim.people.contacts.keys(): 
#         all_contacts.update(sim.people.contacts[layer].find_contacts(all_sx, as_array=False))  # as_array=True will sort indices, which is unnecessary
#     return all_contacts

def sx(sim): 
    '''
    Get people who are currently displaying symptoms.

    Args: 
        sim          (sim) : Simulation object 
        mode         (str) : Choose return option 

    Returns: 
        symptomatic  (set) : Set of people who are symptomatic with COVID or ILI. 
    '''
    all_sx = znu.true(np.logical_or(sim.agents.human.symptomatic, sim.agents.human.symptomatic_ILI))
    return set(all_sx.tolist())

def sx_quar(sim, d=6):
    '''
    Get people who are currently quarantined and the date of symptom onset was more than d days ago.
    Asymptomatic individuals will never invoke this criterion (date_symptomatic and date_symptomatic_ILI are both np.nan).
    '''
    sx_covid = np.logical_and(sim.t - sim.agents.human.date_symptomatic >= d, sim.agents.human.symptomatic)
    sx_ILI = np.logical_and(sim.t - sim.agents.human.date_symptomatic_ILI >= d, sim.agents.human.symptomatic_ILI)
    sx = np.logical_or(sx_covid, sx_ILI)
    sx_quar = znu.true(np.logical_and(sim.agents.human.quarantined == True, sx))
    return set(sx_quar.tolist())

def sev_crit(sim): 
    '''
    Get people who are currently severe or critical.
    '''
    #sev_crit = cvu.true(np.logical_or(sim.people.severe, sim.people.critical))
    sev_crit = znu.true(sim.agents.human.severe)
    return set(sev_crit.tolist())

def RAT_testers(sim, RAT, mode='pos_RAT', d=3): 
    '''
    Get people who received a positive RAT or negative RAT in the past d days.
    '''
    if RAT is None: 
        raise RuntimeError('RAT_testers called but no RAT system was provided.')
    cur = sim.t
    date_positive = RAT.date_positive
    date_negative = RAT.date_negative

    # Get people who had a positive/negative RAT previously (and are not dead) 
    pos_inds = znu.true(np.logical_and(date_positive >= 0, ~sim.agents.human.dead))
    neg_inds = znu.true(np.logical_and(date_negative >= 0, ~sim.agents.human.dead))

    # Get people who were confirmed positive/negative by a RAT test in the last d days
    pos_d_inds = pos_inds[np.logical_and(cur - date_positive[pos_inds] <= d, cur - date_positive[pos_inds] >= 0)]  # 2nd part of "and" statement is needed since reporting date may be in the future
    neg_d_inds = neg_inds[np.logical_and(cur - date_negative[neg_inds] <= d, cur - date_negative[neg_inds] >= 0)]

    if mode == 'pos_RAT': 
        # Return set of individuals who got a positive RAT test in past d days
        return set(pos_d_inds.tolist())

    elif mode == 'neg_RAT': 
        return set(neg_d_inds.tolist())
    
def alerted(sim):
    '''
    Get people who have been alerted by the system to seek a test. 
    '''
    if (not sim.pars['enable_smartwatches']):
        raise RuntimeError('Smartwatches are not enabled, but alerted function was called. Please enable smartwatches in the simulation parameters.')
    
    all_alerted = znu.true(sim.agents.human.alerted)
    return set(all_alerted.tolist())

def other(sim, p=0.05):
    return set(np.random.choice(sim.pars['pop_size_by_type']['human'], int(p * sim.pars['pop_size_by_type']['human']), replace=False).tolist())

def work(sim, base_testobj, p=0.10):
    '''
    Get people who have to do mandatory workplace testing. Randomly select some number of workplaces. Every week, these people will seek a test. 
    To consider: Do we want them to have priority access to testing resources? 

    Args: 
        p    (float) : Proportion of workplaces to undergo mandatory testing
    '''
    if sim.t == 0:
        # Identify which workplaces are chosen to do mandatory testing, and save the workers
        test_fid = np.random.choice(list(base_testobj.fid2uid.keys()), int(base_testobj.n_fids * p), replace=False)
        #print(base_testobj.fid2uid)
        #print(f'Workplaces selected for mandatory testing: {test_fid}')
        worker_uids = np.array([], dtype=int)
        for fid in test_fid:
            worker_uids = np.union1d(worker_uids, base_testobj.fid2uid[fid]).astype(int) 
            # try: # For debugging purposes, in case fid2uid does not contain the fid
            #     worker_uids = np.union1d(worker_uids, base_testobj.fid2uid[fid]).astype(int)
            #     print(f'fid {fid} has {len(base_testobj.fid2uid[fid])} workers.')
            # except KeyError:
            #     print(f'fid {fid} not found in fid2uid mapping. Skipping this workplace.')
            #     print(f'fid2uid keys: {base_testobj.fid2uid.keys()}')
        
        worker_test_dates = np.random.choice(7, len(worker_uids))

        setattr(base_testobj, 'wuid', worker_uids)  # Set of all workers
        setattr(base_testobj, 'wuid_td', worker_test_dates)  # Days of the week they are assigned to test

    # Remove anyone who died
    base_testobj.wuid_td = base_testobj.wuid_td[~sim.agents.human.dead[base_testobj.wuid]]  # Order matters
    base_testobj.wuid = base_testobj.wuid[~sim.agents.human.dead[base_testobj.wuid]]
    #print(f'workers selected for mandatory testing: {base_testobj.wuid}') # For debugging purposes

    # Same workplaces, and thus workers, test each time
    return set(base_testobj.wuid[base_testobj.wuid_td == sim.t % 7].tolist())

# ----------------
# Helper functions
# ----------------
def apply_restrictions(testobj, base_testobj, eligible_dict): 
    '''
    Remove individuals from the eligibility dictionary according to a set of hard restrictions. 
    These individuals will not proceed to seek a test
    '''
    restricted = set()

    if base_testobj.system == 'Diagnostic':

        # Restrict individuals who have tested recently (positive/negative) and are trying to test again
        uids = (testobj.retest_tracker > 0).nonzero()[0]  # negative values are fine
        restricted.update(uids)

        # Iterate through diagnostic criteria and remove anyone who is restricted
        for criterion in testobj.di_criteria:
            eligible_dict[criterion] = eligible_dict[criterion] - eligible_dict[criterion].intersection(restricted)

    elif base_testobj.system == 'Screening':
        # Iterate through screening criteria and remove anyone who is restricted -- currently no restrictions
        for criterion in testobj.sc_criteria:
            eligible_dict[criterion] = eligible_dict[criterion] - eligible_dict[criterion].intersection(restricted)

    return eligible_dict