import numpy as np
from math import ceil

from .. import utils as znu

'''
This file implements a simple model of background ILI.
'''

def infect_ILI(sim): 
    '''
    Void function that infects some proportion of the population with non-COVID ILI, refreshed each week. Meant to be called at each time step. We assume 
    that people cannot be infected with both COVID and ILI at the same time. Unsure if this is realistic. 
    
    Assumptions: 
        - Everyone infected is symptomatic
        - Everyone infected will remain infected for a week
        - Infections occur spontaneously, and does not involve transmission

    Args: 
        sim  (simulation) : Initialized simulation object 
        p_infect  (float) : Proportion of population to infect. For now, assume 8% of population has ILI on any given day
    '''
    if isinstance(sim.pars['bkg_ILI'], (float, np.floating)):
        ILI_float(sim)
    else: 
        ILI_array(sim)
        return


def ILI_float(sim):
    '''
    Args:
        sim          : Simulation object
        COVID_to_ILI : Whether to initialize ILI infections in individuals who have COVID
    '''
    p_infect     = sim.pars['bkg_ILI']
    Avian_to_ILI = sim.pars['Avian_to_ILI']
    cur          = sim.t
    pop_size     = sim.pars['pop_size_by_type']['human']
    n_infected   = int(pop_size * p_infect)  # background level of ILI

    if cur == 0:  # Initialize tracker of whether people have symptomatic ILI or not
        sim.agents.human.unlock()  # Redundancy 
        sim.agents.human['symptomatic_ILI'] = np.full(pop_size, False)  # Whether each individual has symptomatic ILI 
        sim.agents.human['date_symptomatic_ILI'] = np.full(pop_size, np.nan)  # Date that a person last became symptomatic with ILI 
        sim.agents.human['date_recovered_ILI'] = np.full(pop_size, np.nan)  # Date that a person last recovered from ILI
        sim.agents.human['cons_days_ILI'] = np.zeros(pop_size)  # NOTE: Stored to ensure no back-to-back infections
        
        if Avian_to_ILI:
            infect_inds = np.random.choice(sim.agents.human['uid'], n_infected, replace=False).astype(int)
        else:
            no_avian = znu.true(sim.agents.people.exposed == False)
            if len(no_avian) >= n_infected:
                infect_inds = np.random.choice(no_avian, n_infected, replace=False).astype(int)
            elif len(no_avian) == 0:
                infect_inds = np.array([], dtype=int)
            else:
                infect_inds = no_avian  # everyone without covid gets selected

        # Initialize how long until they recover
        recov_times = np.random.choice(np.arange(7, 11), len(infect_inds)).astype(int)  # Recover individuals in 7-10 days

        # Artificially run the infection-recovery process until convergence in distribution of recovery times
        num_days = 1000
        for k in range(num_days):
            recover_today = np.nonzero(recov_times == k)[0]
            new_recov_times = k + np.random.choice(np.arange(7, 11), len(recover_today)).astype(int)
            recov_times[recover_today] = new_recov_times

        # Shift the earliest recovery time back to day one (not day zero! otherwise they will not get checked for recovery)
        recov_times = recov_times - np.min(recov_times) + 1

        # Update trackers
        sim.agents.human['date_symptomatic_ILI'][infect_inds] = cur
        sim.agents.human['date_recovered_ILI'][infect_inds] = cur + recov_times
        sim.agents.human['symptomatic_ILI'][infect_inds] = True

    else: 
        # Determine candidates for new infections
        no_ILI          = znu.true(sim.agents.human['symptomatic_ILI'] == False)
        no_avian_or_ILI = znu.true(np.logical_and(sim.agents.human.exposed == False, sim.agents.human['symptomatic_ILI'] == False))

        # Check who is supposed to recover today, and recover them. This happens after previous step to prevent back-to-back ILI infections.
        recov_today = znu.true(sim.agents.human['date_recovered_ILI'] == cur)
        sim.agents.human['symptomatic_ILI'][recov_today] = False

        num_to_infect = n_infected - np.sum(sim.agents.human['symptomatic_ILI'])  # always positive, size of the gap we need to fill with new ILI infections
        # num_to_infect = len(recov_today)  # NOTE: This creates issues when there are not enough agents to infect (permanently decrements level of background ILI)

        # Infect new people: (Case 1) Allow agents with COVID to be infected with ILI (Case 2) Do not allow agents with COVID to be infected with ILI
        if Avian_to_ILI:
            infect_inds = np.random.choice(no_ILI, num_to_infect, replace=False).astype(int)
        else:
            if len(no_avian_or_ILI) >= num_to_infect:
                infect_inds = np.random.choice(no_avian_or_ILI, num_to_infect, replace=False).astype(int)
            elif len(no_avian_or_ILI) == 0:
                infect_inds = np.array([], dtype=int)
            else:
                infect_inds = no_avian_or_ILI

        # Determine how long until they recover
        recov_times = np.random.choice(np.arange(7, 11), len(infect_inds))  # Recover individuals in 7-10 days
        sim.agents.human['date_symptomatic_ILI'][infect_inds] = cur
        sim.agents.human['date_recovered_ILI'][infect_inds] = cur + recov_times
        sim.agents.human['symptomatic_ILI'][infect_inds] = True

    # Update counter for consecutive days with ILI
    sim.agents.human['cons_days_ILI'][sim.agents.human['symptomatic_ILI']] += 1
    sim.agents.human['cons_days_ILI'][~sim.agents.human['symptomatic_ILI']] = 0


def ILI_array(sim):
    p_infect = sim.pars['bkg_ILI']
    cur = sim.t
    cur_week = ceil( (cur+1) / 7)
    pop_size = sim.pars['pop_size_by_type']['human']  # Population size of humans
    n_days = sim.pars['n_days']
    recov_time = 7

    if len(p_infect) < ceil( (n_days+1) / 7):
        raise RuntimeError('Array length of weekly background ILI rates is not enough to cover the full simulation.')
    else:
        if cur == 0: 
            sim.agents.human.unlock()
            sim.agents.human['symptomatic_ILI'] = np.full(pop_size, False)  # Whether each individual has symptomatic ILI 
            sim.agents.human['date_symptomatic_ILI'] = np.full(pop_size, np.nan)  # Date that a person last became symptomatic with ILI 
            sim.agents.human['date_recovered_ILI'] = np.full(pop_size, np.nan)  # Date that a person last recovered from ILI

            n_infected = int(pop_size * p_infect[0])  # Use rate for first week
            no_avian = znu.true(sim.agents.human.exposed == False)
            infect_inds = np.random.choice(no_avian, n_infected, replace=False)

            sim.agents.human['date_symptomatic_ILI'][infect_inds] = cur 
            sim.agents.human['date_recovered_ILI'][infect_inds] = cur + recov_time
            sim.agents.human['symptomatic_ILI'][infect_inds] = True
        elif cur % recov_time == 0:  # First day of a new week
            no_covid_or_ILI = znu.true(np.logical_and(sim.agents.human.exposed == False, sim.agents.human['symptomatic_ILI'] == False))  

            sim.agents.human['symptomatic_ILI'][:] = False  # Everyone currently infected recovers

            # Number of people to infect based on ILI prevalence of current week
            n_infected = int(pop_size * p_infect[cur_week-1])
            infect_inds = np.random.choice(no_covid_or_ILI, n_infected, replace=False)  # This could break if too many people have covid

            # Cleanup and set recovery dates
            sim.agents.human['date_symptomatic_ILI'][infect_inds] = cur
            sim.agents.human['date_recovered_ILI'][infect_inds] = cur + recov_time
            sim.agents.human['symptomatic_ILI'][infect_inds] = True