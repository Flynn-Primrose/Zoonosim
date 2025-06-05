'''
Defines the behaviour of background influenza like illnesses
'''

import numpy as np
from .. import utils as znu
from math import ceil

# Very simple module of what background ILI looks like

def infect_ILI(sim, pathogen = 0): 
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
    if isinstance(sim.pathogens[pathogen].bkg_ILI, (float, np.floating)):
        ILI_float(sim, pathogen)
    else: 
        ILI_array(sim, pathogen)
        return


def ILI_float(sim, pathogen):
    p_infect = sim.pathogens[pathogen].bkg_ILI
    cur = sim.t
    pop_size = sim.people.pars['pop_size']
    n_infect = int(pop_size * p_infect)

    if cur == 0:  # Initialize tracker of whether people have symptomatic ILI or not
        sim.people.unlock()  # Redundancy 
        sim.people['symptomatic_ILI'] = np.full(pop_size, False)  # Whether each individual has symptomatic ILI 
        sim.people['date_symptomatic_ILI'] = np.full(pop_size, np.nan)  # Date that a person last became symptomatic with ILI 
        sim.people['date_recovered_ILI'] = np.full(pop_size, np.nan)  # Date that a person last recovered from ILI

        # Initialize infections in individuals who do not have COVID
        no_covid = znu.true(sim.people.exposed == False)
        if len(no_covid) >= n_infect:
            infect_inds = np.random.choice(no_covid, n_infect, replace=False).astype(int)
        elif len(no_covid) == 0:
            infect_inds = np.array([], dtype=int)
        else:
            infect_inds = np.random.choice(no_covid, n_infect, replace=True).astype(int)  # everyone without covid gets selected

        # Determine how long until they recover
        recov_times = np.random.choice(np.arange(7, 11), len(infect_inds)).astype(int)  # Recover individuals in 7-10 days
        sim.people['date_symptomatic_ILI'][infect_inds] = cur
        sim.people['date_recovered_ILI'][infect_inds] = cur + recov_times
        sim.people['symptomatic_ILI'][infect_inds] = True

    else: 
        # Find candidates for new infections - not infected with COVID or ILI
        no_covid_or_ILI = znu.true(np.logical_and(sim.people.exposed == False, sim.people['symptomatic_ILI'] == False))  

        # Check who is supposed to recover today, and recover them. This happens after previous step to prevent back-to-back infections.
        recov_today = znu.true(sim.people['date_recovered_ILI'] == cur)
        sim.people['symptomatic_ILI'][recov_today] = False

        # Infect new people, reproduction number is not 1 if n_infect quota was not met in the last infection cycle
        if len(no_covid_or_ILI) >= n_infect:
            infect_inds = np.random.choice(no_covid_or_ILI, n_infect, replace=False).astype(int)
        elif len(no_covid_or_ILI) == 0:
            infect_inds = np.array([], dtype=int)
        else:
            infect_inds = np.random.choice(no_covid_or_ILI, n_infect, replace=True).astype(int)  # everyone without covid or ILI gets selected

        # Determine how long until they recover
        recov_times = np.random.choice(np.arange(7, 11), len(infect_inds))  # Recover individuals in 7-10 days
        sim.people['date_symptomatic_ILI'][infect_inds] = cur
        sim.people['date_recovered_ILI'][infect_inds] = cur + recov_times
        sim.people['symptomatic_ILI'][infect_inds] = True


def ILI_array(sim, pathogen):
    p_infect = sim.pathogens[pathogen].bkg_ILI
    cur = sim.t
    cur_week = ceil( (cur+1) / 7)
    pop_size = sim.people.pars['pop_size']
    n_days = sim.pars['n_days']
    recov_time = 7

    if len(p_infect) < ceil( (n_days+1) / 7):
        raise RuntimeError('Array length of weekly background ILI rates is not enough to cover the full simulation.')
    else:
        if cur == 0: 
            sim.people.unlock()
            sim.people['symptomatic_ILI'] = np.full(pop_size, False)  # Whether each individual has symptomatic ILI 
            sim.people['date_symptomatic_ILI'] = np.full(pop_size, np.nan)  # Date that a person last became symptomatic with ILI 
            sim.people['date_recovered_ILI'] = np.full(pop_size, np.nan)  # Date that a person last recovered from ILI

            n_infect = int(pop_size * p_infect[0])  # Use rate for first week
            no_covid = znu.true(sim.people.exposed == False)
            infect_inds = np.random.choice(no_covid, n_infect, replace=False)

            sim.people['date_symptomatic_ILI'][infect_inds] = cur 
            sim.people['date_recovered_ILI'][infect_inds] = cur + recov_time
            sim.people['symptomatic_ILI'][infect_inds] = True
        elif cur % recov_time == 0:  # First day of a new week
            no_covid_or_ILI = znu.true(np.logical_and(sim.people.exposed == False, sim.people['symptomatic_ILI'] == False))  

            sim.people['symptomatic_ILI'][:] = False  # Everyone currently infected recovers

            # Number of people to infect based on ILI prevalence of current week
            n_infect = int(pop_size * p_infect[cur_week-1])
            infect_inds = np.random.choice(no_covid_or_ILI, n_infect, replace=False)  # This could break if too many people have covid

            # Cleanup and set recovery dates
            sim.people['date_symptomatic_ILI'][infect_inds] = cur
            sim.people['date_recovered_ILI'][infect_inds] = cur + recov_time
            sim.people['symptomatic_ILI'][infect_inds] = True