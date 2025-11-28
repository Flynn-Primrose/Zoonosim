import zoonosim as zn
import numpy as np 
# Define "best" parameters from calibration"
best_pars = dict(
    n_farms = 175,
    start_day = '2022-01-01',
    end_day = '2025-12-31',
    rand_seed = 42,
    beta = dict(
        human = 0.001,
        flock = 0.826,
        barn = 1.390,
        water = 0.962,
    ),
    n_imports = dict(
        human = 0,
        flock = 0,
        barn = 0.536,
        water = 0.229,
    ),
    dur = dict(
        human = {
            # Duration: disease progression
            'exp2inf': dict(dist='lognormal_int', par1=2.0, par2=1.0), # Duration from exposed to infectious
            'inf2sym': dict(dist='lognormal_int', par1=2.0, par2=1.0), # Duration from infectious to symptomatic
            'sym2sev': dict(dist='lognormal_int', par1=5.0, par2=2.0), # Duration from symptomatic to severe symptoms
            # Duration: Recovery
            'asym2rec': dict(dist='lognormal_int', par1=9.0,  par2=4.0), # Duration for asymptomatic people to recover
            'mild2rec': dict(dist='lognormal_int', par1=9.0,  par2=4.0), # Duration for people with mild symptoms to recover
            'sev2rec': dict(dist='lognormal_int', par1=9.0, par2=4.0), # Duration for people with severe symptoms to recover
            'sev2die': dict(dist='lognormal_int', par1=9.0, par2=4.0), # Duration from critical symptoms to death
            # Duration: quarantine
            'quar': 14,   
            # Duration: diagnosis
            'diag': 14
        },
        flock = {
            # Duration: disease progression
            'exp2inf': dict(dist='lognormal_int', par1=2.0, par2=1.0), # Duration from exposed to infectious. 
            'inf2out': dict(dist='lognormal_int', par1=2.0, par2=1.0), # Duration from infectious to recovery/removal. 
            'susp2res': dict(dist='lognormal_int', par1=5.0, par2=1.0), # Duration from suspicion to a definitive test result. 
            # Duration: Quarantine
            'quar': 14
        },
        barn = {
            'contamination': dict(dist='lognormal_int', par1=14, par2=5.0), # Duration of contamination.
            'composting': dict(dist='lognormal_int', par1=7.0, par2=1.0), # Duration of composting. 
            'cleaning': dict(dist='lognormal_int', par1=7.0, par2=1.0), # Duration of cleaning process. 
        },
        water = {
            'contamination': dict(dist='lognormal_int', par1=14, par2=5.0), # Duration of contamination.
        }
    ),
    prognoses = dict(
        human = zn.parameters.relative_human_prognoses(dict(
            age_cutoffs   = np.array([0,       10,      20,      30,      40,      50,      60,      70,      80,     90,]),     # Age cutoffs (lower limits)
            sus_ORs       = np.array([1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.0,     1.00,   1.00]),    # Odds ratios for relative susceptibility 
            trans_ORs     = np.array([1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,   1.00]),    # Odds ratios for relative transmissibility
            comorbidities = np.array([1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,   1.00]),    # Comorbidities by age -- set to 1 by default since already included in disease progression rates
            symp_probs    = np.array([0.66,    0.66,    0.66,    0.66,    0.66,    0.66,    0.66,    0.66,    0.66,   0.66]),    # Overall probability of developing symptoms 
            severe_probs  = np.array([0.33,    0.33,    0.33,    0.33,    0.33,    0.33,    0.33,    0.33,    0.33,   0.33]),     # Overall probability of developing severe symptoms
            death_probs   = np.array([0.33,    0.33,    0.33,    0.33,    0.33,    0.33,    0.33,    0.33,    0.33,   0.33]),    # Overall probability of dying
        )),
        flock = dict(
            breeds = np.array(['duck', 'broiler', 'layer'], dtype=zn.default_str),
            sus_ORs = np.array([2.00, 1.00, 1.00]),
            trans_ORs = np.array([1.00, 1.00, 1.00]),
            baseline_symptomatic_rate = np.array([0.001, 0.001, 0.001]),
            mean_symptomatic_rate_increase = np.array([0.001, 0.0001, 0.0001]),
            baseline_mortality_rate = np.array([0.001, 0.001, 0.001]),
            mean_mortality_rate_increase = np.array([0.002, 0.002, 0.002]),
            baseline_water_rate = np.array([1.00, 1.00, 1.00]),
            mean_water_rate_increase = np.array([1.00, 1.00, 1.00]),
        )
    ),
    production_cycle = dict(
        breeds = np.array(['duck', 'broiler', 'layer'], dtype=zn.default_str),
        cycle_dur = [dict(dist = 'normal_pos', par1 = 75, par2 = 10),
                    dict(dist = 'normal_pos', par1 = 45, par2 = 5),
                    dict(dist = 'normal_pos', par1 = 500, par2=100)],
        flock_size = [dict(dist = 'normal_pos', par1 = 500, par2 = 100),
                    dict(dist = 'normal_pos', par1 = 20000, par2 = 1000),
                    dict(dist = 'normal_pos', par1 = 2000, par2 = 100)]
    ),
    wild = dict(
        human = dict(
            rel_beta = 1.0,
            rel_symp_prob = 0.33,
            rel_severe_prob = 0.25,
            rel_death_prob = 0.01,
            rel_asymp_fact = 0.5

        ),
        flock = dict(
            rel_beta = 1.0,
            rel_symp_delta = 1.0,
            rel_death_delta = 1.0,
            rel_water_delta = 1.0
        ),
        barn = dict(
            rel_beta = 1.0,
            rel_dur_contamination = 1.0
        ),
        water = dict(
            rel_beta = 1.0,
            rel_dur_contamination = 1.0
        )
    ),
)

best_sim = zn.Sim(datafile="zoonosim/data/H5N1_cases_in_QC_poultry.csv", label = "Calibration8_best", pars=best_pars)

msim = zn.MultiSim(best_sim, n_runs = 25, verbose = 0)  # Wrap the simulation in a MultiSim object.

if __name__ == "__main__":
    msim.run()                    # Run the simulations.
    msim.save('msims/calibration1_best.msim')  # Save the multi-simulation object.
    msim.combine()                # Combine the results from all simulations.
    msim.summarize()              # Summarize the combined results.
    msim.plot()                   # Plot the results.
    msim.plot_result(key = 'cum_duck_flock_infectious')            # Plot detailed results.
    msim.plot_result(key = 'cum_layer_flock_infectious')                # Plot detailed results.
    msim.plot_result(key = 'cum_broiler_flock_infectious')                # Plot detailed results.