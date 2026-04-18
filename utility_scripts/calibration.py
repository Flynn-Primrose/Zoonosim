import zoonosim as zn
import numpy as np
import pandas as pd

project_name = "Calibration_single_breed_AOQ_4th_iteration"

# Define new parameters for the simulation
new_pars = dict(
    rand_seed = 92,
    n_farms = 50,
    start_day = '2022-01-01',
    end_day = '2025-12-31',
    record_all_events = False, # Whether to record all events in the simulation (True) or just transmission events.
    pop_pars = dict(
        avg_barns_per_farm = 5.0,
        avg_humans_per_barn = 1.5,
        avg_water_per_farm = 0.75,
        number_of_transients = 3,
        visits_per_day = 3, # Number of farms each transient visits in a day
    ),
    beta = dict(
        human = 0.0,
        ppe = 0.2,
        flock = 0.3,
        barn = 0.1,
        water = 0.1,
    ),
    n_imports = dict(
        human=None,  # Number of imported human cases per day; None = disabled
        flock=None,  # Number of imported flock cases per day; None = disabled
        ppe = None,   # Number of imported PPE contaminations per day; None = disabled
        barn=dict(import_pattern='seasonal', max_import_rate=0.2, peak_day=300),  # Number of imported barn contaminations per day; None = disabled
        water=dict(import_pattern='seasonal', max_import_rate=0.2, peak_day=300),  # Number of imported water contaminations per day; None = disabled
    ), 
    enable_smartwatches = False,
    smartwatch_pars = dict(
        who                       = 'all', # Must be one of 'all', 'permanent', 'transient'; controls who receives smartwatches
        mean_fpr                  =   0.08, # mean false positive rate
        use_variable_fpr          = True, # Whether to use a variable false positive rate
        day_i                     = np.arange(-21, 22, 1), #
        loc                       = 3.25, # Day of max probability of alert, relative to the day of symptom onset.
        alpha                     = 1, # Scales the probability of receiving an alert
        usage_rate                = 1, # Out of people who have smartwatches, the amount who use download the alerting app and stick with it.
        compliance_rate           =   0.05, # probability of quarantining if a smartwatch detects symptoms (only used if testobjs are not available)
        participation_rate        =   0.3,  # proportion of the population that has a smartwatch
    ),
    dynam_layer = dict(
                hp = 0.0, 
                hh = 0.0, 
                hf = 0.0, 
                hb = 0.0, 
                hw = 0.0, 
                pp = 0.0, 
                pf = 0.0, 
                pb = 0.0, 
                pw = 0.0, 
                fb = 0.0, 
                fw = 0.0, 
                bw = 0.0, 
                transient = 1.0
                ),
    beta_layer = dict(
                hp = 0.1, 
                hh = 0.1, 
                hf = 0.1, 
                hb = 0.1, 
                hw = 0.1, 
                pp = 0.1, 
                pf = 0.1, 
                pb = 0.1, 
                pw = 0.1, 
                fb = 0.1, 
                fw = 0.1, 
                bw = 0.1, 
                transient = 0.1  
                ),
    quar_factor = dict(
                hp = 0.0, 
                hh = 0.0, 
                hf = 0.0, 
                hb = 0.0, 
                hw = 0.0, 
                pp = 0.0, 
                pf = 0.0, 
                pb = 0.0, 
                pw = 0.0, 
                fb = 0.0, 
                fw = 0.0, 
                bw = 0.0, 
                transient = 0.0
                ),
    dur = dict(
        human = {
            # Duration: disease progression
            'exp2inf': dict(dist='lognormal_int', par1=5.0, par2=2.0), # Duration from exposed to infectious
            'inf2sym': dict(dist='lognormal_int', par1=5.0, par2=2.0), # Duration from infectious to symptomatic
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
        ppe = {
            'contamination': dict(dist='lognormal_int', par1=2.0, par2=1.0), # Duration of PPE contamination
            'quar': 14, # Duration of PPE quarantine (should match human quarantine duration since PPE is quarantined when associated human case is quarantined)
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
            trans_ORs     = np.array([0.01,    0.01,    0.01,    0.01,    0.01,    0.01,    0.01,    0.01,    0.01,   0.01]),    # Odds ratios for relative transmissibility
            comorbidities = np.array([1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,   1.00]),    # Comorbidities by age -- set to 1 by default since already included in disease progression rates
            symp_probs    = np.array([0.66,    0.66,    0.66,    0.66,    0.66,    0.66,    0.66,    0.66,    0.66,   0.66]),    # Overall probability of developing symptoms 
            severe_probs  = np.array([0.33,    0.33,    0.33,    0.33,    0.33,    0.33,    0.33,    0.33,    0.33,   0.33]),     # Overall probability of developing severe symptoms
            death_probs   = np.array([0.33,    0.33,    0.33,    0.33,    0.33,    0.33,    0.33,    0.33,    0.33,   0.33]),    # Overall probability of dying
        )),
        ppe = dict(
            sus_ORs = np.array([1.00]),
            trans_ORs = np.array([1.00]),
        ),
        flock = dict(
            breeds = np.array(['poultry'], dtype=zn.default_str),
            sus_ORs = np.array([1.00]),
            trans_ORs = np.array([1.00]),
            baseline_symptomatic_rate = np.array([0.001]),
            mean_symptomatic_rate_increase = np.array([0.01]),
            baseline_mortality_rate = np.array([0.001]),
            mean_mortality_rate_increase = np.array([0.01]),
            baseline_water_rate = np.array([1.00]),
            mean_water_rate_increase = np.array([1.00]),
        ),
        barn = dict(
            sus_ORs = np.array([1.00]),
            trans_ORs = np.array([1.00]),
        ),
        water = dict(
            sus_ORs = np.array([0.00]),
            trans_ORs = np.array([1.00]),
        )
    ),
    poultry_pars = dict(
        breeds = np.array(['poultry'], dtype=zn.default_str),
        breed_freqs = np.array([1.0]),
        mortality_suspicion_threshold = [0.0012], # I.E a deviation from the expected mortality rate of 0.01*expected_value will trigger suspicion
        symptomatic_suspicion_threshold = [0.0001], # I.E a deviation from the expected symptomatic rate of 0.01*expected_value will trigger suspicion
        consumption_suspicion_threshold = [0.001], # I.E a deviation from the expected rate of water consumption of 0.01*expected_value will trigger suspicion
        cycle_dur = [dict(dist = 'normal_pos', par1 = 100, par2 = 25)],
        flock_size = [dict(dist = 'normal_pos', par1 = 20000, par2 = 10000)]
    ),
    wild = dict(
        human = dict(
            rel_beta = 1.0,
            rel_symp_prob = 0.33,
            rel_severe_prob = 0.25,
            rel_death_prob = 0.01,
            rel_asymp_fact = 0.5

        ),
        ppe = dict(
            rel_beta = 1.0,
            rel_dur_contamination = 1.0
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

zn.options.set(verbose=0)
zn.options.set(numba_parallel='safe')
zn.options.set(numba_cache=False)


data = pd.read_csv("zoonosim/data/CFIA_monthly_incidence.csv")
data['date'] = pd.to_datetime(data['date'])
data = data.rename(columns={"AOQ_adjusted_50":"monthly_new_poultry_flock_infectious"})
data['monthly_new_human_infectious'] = np.repeat(0, len(data))
data = data[["date", "monthly_new_human_infectious", "monthly_new_poultry_flock_infectious"]]

# Create Simulation
sim = zn.Sim(datafile=data, label = project_name, pars=new_pars)
sim.export_pars(f"saved_pars/{project_name}.json")
# Define calibration parameters
calib_pars = dict(
        beta = dict(
        human = [0.2, 0.0, 0.5],
        ppe = [0.2, 0.0, 0.5],
        flock = [0.2, 0.0, 0.5],
        barn = [0.2, 0.0, 0.5],
        water = [0.2, 0.0, 0.5],
    ),   
    beta_layer = dict(
                hp = [0.2, 0.0, 0.5], 
                hh = [0.2, 0.0, 0.5], 
                hf = [0.2, 0.0, 0.5], 
                hb = [0.2, 0.0, 0.5], 
                hw = [0.2, 0.0, 0.5], 
                pp = [0.2, 0.0, 0.5], 
                pf = [0.2, 0.0, 0.5], 
                pb = [0.2, 0.0, 0.5], 
                pw = [0.2, 0.0, 0.5], 
                fb = [0.2, 0.0, 0.5], 
                fw = [0.2, 0.0, 0.5], 
                bw = [0.2, 0.0, 0.5], 
                transient = [0.2, 0.0, 0.5] 
                ),
)

calib = zn.Calibration(sim, calib_pars, name = project_name, n_reps = 10, total_trials=100, die=True, keep_db=True)

if __name__ == "__main__":
    calib.calibrate()