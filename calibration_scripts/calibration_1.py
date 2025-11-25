import zoonosim as zn
import numpy as np

# Define new parameters for the simulation
new_pars = dict(
    start_day = '2022-01-01',
    end_day = '2025-12-31',
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
        cycle_dur = [dict(dist = 'normal_pos', par1 = 600, par2 = 50),
                    dict(dist = 'normal_pos', par1 = 45, par2 = 5),
                    dict(dist = 'normal_pos', par1 = 150, par2=25)],
        flock_size = [dict(dist = 'normal_pos', par1 = 20000, par2 = 1000),
                    dict(dist = 'normal_pos', par1 = 20000, par2 = 1000),
                    dict(dist = 'normal_pos', par1 = 20000, par2 = 1000)]
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

# Create Agents 
def alt_make_agents():
    return

def alt_make_popdict():
    return

def alt_make_contacts():
    return

# Create Simulation
sim = zn.Sim(datafile="zoonosim/data/H5N1_cases_in_QC_poultry.csv", label = "Calibration1", pars=new_pars, agents = alt_make_agents())

# Define calibration parameters
calib_pars = dict(    
    beta = dict(
        human = [0.5, 0.0, 2.0],
        flock = [0.5, 0.0, 2.0],
        barn = [0.5, 0.0, 2.0],
        water = [0.5, 0.0, 2.0],
    ),
    n_imports = dict(
        barn = [0.1, 0, 1],
        water = [0.1, 0, 1],
    )
)

calib = zn.Calibration(sim, calib_pars, name = "Calibration1", n_trials=25, die=True, keep_db=True)

if __name__ == "__main__":
    calib.calibrate()
