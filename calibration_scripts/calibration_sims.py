import zoonosim as zn
import numpy as np
import json

project_name = 'Calibration3_3'

pars_filename = f"saved_pars/{project_name}.json"

with open(pars_filename, 'r') as file:
    project_pars = json.load(file)

for key, value in project_pars['prognoses']['human'].items():
    project_pars['prognoses']['human'][key] = np.array(value)
for key, value in project_pars['prognoses']['flock'].items():
    if key == 'breed':
        project_pars['prognoses']['flock'][key] = np.array(value, dtype=zn.default_str)
    else:
        project_pars['prognoses']['flock'][key] = np.array(value)
for key, value in project_pars['production_cycle'].items():
    if key == 'breeds':
        project_pars['production_cycle'][key] = np.array(value, dtype=zn.default_str)

project_sim = zn.Sim(pars = project_pars, label = project_name, datafile="zoonosim/data/H5N1_cases_in_QC_poultry.csv", recursive = False)
project_sim.initialize()
best_pars = dict(
    rand_seed = 42,
    # beta = dict(
    #     human = 0.25,
    #     flock = 1.94,
    #     barn = 2.92,
    #     water = 3.24
    # ),
    beta_layer = dict(
        fb = 4.72, 
        bw = 0.77, 
        fw = 0.09, 
        hb = 2.19,  
        hf = 2.30,   
        dw = 0.73, 
    )
)

project_sim.update_pars(best_pars, recursive=True)

project_msim = zn.MultiSim(project_sim, n_runs = 50, verbose = 0)

if __name__ == "__main__":
    project_msim.run()                    # Run the simulations.
    project_msim.save(f'msims/{project_name}_best.msim')  # Save the multi-simulation object.
    project_msim.combine()                # Combine the results from all simulations.
    project_msim.summarize()              # Summarize the combined results.
    project_msim.plot()                   # Plot the results.
    project_msim.plot_result(key = 'cum_duck_flock_infectious')            # Plot detailed results.
    project_msim.plot_result(key = 'cum_layer_flock_infectious')                # Plot detailed results.
    project_msim.plot_result(key = 'cum_broiler_flock_infectious')                # Plot detailed results.