import zoonosim as zn
import numpy as np
import json
import pandas as pd

project_name = 'Calibration_single_breed_AOQ_4th_iteration'

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
for key, value in project_pars['poultry_pars'].items():
    if key == 'breeds':
        project_pars['poultry_pars'][key] = np.array(value, dtype=zn.default_str)

data = pd.read_csv("zoonosim/data/CFIA_monthly_incidence.csv")
data['date'] = pd.to_datetime(data['date'])
data = data.rename(columns={"AOQ_adjusted_50":"monthly_new_poultry_flock_infectious"})
data['monthly_new_human_infectious'] = np.repeat(0, len(data))
data = data[["date", "monthly_new_human_infectious", "monthly_new_poultry_flock_infectious"]]

project_sim = zn.Sim(pars = project_pars, label = project_name, datafile=data, recursive = False)
project_sim.initialize()
best_pars = dict(
    rand_seed = 49,
)

project_sim.update_pars(best_pars, recursive=True)

# project_sim.run()
# project_sim.plot()

project_msim = zn.MultiSim(project_sim, n_runs = 250, verbose = 0)

if __name__ == "__main__":
    project_msim.run(keep_people = True, run_args=dict(auto_finalize=False)) 
    project_msim.finalize()
    project_msim.shrink()                   
    project_msim.save(f'msims/{project_name}_w_data.msim')  # Save the multi-simulation object.
    project_msim.combine()                # Combine the results from all simulations.
    project_msim.summarize()              # Summarize the combined results.
    project_msim.plot()                   # Plot the results.
    project_msim.plot_result(key = 'monthly_new_poultry_flock_infectious')                # Plot detailed results.
    project_msim.plot_result(key = 'monthly_new_human_infectious')            # Plot detailed results.
    #project_msim.plot_result(key = 'cum_layer_flock_infectious')                # Plot detailed results.
    #project_msim.plot_result(key = 'cum_broiler_flock_infectious')                # Plot detailed results.