import zoonosim as zn
import numpy as np
import json
import pandas as pd

project_name = 'Calibration_single_breed_ON'

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

data = pd.read_csv("zoonosim/data/CFIA_cumulative_timeseries.csv")
data['date'] = pd.to_datetime(data['date'])
data = data.rename(columns={"cum_ON_poultry_infectious":"cum_poultry_flock_infectious"})

project_sim = zn.Sim(pars = project_pars, label = project_name, datafile=data, recursive = False)
project_sim.initialize()
best_pars = dict(
    rand_seed = 42,
    n_imports = dict(
        human = None,
        flock = None,
        barn = dict(peak_day = 350, max_import_rate = 0.2),
        water = dict(peak_day = 350, max_import_rate = 1.3),
    )
)

project_sim.update_pars(best_pars, recursive=True)

project_msim = zn.MultiSim(project_sim, n_runs = 250, verbose = 0)

if __name__ == "__main__":
    project_msim.run()                    # Run the simulations.
    project_msim.save(f'msims/{project_name}_peak_at_350_w_data.msim')  # Save the multi-simulation object.
    project_msim.combine()                # Combine the results from all simulations.
    project_msim.summarize()              # Summarize the combined results.
    project_msim.plot()                   # Plot the results.
    project_msim.plot_result(key = 'cum_poultry_flock_infectious')                # Plot detailed results.
    #project_msim.plot_result(key = 'cum_duck_flock_infectious')            # Plot detailed results.
    #project_msim.plot_result(key = 'cum_layer_flock_infectious')                # Plot detailed results.
    #project_msim.plot_result(key = 'cum_broiler_flock_infectious')                # Plot detailed results.