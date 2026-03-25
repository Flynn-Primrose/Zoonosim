import zoonosim as zn
import json
import numpy as np

saved_pars_filename = "saved_pars/single_breed.json"

with open(saved_pars_filename, 'r') as file:
    saved_pars = json.load(file)

for key, value in saved_pars['prognoses']['human'].items():
    saved_pars['prognoses']['human'][key] = np.array(value)
for key, value in saved_pars['prognoses']['flock'].items():
    if key == 'breed':
        saved_pars['prognoses']['flock'][key] = np.array(value, dtype=zn.default_str)
    else:
        saved_pars['prognoses']['flock'][key] = np.array(value)
for key, value in saved_pars['poultry_pars'].items():
    if key == 'breeds':
        saved_pars['poultry_pars'][key] = np.array(value, dtype=zn.default_str)

saved_sim = zn.Sim(pars = saved_pars)

saved_sim.run()                    # Run the simulation.
saved_sim.plot()                   # Plot the results.