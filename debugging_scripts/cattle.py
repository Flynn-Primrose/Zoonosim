import zoonosim as zn
import json
import numpy as np

saved_pars_filename = "saved_pars/default_single_breed_poultry_and_cattle.json"
sim_label = "default single breed poultry and cattle"

saved_pars = zn.pars_from_json(saved_pars_filename)

saved_sim = zn.Sim(pars = saved_pars, label = sim_label)

saved_sim.run()                    # Run the simulations.
saved_sim.plot()                   # Plot the results.
saved_sim.plot('flock')
saved_sim.plot('herd')
saved_sim.plot('human')