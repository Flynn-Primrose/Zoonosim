import zoonosim as zn
import json
import numpy as np

saved_pars_filename = "saved_pars/default_single_breed_poultry_and_cattle.json"
sim_label = "default single breed poultry and cattle"

saved_pars = zn.pars_from_json(saved_pars_filename)

saved_sim = zn.Sim(pars = saved_pars, label = sim_label)
# new_pars = dict(
#     rand_seed = 21,
#     beta = dict(
#         human = saved_pars['beta']['flock']
#     )
# )
# saved_sim.update_pars(new_pars, recursive=True)


msim = zn.MultiSim(saved_sim, label=sim_label, n_runs=100, verbose=0.3)  # Wrap the simulation in a MultiSim object.

if __name__ == "__main__":
    msim.run(keep_people = True, run_args=dict(auto_finalize=False))                    # Run the simulations.
    msim.finalize()                # Finalize the simulations (necessary if auto_finalize=False).
    msim.shrink()                   # Shrink the simulations. This is necessary because we set keep_people to true so we could finalize the sims after running, but means we have to shrink the sims before we can combine them.
    msim.save(f'msims/{sim_label}.msim')  # Save the multi-simulation object.
    msim.reduce(use_mean=True)                # Reduce the results from all simulations.
    msim.summarize()              # Summarize the Reduced results.
    msim.plot()                   # Plot the results.
    msim.plot('flock')
    msim.plot('herd')
    msim.plot('human')
    msim.plot('ppe')

