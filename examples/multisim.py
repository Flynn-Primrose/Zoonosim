# An example using the multisim functionality in Zoonosim

import Zoonosim as zn

sim = zn.Sim()  # Initialize the simulation object.
new_pars = {'rand_seed': 42}
sim.update_pars(new_pars)  # Update parameters if needed.
sim.initialize()  # Initialize the simulation.

msim = zn.MultiSim(sim, n_runs = 3)  # Wrap the simulation in a MultiSim object.

if __name__ == "__main__":
    # Run the multi-simulation
    msim.run()


#print(msim.sims)
# Combine results from all runs
#msim.combine()

    # Plot the results
#msim.plot()