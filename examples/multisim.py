# An example using the multisim functionality in Zoonosim

import Zoonosim as zn

sim = zn.Sim()  # Initialize the simulation object.
msim = zn.MultiSim(sim, n_runs = 5)  # Wrap the simulation in a MultiSim object.

if __name__ == "__main__":
    # Run the multi-simulation
    msim.run()


#print(msim.sims)
# Combine results from all runs
#msim.combine()

    # Plot the results
#msim.plot()