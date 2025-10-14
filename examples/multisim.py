# An example using the multisim functionality in Zoonosim

import src as zn
import sciris as sc
#import warnings

import copy #for debugging
import traceback
import gc #for debugging

#warnings.filterwarnings("ignore", category=UserWarning)  # Ignore some warnings for cleaner output.

sim = zn.Sim()  # Initialize the simulation object.
new_pars = {'rand_seed': 42}
sim.update_pars(new_pars)  # Update parameters if needed.


msim = zn.MultiSim(sim, n_runs = 5)  # Wrap the simulation in a MultiSim object.

if __name__ == "__main__":
    # Run the multi-simulation
    #msim.run(keep_people = True, par_args=dict(parallelizer = 'serial'))  # Use n_cpus > 1 for parallel execution.
    msim.run(verbose = 2, parallel = True)  # Use n_cpus > 1 for parallel execution.


    # msim.combine()  # Combine the results from all simulations.

    # msim.summarize()  # Summarize the combined results.


    # msim.plot()  # Plot the results.


