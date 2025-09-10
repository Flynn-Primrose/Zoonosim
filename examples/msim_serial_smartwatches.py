import Zoonosim as zn
import warnings

warnings.filterwarnings("ignore", category=UserWarning)  # Ignore some warnings for cleaner output.

sim = zn.Sim()  # Initialize the simulation object.
new_pars = {'enable_smartwatches': True, 'rand_seed': 21}
sim.update_pars(new_pars)  # Update parameters if needed.

msim = zn.MultiSim(sim, n_runs = 10)  # Wrap the simulation in a MultiSim object.

if __name__ == "__main__":
    msim.run(verbose = 0.1, parallel = False) 
    msim.save('msims/smartwatches_serial.msim')
    msim.combine()  # Combine the results from all runs.
    msim.plot()  # Plot the combined results.
    msim.plot('smartwatch')  # Plot smartwatch-specific results.