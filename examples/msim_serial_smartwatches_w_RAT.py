import src as zn
import warnings

warnings.filterwarnings("ignore", category=UserWarning)  # Ignore some warnings for cleaner output.

rat = zn.RAT_disc(di_criteria = ['sx'], di_seekprobs = {'sx':0.0}, sc_criteria = ['alerted'], sc_seekprobs = {'alerted': 0.1},)  # Configure the RAT test.

testobj = [rat]  # Create a list of test objects, in this case just the RAT test.
new_pars = {'enable_smartwatches': True ,'enable_testobjs': True, 'testing' : testobj, 'rand_seed': 13}

sim = zn.Sim()  # Initialize the simulation object.
sim.update_pars(new_pars)  # Update the simulation parameters with our new parameters.

msim = zn.MultiSim(sim, n_runs = 10)  # Wrap the simulation in a MultiSim object.

if __name__ == "__main__":
    msim.run(verbose = 0.1, parallel = False) 
    msim.save('msims/smartwatches_w_RAT_serial.msim')
    msim.combine()  # Combine the results from all runs.
    msim.plot('testing')  # Plot the results of the testing.
    msim.plot('smartwatch')  # Plot the smartwatch data.
    msim.plot()  # Plot the overall simulation results.