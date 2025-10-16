import zoonosim as zn
import warnings

warnings.filterwarnings("ignore", category=UserWarning)  # Ignore some warnings for cleaner output.

rat = zn.RAT_disc(di_criteria = ['sx'], di_seekprobs = {'sx':0.0}, sc_criteria = ['work'], sc_seekprobs = {'work': 0.5},) # Configure the RAT test. 
#It appears that there must be arguments for diagnostic testing, here we have used sx (symptomatic) as a criterion with a seek probability of 0.0, in effect this means that the RAT test is not used for diagnostic purposes.
# For Screening, we have set the criteria to 'work' with sc_seekprobs set to 0.5, this means that half of work places will screen workers using RAT tests.

testobj = [rat]  # Create a list of test objects, in this case just the RAT test.
new_pars = {'enable_testobjs': True, 'testing' : testobj, 'rand_seed' : 37}  # create a dictionary of parameters to be modified from the default values, here we enable the test objects and set the testing parameter to our testobj list.

sim = zn.Sim()  # Initialize the simulation object.
sim.update_pars(new_pars)  # Update the simulation parameters with our new parameters.

msim = zn.MultiSim(sim, n_runs = 10)  # Wrap the simulation in a MultiSim object.

if __name__ == "__main__":
    msim.run(verbose = 0.1, parallel = False) 
    msim.save('msims/RAT_only_screening_serial.msim')
    msim.combine()  # Combine the results from all runs.
    msim.plot()  # Plot the combined results.
    msim.plot('testing')
