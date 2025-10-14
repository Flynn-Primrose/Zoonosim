# Example of screening using RAT only
# This code sets up a simulation with a RAT test configured as a screening test.

import Zoonosim as zn


rat = zn.RAT_disc(di_criteria = ['sx'], di_seekprobs = {'sx':0.0}, sc_criteria = ['work'], sc_seekprobs = {'work': 0.5},) # Configure the RAT test. 
#It appears that there must be arguments for diagnostic testing, here we have used sx (symptomatic) as a criterion with a seek probability of 0.0, in effect this means that the RAT test is not used for diagnostic purposes.
# For Screening, we have set the criteria to 'work' with sc_seekprobs set to 0.5, this means that half of work places will screen workers using RAT tests.

testobj = [rat]  # Create a list of test objects, in this case just the RAT test.
new_pars = {'enable_testobjs': True, 'testing' : testobj}  # create a dictionary of parameters to be modified from the default values, here we enable the test objects and set the testing parameter to our testobj list.

sim = zn.Sim()  # Initialize the simulation object.
sim.update_pars(new_pars)  # Update the simulation parameters with our new parameters.
sim.initialize()  # Initialize the simulation with the updated parameters.
sim.run()  # Run the simulation.

sim.plot('testing')  # Plot the results of the testing.
sim.plot()  # Plot the overall simulation results.
