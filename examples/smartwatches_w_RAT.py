# Example of using RAT in combination with smartwatch alerts

import zoonosim as zn

rat = zn.RAT_disc(di_criteria = ['sx'], di_seekprobs = {'sx':0.0}, sc_criteria = ['alerted'], sc_seekprobs = {'alerted': 0.1},)  # Configure the RAT test.

testobj = [rat]  # Create a list of test objects, in this case just the RAT test.
new_pars = {'enable_smartwatches': True ,'enable_testobjs': True, 'testing' : testobj}

sim = zn.Sim()  # Initialize the simulation object.
sim.update_pars(new_pars)  # Update the simulation parameters with our new parameters.
sim.initialize()  # Initialize the simulation with the updated parameters.
sim.run()  # Run the simulation.

sim.plot('testing')  # Plot the results of the testing.
sim.plot('smartwatch')  # Plot the smartwatch data.
sim.plot()  # Plot the overall simulation results.
