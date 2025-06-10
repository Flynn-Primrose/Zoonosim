# This example shows how to run a simulation with the default parameters.

import Zoonosim as zn
sim=zn.Sim()
sim.initialize()
sim.run()

sim.plot()
sim.plot('human')
sim.plot('flock')
sim.plot('barn')
sim.plot('water')