# This example shows how to run a simulation with the default parameters.

import zoonosim as zn
sim=zn.Sim()
sim.initialize()
#print(sim.agents.flock['quarantined'])
#sim.step()
#print(sim.agents.flock['quarantined'])
#print(sim.agents.flock['suspected'])
#print(sim.agents.flock['_pending_quarantine'][sim.agents.flock.t])


sim.run()

sim.plot()
# sim.plot('human')
sim.plot('flock')
# sim.plot('barn')
# sim.plot('water')