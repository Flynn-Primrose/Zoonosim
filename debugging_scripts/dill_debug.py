import dill
import Zoonosim as zn

sim = zn.Sim()
sim.initialize()

dump = dill.dumps(type(sim.agents.meta))
load = dill.loads(dump)

