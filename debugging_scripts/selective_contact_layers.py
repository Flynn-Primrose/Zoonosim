import zoonosim as zn

sim = zn.Sim()
sim.initialize(skip_layers=['hh'])
print(sim.agents.contacts.keys())