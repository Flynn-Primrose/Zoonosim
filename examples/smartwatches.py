import Zoonosim as zn

sim=zn.Sim()
new_pars = {'enable_smartwatches': True}
sim.update_pars(new_pars)
sim.initialize()
sim.run()
sim.plot()
sim.plot('smartwatch')