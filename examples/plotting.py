import zoonosim as zn

sim=zn.Sim()
sim.run()

sim.plot()
sim.plot('human')
sim.plot('flock')
sim.plot('breed')
sim.plot('water')
sim.plot('barn')

sim.plot_result('n_human_alive')
sim.plot_result('n_human_naive')
sim.plot_result('n_human_removed')
sim.plot_result('human_prevalence')
sim.plot_result('human_incidence')