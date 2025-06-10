# This code sets up a simulation with a RAT and PCR test, both of which are configured to be screening tests.

import Zoonosim as zn




rat = zn.RAT_disc(sc_criteria = ['work'], sc_seekprobs = {'work': 0.5},)
pcr = zn.PCR_disc(sc_criteria = ['work'], sc_seekprobs = {'work': 0.5}, RAT_ind = 0)

testobj = [rat, pcr]

new_pars = {'enable_testobjs': True, 'testing' : testobj}

sim=zn.Sim()
sim.update_pars(new_pars)
sim.initialize()
sim.run()

sim.plot()
sim.plot('human')
sim.plot('flock')
sim.plot('barn')
sim.plot('water')