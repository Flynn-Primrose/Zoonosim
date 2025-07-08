# This code sets up a simulation with a RAT and PCR test, both of which are configured to be screening tests.

import Zoonosim as zn




rat = zn.RAT_disc(di_criteria = ['sx'], di_seekprobs = {'sx':0.0}, sc_criteria = ['work'], sc_seekprobs = {'work': 0.5},)
pcr = zn.PCR_disc(di_criteria = ['pos_RAT'], di_seekprobs = {'pos_RAT':1.0}, sc_criteria = ['work'], sc_seekprobs = {'work':0.0}, RAT_ind = 0)

testobj = [rat, pcr]

new_pars = {'enable_testobjs': True, 'testing' : testobj}

sim=zn.Sim()
sim.update_pars(new_pars)
sim.initialize()
sim.run()

sim.plot('testing')
sim.plot()
sim.plot('human')
sim.plot('flock')
sim.plot('barn')
sim.plot('water')