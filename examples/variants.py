# An example showing the use af variants

import Zoonosim as zn

hpai = zn.variant(variant = 'HPAI', days = 0, n_imports = 2 , import_type = ['water', 'barn'])
lpai = zn.variant(variant = 'LPAI', days = 0, n_imports = 10, import_type = ['water', 'barn'])

new_pars = {}

new_pars['n_imports'] = {
    'human': 0,
    'flock': 0,
    'barn': 0,
    'water': 0,
}
new_pars['variants'] = [hpai, lpai]

sim=zn.Sim(pars=new_pars)
sim.initialize()
sim.run()
# Print the results
sim.plot()
sim.plot('human')
sim.plot('flock')
sim.plot('barn')
sim.plot('water')