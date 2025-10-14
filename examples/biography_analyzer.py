import Zoonosim as zn

bio_analyzer = zn.biography(uid = None, agent_type = 'flock', days = range(1, 150))

sim=zn.Sim(analyzers = bio_analyzer)
sim.initialize()
sim.run()

biography = sim['analyzers'][0]

props_to_plot = [
    'headcount',
    'exposed_headcount',
    'infectious_headcount',
    'symptomatic_headcount',
    'total_dead_headcount',
    #'water_consumption',
]

biography.plot(
    props_to_plot,
    do_show=True,
)

#sim.plot('flock')