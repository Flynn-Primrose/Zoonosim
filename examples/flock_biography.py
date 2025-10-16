import zoonosim as zn

bio_analyzer = zn.biography(uid = 266, agent_type = 'flock', days = range(1, 150))

sim=zn.Sim(analyzers = bio_analyzer)

new_pars = { 'rand_seed': 24}

sim.update_pars(new_pars)


sim.run(reset_seed=False)

# print(sim.agents.flock.event_log)

flock_biography = sim['analyzers'][0]

headline_plots = [
    'headcount',
    'total_dead_headcount',
]

SEIR_plots = [
    'exposed_headcount',
    'infectious_headcount',
    'symptomatic_headcount'
]

flock_biography.plot(
    headline_plots,
    do_show=True,
)

flock_biography.plot(
    SEIR_plots,
    do_show=True,
)