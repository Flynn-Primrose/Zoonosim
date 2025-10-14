import src as zn

bio_analyzer = zn.biography(uid = 44, agent_type = 'human', days = range(1, 150))
sim=zn.Sim(analyzers = bio_analyzer)

new_pars = { 'rand_seed': 24}

sim.update_pars(new_pars)

sim.run()

# print(sim.agents.human.infection_log)

human_biography = sim['analyzers'][0]

props_to_plot = [
    'viral_load',
]

human_biography.plot(
    props_to_plot,
    do_show=True,
)