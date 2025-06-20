import Zoonosim as zn

bio_analyzer = zn.biography(uid = None, agent_type = 'human', days = range(1, 100))

sim=zn.Sim(analyzers = bio_analyzer)
sim.initialize()
sim.run()

biography = sim['analyzers'][0]
print(biography.bio)