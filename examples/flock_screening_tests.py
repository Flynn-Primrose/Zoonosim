# An example simulation using regular screening of flocks

import zoonosim as zn

flock_screening = zn.test_flock()

flock_screening_sim = zn.Sim(interventions = flock_screening)

flock_screening_sim.initialize()
flock_screening_sim.run()
flock_screening_sim.plot()