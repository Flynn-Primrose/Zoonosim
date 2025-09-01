import Zoonosim as zn
import dill

sim = zn.Sim()  # Initialize the simulation object.
new_pars = {'rand_seed': 42}
sim.update_pars(new_pars)  # Update parameters if needed.
sim.initialize()  # Initialize the simulation.

run_sim = zn.single_run(sim=sim, keep_people=False)

dump = dill.dumps(run_sim)
print(len(dump))

# for name, value in sim.agents.__dict__.items():
#     try:
#         dill.loads(dill.dumps(value))
#     except Exception as e:
#         print(f"❌ Class attribute {name} failed: {type(value)} -> {e}")
#     else:
#         print(f"✅ Class attribute {name} OK")