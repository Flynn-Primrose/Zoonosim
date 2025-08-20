import Zoonosim as zn
import dill

sim = zn.Sim()
sim.initialize()
# Get the class of your object
meta_cls = type(sim.agents.meta)

# for name, value in meta_cls.__dict__.items():
#     try:
#         dill.loads(dill.dumps(value))
#     except Exception as e:
#         print(f"❌ Class attribute {name} failed: {type(value)} -> {e}")
#     else:
#         print(f"✅ Class attribute {name} OK")


print(meta_cls)
print(meta_cls.__module__)
print(meta_cls.__qualname__)