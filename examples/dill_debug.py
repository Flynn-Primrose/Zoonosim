import dill
from dill import detect
#dill.detect.trace(True)  # turns on tracing
import Zoonosim as zn

sim = zn.Sim()
sim.initialize()

# meta_class = type(sim.agents.meta)
# if hasattr(meta_class, "__reduce__"):
#     print("Reduce info:", meta_class.__reduce__)
# if hasattr(meta_class, "__getstate__"):
#     state = meta_class.__getstate__(sim.agents.meta)
#     print("State keys:", state.keys() if hasattr(state, 'keys') else state)

#import builtins

# meta_class = type(sim.agents.meta)

# original_setattr = meta_class.__setattr__

# def debug_setattr(self, name, value):
#     print(f"Assigning {name} = {value}")
#     original_setattr(self, name, value)

# meta_class.__setattr__ = debug_setattr

# try:
#     dill.loads(dill.dumps(sim.agents.meta))
# except Exception as e:
#     print("Error during unpickling:", e)
# finally:
#     meta_class.__setattr__ = original_setattr

#detect.errors(sim.agents.meta, depth = 3)

# for key, val in sim.agents.flock.__dict__.items():
#     try:
#         dill.dumps(val)
#         dill.loads(dill.dumps(val))
#     except Exception as e:
#         print(f"Problem with attribute {key!r}: {type(val)} -> {e}")

dump = dill.dumps(type(sim.agents.meta))
sim2 = dill.loads(dump)

# bad_attrs = []
# for k,v in sim.agents.__dict__.items():
#     try:
#         dill.loads(dill.dumps(v))
#     except Exception as e:
#         print(f"❌ Cannot unpickle sim.agents.{k}: {e}")
#         bad_attrs.append((k, type(v), e))
#     else:
#         print(f"✅ unPickle OK: sim.agents.{k}")

# attr = sim.agents  # or whatever failed
# for k,v in attr.__dict__.items():
#     try:
#         dill.dumps(v)
#     except Exception as e:
#         print(f"❌ Cannot pickle agents.{k}: {e}")

# dill.dumps(type(sim))

#dill.detect.trace(True)   # turns on tracing
#dill.dumps(sim)