import Zoonosim as zn
import pathosim as inf

import dill

zn_sim = zn.Sim()
zn_sim.initialize()


pop_size = 20000 
  
#Create the circulating pathogen
MyNewPathogen = inf.SARS_COV_2(50)   

#Configure the immunity of that pathogen. 4 parameters have to be passed in this function. See Pathosim Documentation.
MyNewPathogen.configure_generalized_immunity(0.35, 0.95, 180, 14)  #Here, we pass parameters such that the spread of MyNewPathogen ressembles COVID-19's spread.

#Initialize the simulation
inf_sim = inf.Sim(pop_size = pop_size, pathogens = [MyNewPathogen], n_days = 365)
inf_sim.initialize()

# Original meta class (from unmodified pathosim)
meta_orig_cls = inf_sim.people.meta.__class__


# Your modified meta class
meta_mod_cls = zn_sim.agents.meta.__class__

print("inf_sim.people.meta:", type(inf_sim.people.meta))
print("zn_sim.agents.meta:", type(zn_sim.agents.meta))

# print("=== Class info ===")
# print("Original:", meta_orig_cls, "Module:", meta_orig_cls.__module__, "Qualname:", meta_orig_cls.__qualname__)
# print("Modified:", meta_mod_cls, "Module:", meta_mod_cls.__module__, "Qualname:", meta_mod_cls.__qualname__)
# print("Base classes original:", meta_orig_cls.__bases__)
# print("Base classes modified:", meta_mod_cls.__bases__)

# # Compare class-level attributes
# print("\n=== Class-level attributes added ===")
# for name in meta_mod_cls.__dict__:
#     if name not in meta_orig_cls.__dict__:
#         print("Added/changed attribute:", name)

# print("\n=== Class-level attributes missing ===")
# for name in meta_orig_cls.__dict__:
#     if name not in meta_mod_cls.__dict__:
#         print("missing/changed attribute:", name)

# # Compare __init__ specifically
# init_orig = meta_orig_cls.__init__
# init_mod  = meta_mod_cls.__init__

# print("\n=== __init__ comparison ===")
# print("Original __init__ module:", init_orig.__module__, "id:", id(init_orig))
# print("Modified __init__ module:", init_mod.__module__, "id:", id(init_mod))

# if init_orig is not init_mod:
#     print("❌ __init__ differs — likely cause of dill pickling failure")

# # Optional: test if dill would try to pickle it
# try:
#     mod_dump = dill.dumps(init_mod)
# except Exception as e:
#     print("❌ Dill cannot pickle modified __init__:", e)
# else:
#     print("✅ Dill can pickle modified __init__")
#     try:
#         mod_load = dill.loads(mod_dump)
#     except Exception as e:
#         print("❌ Dill cannot unpickle modified __init__:", e)
#     else:
#         print("✅ Dill can unpickle modified __init__")

# # Optional: test if dill would try to pickle it
# try:
#     orig_dump = dill.dumps(init_orig)
# except Exception as e:
#     print("❌ Dill cannot pickle original __init__:", e)
# else:
#     print("✅ Dill can pickle original__init__")
#     try:
#         orig_load = dill.loads(orig_dump)
#     except Exception as e:
#         print("❌ Dill cannot unpickle original __init__:", e)
#     else:
#         print("✅ Dill can unpickle original __init__")

import inspect
# print(inspect.getsource(meta_orig_cls.__init__))
# print(inspect.getsource(meta_mod_cls.__init__))

import difflib
import textwrap

def diff_inits(cls1, cls2, name1="Class1", name2="Class2"):
    # Get source code (dedent to align nicely)
    src1 = textwrap.dedent(inspect.getsource(cls1.__init__)).splitlines()
    src2 = textwrap.dedent(inspect.getsource(cls2.__init__)).splitlines()

    # Run unified diff
    diff = difflib.unified_diff(
        src1, src2,
        fromfile=f"{name1}.__init__",
        tofile=f"{name2}.__init__",
        lineterm=""
    )

    # Print result
    print("\n".join(diff))


# Usage with your classes
# diff_inits(meta_orig_cls, meta_mod_cls, "PeopleMeta", "AgentsMeta")

# print(inspect.getsource(meta_mod_cls))   # show whole class definition
# print(meta_mod_cls.__init__)
# print(type(meta_mod_cls.__init__))

# print(meta_orig_cls.__init__)
# print(type(meta_orig_cls.__init__))

# print(meta_orig_cls.__init__.__qualname__)
# print(meta_mod_cls.__init__.__qualname__)
# print(meta_orig_cls.__init__.__module__)
# print(meta_mod_cls.__init__.__module__)


# print(meta_mod_cls.__init__.__module__)    # where Python thinks it lives
# print(inspect.getsourcefile(meta_mod_cls.__init__))  # which file it comes from
# print(meta_mod_cls.__module__)             # where the class lives

# print("PeopleMeta.__init__:")
# print("  __module__:", PeopleMeta.__init__.__module__)
# print("  __qualname__:", PeopleMeta.__init__.__qualname__)
# print()

# print("AgentsMeta.__init__:")
# print("  __module__:", AgentsMeta.__init__.__module__)
# print("  __qualname__:", AgentsMeta.__init__.__qualname__)
