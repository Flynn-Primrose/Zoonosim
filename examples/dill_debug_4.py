import Zoonosim as zn
#import pathosim as inf

import dill
#dill.detect.trace(True)
dill.settings['recurse'] = False
# import faulthandler
# faulthandler.enable()

zn_sim = zn.Sim()
zn_sim.initialize()



# pop_size = 20000 
  
# #Create the circulating pathogen
# MyNewPathogen = inf.SARS_COV_2(50)   

# #Configure the immunity of that pathogen. 4 parameters have to be passed in this function. See Pathosim Documentation.
# MyNewPathogen.configure_generalized_immunity(0.35, 0.95, 180, 14)  #Here, we pass parameters such that the spread of MyNewPathogen ressembles COVID-19's spread.

# #Initialize the simulation
# inf_sim = inf.Sim(pop_size = pop_size, pathogens = [MyNewPathogen], n_days = 365)
# inf_sim.initialize()




# inf_meta_cls = inf_sim.people.meta
zn_cls = zn_sim.agents.human

for name, value in zn_cls.__dict__.items():
    try:
        # dill.loads(dill.dumps(value))
        dill.dumps(value)
    except Exception as e:
        print(f"❌ Class attribute {name} failed: {type(value)} -> {e}")
    else:
        print(f"✅ Class attribute {name} OK: {type(value)}")

# for name, value in inf_meta_cls.__dict__.items():
#     try:
#         dill.loads(dill.dumps(value))
#     except Exception as e:
#         print(f"❌ Class attribute {name} failed: {type(value)} -> {e}")
#     else:
#         print(f"✅ Class attribute {name} OK: {type(value)}")