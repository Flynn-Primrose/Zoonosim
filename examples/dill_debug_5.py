import Zoonosim as zn
import pathosim as inf

import dill

zn_sim = zn.Sim()
zn_sim.initialize()

zn_cls = zn_sim.agents

pop_size = 20000 
  
#Create the circulating pathogen
MyNewPathogen = inf.SARS_COV_2(50)   

#Configure the immunity of that pathogen. 4 parameters have to be passed in this function. See Pathosim Documentation.
MyNewPathogen.configure_generalized_immunity(0.35, 0.95, 180, 14)  #Here, we pass parameters such that the spread of MyNewPathogen ressembles COVID-19's spread.

#Initialize the simulation
inf_sim = inf.Sim(pop_size = pop_size, pathogens = [MyNewPathogen], n_days = 365)
inf_sim.initialize()

inf_cls = inf_sim.people

def test_class(cls):
    for name, value in cls.__dict__.items():
        try:
            # dill.loads(dill.dumps(value))
            dill.loads(dill.dumps(value))
        except Exception as e:
            print(f"❌ Class attribute {name} failed: {type(value)} -> {e}")
        else:
            print(f"✅ Class attribute {name} OK: {type(value)}")

print("Testing Zoonosim class attributes:")
print("--------------------------------")
test_class(zn_cls)
print("\nTesting Pathosim class attributes:")
print("--------------------------------")
test_class(inf_cls)