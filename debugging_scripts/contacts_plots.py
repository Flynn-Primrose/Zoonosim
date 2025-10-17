import zoonosim as zn
import networkx as nx
import matplotlib.pyplot as plt

sim = zn.Sim()
sim.initialize()
contact_graph = sim.agents.contacts.to_graph()
nx.draw(contact_graph)
plt.show()