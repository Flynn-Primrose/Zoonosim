import zoonosim as zn
import pandas as pd
import numpy as np

AOQ_4th_iteration = zn.MultiSim.load("./msims/AOQ_4th_iteration_RHB.msim")


AOQ_4th_iteration.reduce(use_mean=True)



# data = pd.read_csv("zoonosim/data/CFIA_monthly_incidence.csv")
# data['date'] = pd.to_datetime(data['date'])
# data = data.rename(columns={"AOQ_adjusted_50":"monthly_new_poultry_flock_infectious"})
# data['monthly_new_human_infectious'] = np.repeat(0, len(data))
# data = data[["date", "monthly_new_human_infectious", "monthly_new_poultry_flock_infectious"]]

# AOQ_4th_iteration.base_sim.load_data(data)

AOQ_4th_iteration.plot()
AOQ_4th_iteration.plot_result(key = 'monthly_new_human_infectious')
AOQ_4th_iteration.plot_result(key = "monthly_new_poultry_flock_infectious")
AOQ_4th_iteration.plot('human')
AOQ_4th_iteration.plot_transmission_vectors('human')
AOQ_4th_iteration.plot('ppe')
AOQ_4th_iteration.plot_transmission_vectors('ppe')
AOQ_4th_iteration.plot('flock')
AOQ_4th_iteration.plot_transmission_vectors('flock')
AOQ_4th_iteration.plot('barn')
AOQ_4th_iteration.plot_transmission_vectors('barn')
