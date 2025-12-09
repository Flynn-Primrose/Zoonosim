import zoonosim as zn

msim = zn.MultiSim.load("./msims/calibration2_1_best2.msim")

msim.combine()
msim.summarize()
msim.plot()
msim.plot_result(key='cum_duck_flock_infectious')
msim.plot_result(key='cum_layer_flock_infectious')
msim.plot_result(key='cum_broiler_flock_infectious')