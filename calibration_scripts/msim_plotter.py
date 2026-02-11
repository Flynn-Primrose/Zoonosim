import zoonosim as zn

qc_200 = zn.MultiSim.load("./msims/Calibration_single_breed_2_peak_at_200_w_data.msim")
qc_250 = zn.MultiSim.load("./msims/Calibration_single_breed_2_peak_at_250_w_data.msim")
qc_300 = zn.MultiSim.load("./msims/Calibration_single_breed_2_peak_at_300_w_data.msim")
qc_350 = zn.MultiSim.load("./msims/Calibration_single_breed_2_peak_at_350_w_data.msim")

on_200 = zn.MultiSim.load("./msims/Calibration_single_breed_ON_peak_at_200_w_data.msim")
on_250 = zn.MultiSim.load("./msims/Calibration_single_breed_ON_peak_at_250_w_data.msim")
on_300 = zn.MultiSim.load("./msims/Calibration_single_breed_ON_peak_at_300_w_data.msim")
on_350 = zn.MultiSim.load("./msims/Calibration_single_breed_ON_peak_at_350_w_data.msim")

qc_200.reduce()
qc_250.reduce()
qc_300.reduce()
qc_350.reduce()

# qc_all = zn.MultiSim.merge([qc_200, qc_250, qc_300, qc_350])

on_200.reduce()
on_250.reduce()
on_300.reduce()
on_350.reduce()

# on_all = zn.MultiSim.merge([on_200, on_250, on_300, on_350])


# qc_200.plot('seasonality')
# qc_250.plot('seasonality')
# qc_300.plot('seasonality')
# qc_350.plot('seasonality')
# qc_all.plot('seasonality')

# on_200.plot('seasonality')
# on_250.plot('seasonality')
# on_300.plot('seasonality')
# on_350.plot('seasonality')
# on_all.plot('seasonality')





qc_200.plot_result(key = 'cum_poultry_flock_infectious')
qc_250.plot_result(key = 'cum_poultry_flock_infectious')
qc_300.plot_result(key = 'cum_poultry_flock_infectious')
qc_350.plot_result(key = 'cum_poultry_flock_infectious')
# qc_all.plot_result(key = 'cum_poultry_flock_infectious')

on_200.plot_result(key = 'cum_poultry_flock_infectious')
on_250.plot_result(key = 'cum_poultry_flock_infectious')
on_300.plot_result(key = 'cum_poultry_flock_infectious')
on_350.plot_result(key = 'cum_poultry_flock_infectious')
# on_all.plot_result(key = 'cum_poultry_flock_infectious')


qc_200.plot('flock')
qc_250.plot('flock')
qc_300.plot('flock')
qc_350.plot('flock')
# qc_all.plot('flock')

on_200.plot('flock')
on_250.plot('flock')
on_300.plot('flock')
on_350.plot('flock')
# on_all.plot('flock')