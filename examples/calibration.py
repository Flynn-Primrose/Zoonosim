import zoonosim as zn

sim = zn.Sim(datafile="zoonosim/data/H5N1_cases_in_QC_poultry.csv")
calib_pars = dict(
    beta = dict(
        human = [0.5, 0.0, 2.0],
        flock = [0.5, 0.0, 2.0],
        barn = [0.5, 0.0, 2.0],
        water = [0.5, 0.0, 2.0],
    )
)

calib = zn.Calibration(sim, calib_pars, total_trials=10, die=True)

if __name__ == "__main__":
    calib.calibrate()
#calib.plot()