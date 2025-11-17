import zoonosim as zn

# import optuna as op

# storage = 'sqlite:///zoonosim_calibration.db'

# study = op.load_study(storage=storage, study_name='zoonosim_calibration')

# study_dataframe = study.trials_dataframe()
# study_dataframe.to_csv('calibration_study_results.csv', index=False)

new_pars = dict(
    start_day = '2022-01-01',
    end_day = '2025-12-31'
)

sim = zn.Sim(datafile="zoonosim/data/H5N1_cases_in_QC_poultry.csv", pars=new_pars)
calib_pars = dict(
    beta = dict(
        human = [0.5, 0.0, 2.0],
        flock = [0.5, 0.0, 2.0],
        barn = [0.5, 0.0, 2.0],
        water = [0.5, 0.0, 2.0],
    ),
    n_imports = dict(
        barn = [0.1, 0, 1],
        water = [0.1, 0, 1],
    )
)

calib = zn.Calibration(sim, calib_pars, n_trials=5, die=True)

if __name__ == "__main__":
    calib.calibrate()
    calib.plot_trend()