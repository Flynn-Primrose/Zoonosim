import optuna as op

project_name = "Calibration_single_breed_AOQ_4th_iteration"
filename = f'./studies/{project_name}.db'
storage   = op.storages.JournalStorage(op.storages.journal.JournalFileBackend(filename, op.storages.journal.JournalFileOpenLock(filename)))
study = op.load_study(study_name= project_name, storage = storage)
# print(study.best_params)
trails = study.trials_dataframe()
good_trials = trails[trails['value'] < 6.0]
good_avg_values = good_trials.filter(like='params_').mean().to_dict()
print(good_avg_values)

