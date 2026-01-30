import optuna as op

name = "Calibration_single_breed_2"
db_name = f'./studies/{name}.db'
study = op.load_study(study_name=name, storage=op.storages.JournalStorage(op.storages.journal.JournalFileBackend(db_name, op.storages.journal.JournalFileOpenLock(db_name))))
op.visualization.plot_optimization_history(study).show()
print(study.best_params)