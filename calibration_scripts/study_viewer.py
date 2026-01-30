import optuna as op

project_name = "Calibration_single_breed_2"
filename = f'./studies/{project_name}.db'
storage   = op.storages.JournalStorage(op.storages.journal.JournalFileBackend(filename, op.storages.journal.JournalFileOpenLock(filename)))
study = op.load_study(study_name= project_name, storage = storage)
print(study.best_params)