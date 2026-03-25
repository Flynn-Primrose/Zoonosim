import optuna as op

name1 = "Calibration_single_breed_ON"
db_name1 = f'./studies/{name1}.db'
study = op.load_study(study_name=name1, storage=op.storages.JournalStorage(op.storages.journal.JournalFileBackend(db_name1, op.storages.journal.JournalFileOpenLock(db_name1))))

name2 = "Calibration_single_breed"
db_name2 = f'./studies/{name2}.db'
study2 = op.load_study(study_name=name2, storage=op.storages.JournalStorage(op.storages.journal.JournalFileBackend(db_name2, op.storages.journal.JournalFileOpenLock(db_name2))))

name3 = "Calibration_single_breed_2"
db_name3 = f'./studies/{name3}.db'
study3 = op.load_study(study_name=name3, storage=op.storages.JournalStorage(op.storages.journal.JournalFileBackend(db_name3, op.storages.journal.JournalFileOpenLock(db_name3))))


op.visualization.plot_edf([study, study2, study3]).show()

#op.visualization.plot_optimization_history(study).show()
#op.visualization.plot_param_importances(study).show()
#op.visualization.plot_contour(study).show()
#print(study.best_params)