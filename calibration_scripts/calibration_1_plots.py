import zoonosim as zn
import optuna as op

name = "Calibration8"
db_name = f'./studies/{name}.db'
study8 = op.load_study(study_name="Calibration8", storage=op.storages.JournalStorage(op.storages.journal.JournalFileBackend(db_name, op.storages.journal.JournalFileOpenLock(db_name))))
op.visualization.plot_optimization_history(study8).show()
print(study8.best_params)