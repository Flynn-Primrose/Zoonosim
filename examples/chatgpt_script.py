import sys
import importlib.util
import multiprocessing as mp
import dill

# ---------- CONFIGURE THIS ----------
PACKAGE_NAME = "my_package"   # top-level package
SUBMODULES = ["A", "B"]       # submodules containing classes to test
CLASS_NAMES = ["classA", "classB"]
# -----------------------------------

def check_import(module_name):
    spec = importlib.util.find_spec(module_name)
    if spec is None:
        return f"FAIL: cannot find {module_name}"
    try:
        module = importlib.import_module(module_name)
        return f"OK: {module_name} imported ({module})"
    except Exception as e:
        return f"FAIL: error importing {module_name}: {e}"

def test_dill_pickle(module_name, class_name):
    """Try to pickle and unpickle the class by reference using dill"""
    try:
        module = importlib.import_module(module_name)
        cls = getattr(module, class_name)
        # dill.dumps should pickle by reference if importable
        pickled = dill.dumps(cls)
        unpickled = dill.loads(pickled)
        # check if unpickled class is the same object as original
        same_reference = cls is unpickled
        return f"{class_name} pickled by reference? {same_reference}"
    except Exception as e:
        return f"{class_name} pickling failed: {e}"

def worker_test(_):
    results = []
    # Check top-level package
    results.append(check_import(PACKAGE_NAME))
    # Check submodules
    for submod, cls_name in zip(SUBMODULES, CLASS_NAMES):
        full_mod = f"{PACKAGE_NAME}.{submod}"
        results.append(check_import(full_mod))
        results.append(test_dill_pickle(full_mod, cls_name))
    return results

if __name__ == "__main__":
    print("=== Checking imports and dill pickling in main process ===")
    for r in worker_test(None):
        print(r)

    print("\n=== Checking in a separate worker process ===")
    with mp.Pool(1) as pool:
        worker_results = pool.map(worker_test, [0])
    for r in worker_results[0]:
        print(r)
