
import dill
import io

def dill_pickle_info(obj):
    """
    Inspect how dill pickles an object.
    
    Returns:
        dict with:
            - module:  obj.__module__
            - qualname: obj.__qualname__ (if available)
            - by_value: True if dill pickled inline, False if by reference
            - snippet:  first 200 bytes of dill dump for inspection
    """
    info = {}
    
    # Get basic module + qualname if available
    info["module"] = getattr(obj, "__module__", None)
    info["qualname"] = getattr(obj, "__qualname__", None)
    
    # Serialize with dill
    data = dill.dumps(obj)
    snippet = data[:200]
    info["snippet"] = snippet
    
    # Heuristic: if dill stored a dotted path, it's by reference.
    # If you see FUNCTION/CLASS code objects, it's by value.
    info["by_value"] = b'code' in snippet or b'c__builtin__' in snippet
    
    return info


# Example usage:
if __name__ == "__main__":
    import Zoonosim as zn

    sim =zn.Sim()
    sim.initialize()
    obj = sim.agents.human.meta
    
    print(dill_pickle_info(obj))
    
