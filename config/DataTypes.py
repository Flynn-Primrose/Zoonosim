import numpy as np
import numba as nb

from .Options import options as zno

__all__ = ['default_float', 'default_int', 'nbfloat', 'nbint', 'result_float']

result_float = np.float64 # Always use float64 for results, for simplicity
if zno.precision == 32:
    default_float = np.float32
    default_int   = np.int32
    nbfloat       = nb.float32
    nbint         = nb.int32
elif zno.precision == 64: # pragma: no cover
    default_float = np.float64
    default_int   = np.int64
    nbfloat       = nb.float64
    nbint         = nb.int64
else:
    raise NotImplementedError(f'Precision must be either 32 bit or 64 bit, not {zno.precision}')