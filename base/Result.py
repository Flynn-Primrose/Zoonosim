import sciris as sc
import numpy as np

from ... import defaults as znd 

__all__ = ['Result']

class Result(object):
    '''
    Stores a single result -- by default, acts like an array.

    Args:
        name (str): name of this result, e.g. new_infections
        npts (int): if values is None, precreate it to be of this length
        scale (bool): whether or not the value scales by population scale factor
        color (str/arr): default color for plotting (hex or RGB notation)
        n_variants (int): the number of variants the result is for (0 for results not by variant)

    **Example**::

        import covasim as cv
        r1 = cv.Result(name='test1', npts=10)
        r1[:5] = 20
        print(r1.values)
    '''

    def __init__(self, name=None, npts=None, scale=True, color=None, n_variants=0):
        self.name =  name  # Name of this result
        self.scale = scale # Whether or not to scale the result by the scale factor
        if color is None:
            color = znd.get_default_colors()['default'] 
        self.color = color # Default color
        if npts is None:
            npts = 0
        npts = int(npts)

        if n_variants>0:
            self.values = np.zeros((n_variants, npts), dtype=znd.result_float)
        else:
            self.values = np.zeros(npts, dtype=znd.result_float)

        self.low  = None
        self.high = None
        return


    def __repr__(self):
        ''' Use pretty repr, like sc.prettyobj, but displaying full values '''
        output  = sc.prepr(self, skip=['values', 'low', 'high'], use_repr=False)
        output += 'values:\n' + repr(self.values)
        if self.low is not None:
            output += '\nlow:\n' + repr(self.low)
        if self.high is not None:
            output += '\nhigh:\n' + repr(self.high)
        return output


    def __getitem__(self, key):
        ''' To allow e.g. result['high'] instead of result.high, and result[5] instead of result.values[5] '''
        if isinstance(key, str):
            output = getattr(self, key)
        else:
            output = self.values.__getitem__(key)
        return output


    def __setitem__(self, key, value):
        ''' To allow e.g. result[:] = 1 instead of result.values[:] = 1 '''
        if isinstance(key, str):
            setattr(self, key, value)
        else:
            self.values.__setitem__(key, value)
        return

    def __eq__(self, other):
        #if we are comparing to a Result object
        if isinstance(other, Result):
            return np.array_equal(self.values, other.values)

        #if we are comparing to a regular array
        else: 
            if len(self.values) != len(other):
                return False
            for i in range(len(self.values)):
                if(self.values[i] != other[i]):
                    return False
            return True

    def __len__(self):
        ''' To allow len(result) instead of len(result.values) '''
        return len(self.values)


    @property
    def npts(self):
        return len(self.values)