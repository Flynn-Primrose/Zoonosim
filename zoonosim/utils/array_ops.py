'''
 Simple array utilities for running Zoonosim.


'''

import numba as nb # For faster computations
import numpy as np # For numerics

__all__ = ['true',   'false',   'defined',   'undefined',
            'itrue',  'ifalse',  'idefined',  'iundefined',
            'itruei', 'ifalsei', 'idefinedi', 'iundefinedi', 'custom_np_fmin']

@nb.njit()
def custom_np_fmin(arr1, arr2, size):
    for j in range(size):
        if not (np.isnan(arr1[j]) and np.isnan(arr2[j])):
            if np.isnan(arr2[j]):
                arr2[j] = arr1[j]
            elif not np.isnan(arr1[j]):
                arr2[j] = min(arr1[j], arr2[j]) 
            else:
                arr2[j] = np.nan
        else:
            arr2[j] = np.nan
    return arr2

@nb.njit()
def custom_np_fmax(arr1, arr2, size):
    for j in range(size):
        if not (np.isnan(arr1[j]) and np.isnan(arr2[j])):
            if np.isnan(arr2[j]):
                arr2[j] = arr1[j]
            elif not np.isnan(arr1[j]):
                arr2[j] = max(arr1[j], arr2[j]) 
            else:
                arr2[j] = np.nan
        else:
            arr2[j] = np.nan
    return arr2

@nb.njit()
def clamp_np_arr(arr, minimum, maximum):
    for j in range(len(arr)):
        if arr[j] < minimum:
            arr[j] = minimum
        if arr[j] > maximum:
            arr[j] = maximum

    return arr

def np_arr_equal_with_nan(arr1,arr2):
    if(len(arr1) != len (arr2)): 
        return False

    for i in range(len(arr1)):
            if arr1[i] != arr2[i]:
                if(not np.isnan(arr1[i]) or not  np.isnan(arr2[i])): 
                    return False

    return True

def true(arr):
    '''
    Returns the indices of the values of the array that are true: just an alias
    for arr.nonzero()[0].

    Args:
        arr (array): any array

    **Example**::

        inds = cv.true(np.array([1,0,0,1,1,0,1])) # Returns array([0, 3, 4, 6])
    '''
    return arr.nonzero()[0]


def false(arr):
    '''
    Returns the indices of the values of the array that are false.

    Args:
        arr (array): any array

    **Example**::

        inds = cv.false(np.array([1,0,0,1,1,0,1]))
    '''
    return np.logical_not(arr).nonzero()[0]


def defined(arr):
    '''
    Returns the indices of the values of the array that are not-nan.

    Args:
        arr (array): any array

    **Example**::

        inds = cv.defined(np.array([1,np.nan,0,np.nan,1,0,1]))
    '''
    return (~np.isnan(arr)).nonzero()[0]


def undefined(arr):
    '''
    Returns the indices of the values of the array that are not-nan.

    Args:
        arr (array): any array

    **Example**::

        inds = cv.defined(np.array([1,np.nan,0,np.nan,1,0,1]))
    '''
    return np.isnan(arr).nonzero()[0]


def itrue(arr, inds):
    '''
    Returns the indices that are true in the array -- name is short for indices[true]

    Args:
        arr (array): a Boolean array, used as a filter
        inds (array): any other array (usually, an array of indices) of the same size

    **Example**::

        inds = cv.itrue(np.array([True,False,True,True]), inds=np.array([5,22,47,93]))
    '''
    return inds[arr]


def ifalse(arr, inds):
    '''
    Returns the indices that are true in the array -- name is short for indices[false]

    Args:
        arr (array): a Boolean array, used as a filter
        inds (array): any other array (usually, an array of indices) of the same size

    **Example**::

        inds = cv.ifalse(np.array([True,False,True,True]), inds=np.array([5,22,47,93]))
    '''
    return inds[np.logical_not(arr)]


def idefined(arr, inds):
    '''
    Returns the indices that are defined in the array -- name is short for indices[defined]

    Args:
        arr (array): any array, used as a filter
        inds (array): any other array (usually, an array of indices) of the same size

    **Example**::

        inds = cv.idefined(np.array([3,np.nan,np.nan,4]), inds=np.array([5,22,47,93]))
    '''
    return inds[~np.isnan(arr)]


def iundefined(arr, inds):
    '''
    Returns the indices that are undefined in the array -- name is short for indices[undefined]

    Args:
        arr (array): any array, used as a filter
        inds (array): any other array (usually, an array of indices) of the same size

    **Example**::

        inds = cv.iundefined(np.array([3,np.nan,np.nan,4]), inds=np.array([5,22,47,93]))
    '''
    return inds[np.isnan(arr)]



def itruei(arr, inds):
    '''
    Returns the indices that are true in the array -- name is short for indices[true[indices]]

    Args:
        arr (array): a Boolean array, used as a filter
        inds (array): an array of indices for the original array

    **Example**::

        inds = cv.itruei(np.array([True,False,True,True,False,False,True,False]), inds=np.array([0,1,3,5]))
    '''
    return inds[arr[inds]]


def ifalsei(arr, inds):
    '''
    Returns the indices that are false in the array -- name is short for indices[false[indices]]

    Args:
        arr (array): a Boolean array, used as a filter
        inds (array): an array of indices for the original array

    **Example**::

        inds = cv.ifalsei(np.array([True,False,True,True,False,False,True,False]), inds=np.array([0,1,3,5]))
    '''
    return inds[np.logical_not(arr[inds])]


def idefinedi(arr, inds):
    '''
    Returns the indices that are defined in the array -- name is short for indices[defined[indices]]

    Args:
        arr (array): any array, used as a filter
        inds (array): an array of indices for the original array

    **Example**::

        inds = cv.idefinedi(np.array([4,np.nan,0,np.nan,np.nan,4,7,4,np.nan]), inds=np.array([0,1,3,5]))
    '''
    return inds[~np.isnan(arr[inds])]


def iundefinedi(arr, inds):
    '''
    Returns the indices that are undefined in the array -- name is short for indices[defined[indices]]

    Args:
        arr (array): any array, used as a filter
        inds (array): an array of indices for the original array

    **Example**::

        inds = cv.iundefinedi(np.array([4,np.nan,0,np.nan,np.nan,4,7,4,np.nan]), inds=np.array([0,1,3,5]))
    '''
    return inds[np.isnan(arr[inds])]