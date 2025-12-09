import zoonosim as zn
import numpy as np

def find_cutoff(age_cutoffs, age):
    '''
    Find which age bin each person belongs to -- e.g. with standard
    age bins 0, 10, 20, etc., ages [5, 12, 4, 58] would be mapped to
    indices [0, 1, 0, 5]. Age bins are not guaranteed to be uniform
    width, which is why this can't be done as an array operation.
    '''
    #print(age_cutoffs) #debug
    #print(age) #debug
    return np.nonzero(age_cutoffs <= age)[0][-1]  # Index of the age bin to use

sim=zn.Sim()
progs = sim.pars['prognoses']['human']
inds = np.fromiter((find_cutoff(progs['age_cutoffs'], this_age) for this_age in [30, 31, 60]), dtype=np.int64, count=3)
print(type(inds))
print(type(progs['symp_probs']))
print(progs['symp_probs'][inds])