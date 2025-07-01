import os
import sys
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

import numpy as np
import symplectic as sy
import symplectic_sampler as symsa


##test for symplectic_sampler()
#
def test_symplectic_sampler():
    n = 8
    it = 30
    J = sy.symplectic_J_product_coord(n)
    for i in range(it):
        S = symsa.symplectic_sampler(n)
        assert np.array_equal(S.transpose() @ J @ S % 2, J), "Sampler output is not a symplectic matrix"