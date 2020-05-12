## Simple test for MinDivLP code
import numpy as np
from numpy.linalg import norm
from PythonCode.src.MinDivLP import MinDivLP
from numpy.linalg import norm

A_k_small = np.array(
    [(.5, 0, 0),
     (0, 1 / 3, 1 / 5),
     (.5, 2 / 3, 4 / 5)])

A_k_large = np.array(
    [(0, 0, 1 / 6),
     (1 / 4, 1 / 3, 0),
     (1 / 4, 0, 2 / 6),
     (1 / 4, 1 / 3, 2 / 6),
     (1 / 4, 1 / 3, 1 / 6)])

x_true = np.array([1, 0, 0]).reshape(3, 1)  # to test against

y_small = A_k_small @ x_true
y_large = A_k_large @ x_true

const = 10000
q = 0.1

# the algorithm
x_star = MinDivLP(A_k_small, A_k_large, y_small, y_large, const, q).reshape((3, 1))

# Then do the tests

# a realistic test
thresh = .00001
L1_error = norm(x_star - x_true, 1)
assert L1_error < thresh, 'L1 reconstruction error is above %f' % thresh

# a test for EXACT recovery (which will probably never pass on realistic data)
assert np.sum(x_star != x_true) == 0, 'Reconstructed and true vectors are not exactly the same'
