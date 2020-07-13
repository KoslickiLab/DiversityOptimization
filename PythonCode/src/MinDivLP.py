import numpy as np
from .sparse_nnls import sparse_nnls
from scipy.sparse import vstack


def MinDivLP(A_k_small, A_k_large, y_small, y_large, const, q, thresh=0.01):
    """ MinDivLP
    A basic, regularized version of the MinDivLP algorithm.
    Call via:
    x_star = MinDivLP(A_k_small, A_k_large, y_small, y_large, lambda, q)

    Parameters are:
    A_k_small is the[m_small, N] - sized sensing matrix
    A_k_large is the[m_large, N] - sized sensing matrix
    y_small is the data vector of size[m_small, 1]
    y_large is the data vector of size[m_large, 1]
    lambda is the regularization paramater (larger values indicated better
    fit to constraints, at the cost potentially higher execution time and
    may lead to over - fitting if set too large. Typical value is 10000 or
    1000
    q is the parameter used in the MinDivLP algorithm. Must have 0 < q < 1,
    typically, q is set to something like q = 0.1

    Returns:
    x_star: an [N, 1] vector
    """

    B = A_k_large > 0

    epsilon = 0.0001
    denom = np.power(B.T @ y_large, 1 - q) + epsilon
    f = 1/denom

    x_star = sparse_nnls(vstack((f.T, const * A_k_small)), np.append(0, const * y_small))
    x_star = x_star / sum(x_star)
    x_star[np.where(x_star < thresh)] = 0  # Set threshold
    return x_star
