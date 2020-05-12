import numpy as np
from scipy.optimize import nnls

def MinDivLP(A_k_small, A_k_large, y_small, y_large, const, q):
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
    f = 1 / np.power(B.T @ y_large, 1-q)

    x_star = nnls(np.vstack( (f.T, const*A_k_small.toarray()) ), np.append(0, const*y_small))[0]
    return x_star / sum(x_star)