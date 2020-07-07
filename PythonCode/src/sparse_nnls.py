# /usr/bin/env python
from numpy import (zeros, ones, finfo, inf, argmax)
from scipy.sparse.linalg import (norm, lsqr)
from scipy.sparse._sparsetools import (csr_matvec, csc_matvec)


def sparse_nnls(C, d, tol=-1, itmax_factor=3):
    """ Calculate argmin ||Cx - d||_2 subject to x >= 0 when C is sparse

    Parameters are:
    C is a scipy.sparse matrix of size m by n
    d is an ndarray of size m or scipy.sparse matrix of size m by 1
    tol: tolerance (optional)
    itmax_factor: factor to determine maximum iterations allowed (optional)

    Returns:
    x: an ndarray that minimizes ||Cx - d||_2 subject to x >= 0
    """

    C = C.tocsc()

    # Set the tolerance
    m, n = C.shape
    tol = 10 * finfo(float).eps * norm(C, 1) * (max(C.shape) + 1) if tol == -1 else tol
    itmax = itmax_factor * n

    # Initialize vector of n zeros and Infs (to be used later)
    wz = zeros(n)

    # Initialize set of non-active columns to null
    P = zeros(n, dtype=bool)

    # Initialize set of active columns to all and the initial point to zeros
    Z = ones(n, dtype=bool)
    x = zeros(n)

    Ctrans = C.T  # transpose c
    dtemp = d  # copy of d

    # resid = d - C*x
    resid = -dtemp
    csc_matvec(m, n, C.indptr, C.indices, C.data, x, resid)
    resid = -resid
    # w = Ctrans*resid
    w = zeros(n)
    csr_matvec(n, m, Ctrans.indptr, Ctrans.indices, Ctrans.data, resid, w)

    # Set up iteration criteria
    outeriter = 0
    i = 0

    # Outer loop to put variables into set to hold positive coefficients
    while any(Z) and any(w[Z] > tol):
        # print(f"On iteration {outeriter}\n")

        outeriter += 1

        # Reset intermediate solution z
        z = zeros(n)

        # Create wz, a Lagrange multiplier vector of variables in the zero set.
        # wz must have the same size as w to preserve the correct indices, so
        # set multipliers to -Inf for variables outside of the zero set.
        wz[P] = -inf
        wz[Z] = w[Z]

        # Find variable with largest Lagrange multiplier
        t = argmax(wz)

        # Move variable t from zero set to positive set
        P[t] = True
        Z[t] = False

        # Compute intermediate solution using only variables in positive set
        z[P] = lsqr(C[:, [i for i, e in enumerate(P) if e]], d)[0]

        # inner loop to remove elements from the positive set which no longer belong
        while any(z[P] <= 0):
            # print("Entering inner loop\n")

            i += 1

            if i > itmax:
                print("sparse_nnls:IterationCountExceeded")
                x = z
                return x

            # Find indices where intermediate solution z is approximately negative
            Q = (z <= 0) & P

            # Choose new x subject to keeping new x nonnegative
            alpha = min(x[Q] / (x[Q] - z[Q]))
            x = x + alpha * (z - x)

            # Reset Z and P given intermediate values of x
            Z = ((abs(x) < tol) & P) | Z
            P = ~Z
            z = zeros(n)  # Reset z
            z[P] = lsqr(C[:, [i for i, e in enumerate(P) if e]], d)[0]  # Re-solve for z

        x = z

        # resid = d - C*x
        resid = -dtemp
        csc_matvec(m, n, C.indptr, C.indices, C.data, x, resid)
        resid = -resid
        # w = Ctrans*resid
        w = zeros(n)
        csr_matvec(n, m, Ctrans.indptr, Ctrans.indices, Ctrans.data, resid, w)

    return x
