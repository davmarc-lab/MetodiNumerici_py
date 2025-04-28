import numpy as np
import numpy.linalg as npl


def norm1(x):
    return np.max(np.sum(np.abs(x), axis=0))


def norm2(A):
    M = np.dot(A.T, A)
    autoval = np.linalg.eigvals(M)
    return np.sqrt(np.max(np.abs(autoval)))


def norminf(x):
    return np.max(np.sum(np.abs(x), axis=1))


def cond(A):
    return norminf(A) * norminf(np.linalg.inv(A))


def errrel(x, px):
    return npl.norm(px - x, np.inf) / npl.norm(x, np.inf)
