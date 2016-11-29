
import numpy as np
import scipy.sparse as sps

def kron3(A, B, C):
    return sps.kron(A, sps.kron(B, C))

def kron_IID(D, A):

    n = D.shape[1]
    A = A.reshape((n,n,n))
    return np.einsum('kj,aij', D, A).ravel()

def kron_IDI(D, A):

    n = D.shape[1]
    A = A.reshape((n,n,n))
    return np.einsum('ij,ajk', D, A).ravel()

def kron_DII(D, A):

    n = D.shape[1]
    A = A.reshape((n,n,n))
    return np.einsum('ij,jkl', D, A).ravel()
