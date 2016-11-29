
import unittest

import numpy as np
import scipy.sparse as sps

from sempy.tensortools import kron3, kron_IID, kron_IDI, kron_DII

mnorm = lambda x: np.max(np.abs(x))


class TestTensorProductForms(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        n = 32
        cls.I   = sps.eye(n)
        cls.D   = np.random.rand(n, n)
        cls.A   = np.random.rand(n,n,n).ravel()
        cls.tol = 1e-12

    def test_kron_IID(self):

        I, D, A = self.I, self.D, self.A
        D1 = kron3(I, I, D)
        err = mnorm(kron_IID(D, A)-D1.dot(A))

        self.assertTrue(err<self.tol)

    def test_kron_IDI(self):

        I, D, A = self.I, self.D, self.A
        D2 = kron3(I, D, I)
        err = mnorm(kron_IDI(D, A)-D2.dot(A))

        self.assertTrue(err<self.tol)

    def test_kron_DII(self):

        I, D, A = self.I, self.D, self.A
        D3 = kron3(D, I, I)
        err = mnorm(kron_DII(D, A)-D3.dot(A))

        self.assertTrue(err<self.tol)
