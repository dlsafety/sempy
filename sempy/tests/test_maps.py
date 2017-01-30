
import unittest

import numpy as np

from sempy.maps import LinearIsopMap

mnorm = lambda x: np.max(np.abs(x))


class TestLinearIsopMap(unittest.TestCase):

    tol = 1e-12

    def test_phys_to_ref(self):

        lmap = LinearIsopMap()

        scale = 10.
        nodes = lmap.node_ref/scale+5.0
        # random pertubation of [-0.25,0.25)/scale
        np.random.seed(7654)
        peturb = (np.random.rand(*nodes.shape)-0.5)/2.0
        nodes += peturb/scale

        ref = np.array([[0.0,0.0,0.0],
                        [0.5,0.0,0.0],
                        [-0.5,0.0,0.0],
                        [1.,-1.,0.99],
                        [0.75,-0.75,1.5],
                        [1.5,1.5,1.5],
                        [-2.5,2.,-.5]])

        Y = lmap.ref_to_phys(ref, nodes)
        Z = lmap.phys_to_ref(Y, nodes)
        err = mnorm(Z-ref)

        self.assertTrue(err<self.tol)

        is_interior = np.array([True,True,True,True,
                                False,False,False])
        self.assertTrue(np.all(is_interior==lmap.is_interior(Z)))
