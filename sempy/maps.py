
import numpy as np
na = np.newaxis
import scipy.optimize

def cube_points(x):

    n  = len(x)
    X = np.zeros((n,n,n,3))
    X[:,:,:,0] = x[na,na,:]
    X[:,:,:,1] = x[na,:,na]
    X[:,:,:,2] = x[:,na,na]

    return X.reshape((-1,3))


class LinearIsopMap(object):


    vertex_ref = np.array([[-1,-1,-1],
                           [ 1,-1,-1],
                           [-1, 1,-1],
                           [ 1, 1,-1],
                           [-1,-1, 1],
                           [ 1,-1, 1],
                           [-1, 1, 1],
                           [ 1, 1, 1]],
                          dtype=np.double)

    node_ref = vertex_ref

    def __init__(self):
        pass

    def _l1(self, X):
        return (1.0-X)/2.0
    def _dl1(self, X):
        return -np.ones_like(X)/2.0
    def _l2(self, X):
        return (1.0+X)/2.0
    def _dl2(self, X):
        return np.ones_like(X)/2.0

    def ref_to_phys(self, X, nodes):

        v1 = self._l1(X)
        v2 = self._l2(X)
        v1x, v1y, v1z = v1.T
        v2x, v2y, v2z = v2.T

        P = (v1x*v1y*v1z)[:,na]*nodes[0,:]+\
            (v2x*v1y*v1z)[:,na]*nodes[1,:]+\
            (v1x*v2y*v1z)[:,na]*nodes[2,:]+\
            (v2x*v2y*v1z)[:,na]*nodes[3,:]+\
            (v1x*v1y*v2z)[:,na]*nodes[4,:]+\
            (v2x*v1y*v2z)[:,na]*nodes[5,:]+\
            (v1x*v2y*v2z)[:,na]*nodes[6,:]+\
            (v2x*v2y*v2z)[:,na]*nodes[7,:]

        return P

    # Ref points used to find init guess for non-linear solver
    _Xgref  = cube_points(np.linspace(-1,1,20))
    # Non-linear solver tolerance
    _tol_phys_to_ref = 1e-14
    def phys_to_ref(self, Y, nodes):

        solve = scipy.optimize.broyden1
        rtop   = self.ref_to_phys

        # Build initial guess
        #####################
        center = rtop(np.array([[0.,0.,0.]]), nodes).ravel()
        mins, maxs = np.sort(nodes, axis=0)[[0,-1]]

        # Find the closest point to Y to use as our initial guess
        Xgref = self._Xgref
        Yg    = rtop(Xgref, nodes)
        x0 = np.zeros_like(Y)
        for i in range(len(Y)):
            dist2 = np.sum((Yg-Y[i])**2, axis=-1)
            x0[i] = Xgref[np.argmin(dist2)]

        # Solve the non-linear system
        #############################
        def F(x):
            return Y-rtop(x, nodes)

        f_tol = self._tol_phys_to_ref
        if np.max(np.abs(F(x0)))<f_tol:
            # For some godforsaken reason the non-linear solver fails
            # when the initial guess is too good
            return x0
        else:
            # Same thing for each of the dimensions
            x0 += 0.1*(np.abs(F(x0))<1e-6)
            xref = solve(F, x0, f_tol=f_tol)
            return xref

    def calc_jacb(self, X, nodes):

        v1 = self._l1(X)
        v2 = self._l2(X)
        v1x, v1y, v1z = v1.T
        v2x, v2y, v2z = v2.T

        dv1 = self._dl1(X)
        dv2 = self._dl2(X)
        dv1x, dv1y, dv1z = dv1.T
        dv2x, dv2y, dv2z = dv2.T

        J = np.zeros((len(X),3,3))

        t1x, t2x, t3x = ((dv1x*v1y*v1z)[:,na]*nodes[0,:]+\
                         (dv2x*v1y*v1z)[:,na]*nodes[1,:]+\
                         (dv1x*v2y*v1z)[:,na]*nodes[2,:]+\
                         (dv2x*v2y*v1z)[:,na]*nodes[3,:]+\
                         (dv1x*v1y*v2z)[:,na]*nodes[4,:]+\
                         (dv2x*v1y*v2z)[:,na]*nodes[5,:]+\
                         (dv1x*v2y*v2z)[:,na]*nodes[6,:]+\
                         (dv2x*v2y*v2z)[:,na]*nodes[7,:]).T

        t1y, t2y, t3y = ((v1x*dv1y*v1z)[:,na]*nodes[0,:]+\
                         (v2x*dv1y*v1z)[:,na]*nodes[1,:]+\
                         (v1x*dv2y*v1z)[:,na]*nodes[2,:]+\
                         (v2x*dv2y*v1z)[:,na]*nodes[3,:]+\
                         (v1x*dv1y*v2z)[:,na]*nodes[4,:]+\
                         (v2x*dv1y*v2z)[:,na]*nodes[5,:]+\
                         (v1x*dv2y*v2z)[:,na]*nodes[6,:]+\
                         (v2x*dv2y*v2z)[:,na]*nodes[7,:]).T

        t1z, t2z, t3z = ((v1x*v1y*dv1z)[:,na]*nodes[0,:]+\
                         (v2x*v1y*dv1z)[:,na]*nodes[1,:]+\
                         (v1x*v2y*dv1z)[:,na]*nodes[2,:]+\
                         (v2x*v2y*dv1z)[:,na]*nodes[3,:]+\
                         (v1x*v1y*dv2z)[:,na]*nodes[4,:]+\
                         (v2x*v1y*dv2z)[:,na]*nodes[5,:]+\
                         (v1x*v2y*dv2z)[:,na]*nodes[6,:]+\
                         (v2x*v2y*dv2z)[:,na]*nodes[7,:]).T

        J[:,0,0] = t1x
        J[:,0,1] = t1y
        J[:,0,2] = t1z
        J[:,1,0] = t2x
        J[:,1,1] = t2y
        J[:,1,2] = t2z
        J[:,2,0] = t3x
        J[:,2,1] = t3y
        J[:,2,2] = t3z

        return J

    def is_interior(self, X, tol=1e-12):
        """ Checks if points are in the interior of the reference element
        """

        return np.max(np.abs(X), axis=-1)<=1.0+tol
