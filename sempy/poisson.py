
import numpy as np
na = np.newaxis
import scipy.sparse as sps
from maps import NoConvergence
from sem import SEMhat

from tensortools import kron3, kron_DII, kron_IDI, kron_IID


class PoissonProblem(object):

    def __init__(self, topo, lmap):
        self.topo = topo
        self.lmap = lmap
        self.N    = topo.N
        self.n_elem  = topo.n_elem
        self.n_dofs  = topo.n_dofs

        self.semh = SEMhat(self.N)

    def build(self, node_phys):

        self.node_phys = node_phys
        topo, lmap = self.topo, self.lmap
        N, semh = self.N, self.semh

        n_elem = topo.n_elem
        n_dofs = topo.n_dofs

        wgll = semh.wgll
        wv   = (wgll[:,na,na]*wgll[na,:,na]*wgll[na,na,:]).ravel()

        xgll = semh.xgll
        n = len(xgll)
        quad_ref = np.zeros((n,n,n,3))
        quad_ref[:,:,:,0] = xgll[na,na,:]
        quad_ref[:,:,:,1] = xgll[na,:,na]
        quad_ref[:,:,:,2] = xgll[:,na,na]
        quad_ref = quad_ref.reshape((-1,3))

        nn = (self.N+1)**3
        dof_phys = np.zeros((n_dofs, 3))
        wvals    = np.zeros(n_dofs)

        # build Gij
        G11 = np.zeros((n_elem, nn))
        G12 = np.zeros((n_elem, nn))
        G13 = np.zeros((n_elem, nn))
        G21 = np.zeros((n_elem, nn))
        G22 = np.zeros((n_elem, nn))
        G23 = np.zeros((n_elem, nn))
        G31 = np.zeros((n_elem, nn))
        G32 = np.zeros((n_elem, nn))
        G33 = np.zeros((n_elem, nn))

        etv = topo.elem_to_vertex
        etd = topo.elem_to_dof
        for i in range(n_elem):

            ver = node_phys[etv[i]]
            J   = lmap.calc_jacb(quad_ref, ver)
            Ji  = np.linalg.inv(J)
            j   = np.linalg.det(J).ravel()
            dof_phys[etd[i],:] = lmap.ref_to_phys(quad_ref, ver)

            G0 = np.matmul(Ji, np.transpose(Ji, (0,2,1)))
            G0 *= (wv*j)[:,na,na]
            wvals[etd[i]] += (wv*j)

            G11[i,:] = G0[:,0,0]
            G12[i,:] = G0[:,0,1]
            G13[i,:] = G0[:,0,2]

            G21[i,:] = G0[:,1,0]
            G22[i,:] = G0[:,1,1]
            G23[i,:] = G0[:,1,2]

            G31[i,:] = G0[:,2,0]
            G32[i,:] = G0[:,2,1]
            G33[i,:] = G0[:,2,2]

        self.G11, self.G12, self.G13 = G11, G12, G13
        self.G21, self.G22, self.G23 = G21, G22, G23
        self.G31, self.G32, self.G33 = G31, G32, G33
        self.dof_phys = dof_phys

        # Differentiation operators
        D = semh.Dh
        self.D1  = lambda A: kron_IID(D, A)
        self.D2  = lambda A: kron_IDI(D, A)
        self.D3  = lambda A: kron_DII(D, A)
        DT = D.T
        self.D1T = lambda A: kron_IID(DT, A)
        self.D2T = lambda A: kron_IDI(DT, A)
        self.D3T = lambda A: kron_DII(DT, A)
        self.D = D

        # Build diagonal mass matrix B
        B = sps.dia_matrix((wvals, 0),
                           shape=(n_dofs,n_dofs))
        self.B = B

        # Compute all element centers
        rtop = lmap.ref_to_phys
        zero = np.array([[0.,0.,0.]])
        centers = np.zeros((n_elem, 3))
        for ielem in range(n_elem):
            nodes = node_phys[etv[ielem]]
            centers[ielem] = rtop(zero, nodes).ravel()
        self.centers = centers


    def apply_A(self, x, apply_R=True, apply_Q=True):

        R, Q = self.topo.R, self.topo.Q
        n_elem = self.n_elem

        G11, G12, G13 = self.G11, self.G12, self.G13
        G21, G22, G23 = self.G21, self.G22, self.G23
        G31, G32, G33 = self.G31, self.G32, self.G33

        D1,  D2,  D3  = self.D1,  self.D2,  self.D3
        D1T, D2T, D3T = self.D1T, self.D2T, self.D3T

        if apply_R:
            x = R.T.dot(x)
        if apply_Q:
            x = Q.dot(x)
        x = x.reshape((n_elem, -1))
        y = np.zeros_like(x)
        for i in xrange(n_elem):

            # Keep the steps spit to facilitate profiling
            Dx = D1(x[i])
            a = G11[i,na]*(Dx)
            y[i] += D1T(a)
            a = G21[i,na]*(Dx)
            y[i] += D2T(a)
            a = G31[i,na]*(Dx)
            y[i] += D3T(a)

            Dx = D2(x[i])
            a = G12[i,na]*(Dx)
            y[i] += D1T(a)
            a = G22[i,na]*(Dx)
            y[i] += D2T(a)
            a = G32[i,na]*(Dx)
            y[i] += D3T(a)

            Dx = D3(x[i])
            a = G13[i,na]*(Dx)
            y[i] += D1T(a)
            a = G23[i,na]*(Dx)
            y[i] += D2T(a)
            a = G33[i,na]*(Dx)
            y[i] += D3T(a)

        y = y.ravel()
        if apply_Q:
            y = Q.T.dot(y)
        if apply_R:
            y = R.dot(y)

        return y

    def find_elem_ref(self, X, max_elem_tries=27, check_tol=1e-6):

        assert X.ndim==2 and X.shape[1]==3

        # Return values
        relem = np.zeros(len(X),  dtype=np.int)
        rref  = np.zeros(X.shape, dtype=np.double)

        node_phys = self.node_phys
        etn       = self.topo.elem_to_vertex
        etd       = self.topo.elem_to_dof
        lmap      = self.lmap
        semh      = self.semh
        rtop      = lmap.ref_to_phys

        max_elem_tries = min(self.n_elem, max_elem_tries)

        for ix in range(len(X)):

            x = X[ix,:]

            r2 = np.sum((self.centers-x)**2, axis=-1)
            elem_inds = np.argsort(r2)
            elem_id = -1
            for ielem in elem_inds[:max_elem_tries]:
                nodes = node_phys[etn[ielem]]
                try:
                    xref = lmap.phys_to_ref(x[na,:], nodes)
                    if np.all(lmap.is_interior(xref)):
                        elem_id = ielem
                        ref     = xref.copy()
                        break
                except NoConvergence:
                    # Ignore when the scipy nonlinear solver fails to
                    # converge. The nonlinear solver will often crap
                    # out if the point is far outside the element we
                    # are trying to check.
                    pass

            if elem_id==-1:
                raise LookupError("Element not found: "+str(x))

            # Check that the result makes sense
            d = rtop(ref[na,:], nodes)-x
            if np.max(np.abs(d))>check_tol:
                raise LookupError("Bad (elem, ref) pair: "+str(x))

            relem[ix]  = ielem
            rref[ix,:] = ref

        return (relem, rref)
