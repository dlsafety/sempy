
import numpy as np
from numpy import newaxis
import scipy
import scipy.sparse as sps
from scipy.sparse.linalg import spsolve
from scipy.special import p_roots

from pyfem.sem  import SEMhat
from pyfem.topo import Interval
from pyfem.poly import gll_points
from pyfem.poly import eval_lagrange_d0 as eval_phi1d

from tensortools import kron3, kron_DII, kron_IDI, kron_IID


class HexCubePoisson(object):

    def __init__(self, N, E, L=1.0, periodic=False):
        self.N, self.E = N, E
        self.L = L
        self.periodic = periodic

    def build_mesh(self):

        N, E = self.N, self.E
        L = self.L
        periodic = self.periodic

        n      = N+1
        n_dofs = N*E

        if not periodic:
            n_dofs += 1
        self.n_dofs = n_dofs

        semh = SEMhat(N)

        vertices = np.linspace(0, L, E+1)
        etv      = np.zeros((E, 2), dtype=np.int)
        etv[:,0] = np.arange(E)
        etv[:,1] = np.arange(E)+1
        topo  = Interval()
        xq = topo.ref_to_phys(vertices[etv], semh.xgll)
        if periodic:
            xq = xq[:,:-1]
        self.xq = xq
        jacb_det = topo.calc_jacb(vertices[etv])[0]
        self.jacb_det = jacb_det

        if periodic:
            etv[0,0] = etv[-1,-1]

        # Make 1D elem to dof map
        etd = np.arange(E*(N+1))
        etd = etd.reshape((E, -1))
        etd -= np.arange(E)[:,newaxis]
        if periodic:
            etd[-1,-1] = etd[0,0]

        # Make Q
        cols = etd.ravel()
        rows = np.arange(len(cols))
        vals = np.ones(len(cols))
        Q0 = sps.coo_matrix((vals,(rows,cols)),
                            dtype=np.int)
        Q = kron3(Q0, Q0, Q0)
        self.Q = Q
        self.Q0 = Q0

        # Build etd for full mesh
        nd = E*(N+1)
        dofs = Q.dot(np.arange(Q.shape[1])).reshape((nd,nd,nd))
        etd = np.zeros((E**3, n**3), dtype=np.int)
        ind = 0
        for iz in range(E):
            for iy in range(E):
                for ix in range(E):
                    a = dofs[iz*n:iz*n+n,iy*n:iy*n+n,ix*n:ix*n+n]
                    etd[ind,:] = a.ravel()
                    ind += 1
        self.etd = etd
        dofs = None

        # Q for assembly. NOT the tensor product Q.
        cols = etd.ravel()
        rows = np.arange(len(cols))
        vals = np.ones(len(cols))
        QA   = sps.coo_matrix((vals,(rows,cols)),
                               dtype=np.int)
        self.QA = QA

        # Make R
        if periodic:
            R0 = sps.eye(n_dofs)
        else:
            R0 = sps.dia_matrix((np.ones(n_dofs),1),
                                 shape=(n_dofs-2,n_dofs))
        R = kron3(R0, R0, R0)
        self.R = R

        if (not periodic):
            dofs = np.arange(n_dofs**3)
            bndy_dofs = list(set(dofs)-set(R.dot(dofs)))
            bndy_dofs = np.array(bndy_dofs, dtype=np.int)
            self.bndy_dofs = bndy_dofs
            dofs = None

        ## Build A and B
        Al  = sps.kron(sps.eye(E), semh.Ah/jacb_det)
        A0  = Q0.T.dot(Al.dot(Q0))
        A0  = R0.dot(A0.dot(R0.T))
        Ah0 = semh.Ah/jacb_det

        # Full local mass matrix
        xgll, wgll = gll_points(n)
        xg, wg     = p_roots(n)
        L = eval_phi1d(xgll, xg).T
        Bf = sps.dia_matrix((wg,0),shape=(n,n))
        Bf = L.T.dot(Bf.dot(L))

        # Bl0 = sps.kron(sps.eye(E), semh.Bh*jacb_det)
        Bl0 = sps.kron(sps.eye(E), Bf*jacb_det)
        B0  = Q0.T.dot(Bl0.dot(Q0))
        B0  = R0.dot(B0.dot(R0.T))
        Bh0 = Bf*jacb_det
        Bh  = kron3(Bh0, Bh0, Bh0).toarray()
        self.Bh = Bh

        Ah  = kron3(Bh0,Bh0,Ah0)
        Ah += kron3(Bh0,Ah0,Bh0)
        Ah += kron3(Ah0,Bh0,Bh0)
        self.Ah = Ah.toarray()
        # Bh = kron3(Bh0, Bh0, Bh0)
        # self.Bh = Bh

        eigvals, eigvecs = scipy.linalg.eigh(A0.toarray(), B0.toarray())
        self.eigvals, self.eigvecs = eigvals, eigvecs
        n_eigs = len(eigvals)

        L = sps.dia_matrix((eigvals, 0),
                           shape=(n_eigs,n_eigs))
        eye = sps.eye(n_eigs)
        C1 = kron3(eye,eye,L)+kron3(eye,L,eye)+kron3(L,eye,eye)
        is_zero = np.abs(C1.data)<1e-12
        C1.data[is_zero] = 1.0
        C1.data = 1.0/C1.data
        C1.data[is_zero] = 0.0
        self.C1 = C1

    def solve(self, rhs, u=None,
              assemble=False):
        """Tensor solve

        Assumes that rhs includes boundary dofs if periodic=False.

        :param rhs:
        :param assemble: Perform assembly
        :returns:
        :rtype:

        """

        Q, R = self.QA, self.R
        N, E = self.N, self.E
        periodic = self.periodic
        Ah, Bh = self.Ah, self.Bh

        if assemble:
            rhs = Q.dot(rhs)
            rhs = rhs.reshape((E**3,-1))
            rhs = rhs.dot(Bh.T).ravel()
            rhs = Q.T.dot(rhs)

        if periodic:
            rhs -= np.mean(rhs)

        # Set dirch BCs
        if not (u is None):
            assert (not periodic)
            bdofs = self.bndy_dofs
            v = Q.dot(u)
            v = v.reshape((E**3,-1))
            v = v.dot(Ah.T).ravel()
            v = Q.T.dot(v)
            rhs -= v

        # rhs = R.dot(rhs)

        eigvals, eigvecs = self.eigvals, self.eigvecs
        n_eigs = len(eigvals)
        C1 = self.C1

        U = rhs.reshape((n_eigs, n_eigs, n_eigs))

        eigvecsT = eigvecs.T
        W1 = kron_IID(eigvecsT, U)
        W1 = kron_IDI(eigvecsT, W1)
        W1 = kron_DII(eigvecsT, W1)

        W2 = C1.dot(W1.ravel()).reshape((n_eigs, n_eigs, n_eigs))

        W3 = kron_IID(eigvecs, W2)
        W3 = kron_IDI(eigvecs, W3)
        W3 = kron_DII(eigvecs, W3)

        # sol  = R.T.dot(W3.ravel())
        sol = W3.ravel()
        
        if periodic:
            sol -= np.mean(sol)

        if not (u is None):
            bdofs = self.bndy_dofs
            sol[bdofs] = u[bdofs]

        return sol

    def get_dof_phys(self):

        n_dofs = self.n_dofs
        d = np.unique(self.xq)
        dof_phys = np.zeros((n_dofs**3, 3))
        ind = 0
        for iz in range(n_dofs):
            for iy in range(n_dofs):
                for ix in range(n_dofs):
                    dof_phys[ind, 0] = d[ix]
                    dof_phys[ind, 1] = d[iy]
                    dof_phys[ind, 2] = d[iz]
                    ind += 1

        return dof_phys
