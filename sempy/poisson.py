
import numpy as np
na = np.newaxis
import scipy.sparse as sps

from pyfem.sem import SEMhat

kron3 = lambda x,y,z: sps.kron(x,sps.kron(y,z))

class PoissonProblem(object):

    def __init__(self, topo, lmap):
        self.topo = topo
        self.lmap = lmap
        self.N    = topo.N
        self.n_elem  = topo.n_elem
        self.n_dofs  = topo.n_dofs
        self.nx_dofs = topo.nx_dofs
        self.ny_dofs = topo.ny_dofs
        self.nz_dofs = topo.nz_dofs

        self.semh = SEMhat(self.N)

    def build(self, vertex_phys):

        topo, lmap = self.topo, self.lmap
        N, semh = self.N, self.semh

        nx, ny, nz = topo.nx, topo.ny, topo.nz
        nx_dofs, ny_dofs = topo.nx_dofs, topo.ny_dofs
        nz_dofs, n_elem  = topo.nz_dofs, topo.n_elem
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

        # build Gij
        G11 = []
        G12 = []
        G13 = []
        G21 = []
        G22 = []
        G23 = []
        G31 = []
        G32 = []
        G33 = []

        nn = nx*ny*nz
        s  = (nn, nn)
        dof_phys = np.zeros((nx_dofs*ny_dofs*nz_dofs, 3))
        wvals    = np.zeros(nx_dofs*ny_dofs*nz_dofs)

        etv = topo.elem_to_vertex
        etd = topo.elem_to_dof
        for i in range(n_elem):

            ver = vertex_phys[etv[i]]
            J   = lmap.calc_jacb(quad_ref, ver)
            Ji  = np.linalg.inv(J)
            j   = np.linalg.det(J).ravel()
            dof_phys[etd[i],:] = lmap.ref_to_phys(quad_ref, ver)

            G0 = np.matmul(Ji, np.transpose(Ji, (0,2,1)))
            G0 *= (wv*j)[:,na,na]
            wvals[etd[i]] += (wv*j)

            G11 += [sps.dia_matrix((G0[:,0,0], 0), shape=s)]
            G12 += [sps.dia_matrix((G0[:,0,1], 0), shape=s)]
            G13 += [sps.dia_matrix((G0[:,0,2], 0), shape=s)]

            G21 += [sps.dia_matrix((G0[:,1,0], 0), shape=s)]
            G22 += [sps.dia_matrix((G0[:,1,1], 0), shape=s)]
            G23 += [sps.dia_matrix((G0[:,1,2], 0), shape=s)]

            G31 += [sps.dia_matrix((G0[:,2,0], 0), shape=s)]
            G32 += [sps.dia_matrix((G0[:,2,1], 0), shape=s)]
            G33 += [sps.dia_matrix((G0[:,2,2], 0), shape=s)]

        self.G11, self.G12, self.G13 = G11, G12, G13
        self.G21, self.G22, self.G23 = G21, G22, G23
        self.G31, self.G32, self.G33 = G31, G32, G33
        self.dof_phys = dof_phys

        # Differentiation operators
        D1 = kron3(sps.eye(nz), sps.eye(ny), semh.Dh)
        D2 = kron3(sps.eye(nz), semh.Dh,     sps.eye(nx))
        D3 = kron3(semh.Dh,     sps.eye(ny), sps.eye(nx))
        self.D1, self.D2, self.D3 = D1, D2, D3

        # Build mass matrix B
        B = sps.dia_matrix((wvals, 0),
                           shape=(n_dofs,n_dofs))
        self.B = B

    def apply_A(self, x, apply_R=True, apply_Q=True):

        R, Q = self.topo.R, self.topo.Q
        n_elem = self.n_elem

        G11, G12, G13 = self.G11, self.G12, self.G13
        G21, G22, G23 = self.G21, self.G22, self.G23
        G31, G32, G33 = self.G31, self.G32, self.G33
        D1, D2, D3 = self.D1, self.D2, self.D3

        if apply_R:
            x = R.T.dot(x)
        if apply_Q:
            x = Q.dot(x)
        x = x.reshape((n_elem, -1))
        y = np.zeros_like(x)
        for i in xrange(n_elem):
            Dx = D1.dot(x[i])
            y[i] += D1.T.dot(G11[i].dot(Dx))
            y[i] += D2.T.dot(G21[i].dot(Dx))
            y[i] += D3.T.dot(G31[i].dot(Dx))
            Dx = D2.dot(x[i])
            y[i] += D1.T.dot(G12[i].dot(Dx))
            y[i] += D2.T.dot(G22[i].dot(Dx))
            y[i] += D3.T.dot(G32[i].dot(Dx))
            Dx = D3.dot(x[i])
            y[i] += D1.T.dot(G13[i].dot(Dx))
            y[i] += D2.T.dot(G23[i].dot(Dx))
            y[i] += D3.T.dot(G33[i].dot(Dx))

        y = y.ravel()
        if apply_Q:
            y = Q.T.dot(y)
        if apply_R:
            y = R.dot(y)

        return y

