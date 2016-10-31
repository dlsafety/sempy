
import numpy as np
na = np.newaxis
import scipy
import scipy.sparse as sps

kron3 = lambda x,y,z: sps.kron(x,sps.kron(y,z))

class CubicTopology(object):

    def __init__(self, N, E, periodic=False):

        self.N, self.E = N, E
        self.periodic = periodic

        # Unpack E
        if hasattr(E, "__len__"):
            self.Ex, self.Ey, self.Ez = E
        else:
            self.Ex, self.Ey, self.Ez = E, E, E

    def build(self):

        Ex, Ey, Ez = self.Ex, self.Ey, self.Ez
        N, periodic = self.N, self.periodic

        self.nx = nx = N+1
        self.ny = ny = N+1
        self.nz = nz = N+1
        nx_dofs = N*Ex+1
        ny_dofs = N*Ey+1
        nz_dofs = N*Ez+1
        self.n_elem = n_elem = Ex*Ey*Ez

        if periodic:
            nx_dofs -= 1
            ny_dofs -= 1
            nz_dofs -= 1
        n_dofs = nz_dofs*ny_dofs*nx_dofs
        self.n_dofs = n_dofs
        self.nx_dofs = nx_dofs
        self.ny_dofs = ny_dofs
        self.nz_dofs = nz_dofs

        # Build elem to mesh vertex map
        etv = np.zeros((n_elem, 8), dtype=np.int)
        ind = 0
        for iz in range(Ez):
            for iy in range(Ey):
                for ix in range(Ex):
                    etv[ind, 0] = ix+iy*(Ex+1)
                    etv[ind, 1] = ix+iy*(Ex+1)+1
                    etv[ind, 2] = ix+(iy+1)*(Ex+1)
                    etv[ind, 3] = ix+(iy+1)*(Ex+1)+1
                    etv[ind,:4] += iz*(Ex+1)*(Ey+1)
                    etv[ind,4:] = etv[ind,:4]+(Ex+1)*(Ey+1)
                    ind += 1

        self.elem_to_vertex = etv

        # Build restriction operator
        if periodic:
            R0x = sps.eye(nx_dofs)
            R0y = sps.eye(ny_dofs)
            R0z = sps.eye(nz_dofs)
        else:
            R0x = sps.dia_matrix((np.ones(nx_dofs),1),
                                  shape=(nx_dofs-2,nx_dofs))
            R0y = sps.dia_matrix((np.ones(ny_dofs),1),
                                  shape=(ny_dofs-2,ny_dofs))
            R0z = sps.dia_matrix((np.ones(nz_dofs),1),
                                  shape=(nz_dofs-2,nz_dofs))

        R = kron3(R0z, R0y, R0x)
        self.R = R

        # Compute boundary dofs
        if not periodic:
            bd = set(np.arange(n_dofs))-set(R.dot(np.arange(n_dofs)))
            boundary_dofs = np.sort(np.array(list(bd)))
        else:
            boundary_dofs = np.array([],dtype=np.int)

        boundary_dofs.sort()
        self.boundary_dofs = boundary_dofs

        # Build elem to dof maps
        etd = np.zeros((n_elem, nx*ny*nz), dtype=np.int)
        rngx = np.arange(nx)
        rngy = np.arange(ny)
        rngz = np.arange(nz)
        nxy_dofs = nx_dofs*ny_dofs

        ind = 0
        for iz in range(Ez):
            for iy in range(Ey):
                for ix in range(Ex):
                    indz = iz*N
                    indy = iy*N
                    indx = ix*N

                    e  = (rngx[na,na,:]+indx)%nx_dofs+\
                         ((rngy[na,:,na]+indy)*nx_dofs)%nxy_dofs+\
                         (rngz[:,na,na]+indz)*nxy_dofs
                    e = e%n_dofs
                    etd[ind,:] = e.ravel()

                    ind += 1

        self.elem_to_dof = etd

        cols = etd.ravel()
        rows = np.arange(len(cols))
        vals = np.ones(len(cols))
        Q = sps.coo_matrix((vals,(rows,cols))).tocsr()
        self.Q = Q

    def get_vertex_ref(self):
        """ Convenience function
        """
        Ex, Ey, Ez = self.Ex, self.Ey, self.Ez

        vertx = np.linspace(-1,1,Ex+1)
        verty = np.linspace(-1,1,Ey+1)
        vertz = np.linspace(-1,1,Ez+1)

        XYZ = np.zeros((Ez+1,Ey+1,Ex+1,3))
        XYZ[:,:,:,0] = vertx[na,na,:]
        XYZ[:,:,:,1] = verty[na,:,na]
        XYZ[:,:,:,2] = vertz[:,na,na]
        vertex_ref = XYZ.reshape((-1,3))

        return vertex_ref
