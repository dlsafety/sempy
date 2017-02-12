
import numpy as np


class TensorBasis3D(object):

    vertices = np.array([[-1,-1,-1],
                         [ 1,-1,-1],
                         [-1, 1,-1],
                         [ 1, 1,-1],
                         [-1,-1, 1],
                         [ 1,-1, 1],
                         [-1, 1, 1],
                         [ 1, 1, 1]],
                        dtype=np.int)

    edge_to_vertex = np.array([[0,1],
                               [0,2],
                               [1,3],
                               [2,3],
                               [0,4],
                               [1,5],
                               [2,6],
                               [3,7],
                               [4,5],
                               [4,6],
                               [5,7],
                               [6,7]],
                              dtype=np.int)

    face_to_vertex = np.array([[0,1,2,3],
                               [0,1,4,5],
                               [0,2,4,6],
                               [1,3,5,7],
                               [2,3,6,7],
                               [4,5,6,7]],
                              dtype=np.int)

    n_vertices        = 8
    n_edges           = 12
    n_vertex_per_edge = 2
    n_faces           = 6
    n_edge_per_face   = 4
    n_vertex_per_face = 4


class LagrangeBasisHex(TensorBasis3D):

    is_nodal = True

    def __init__(self, order):
        self.order = order
        self.q = q = order+1
        self.n_dofs = (order+1)**3

        roots = np.linspace(-1, 1, order+1)
        self.roots = roots

        n_dof_per_vertex = 1
        n_dof_per_edge   = order-1
        n_dof_per_face   = (order-1)**2
        n_dof_per_bubble = (order-1)**3

        assert (n_dof_per_vertex*self.n_vertices+\
                n_dof_per_edge*self.n_edges+\
                n_dof_per_face*self.n_faces+\
                n_dof_per_bubble)==self.n_dofs

        vertex_to_dof = np.zeros((self.n_vertices, n_dof_per_vertex),
                                 dtype=np.int)
        edge_to_dof   = np.zeros((self.n_edges, n_dof_per_edge),
                                 dtype=np.int)
        face_to_dof   = np.zeros((self.n_faces, n_dof_per_face),
                                 dtype=np.int)
        bubble_to_dof = np.zeros(n_dof_per_bubble,
                                 dtype=np.int)

        # Assign DOF mappings
        nd = order+1
        dofs = np.arange(nd**3, dtype=np.int).reshape((nd,nd,nd))

        print dofs
        vertex_to_dof[0,0] = dofs[0,0,0]
        vertex_to_dof[1,0] = dofs[0,0,-1]
        vertex_to_dof[2,0] = dofs[0,-1,0]
        vertex_to_dof[3,0] = dofs[0,-1,-1]
        vertex_to_dof[4,0] = dofs[-1,0,0]
        vertex_to_dof[5,0] = dofs[-1,0,-1]
        vertex_to_dof[6,0] = dofs[-1,-1,0]
        vertex_to_dof[7,0] = dofs[-1,-1,-1]

        if n_dof_per_edge>0:
            edge_to_dof[0,:]  = dofs[0,0,1:-1]
            edge_to_dof[1,:]  = dofs[0,1:-1,0]
            edge_to_dof[2,:]  = dofs[0,1:-1,-1]
            edge_to_dof[3,:]  = dofs[0,-1,1:-1]
            edge_to_dof[4,:]  = dofs[1:-1,0,0]
            edge_to_dof[5,:]  = dofs[1:-1,0,-1]
            edge_to_dof[6,:]  = dofs[1:-1,-1,0]
            edge_to_dof[7,:]  = dofs[1:-1,-1,-1]
            edge_to_dof[8,:]  = dofs[-1,0,1:-1]
            edge_to_dof[9,:]  = dofs[-1,1:-1,0]
            edge_to_dof[10,:] = dofs[-1,1:-1,-1]
            edge_to_dof[11,:] = dofs[-1,-1,1:-1]

        if n_dof_per_face>0:
            face_to_dof[0,:] = dofs[0,1:-1,1:-1].ravel()
            face_to_dof[1,:] = dofs[1:-1,0,1:-1].ravel()
            face_to_dof[2,:] = dofs[1:-1,1:-1,0].ravel()
            face_to_dof[3,:] = dofs[1:-1,1:-1,-1].ravel()
            face_to_dof[4,:] = dofs[1:-1,-1,1:-1].ravel()
            face_to_dof[5,:] = dofs[-1,1:-1,1:-1].ravel()

        if n_dof_per_bubble>0:
            bubble_to_dof[:] = dofs[1:-1,1:-1,1:-1].ravel()

        self.vertex_to_dof = vertex_to_dof
        self.edge_to_dof   = edge_to_dof
        self.face_to_dof   = face_to_dof
        self.bubble_to_dof = bubble_to_dof

        dof_check = np.hstack([vertex_to_dof.ravel(),
                               edge_to_dof.ravel(),
                               face_to_dof.ravel(),
                               bubble_to_dof.ravel()])
        dofs = dofs.ravel()
        assert len(dof_check)==len(dofs)
        assert np.all(np.sort(dof_check)==dofs)

        self.n_dof_per_vertex = n_dof_per_vertex
        self.n_dof_per_edge   = n_dof_per_edge
        self.n_dof_per_face   = n_dof_per_face
        self.n_dof_per_bubble = n_dof_per_bubble
