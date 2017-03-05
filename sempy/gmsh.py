
import numpy as np
import scipy.sparse as sps

import meshpy
from meshpy.gmsh_reader import read_gmsh
from meshpy.gmsh_reader import (GmshMeshReceiverNumPy,
                                GmshHexahedralElement,
                                GmshQuadrilateralElement,
                                GmshPoint,
                                GmshIntervalElement)

from basis import HexElement


class InvalidGmshMesh(Exception):
    pass


class MeshGmsh(object):

    def __init__(self):
        pass

    def build(self, file_name):

        ### Parse gmsh file
        ###################
        mr = GmshMeshReceiverNumPy()
        read_gmsh(mr, file_name)

        elem_type_inds = {}

        # Find the element types
        for i in range(len(mr.element_types)):

            e = mr.element_types[i]
            if e in elem_type_inds:
                elem_type_inds[e] += [i]
            else:
                elem_type_inds[e]  = [i]

        # Split out the ones we care about
        hex_type  = None
        quad_type = None
        for t in elem_type_inds.keys():

            if isinstance(t, GmshHexahedralElement):
                hex_type  = t
            if isinstance(t, GmshQuadrilateralElement):
                quad_type = t

        if hex_type is None:
            raise InvalidGmshMesh("No hex elements found.")
        if quad_type is None:
            raise InvalidGmshMesh("No quad face elements found.")

        assert hex_type.order==1, "Only linear maps currently supported"

        hex_inds = elem_type_inds[hex_type]
        hex_inds = np.sort(hex_inds)

        quad_inds = elem_type_inds[quad_type]
        quad_inds = np.sort(quad_inds)

        # Build connectivity maps
        elem_to_node = np.zeros((len(hex_inds),
                                 hex_type.node_count()),
                                dtype=np.int)
        for i in range(len(hex_inds)):
            ind = hex_inds[i]
            elem_to_node[i,:] = mr.elements[ind]

        bndy_face_to_node = np.zeros((len(quad_inds),
                                      quad_type.node_count()),
                                     dtype=np.int)
        bf2n = bndy_face_to_node
        for i in range(len(quad_inds)):
            ind = quad_inds[i]
            bf2n[i,:] = mr.elements[ind]

        # Nodes array
        nodes = np.array(mr.points)
        self.nodes = nodes

        vertices = nodes
        self.vertices = vertices

        # Switch nodes to lex ordering
        inds = hex_type.get_lexicographic_gmsh_node_indices()
        elem_to_node = elem_to_node[:,inds]
        self.elem_to_node = elem_to_node

        inds = quad_type.get_lexicographic_gmsh_node_indices()
        bndy_face_to_node = bf2n[:,inds]
        self.bndy_face_to_node = bndy_face_to_node

        elem_to_vertex = elem_to_node
        self.elem_to_vertex = elem_to_vertex

        boundary_vertices = np.unique(bndy_face_to_node.ravel())
        self.boundary_vertices = boundary_vertices

        # Spot check the numberings
        assert np.all(np.unique(elem_to_node)==np.arange(len(nodes)))


        ### Build connectivity arrays
        #############################
        geom = HexElement
        self.geom = geom

        # Build set of edges and edge maps
        elem_to_edge = np.zeros((len(elem_to_vertex),
                                geom.n_edges), dtype=np.int)
        elem_to_edge_dir = np.zeros((len(elem_to_vertex),
                                     geom.n_edges), dtype=np.int)

        edge_id = {}
        eid = 0
        for ielem in range(len(elem_to_vertex)):
            etv = elem_to_vertex[ielem]
            elem_edges = etv[geom.edge_to_vertex]
            for iedge in range(geom.n_edges):
                edge = elem_edges[iedge]
                # Edge direction
                elem_to_edge_dir[ielem, iedge] = 1 if edge[0]<edge[1] else -1
                edge.sort()
                t = tuple(edge)
                if not t in edge_id:
                    edge_id[t] = eid
                    eid += 1
                elem_to_edge[ielem, iedge] = edge_id[t]

        assert len(edge_id)==eid
        edge_to_vertex = np.zeros((len(edge_id), geom.n_vertex_per_edge),
                                  dtype=np.int)
        for k, v in edge_id.iteritems():
            edge_to_vertex[v, :] = k
        assert np.all(edge_to_vertex[:,0]<edge_to_vertex[:,1])
        self.elem_to_edge = elem_to_edge
        self.elem_to_edge_dir = elem_to_edge_dir
        self.edge_to_vertex = edge_to_vertex
        self.edge_id = edge_id

        # Build set of faces and face maps
        elem_to_face = np.zeros((len(elem_to_vertex),
                                 geom.n_faces), dtype=np.int)
        elem_to_face_dir = np.zeros((len(elem_to_vertex),
                                     geom.n_faces), dtype=np.int)
        face_id = {}
        fid = 0
        for ielem in range(len(elem_to_vertex)):
            etv = elem_to_vertex[ielem]
            elem_faces = etv[geom.face_to_vertex]
            for iface in range(geom.n_faces):
                face = elem_faces[iface]
                # How far to rotate the face
                elem_to_face_dir[ielem,iface] = np.argmin(face)
                face.sort()
                t = tuple(face)
                if not t in face_id:
                    face_id[t] = fid
                    fid += 1
                elem_to_face[ielem, iface] = face_id[t]

        assert len(face_id)==fid
        face_to_vertex = np.zeros((len(face_id), geom.n_vertex_per_face),
                                  dtype=np.int)
        for k, v in face_id.iteritems():
            face_to_vertex[v, :] = k

        assert np.all(face_to_vertex[:,0]<face_to_vertex[:,1])
        assert np.all(face_to_vertex[:,1]<face_to_vertex[:,2])
        assert np.all(face_to_vertex[:,2]<face_to_vertex[:,3])
        self.elem_to_face = elem_to_face
        self.elem_to_face_dir = elem_to_face_dir
        self.face_to_vertex = face_to_vertex
        self.face_id = face_id

        # Add component counts
        self.n_elems    = len(elem_to_vertex)
        self.n_vertices = len(vertices)
        self.n_edges    = len(edge_id)
        self.n_faces    = len(face_id)

        # Find edges on the boundary
        boundary_edges = []
        for iedge in range(self.n_edges):
            edge = self.edge_to_vertex[iedge]
            if (edge[0] in boundary_vertices) and\
               (edge[1] in boundary_vertices):
                boundary_edges.append(iedge)
        boundary_edges = np.array(boundary_edges, dtype=np.int)
        self.boundary_edges = boundary_edges

        # Find faces on the boundary
        boundary_faces = []
        for iface in range(self.n_faces):
            face = self.face_to_vertex[iface]
            if (face[0] in boundary_vertices) and\
               (face[1] in boundary_vertices) and\
               (face[2] in boundary_vertices) and\
               (face[3] in boundary_vertices):
                boundary_faces.append(iface)
        boundary_faces = np.array(boundary_faces, dtype=np.int)
        self.boundary_faces = boundary_faces


    def add_basis(self, basis):

        self.basis = basis
        self.N = basis.N

        # DOF counts
        self.n_dofs = self.n_vertices*basis.n_dof_per_vertex+\
                      self.n_edges*basis.n_dof_per_edge+\
                      self.n_faces*basis.n_dof_per_face+\
                      self.n_elems*basis.n_dof_per_bubble

        # Vertex DOFs
        n_vertex_dofs = self.n_vertices*basis.n_dof_per_vertex
        vertex_to_dof = np.arange(n_vertex_dofs)
        self.vertex_to_dof = vertex_to_dof.reshape((self.n_vertices,-1))

        # Edge DOFs
        n_edge_dofs = self.n_edges*basis.n_dof_per_edge
        edge_dofs   = np.arange(n_edge_dofs)+n_vertex_dofs
        edge_to_dof = edge_dofs.reshape((self.n_edges, -1))
        self.edge_to_dof = edge_to_dof

        # Face DOFs
        n_face_dofs = self.n_faces*basis.n_dof_per_face
        face_dofs   = np.arange(n_face_dofs)+n_vertex_dofs+n_edge_dofs
        face_to_dof = face_dofs.reshape((self.n_faces, -1))
        self.face_to_dof = face_to_dof

        # Bubble DOFs
        n_bubble_dofs = self.n_elems*basis.n_dof_per_bubble
        bubble_dofs = np.arange(n_bubble_dofs)+n_vertex_dofs\
                      +n_edge_dofs+n_face_dofs
        bubble_to_dof = bubble_dofs.reshape((self.n_elems, -1))
        self.bubble_to_dof = bubble_to_dof

        # Build elem to dof map
        self._build_elem_to_dof()
        elem_to_dof = self.elem_to_dof

        bdofs = [vertex_to_dof[self.boundary_vertices].ravel(),
                 edge_to_dof[self.boundary_edges].ravel()]
        bdofs.append(face_to_dof[self.boundary_faces].ravel())
        self.boundary_dofs = np.hstack(bdofs)
        n_boundary_dofs = len(self.boundary_dofs)
        self.n_boundary_dofs = n_boundary_dofs

        assert np.all(elem_to_dof>=0)
        assert np.max(elem_to_dof)==(self.n_dofs-1)

        # Build operators
        cols = elem_to_dof.ravel()
        rows = np.arange(len(cols))
        vals = np.ones(len(cols))
        Q = sps.coo_matrix((vals,(rows,cols))).tocsr()
        self.Q = Q

        n_dofs = self.n_dofs
        mask = np.ones(n_dofs, dtype=bool)
        mask[self.boundary_dofs] = False
        cols = np.arange(n_dofs)[mask]
        rows = np.arange(n_dofs-n_boundary_dofs)
        vals = np.ones(len(cols), dtype=np.int)
        R = sps.coo_matrix((vals, (rows, cols)),
                           shape=(n_dofs-n_boundary_dofs, n_dofs)).tocsr()
        self.R = R



    def _build_elem_to_dof(self):
        """ Assumes vertex,edge and face maps exist
        """

        basis = self.basis
        elem_to_dof = np.zeros((self.n_elems, basis.n_dofs),
                                dtype=np.int)-1

        vertex_to_dof = self.vertex_to_dof
        elem_to_vertex = self.elem_to_vertex
        bvtd = basis.vertex_to_dof.ravel()
        elem_to_dof[:,bvtd] = vertex_to_dof[elem_to_vertex].reshape((self.n_elems, -1))

        edge_to_dof = self.edge_to_dof
        elem_to_edge = self.elem_to_edge
        elem_to_edge_dir = self.elem_to_edge_dir
        betd = basis.edge_to_dof
        for ielem in range(self.n_elems):
            for iedge in range(basis.n_edges):
                dofs = edge_to_dof[elem_to_edge[ielem,iedge]]
                if elem_to_edge_dir[ielem,iedge]==-1:
                    dofs = dofs[::-1]
                elem_to_dof[ielem,betd[iedge]] = dofs

        face_to_dof = self.face_to_dof
        elem_to_face = self.elem_to_face
        elem_to_face_dir = self.elem_to_face_dir
        bftd = basis.face_to_dof
        ne_dofs = basis.n_dof_per_edge
        for ielem in range(self.n_elems):
            for iface in range(basis.n_faces):
                dofs = face_to_dof[elem_to_face[ielem,iface]]
                dofs = dofs.reshape((ne_dofs,ne_dofs))
                d = elem_to_face_dir[ielem, iface]
                # Handle rotation in element faces. I belive that they
                # will never be transposed, but should add a check
                # just in case
                if d==1:
                    dofs = dofs.T
                    dofs = dofs[:,::-1]
                elif d==2:
                    dofs = dofs.T
                    dofs = dofs[::-1,:]
                elif d==3:
                    dofs = dofs[::-1,:]
                    dofs = dofs[:,::-1]

                elem_to_dof[ielem,bftd[iface]] = dofs.ravel()

        bubble_to_dof = self.bubble_to_dof
        bbtd = basis.bubble_to_dof.ravel()
        elem_to_dof[:,bbtd] = bubble_to_dof

        self.elem_to_dof = elem_to_dof
