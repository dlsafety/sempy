
import numpy as np

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

        # Build set of edges and edge maps
        elem_to_edge = np.zeros((len(elem_to_vertex),
                                geom.n_edges), dtype=np.int)
        edge_id = {}
        eid = 0
        for ielem in range(len(elem_to_vertex)):
            etv = elem_to_vertex[ielem]
            elem_edges = etv[geom.edge_to_vertex]
            for iedge in range(geom.n_edges):
                edge = elem_edges[iedge]
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
        self.edge_to_vertex = edge_to_vertex
        self.edge_id = edge_id

        # Build set of faces and face maps
        elem_to_face = np.zeros((len(elem_to_vertex),
                                 geom.n_faces), dtype=np.int)
        face_id = {}
        fid = 0
        for ielem in range(len(elem_to_vertex)):
            etv = elem_to_vertex[ielem]
            elem_faces = etv[geom.face_to_vertex]
            for iface in range(geom.n_faces):
                face = elem_faces[iface]
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
        self.face_to_vertex = face_to_vertex
        self.face_id = face_id

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
