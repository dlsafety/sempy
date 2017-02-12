
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

        # Switch nodes to lex ordering
        inds = hex_type.get_lexicographic_gmsh_node_indices()
        elem_to_node = elem_to_node[:,inds]
        self.elem_to_node = elem_to_node

        inds = quad_type.get_lexicographic_gmsh_node_indices()
        bndy_face_to_node = bf2n[:,inds]
        self.bndy_face_to_node = bndy_face_to_node

        elem_to_vertex = elem_to_node
        self.elem_to_vertex = elem_to_vertex

        # Spot check the numberings
        assert np.all(np.unique(elem_to_node)==np.arange(len(nodes)))
