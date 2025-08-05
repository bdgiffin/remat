# Module for mesh pre-processing

# Import GMSH for meshing
import pygmsh

# Import various needed Python packages
import sys
from math import *
import numpy as np

# --------------------------------------------------------------------------



# --------------------------------------------------------------------------

def find_mesh_boundary(elements):
    """
    Finds the boundary edges of a 2D mixed quad/tri mesh.

    Args:
        elements (np.ndarray): Nelems x 4 array of element connectivity (node indices).

    Returns:
        boundary_nodes (np.ndarray): Nboundary_nodes x 1 array of boundary nodes
        boundary_edges (np.ndarray): Nedges x 2A array of boundary edges (node_idx1, node_idx2).
    """
    edge_counts = {}
    for element in elements:
        # Extract edges of the element
        if (element[2] == element[3]):
            edges = [(element[0], element[1]), (element[1], element[2]), (element[2], element[0])]
        else:
            edges = [(element[0], element[1]), (element[1], element[2]), (element[2], element[3]), (element[3], element[0])]
        for edge in edges:
            # Ensure consistent order for edge (e.g., (min_idx, max_idx))
            orientation = +1
            if (edge[0] > edge[1]):
                orientation = -1
            sorted_edge = tuple(sorted(edge))
            edge_counts[sorted_edge] = (edge_counts.get(sorted_edge, 0) + 1, orientation)

    boundary_edges = []
    for edge, count_and_orientation in edge_counts.items():
        if count_and_orientation[0] == 1:
            if count_and_orientation[1] == +1:
                boundary_edges.append(edge)
            elif count_and_orientation[1] == -1:
                boundary_edges.append((edge[1], edge[0]))
    boundary_edges = np.array(boundary_edges)
    boundary_nodes = np.unique(boundary_edges.flatten())
                
    return boundary_nodes, boundary_edges

# --------------------------------------------------------------------------

def merge_meshes(mesh1,mesh2):

# --------------------------------------------------------------------------
