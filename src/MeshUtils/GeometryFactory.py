import numpy as np
import sys, platform
import os
from Geometry import *

# Packages for reading/writing mesh files
if (not sys.platform == "emscripten"):
    import pyexodus
    from frame_maker_file_io import read_dxf

class GeometryFactory:
    
    def cartesian_grid(self,bottom_corner,top_corner,subdivisions,perturbation=0.0):
        Ndims = min(len(bottom_corner),len(top_corner),len(subdivisions))
        Nnodes = 1
        Nelems = 1
        edge_lengths = [0.0]*Ndims
        for i in range(Ndims):
            Nnodes = Nnodes*(subdivisions[i]+1)
            Nelems = Nelems*(subdivisions[i])
            edge_lengths[i] = (top_corner[i]-bottom_corner[i])/subdivisions[i]

        vertices = np.zeros((Nnodes,Ndims))
        connectivity = np.zeros((Nelems,2**Ndims),dtype=int)
        if (Ndims == 2):
            for j in range(subdivisions[1]+1):
                y = bottom_corner[1] + edge_lengths[1]*j
                for i in range(subdivisions[0]+1):
                    x = bottom_corner[0] + edge_lengths[0]*i
                    index = (subdivisions[0]+1)*j+i
                    vertices[index,0] = x
                    if ((i != 0) and (i != subdivisions[0])):
                        vertices[index,0] = vertices[index,0] + perturbation*edge_lengths[0]*(-1)**(i+j)
                    vertices[index,1] = y
                    if ((j != 0) and (j != subdivisions[1])):
                        vertices[index,1] = vertices[index,1] + perturbation*edge_lengths[1]*(-1)**(i+j)
            for j in range(subdivisions[1]):
                for i in range(subdivisions[0]):
                    elemID = (subdivisions[0])*j+i
                    index  = (subdivisions[0]+1)*j+i
                    connectivity[elemID,0] = index
                    connectivity[elemID,1] = index+1
                    connectivity[elemID,2] = index+1+(subdivisions[0]+1)
                    connectivity[elemID,3] = index  +(subdivisions[0]+1)
            return QuadGeometry(vertices,connectivity)
        elif (Ndims == 3):
            for k in range(subdivisions[2]+1):
                z = bottom_corner[2] + edge_lengths[2]*k
                for j in range(subdivisions[1]+1):
                    y = bottom_corner[1] + edge_lengths[1]*j
                    for i in range(subdivisions[0]+1):
                        x = bottom_corner[0] + edge_lengths[0]*i
                        index = (subdivisions[1]+1)*(subdivisions[0]+1)*k+(subdivisions[0]+1)*j+i
                        vertices[index,0] = x
                        if ((i != 0) and (i != subdivisions[0])):
                            vertices[index,0] = vertices[index,0] + perturbation*edge_lengths[0]*(-1)**(i+j+k)
                        vertices[index,1] = y
                        if ((j != 0) and (j != subdivisions[1])):
                            vertices[index,1] = vertices[index,1] + perturbation*edge_lengths[1]*(-1)**(i+j+k)
                        vertices[index,2] = z
                        if ((k != 0) and (k != subdivisions[2])):
                            vertices[index,2] = vertices[index,2] + perturbation*edge_lengths[2]*(-1)**(i+j+k)
            for k in range(subdivisions[2]):
                for j in range(subdivisions[1]):
                    for i in range(subdivisions[0]):
                        elemID = (subdivisions[1])*(subdivisions[0])*k+(subdivisions[0])*j+i
                        index  = (subdivisions[1]+1)*(subdivisions[0]+1)*k+(subdivisions[0]+1)*j+i
                        connectivity[elemID,0] = index
                        connectivity[elemID,1] = index+1
                        connectivity[elemID,2] = index+1+(subdivisions[0]+1)
                        connectivity[elemID,3] = index  +(subdivisions[0]+1)
                        connectivity[elemID,4] = index                      +(subdivisions[1]+1)*(subdivisions[0]+1)
                        connectivity[elemID,5] = index+1                    +(subdivisions[1]+1)*(subdivisions[0]+1)
                        connectivity[elemID,6] = index+1+(subdivisions[0]+1)+(subdivisions[1]+1)*(subdivisions[0]+1)
                        connectivity[elemID,7] = index  +(subdivisions[0]+1)+(subdivisions[1]+1)*(subdivisions[0]+1)
            return HexGeometry(vertices,connectivity)

    def from_dxf(self,filename,layer_names):
        # read the DXF file, and return:
        # 1.) the list of nodal coordinates
        # 2.) element connectivities (grouped by layer)
        # 3.) the list of layer names
        points, layer_connectivity, layer_names = read_dxf(filename,layer_names)

        # get mesh totals
        Nnodes = points.shape[0]
        Nelems = layer_connectivity[0].shape[0]

        # get node coordinates
        coordinates = np.zeros((Nnodes,2))
        for i in range(0,Nnodes):
            coordinates[i,0] = points[i,0]
            coordinates[i,1] = points[i,1]

        # get element connectivity
        first_layer = layer_connectivity[0]
        connectivity = np.zeros((Nelems,2),dtype=int)
        for i in range(0,Nelems):
            connectivity[i,0] = first_layer[i,0]
            connectivity[i,1] = first_layer[i,1]

        # assume that the first layer is the only layer
        return TrussGeometry(coordinates,connectivity)

    def from_exodus(self,filename,block_id):
        # This assumes each Exodus file contains only one element block
        exo = pyexodus.exodus(file=filename, mode='r', array_type='numpy')

        # get node coordinates
        coordinates = exo.get_coords()
        Ndims  = len(coordinates)
        Nnodes = len(coordinates[0])
        vertices = np.zeros((Nnodes,Ndims))
        for i in range(Ndims):
            vertices[:,i] = coordinates[i]

        # get element connectivity
        elem_connectivity, Nelements, Nnodes_per_elem = exo.get_elem_connectivity(id=block_id)
        connectivity = np.zeros((Nelements,Nnodes_per_elem),dtype=int)
        for i in range(Nelements):
            for j in range(Nnodes_per_elem):
                connectivity[i,j] = elem_connectivity[i,j]-1

        if ((Ndims == 3) and (Nnodes_per_elem == 8)):
            return HexGeometry(vertices,connectivity)
        elif ((Ndims == 2) and (Nnodes_per_elem == 4)):
            return QuadGeometry(vertices,connectivity)

    def from_meshio(self,mesh):
        # Get the total number of elements
        Nelems = 0
        for cell_block in mesh.cells:
            if ((cell_block.type == "triangle") or (cell_block.type == "quad")):
                Nelems = Nelems + cell_block.data.shape[0]

        # Consolidate all (2d) elements into possibily degenerated quadrilaterals
        connectivity = np.zeros((Nelems,4), dtype=int)
        Nelems = 0
        for cell_block in mesh.cells:
            if (cell_block.type == "triangle"):
                for i, cell_connectivity in enumerate(cell_block.data):
                    connectivity[Nelems+i,:] = [ cell_connectivity[0], cell_connectivity[1], cell_connectivity[2], cell_connectivity[2] ]
                Nelems = Nelems + cell_block.data.shape[0]
            elif (cell_block.type == "quad"):
                for i, cell_connectivity in enumerate(cell_block.data):
                    connectivity[Nelems+i,:] = [ cell_connectivity[0], cell_connectivity[1], cell_connectivity[2], cell_connectivity[3] ]
                Nelems = Nelems + cell_block.data.shape[0]  

        # Get coordinates of all points
        coordinates = mesh.points[:,0:2]
        Nnodes = coordinates.shape[0]

        # Get the list of unique node IDs
        unique_nodes, unique_indices = np.unique(connectivity.flatten(), return_inverse=True)

        # Remove any nodes not appearing in the connectivity array, and renumber unique IDs
        coordinates  = coordinates[unique_nodes,:]
        connectivity = unique_indices.reshape((Nelems,4))
                
        return QuadGeometry(coordinates,connectivity)

        
