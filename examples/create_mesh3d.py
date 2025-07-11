# Append the location of the locally installed REMAT package to sys.path
import sys
sys.path.append("../install/package/")

# Load MeshUtils modules for geometry processing
from GeometryFactory import *
from Model import *

# Other needed Python packages
import os
import pyexodus
import numpy as np
import math

# create a simple mesh
points = np.array([[0,0,0],[0,1,0],[1,1,0],[1,0,0],[0,0,1],[0,1,1],[1,1,1],[1,0,1],[0,0,2],[0,1,2],[1,1,2],[1,0,2]])
connectivity = np.array([[0,1,2,3,4,5,6,7],[4,5,6,7,8,9,10,11]])

# create a new Exodus file
filename = 'cube.exo'
try:
    os.remove(filename)
except OSError:
    pass
exo = pyexodus.exodus(file=filename, mode='w', array_type='numpy', title=None, numDims=3, numNodes=points.shape[0], numElems=connectivity.shape[0], numBlocks=1, numNodeSets=0, numSideSets=0, io_size=8, compression=None)

# put node coordinates
exo.put_coords(xCoords=points[:,0],yCoords=points[:,1],zCoords=points[:,2])

# put element block 1 data
exo.put_elem_blk_info(id=1, elemType='HEX8', numElems=connectivity.shape[0], numNodesPerElem=connectivity.shape[1], numAttrsPerElem=0)

# put element block 1 connectivity (shift indicies by 1 to ensure compatibility with 1-based indexing in Exodus files, whereas Python uses 0-based indexing)
exo.put_elem_connectivity(id=1, connectivity=connectivity, shift_indices=1, chunk_size_in_mb=128)

# close the Exodus object (avoids error message)
exo.close()

# Read in data to geometry factory
geom_factory = GeometryFactory()
grid = geom_factory.cartesian_grid((0.0,0.0,0.0),(1.0,1.0,1.0),(2,2,2),0.1)
print(grid.geometry_id)

print(grid.vertices)
print(grid.connectivity)

cube = geom_factory.from_exodus('cube.exo')
print(cube.geometry_id)

bottom_nodes = cube.select_nodes(Select_Z_eq(0.0))
print(bottom_nodes.geometry_id)
print(bottom_nodes.vertex_ids)

top_nodes = cube.select_nodes(Select_Z_eq(2.0))
print(top_nodes.geometry_id)
print(top_nodes.vertex_ids)

side_faces = cube.select_faces(Select_X_eq(0.0))
print(side_faces.geometry_id)
print(side_faces.connectivity)

cube.mirror(np.array([0,0,0]),np.array([1.0,0.0,0.0]),copy=True)
cube.rotate(np.array([0,0,0]),np.array([0.0,0.0,1.0]),np.pi/4.0)
cube.translate(np.array([2,0,0]))
print(cube.geometry_id)
print(cube.vertices)
print(cube.connectivity)

# Export mesh for visualization
model = Model()
model.add_part(Part(grid,Material(None,None)))
model.add_part(Part(cube,Material(None,None)))
model.write_exodus('output.exo')
