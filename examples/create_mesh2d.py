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
