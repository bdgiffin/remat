# Append the location of the locally installed REMAT package to sys.path
import sys
sys.path.append("../install/package/")

# Load the REMAT and MeshUtils packages
import REMAT
from GeometryFactory import *
from Model import *

# Python package for reading/writing data in the Exodus mesh database format
# NOTE: PYEXODUS V0.1.5 NEEDS TO BE MODIFIED TO WORK CORRECTLY WITH PYTHON 3.12
# https://pypi.org/project/pyexodus/
import pyexodus

# Other needed Python packages
from math import *
import numpy as np
import time as timer

# Import GMSH for meshing
import pygmsh

# --------------------------------------------------------------------------

# Call C/C++ library API functions from Python:

# Define global parameters
REMAT.API.define_parameter(b"body_force_y",       -0.0e-1)
REMAT.API.define_parameter(b"initial_velocity_x", +0.0e-1)
REMAT.API.define_parameter(b"initial_velocity_y", -0.0e-1)
REMAT.API.define_parameter(b"mass_damping_factor", 1.0e-0)
REMAT.API.define_parameter(b"contact_stiffness",   0.0e+0)
REMAT.API.define_parameter(b"search_radius",       1.0e+0)
REMAT.API.define_parameter(b"overflow_limit",      1001.0) # steps

# Define material parameters
REMAT.API.define_parameter(b"density",        1.0)
REMAT.API.define_parameter(b"youngs_modulus", 5.0)
REMAT.API.define_parameter(b"poissons_ratio", 0.28)

# Set the integrator type: "float" (default), or "fixed"
REMAT.API.set_integrator_type(b"float")

# Pre-process mesh/geometry ------------------------------------------------

# Create geometry factory
geom_factory = GeometryFactory()

# create half-space
grid = geom_factory.cartesian_grid([0.0,0.0],[10.0,3.0],[800,240])
window_size  = 0.1
top_nodes    = grid.select_faces(Select_XY_window([5.0-window_size,3.0-window_size],[5.0+window_size,3.0+window_size]))
bottom_nodes = grid.select_nodes(Select_Y_eq(0.0))
left_nodes   = grid.select_nodes(Select_X_eq(0.0))
right_nodes  = grid.select_nodes(Select_X_eq(10.0))

# Create Model object
model = Model()
model.add_part(Part(grid,Material(None,None)))
model.add_initial_condition(top_nodes,[0.0,-4.0e-1])
model.add_boundary_condition(bottom_nodes,[True,True]) # fully fix the bottom surface ndoes
model.add_boundary_condition(left_nodes,[True,True]) # fully fix the left surface ndoes
model.add_boundary_condition(right_nodes,[True,True]) # fully fix the right surface ndoes

# Generate the problem data from the Model object
coordinates, velocities, fixity, connectivity, contacts, _ = model.generate_problem()
Nnodes = coordinates.shape[0]
Nelems = connectivity.shape[0]
        
# Define the problem geometry
filename = "seismic_wave_float.exo"
REMAT.create_geometry(coordinates,velocities,fixity,connectivity,contacts,filename)

# Define variable material properties
REMAT.define_variable_properties(lambda x, y: 1.0 - 0.5*(y > 0.5 + 0.07*(x-5.0)) + 0.4*(y > 1.0 - 0.07**(x-5.0)) - 0.5*(y > 1.5 + 0.1*(x-5.0)) + 0.25*(y > 2.0 + 0.05*(x-5.0)) + 0.5*(y > 2.5 - 0.02*(x-5.0)))

# Run analysis -------------------------------------------------------------

# initialization
time = 0.0 # [s] starting time
REMAT.API.initialize()
REMAT.output_state()

# set analysis time-stepping parameters
dt = 0.25e-2 # [s] time increment
step_id = 0
Nsteps = 100
Nsub_steps = 10

start_time = timer.time()

# perform the reversible analysis
while (step_id < Nsteps):
    step_id = step_id + 1
    time = REMAT.API.update_state(+dt,Nsub_steps)
    REMAT.output_state()
    
while (step_id > 0):
    step_id = step_id - 1
    time = REMAT.API.update_state(-dt,Nsub_steps)
    REMAT.output_state()

end_time = timer.time()
elapsed_time = end_time - start_time
print(f"Elapsed time: {elapsed_time:.4f} seconds")

REMAT.finalize()

# --------------------------------------------------------------------------

