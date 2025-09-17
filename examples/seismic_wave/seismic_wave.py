# /// script
# dependencies = [
#  "pygame_widgets"
# ]
# ///
import sys
sys.path.append("../../install/package")

# Other needed Python packages
from math import *
import numpy as np
import pickle
import time as timer

# Import the "asyncio" Python package
import asyncio

# Load MeshUtils packages
from GeometryFactory import *
from Model import *

# Append the location of the locally installed REMAT package to sys.path
import sys, platform
if (not sys.platform == "emscripten"):
    sys.path.append("../../install/package/")

    # Import ExodusIO for output
    from ExodusIO import *

# Load the REMAT and MeshUtils packages
import REMAT
from Animation import *

# --------------------------------------------------------------------------

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

# Define viscous parameters
REMAT.API.define_parameter(b"relaxation_time", 0.1)

# Set the integrator type: "float" (default), or "fixed"
REMAT.API.set_integrator_type(b"float")

# Pre-process mesh/geometry ------------------------------------------------

Nx = 800
Ny = 240
if (sys.platform == "emscripten"):
    Nx = int(Nx/4)
    Ny = int(Ny/4)

# Create geometry factory
geom_factory = GeometryFactory()

# create half-space
grid = geom_factory.cartesian_grid([0.0,0.0],[10.0,3.0],[Nx,Ny])
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
coordinates, velocities, fixity, connectivity, contacts, truss_connectivity = model.generate_problem()
        
# Define the problem geometry
REMAT.create_geometry(coordinates,velocities,fixity,connectivity,contacts,truss_connectivity)

# Define variable material properties
REMAT.define_variable_properties(lambda x, y: 1.0 - 0.5*(y > 0.5 + 0.07*(x-5.0)) + 0.4*(y > 1.0 - 0.07**(x-5.0)) - 0.5*(y > 1.5 + 0.1*(x-5.0)) + 0.25*(y > 2.0 + 0.05*(x-5.0)) + 0.5*(y > 2.5 - 0.02*(x-5.0)))

# Run analysis -------------------------------------------------------------

# set analysis time-stepping parameters
dt = 0.25e-2 # [s] time increment
step_id = 0
Nsteps = 100
Nsub_steps = 10

if (not sys.platform == "emscripten"):
    # Create output Exodus file
    exo = ExodusIO()
    exo.create("seismic_wave_test.exo")

    # initialization
    time = 0.0 # [s] starting time
    REMAT.API.initialize()
    exo.output_state()

    start_time = timer.time()

    # perform the forward-in-time analysis
    while (step_id < Nsteps):
        step_id = step_id + 1
        time = REMAT.API.update_state(+dt,Nsub_steps)
        exo.output_state()

    # perform the time-reversed analysis
    while (step_id > 0):
        step_id = step_id - 1
        time = REMAT.API.update_state(-dt,Nsub_steps)
        exo.output_state()

    exo.finalize()

    end_time = timer.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.4f} seconds")

else:

    anim = Animation(Nsteps,Nsub_steps,dt,element_color="pressure",element_field_max=0.1)
    asyncio.run(anim.start())

# --------------------------------------------------------------------------

