# /// script
# dependencies = [
#  "pygame_widgets"
# ]
# ///

# Other needed Python packages
from math import *
import numpy as np
import pickle

# Import the "asyncio" Python package
import asyncio

# Append the location of the locally installed REMAT package to sys.path
import sys, platform
if (not sys.platform == "emscripten"):
    sys.path.append("../../install/package/")

    # Load MeshUtils packages
    from GeometryFactory import *
    from Model import *

    # Import GMSH for meshing
    import pygmsh

# Load the REMAT and MeshUtils packages
import REMAT
from Animation import *

# --------------------------------------------------------------------------

# Define global parameters
REMAT.API.define_parameter(b"body_force_y",       -1.0e-2)
REMAT.API.define_parameter(b"initial_velocity_x", +1.5e-1)
REMAT.API.define_parameter(b"initial_velocity_y", -0.0e-1)
REMAT.API.define_parameter(b"mass_damping_factor", 2.0e-2)
REMAT.API.define_parameter(b"contact_stiffness",   5.0e+0)
REMAT.API.define_parameter(b"search_radius",       1.0e+0)
REMAT.API.define_parameter(b"overflow_limit",      50000.0) # steps

# Define material parameters
REMAT.API.define_parameter(b"density",        1.0)
REMAT.API.define_parameter(b"youngs_modulus", 5.0)
REMAT.API.define_parameter(b"poissons_ratio", 0.28)

# Set the integrator type: "float" (default), or "fixed"
REMAT.API.set_integrator_type(b"fixed")

# Pre-process mesh/geometry ------------------------------------------------

# Manually define point
Nnodes = 1
Nelems = 0
Ntruss = 0
coordinates  = np.zeros((Nnodes,2))
velocities   = np.zeros((Nnodes,2))
fixity       = np.zeros((Nnodes,2),dtype=np.bool_)
for i in range(0,Nnodes):
    coordinates[i,:] = [1.0, 4.0]

# Define the problem geometry
connectivity = np.empty((0,0),dtype=np.int32)
truss_connectivity = np.empty((0,0),dtype=np.int32)
contacts = []
REMAT.create_geometry(coordinates,velocities,fixity,connectivity,contacts,truss_connectivity)

# Define point masses
Npoints = 1
point_ids = np.array([Nnodes-1],dtype=np.int32)
point_mass = np.array([3.0])
REMAT.API.define_point_mass(point_ids,point_mass,Npoints)

# Run analysis -------------------------------------------------------------

# Define "main" animation loop
Nsteps = 300
Nsub_steps = 100
dt = 1.0e-2 # [s] time increment
anim = Animation(Nsteps,Nsub_steps,dt,element_color="Dark Red")
recording = "output_animation"
asyncio.run(anim.start(recording))

# --------------------------------------------------------------------------

