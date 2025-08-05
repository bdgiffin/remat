# /// script
# dependencies = [
#  "pygame_widgets"
# ]
# ///

# Other needed Python packages
import numpy as np
from math import *
import random

# Import the "asyncio" Python package
import asyncio

# Append the location of the locally installed REMAT package to sys.path
import sys, platform
if (not sys.platform == "emscripten"):
    sys.path.append("../../install/package/")

# Load the REMAT and MeshUtils packages
import REMAT
from Animation import *

# --------------------------------------------------------------------------

# Define global parameters
REMAT.API.define_parameter(b"body_force_y",       -2.0e-2)
REMAT.API.define_parameter(b"initial_velocity_x", +0.5e-1)
REMAT.API.define_parameter(b"initial_velocity_y", -0.0e-1)
REMAT.API.define_parameter(b"mass_damping_factor", 0.0e-2)
REMAT.API.define_parameter(b"contact_stiffness",   0.0e+0)

# Define material parameters
REMAT.API.define_parameter(b"density",         1.0)
REMAT.API.define_parameter(b"youngs_modulus", 1000.0)
REMAT.API.define_parameter(b"yield_stress",   1.0e+9)
REMAT.API.define_parameter(b"viscosity",       1.0)
REMAT.API.define_parameter(b"area",            1.0)
REMAT.API.define_parameter(b"poissons_ratio", 0.28)

# Set the integrator type: "float" (default), "fixed", or "mixed"
REMAT.API.set_integrator_type(b"fixed")

# Pre-process mesh/geometry ------------------------------------------------
    
# Manually define nodal coordinates and element geometry
length = 3.5
Nnodes = 4
Nelems = 0
Ntruss = Nnodes-1
coordinates  = np.zeros((Nnodes,2))
velocities   = np.zeros((Nnodes,2))
fixity       = np.zeros((Nnodes,2),dtype=np.bool_)
truss_connectivity = np.zeros((Ntruss,2),dtype=np.int32)
random.seed(8)
for i in range(0,Nnodes):
    coordinates[i,:] = [5.0+0.001*random.random(), 3.0+(length/Nnodes)*i]
    if (i == 0):
        fixity[i,:] = [True,True]
    else:
        truss_connectivity[i-1,:] = [i-1,i]

# Define the problem geometry
connectivity = np.empty((0,0),dtype=np.int32)
contacts = []
REMAT.create_geometry(coordinates,velocities,fixity,connectivity,contacts,truss_connectivity)

# Define point masses
Npoints = 1
point_ids = np.array([Nnodes-1],dtype=np.int32)
point_mass = np.array([3.0])
REMAT.API.define_point_mass(point_ids,point_mass,Npoints)

# Run analysis -------------------------------------------------------------

# Define "main" animation loop
Nsteps = 400
Nsub_steps = 50
dt = 1.0e-2 # [s] time increment
anim = Animation(Nsteps,Nsub_steps,dt)
recording = ""
asyncio.run(anim.start(recording))

# --------------------------------------------------------------------------
