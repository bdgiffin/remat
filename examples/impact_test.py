# Append the location of the locally installed REMAT package to sys.path
import sys
sys.path.append("../install/package/")

# Load the REMAT package
import REMAT

# Python package for reading/writing data in the Exodus mesh database format
# NOTE: PYEXODUS V0.1.5 NEEDS TO BE MODIFIED TO WORK CORRECTLY WITH PYTHON 3.12
# https://pypi.org/project/pyexodus/
import pyexodus

# Other needed Python packages
from math import *
import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------

# Call C/C++ library API functions from Python:

# Define global parameters
REMAT.API.define_parameter(b"dt_scale_factor",     1.0e-5)
REMAT.API.define_parameter(b"body_force_y",       -1.0e-5)
REMAT.API.define_parameter(b"initial_velocity_y", -1.0e-5)
REMAT.API.define_parameter(b"mass_damping_factor", 1.0e-9)

# Define material parameters
REMAT.API.define_parameter(b"density",        0.5)
REMAT.API.define_parameter(b"youngs_modulus", 215.0e+1)
REMAT.API.define_parameter(b"poissons_ratio", 0.28)

# Pre-process mesh/geometry ------------------------------------------------

Nnodes = 9
Nelems = 4
coordinates = np.array([ 0.0, 0.0,
                         1.0, 0.0,
                         2.0, 0.0,
		         0.0, 1.0,
                         1.0, 1.0,
                         2.0, 1.0,
		 	 0.0, 2.0,
                         1.0, 2.0,
                         2.0, 2.0 ])
connectivity = np.array([[ 0, 1, 4, 3 ],
                         [ 1, 2, 5, 4 ],
                         [ 3, 4, 7, 6 ],
                         [ 4, 5, 8, 7 ]], dtype=np.int32);

# Define the problem geometry
filename = "output.exo"
REMAT.create_geometry(coordinates,connectivity,Nnodes,Nelems,filename)

# Run analysis -------------------------------------------------------------

# initialization
time = 0.0 # [s] starting time
REMAT.API.initialize(time)
REMAT.output_state()

# perform the analysis
dt = 1.0e-7 # [s] time increment
Nsteps = 10
for step_id in range(1,Nsteps):
    time = time + dt
    REMAT.API.update_state(time,dt)
    REMAT.output_state()

# finalize the REMAT module (close the Exodus files)
REMAT.finalize()

# --------------------------------------------------------------------------
