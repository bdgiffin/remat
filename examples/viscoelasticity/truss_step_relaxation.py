# A unit-length bar is fixed on the left and receives a 1% step strain on the right node.
# Stress relaxes as the Maxwell dashpot flows, and reversing time recovers the
# initial configuration bit-for-bit.


from math import *
import sys
import time as timer

import numpy as np

# Append the location of the locally installed REMAT package to sys.path
exodus_available = False
if not sys.platform == "emscripten":
    sys.path.append("../../install/package/")
    from ExodusIO import ExodusIO
    exodus_available = True
else:
    exodus_available = False

import REMAT

# --------------------------------------------------------------------------

# Problem constants
gauge_length = 1.0
epsilon0 = 1.0e-2

# Define global parameters
REMAT.API.define_parameter(b"body_force_x", 0.0)
REMAT.API.define_parameter(b"body_force_y", 0.0)
REMAT.API.define_parameter(b"mass_damping_factor", 0.0)
REMAT.API.define_parameter(b"dt_scale_factor", 1.0)

# Define material parameters
REMAT.API.define_parameter(b"density", 1.0)
REMAT.API.define_parameter(b"youngs_modulus", 1.0)
REMAT.API.define_parameter(b"poissons_ratio", 0.25)
REMAT.API.define_parameter(b"truss_density", 1.0)
REMAT.API.define_parameter(b"truss_youngs_modulus", 1.0)
REMAT.API.define_parameter(b"area", 1.0)
REMAT.API.define_parameter(b"viscosity", 5.0e-1)
REMAT.API.define_parameter(b"mat_overflow_limit", 100.0)

# Set the integrator type 
REMAT.API.set_integrator_type(b"fixed_truss_visco")

# --------------------------------------------------------------------------

# Two-node truss geometry (aligned with global x-axis)
coordinates =   np.array([[0.0, 0.0],[gauge_length, 0.0]], dtype=np.double)
velocities =    np.array([[0.0, 0.0],[0.0, 0.0]], dtype=np.double)
fixity =        np.array([[0.0, 0.0],[0.0, 0.0]], dtype=np.bool_)
fixity[0, :] = True          # left node fully fixed
fixity[1, 1] = True          # right node vertical DOF fixed

connectivity = np.zeros((0, 4), dtype=np.int32)  
contacts = []                                    
truss_connectivity = np.array([[0, 1]], dtype=np.int32)

# Define the problem geometry
REMAT.create_geometry(coordinates,
                      velocities,
                      fixity,
                      connectivity,
                      contacts,
                      truss_connectivity)

# --------------------------------------------------------------------------

# Step strain displacement boundary condition for the right node
left_x = coordinates[0, 0]


def right_node_step(time, x, y):
    if time == 0.0:
        return 0.0
    else:
        return epsilon0 * (x - left_x)


REMAT.define_displacement_bc(np.array([1], dtype=np.int32), 0, right_node_step)

REMAT.API.initialize()

exo = None
if exodus_available:
    exo = ExodusIO()
    exo.create("truss_step_relaxation.exo")
    exo.output_state()





# Run analysis -------------------------------------------------------------

dt = 1.0e-3
Nsteps = 40
Nsub_steps = 100

# Forward integration
for step in range(1, Nsteps + 1):
    REMAT.API.update_state(+dt, Nsub_steps)
    if exo:
        exo.output_state()

# Backward integration
for step in range(Nsteps, 0, -1):
    REMAT.API.update_state(-dt, Nsub_steps)
    if exo:
        exo.output_state()

if exo:
    exo.output_state()
    exo.finalize()


# --------------------------------------------------------------------------
# Quick check of node and truss state arrays

# num_nodes = REMAT.API.get_num_entities(b"node")
# node_fields = REMAT.API.get_num_fields(b"node")
# node_state = np.zeros((num_nodes, node_fields), dtype=np.double)
# REMAT.API.get_fields(b"node", node_state)
# node_state_copy = node_state.copy()
# print("node state:")
# print(node_state_copy)

# num_truss = REMAT.API.get_num_entities(b"truss")
# truss_fields = REMAT.API.get_num_fields(b"truss")
# truss_state = np.zeros((num_truss, truss_fields), dtype=np.double)
# REMAT.API.get_fields(b"truss", truss_state)
# truss_state_copy = truss_state.copy()
# print("truss state:")
# print(truss_state_copy)
# --------------------------------------------------------------------------