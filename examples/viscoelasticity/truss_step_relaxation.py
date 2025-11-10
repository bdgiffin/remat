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
epsilon0 = 1.0e-2 # Suggested for right_node_step and right_node_constant_rate
epsilon0 = 0.05     # Suggested for right_node_sinusoidal and right_node_clipped_sinusoid


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
REMAT.API.define_parameter(b"viscosity", 10e-1)
REMAT.API.define_parameter(b"mat_overflow_limit", 30.0)
# REMAT.API.define_parameter(b"mat_overflow_limit", 3000000000000.0)

# Set the integrator type 
REMAT.API.set_integrator_type(b"fixed_truss_visco")
# REMAT.API.set_integrator_type(b"float_truss_visco")

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
# Define displacement boundary conditions
left_x = coordinates[0, 0]

# Step strain
def right_node_step(time, x, y):
    if time == 0.0:
        return 0.0
    else:
        eps_t = epsilon0
        return eps_t * (x - left_x)


# Constant strain rate
def right_node_constant_rate(time, x, y):
    R = epsilon0
    eps_t = R * time
    return eps_t * (x - left_x)


# Sinusoidal cyclic loading (DMA)
def right_node_sinusoidal(time, x, y):
    A = epsilon0
    omega = 2.0 * pi  # 1 Hz by default
    eps_t = A * sin(omega * time)
    return eps_t * (x - left_x)


# Clipped sinusoid
def right_node_clipped_sinusoid(time, x, y):
    C = 0.0
    A = epsilon0
    omega = 2.0 * pi
    Amin = -0.5 * epsilon0
    Amax = +0.5 * epsilon0
    raw = A * sin(omega * time)
    clipped = raw
    if clipped < Amin:
        clipped = Amin
    elif clipped > Amax:
        clipped = Amax
    eps_t = C + clipped
    return eps_t * (x - left_x)


REMAT.define_displacement_bc(np.array([1], dtype=np.int32), 0, right_node_sinusoidal)

REMAT.API.initialize()

exo = None
if exodus_available:
    exo = ExodusIO()
    exo.create("truss_step_relaxation.exo")
    exo.output_state()


# Run analysis -------------------------------------------------------------

dt = 1.0e-3
Nsteps = 550
Nsub_steps = 10

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