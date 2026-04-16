# /// script
# dependencies = [
#  "pygame_widgets"
# ]
# ///
import sys
sys.path.append("../../install/package")

from math import *
import numpy as np
import time as timer

import asyncio

from GeometryFactory import *
from Model import *

import sys, platform
if (not sys.platform == "emscripten"):
    sys.path.append("../../install/package/")
    from ExodusIO import *

import REMAT
from Animation import *

# --------------------------------------------------------------------------
# Case-3-style forward/backward replay only (no adjoint pass)
# --------------------------------------------------------------------------

# Global parameters
REMAT.API.define_parameter(b"body_force_y",       -0.0)
REMAT.API.define_parameter(b"initial_velocity_x", +0.0)
REMAT.API.define_parameter(b"initial_velocity_y", -0.0)
REMAT.API.define_parameter(b"contact_stiffness",   0.0)
REMAT.API.define_parameter(b"search_radius",       1.0)
REMAT.API.define_parameter(b"overflow_limit",      10.0)

# Material parameters (case-3 baseline)
REMAT.API.define_parameter(b"density",        1.0)
REMAT.API.define_parameter(b"youngs_modulus", 5.0)
REMAT.API.define_parameter(b"poissons_ratio", 0.28)

# Viscous parameters (case-3 true tau)
REMAT.API.define_parameter(b"relaxation_time", 0.35)
REMAT.API.define_parameter(b"shear_modulus_Maxwell_element", 2.0)

# Integrator and material overflow settings requested
REMAT.API.set_integrator_type(b"fixed_visco")
REMAT.API.define_parameter(b"mat_overflow_limit", 10.0)

# --------------------------------------------------------------------------
# Geometry / mesh (case-3 settings)
# --------------------------------------------------------------------------

width = 10.0
height = 3.0

Nx = 88
Ny = 32
if (sys.platform == "emscripten"):
    Nx = int(Nx/2)
    Ny = int(Ny/2)

# Case-3 source setup
source_window_fraction = 0.09
impact_velocity = 1.05
source_half_width = 0.5*source_window_fraction*width
source_window_y = 0.1
xmid = 0.5*width

geom_factory = GeometryFactory()
grid = geom_factory.cartesian_grid([0.0,0.0],[width,height],[Nx,Ny])

# Source nodes near the top-center window
source_nodes = grid.select_nodes(
    Select_XY_window(
        [xmid-source_half_width, height-source_window_y],
        [xmid+source_half_width, height+source_window_y]
    )
)
source_node_ids = np.asarray(source_nodes.global_node_ids(),dtype=np.int32).reshape(-1)
if source_node_ids.size == 0:
    raise RuntimeError("No source nodes selected for excitation; broaden source window.")

bottom_nodes = grid.select_nodes(Select_Y_eq(0.0))
left_nodes   = grid.select_nodes(Select_X_eq(0.0))
right_nodes  = grid.select_nodes(Select_X_eq(width))

model = Model()
model.add_part(Part(grid,Material(None,None)))
model.add_initial_condition(source_nodes,[0.0,-impact_velocity])
model.add_boundary_condition(bottom_nodes,[True,True])
model.add_boundary_condition(left_nodes,[True,True])
model.add_boundary_condition(right_nodes,[True,True])

coordinates, velocities, fixity, connectivity, contacts, truss_connectivity = model.generate_problem()

# Debug sanity checks for excitation plumbing
source_speed = np.linalg.norm(velocities[source_node_ids,:],axis=1)
if np.max(source_speed) <= 0.0:
    raise RuntimeError("Source excitation is zero after model.generate_problem().")
if np.all(fixity[source_node_ids,1]):
    raise RuntimeError("All source nodes are fixed in Y; wave cannot be excited.")

REMAT.create_geometry(coordinates,velocities,fixity,connectivity,contacts,truss_connectivity)

# Case-3 layered stiffness profile (bottom -> top)
layer_bounds = np.linspace(0.0,height,5)  # 4 layers
layer_values = np.array([0.52,1.82,0.60,1.70],dtype=np.double)

def layered_stiffness_scaling(_, y):
    if y < layer_bounds[1]:
        return float(layer_values[0])
    elif y < layer_bounds[2]:
        return float(layer_values[1])
    elif y < layer_bounds[3]:
        return float(layer_values[2])
    return float(layer_values[3])

REMAT.define_variable_properties(layered_stiffness_scaling)

# --------------------------------------------------------------------------
# Time integration (forward + backward only)
# --------------------------------------------------------------------------

dt = 9e-3
Nsteps = 50
Nsub_steps = 10
step_id = 0

if (not sys.platform == "emscripten"):
    print(f"Selected source nodes: {source_node_ids.size}")
    print(f"Source speed max: {float(np.max(source_speed)):.6e}")
    print(f"Mesh resolution: Nx={Nx}, Ny={Ny}, dt={dt}, Nsteps={Nsteps}, Nsub_steps={Nsub_steps}")

    exo = ExodusIO()
    exo.create("dissipative_wave_forward_backward_only.exo")

    REMAT.API.initialize()
    exo.output_state()

    start_time = timer.time()

    # Forward pass
    while (step_id < Nsteps):
        step_id += 1
        REMAT.API.update_state(dt,Nsub_steps,REMAT.PASS_FORWARD)
        exo.output_state()

    # Backward replay pass (no adjoint)
    while (step_id > 0):
        step_id -= 1
        REMAT.API.update_state(dt,Nsub_steps,REMAT.PASS_BACKWARD)
        exo.output_state()

    exo.finalize()

    end_time = timer.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.4f} seconds")
    print("Wrote: dissipative_wave_forward_backward_only.exo")

else:
    # Optional browser animation branch (forward-only visualization)
    anim = Animation(Nsteps,Nsub_steps,dt,element_color="viscous_strain_yy",element_field_max=0.01)
    asyncio.run(anim.start())

# --------------------------------------------------------------------------
