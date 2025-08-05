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

# Call C/C++ library API functions from Python:

# Define global parameters
REMAT.API.define_parameter(b"body_force_y",       -1.0e-1)
REMAT.API.define_parameter(b"initial_velocity_x", +0.0e-1)
REMAT.API.define_parameter(b"initial_velocity_y", -0.0e-1)
REMAT.API.define_parameter(b"mass_damping_factor", 4.0e-2)
REMAT.API.define_parameter(b"contact_stiffness",   5.0e+0)
REMAT.API.define_parameter(b"search_radius",       1.0e+0)

# Define material parameters
REMAT.API.define_parameter(b"density",        1.0)
REMAT.API.define_parameter(b"youngs_modulus", 50.0)
REMAT.API.define_parameter(b"poissons_ratio", 0.28)

# Set the integrator type: "float" (default), or "fixed"
REMAT.API.set_integrator_type(b"fixed")

# Pre-process mesh/geometry ------------------------------------------------

# Generate and save problem geometry
if (not sys.platform == "emscripten"):
    # Create geometry factory
    geom_factory = GeometryFactory()

    # create block
    block = geom_factory.cartesian_grid([1.0,4.5],[2.0,5.5],[4,4],0.0)
    boundary_nodes = block.select_boundary_nodes()

    with pygmsh.geo.Geometry() as geom:
        lcar = 0.5
        p1 = geom.add_point([0.0,  5.0], lcar)
        c2 = geom.add_point([4.0,  0.25], lcar)
        p2 = geom.add_point([0.0,  0.0], lcar)
        p3 = geom.add_point([10.0, 0.0], lcar)
        c3 = geom.add_point([6.0,  0.25], lcar)
        p4 = geom.add_point([10.0, 5.0], lcar)
        s1 = geom.add_bspline([p4, c3, c2, p1])

        s2 = geom.add_line(p3, p4)
        s3 = geom.add_line(p2, p3)
        s4 = geom.add_line(p1, p2)

        ll = geom.add_curve_loop([s4, s3, s2, s1])
        pl = geom.add_plane_surface(ll)

        geom.set_recombined_surfaces([pl])
        well = geom_factory.from_meshio(geom.generate_mesh(dim=2))

    # create half-space
    well_faces   = well.select_boundary_faces()
    bottom_nodes = well.select_nodes(Select_Y_eq(0.0))
    left_nodes   = well.select_nodes(Select_X_eq(0.0))
    right_nodes  = well.select_nodes(Select_X_eq(10.0))

    # Create Model object
    model = Model()
    model.add_part(Part(block,Material(None,None)))
    model.add_part(Part(well,Material(None,None)))
    model.add_boundary_condition(bottom_nodes,[True,True]) # fully fix the bottom surface ndoes
    model.add_boundary_condition(left_nodes,[True,True]) # fully fix the left surface ndoes
    model.add_boundary_condition(right_nodes,[True,True]) # fully fix the right surface ndoes
    model.add_contact_interaction(boundary_nodes,well_faces)

    # Generate the problem data from the Model object and save in a pickle file
    pickle.dump(model.generate_problem(),open("problem.pkl",'wb'))

# Load the problem data from pickle file
coordinates, velocities, fixity, connectivity, contacts, truss_connectivity = pickle.load(open("problem.pkl",'rb'))

# Define the problem geometry
REMAT.create_geometry(coordinates,velocities,fixity,connectivity,contacts,truss_connectivity)

# Run analysis -------------------------------------------------------------

# Define "main" animation loop
Nsteps = 220
Nsub_steps = 80
dt = 1.0e-2 # [s] time increment
anim = Animation(Nsteps,Nsub_steps,dt,display_energy=True,element_color="Grey",edge_color="Dark Grey")
recording = "" #"output_animation"
asyncio.run(anim.start(recording))

# --------------------------------------------------------------------------

