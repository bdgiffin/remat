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
REMAT.API.define_parameter(b"body_force_y",       -0.0e-1)
REMAT.API.define_parameter(b"initial_velocity_x", +0.0e-1)
REMAT.API.define_parameter(b"initial_velocity_y", -0.0e-1)
REMAT.API.define_parameter(b"mass_damping_factor", 0.0e-1)
REMAT.API.define_parameter(b"contact_stiffness",   5.0e+0)
REMAT.API.define_parameter(b"search_radius",       1.0e+0)

# Define material parameters
REMAT.API.define_parameter(b"density",        1.0)
REMAT.API.define_parameter(b"youngs_modulus", 5.0)
REMAT.API.define_parameter(b"poissons_ratio", 0.28)

# Set the integrator type: "float" (default), or "fixed"
REMAT.API.set_integrator_type(b"float")

# Pre-process mesh/geometry ------------------------------------------------

# Generate and save problem geometry
if (not sys.platform == "emscripten"):
    # Create geometry factory
    geom_factory = GeometryFactory()

    # create impactor
    with pygmsh.geo.Geometry() as geom:
        obj = geom.add_circle([4.0,4.0],0.75,mesh_size=0.25,compound=True,num_sections=4,make_surface=True)
        print(obj)
        #geom.set_recombined_surfaces([obj.surface])
        impactor = geom_factory.from_meshio(geom.generate_mesh(dim=2))
        boundary_nodes = impactor.select_boundary_nodes()
        impactor_nodes = impactor.select_nodes(Select_all())

    # create half-space
    grid = geom_factory.cartesian_grid([0.0,0.0],[10.0,3.0],[200,60])
    top_faces = grid.select_faces(Select_Y_eq(3.0))
    bottom_nodes = grid.select_nodes(Select_Y_eq(0.0))
    left_nodes = grid.select_nodes(Select_X_eq(0.0))
    right_nodes = grid.select_nodes(Select_X_eq(10.0))

    # Create Model object
    model = Model()
    model.add_part(Part(impactor,Material(None,None)))
    model.add_part(Part(grid,Material(None,None)))
    model.add_initial_condition(impactor_nodes,[+2.0e-1,-4.0e-1])
    model.add_boundary_condition(bottom_nodes,[True,True]) # fully fix the bottom surface ndoes
    model.add_boundary_condition(left_nodes,[True,True]) # fully fix the left surface ndoes
    model.add_boundary_condition(right_nodes,[True,True]) # fully fix the right surface ndoes
    model.add_contact_interaction(boundary_nodes,top_faces)

    # Generate the problem data from the Model object and save in a pickle file
    pickle.dump(model.generate_problem(),open("problem.pkl",'wb'))

# Load the problem data from pickle file
coordinates, velocities, fixity, connectivity, contacts, truss_connectivity = pickle.load(open("problem.pkl",'rb'))
        
# Define the problem geometry
REMAT.create_geometry(coordinates,velocities,fixity,connectivity,contacts,truss_connectivity)

# Run analysis -------------------------------------------------------------

# Define "main" animation loop
Nsteps = 1000
Nsub_steps = 10
dt = 1.0e-2 # [s] time increment
anim = Animation(Nsteps,Nsub_steps,dt,element_color="pressure",edge_color="Black")
recording = "" #"output_animation"
asyncio.run(anim.start(recording))

# --------------------------------------------------------------------------

