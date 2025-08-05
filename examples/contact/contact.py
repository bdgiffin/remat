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
REMAT.API.define_parameter(b"initial_velocity_y", -2.0e-1)
REMAT.API.define_parameter(b"mass_damping_factor", 2.0e-1)
REMAT.API.define_parameter(b"contact_stiffness",   5.0e+0)
REMAT.API.define_parameter(b"search_radius",       1.0e+0)
REMAT.API.define_parameter(b"overflow_limit",      10000.0) # steps

# Define material parameters
REMAT.API.define_parameter(b"density",        1.0)
REMAT.API.define_parameter(b"youngs_modulus", 5.0)
REMAT.API.define_parameter(b"poissons_ratio", 0.28)

# Set the integrator type: "float" (default), or "fixed"
REMAT.API.set_integrator_type(b"fixed")

# Pre-process mesh/geometry ------------------------------------------------

# Generate and save problem geometry
if (not sys.platform == "emscripten"):
    # Create geometry factory
    geom_factory = GeometryFactory()

    with pygmsh.geo.Geometry() as geom:
        # create rectangular object
        obj = geom.add_rectangle(1.6+2, 2.6+2, 4.0, 5.0, 1.0, 0.45)
        #obj = geom.add_circle([2.5,4.0],1.0,mesh_size=2.0)
        geom.set_recombined_surfaces([obj.surface])
        rect = geom_factory.from_meshio(geom.generate_mesh(dim=2))
        #bottom_nodes = rect.select_nodes(Select_Y_eq(4.0))
        bottom_nodes = rect.select_boundary_nodes()
        all_nodes    = rect.select_nodes(Select_all())

    with pygmsh.geo.Geometry() as geom:
        # create polygonal object
        obj = geom.add_polygon(
            [
                [1.0+2, 1.0],
                [3.0+2, 1.0],
                [3.0+2, 3.0],
                [1.0+2, 3.0],
            ],
            mesh_size=0.35,
        )
        geom.set_recombined_surfaces([obj.surface])
        poly = geom_factory.from_meshio(geom.generate_mesh(dim=2))
        #top_faces = poly.select_faces(Select_Y_eq(3.0))
        top_faces = poly.select_boundary_faces()

    # Create Model object
    model = Model()
    model.add_part(Part(rect,Material(None,None)))
    model.add_part(Part(poly,Material(None,None)))
    model.add_initial_condition(all_nodes,[+3.0e-1,-2.0e-1])
    model.add_contact_interaction(bottom_nodes,top_faces)
    
    # Generate the problem data from the Model object and save in a pickle file
    pickle.dump(model.generate_problem(),open("problem.pkl",'wb'))

# Load the problem data from pickle file
coordinates, velocities, fixity, connectivity, contacts, truss_connectivity = pickle.load(open("problem.pkl",'rb'))
        
# Define the problem geometry
REMAT.create_geometry(coordinates,velocities,fixity,connectivity,contacts,truss_connectivity)

# Run analysis -------------------------------------------------------------

# Define "main" animation loop
Nsteps = 300
Nsub_steps = 100
dt = 0.25e-2 # [s] time increment
anim = Animation(Nsteps,Nsub_steps,dt,element_color="pressure",edge_color="Black")
recording = "" # "output_animation"
asyncio.run(anim.start(recording))

# --------------------------------------------------------------------------

