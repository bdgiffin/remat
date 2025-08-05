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
REMAT.API.define_parameter(b"body_force_y",       -1.0e-1)
REMAT.API.define_parameter(b"initial_velocity_x", +3.0e-1)
REMAT.API.define_parameter(b"initial_velocity_y", -3.0e-1)
REMAT.API.define_parameter(b"mass_damping_factor", 1.0e-1)
REMAT.API.define_parameter(b"contact_stiffness",   2.0e+0)

# Define material parameters
REMAT.API.define_parameter(b"density",        1.0)
REMAT.API.define_parameter(b"youngs_modulus", 1.0)
REMAT.API.define_parameter(b"poissons_ratio", 0.28)

# Pre-process mesh/geometry ------------------------------------------------

# Generate and save problem geometry
if (not sys.platform == "emscripten"):
    # Create geometry factory
    geom_factory = GeometryFactory()

    with pygmsh.geo.Geometry() as geom:
            #obj = geom.add_rectangle(1.0, 2.0, 1.0, 2.0, 1.0, 0.1)
            obj = geom.add_polygon(
                [
                    [1.0, 1.0],
                    [5.0, 1.2],
                    [5.1, 5.2],
                    [1.1, 4.0],
                ],
                mesh_size=0.5,
            )
            geom.set_recombined_surfaces([obj.surface])
            mesh = geom_factory.from_meshio(geom.generate_mesh(dim=2))

    # Create Model object
    model = Model()
    model.add_part(Part(mesh,Material(None,None)))

    # Generate the problem data from the Model object and save in a pickle file
    pickle.dump(model.generate_problem(),open("problem.pkl",'wb'))

# Load the problem data from pickle file
coordinates, velocities, fixity, connectivity, contacts, truss_connectivity = pickle.load(open("problem.pkl",'rb'))

# Define the problem geometry
REMAT.create_geometry(coordinates,velocities,fixity,connectivity,contacts,truss_connectivity)

# Run analysis -------------------------------------------------------------

# Define "main" animation loop
Nsteps = 400
Nsub_steps = 10
dt = 2.0e-2 # [s] time increment
anim = Animation(Nsteps,Nsub_steps,dt,element_color="pressure",edge_color="Black")
recording = "" #"output_animation"
asyncio.run(anim.start(recording))

# --------------------------------------------------------------------------

