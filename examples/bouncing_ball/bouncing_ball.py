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
REMAT.API.define_parameter(b"overflow_limit",      5000.0) # steps

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

    # create impactor
    with pygmsh.geo.Geometry() as geom:
        cx = 1.5
        cy = 5.4
        cr = 0.06
        obj = geom.add_circle([cx,cy],0.5,mesh_size=0.15,compound=True,num_sections=4,make_surface=True)
        print(obj)
        impactor = geom_factory.from_meshio(geom.generate_mesh(dim=2))
        tracer_nodes = impactor.select_nodes(Select_XY_window([cx-cr,cy-cr],[cx+cr,cy+cr]))

    # Create Model object
    model = Model()
    model.add_part(Part(impactor,Material(None,None)))

    # Generate the problem data from the Model object and save in a pickle file
    pickle.dump((model.generate_problem(), tracer_nodes),open("problem.pkl",'wb'))

# Load the problem data from pickle file
(coordinates, velocities, fixity, connectivity, contacts, truss_connectivity), tracer_nodes = pickle.load(open("problem.pkl",'rb'))
        
# Define the problem geometry
REMAT.create_geometry(coordinates,velocities,fixity,connectivity,contacts,truss_connectivity)

# Define point masses
point_ids  = tracer_nodes.vertex_ids
Npoints = point_ids.shape[0]
point_mass = np.ones(Npoints) * 0.0 # massless tracer point
REMAT.API.define_point_mass(point_ids.astype(np.int32),point_mass,Npoints)

# Run analysis -------------------------------------------------------------

# Define "main" animation loop
Nsteps = 300
Nsub_steps = 100
dt = 1.0e-2 # [s] time increment
anim = Animation(Nsteps,Nsub_steps,dt,element_color="Dark Red")
recording = "output_animation"
asyncio.run(anim.start(recording))

# --------------------------------------------------------------------------

