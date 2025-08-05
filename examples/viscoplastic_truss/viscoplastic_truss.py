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

# Load the REMAT and Animation packages
import REMAT
from Animation import *

# --------------------------------------------------------------------------

# Define global parameters
REMAT.API.define_parameter(b"body_force_y",       -0.0e-2)
REMAT.API.define_parameter(b"initial_velocity_x", +0.0e-1)
REMAT.API.define_parameter(b"initial_velocity_y",  0.0e+0)
#REMAT.API.define_parameter(b"initial_velocity_y", -4.0e+0)
REMAT.API.define_parameter(b"mass_damping_factor", 1.0e-2)
REMAT.API.define_parameter(b"contact_stiffness",   2.0e+3)

# Define material parameters
eps_fail = 0.3
REMAT.API.define_parameter(b"truss_density",   1.0)
REMAT.API.define_parameter(b"truss_youngs_modulus", 2000.0)
REMAT.API.define_parameter(b"density",        3.0)
REMAT.API.define_parameter(b"youngs_modulus", 2000.0)
REMAT.API.define_parameter(b"yield_stress",   30.0)
REMAT.API.define_parameter(b"viscosity",      1.0e+1)
REMAT.API.define_parameter(b"area",            1.0)
REMAT.API.define_parameter(b"poissons_ratio", 0.28)
REMAT.API.define_parameter(b"eps_fail",   eps_fail)
#REMAT.API.define_parameter(b"overflow_limit", 1000.0)
REMAT.API.define_parameter(b"mat_overflow_limit", 50.0)

# Set the integrator type: "float" (default), "fixed", or "mixed"
REMAT.API.set_integrator_type(b"fixed")

# Pre-process mesh/geometry ------------------------------------------------

# Generate and save problem geometry
if (not sys.platform == "emscripten"):
    # Create geometry factory
    geom_factory = GeometryFactory()

    # create impactor
    with pygmsh.geo.Geometry() as geom:
        #obj = geom.add_rectangle(2.5, 3.5, 3.2, 4.2, 1.0, 0.5)
        obj = geom.add_circle([6.5,5.0],1.0,mesh_size=0.5,compound=True,num_sections=4,make_surface=True)
        #geom.set_recombined_surfaces([obj.surface])
        impactor = geom_factory.from_meshio(geom.generate_mesh(dim=2))
        boundary_nodes = impactor.select_boundary_nodes()
        impactor_nodes = impactor.select_nodes(Select_all())

    # read truss geometry from dxf file
    dxf_filename = "platform.dxf" # "arch.dxf"
    layer_names  = ["truss"]
    truss = geom_factory.from_dxf(dxf_filename,layer_names)

    # select faces and nodes on boundary of truss
    top_faces    = truss.select_faces(Select_Y_eq(3.0))
    bottom_nodes = truss.select_nodes(Select_Y_eq(0.0))
    left_nodes   = truss.select_nodes(Select_X_eq(0.0))
    right_nodes  = truss.select_nodes(Select_X_eq(10.0))

    # Create Model object
    model = Model()
    model.add_part(Part(impactor,Material(None,None)))
    model.add_part(Part(truss,Material(None,None)))
    model.add_initial_condition(impactor_nodes,[+5.0e-1,-4.0e-0])
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
Nsteps = 150
Nsub_steps = 10
dt = 1.0e-3 # [s] time increment
anim = Animation(Nsteps,Nsub_steps,dt,
                 element_color="Black",truss_color="equivalent_plastic_strain",truss_field_max=eps_fail)
recording = "" # "output_animation"
asyncio.run(anim.start(recording))

# --------------------------------------------------------------------------
