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
REMAT.API.define_parameter(b"body_force_y",       -0.0e-3)
REMAT.API.define_parameter(b"initial_velocity_x", +0.0e-1)
REMAT.API.define_parameter(b"initial_velocity_y", -0.0e-1)
REMAT.API.define_parameter(b"mass_damping_factor", 5.0e-2)
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
        obj = geom.add_circle([5.0,5.0],0.5,mesh_size=0.3,compound=True,num_sections=4,make_surface=True)
        print(obj)
        #geom.set_recombined_surfaces([obj.surface])
        impactor = geom_factory.from_meshio(geom.generate_mesh(dim=2))
        boundary_nodes = impactor.select_boundary_nodes()
        impactor_nodes = impactor.select_nodes(Select_all())

    with pygmsh.geo.Geometry() as geom:
        elev   = 2.0
        height = 3.5
        thick  = 0.25
        lcar = 0.125
        p1  = geom.add_point([0.0,  elev], lcar)
        c1  = geom.add_point([2.0,  elev], lcar)
        c13 = geom.add_point([5.0,  elev+height], lcar)
        c3  = geom.add_point([8.0,  elev], lcar)
        p3  = geom.add_point([10.0, elev], lcar)

        p2  = geom.add_point([0.0,  elev+thick], lcar)
        c2  = geom.add_point([2.0,  elev+thick], lcar)
        c24 = geom.add_point([5.0,  elev+thick+height], lcar)
        c4  = geom.add_point([8.0,  elev+thick], lcar)
        p4  = geom.add_point([10.0, elev+thick], lcar)

        s1 = geom.add_bspline([p1, c1, c13, c3, p3])
        s2 = geom.add_line(p3, p4)
        s3 = geom.add_bspline([p4, c4, c24, c2, p2])
        s4 = geom.add_line(p2, p1)

        ll = geom.add_curve_loop([s1, s2, s3, s4])
        pl = geom.add_surface(ll)

        geom.set_recombined_surfaces([pl])
        arch = geom_factory.from_meshio(geom.generate_mesh(dim=2))

    # create half-space
    arch_faces  = arch.select_boundary_faces()
    left_nodes  = arch.select_nodes(Select_X_eq(0.0))
    right_nodes = arch.select_nodes(Select_X_eq(10.0))

    # Create Model object
    model = Model()
    model.add_part(Part(impactor,Material(None,None)))
    model.add_part(Part(arch,Material(None,None)))
    #model.add_initial_condition(impactor_nodes,[+0.0e-1,-5.0e-1])
    model.add_initial_condition(impactor_nodes,[+0.0e-1,-4.2e-1])
    model.add_boundary_condition(left_nodes,[True,True]) # fully fix the left surface ndoes
    model.add_boundary_condition(right_nodes,[True,True]) # fully fix the right surface ndoes
    model.add_contact_interaction(boundary_nodes,arch_faces)
    
    # Generate the problem data from the Model object and save in a pickle file
    pickle.dump(model.generate_problem(),open("problem.pkl",'wb'))

# Load the problem data from pickle file
coordinates, velocities, fixity, connectivity, contacts, truss_connectivity = pickle.load(open("problem.pkl",'rb'))
        
# Define the problem geometry
REMAT.create_geometry(coordinates,velocities,fixity,connectivity,contacts,truss_connectivity)

# Run analysis -------------------------------------------------------------

# Define "main" animation loop
Nsteps = 180
Nsub_steps = 100
dt = 1.0e-2 # [s] time increment
anim = Animation(Nsteps,Nsub_steps,dt,element_color="pressure",edge_color="Black")
recording = "" # "output_animation"
asyncio.run(anim.start(recording))

# --------------------------------------------------------------------------

