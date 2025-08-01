# Append the location of the locally installed REMAT package to sys.path
import sys, platform

# Import the "asyncio" Python package
import asyncio

# Append location of REMAT install package directory
if (not sys.platform == "emscripten"):
    sys.path.append("../install/package/")

# Load the REMAT and MeshUtils packages
#import package
import REMAT
from GeometryFactory import *
from Model import *

# Python package for reading/writing data in the Exodus mesh database format
#import pyexodus

# Other needed Python packages
from math import *
import numpy as np

# (Optional) use PyGame to visualize the problem as it is running
import pygame

# Import GMSH for meshing
import pygmsh

# --------------------------------------------------------------------------

# (include a "progress bar" or a clock to show advancement of time forward or backward)

# Initialize pygame window for interactive visualization
pygame.init()
window = pygame.display.set_mode((1000,600))

# Set the window title and icon
pygame.display.set_caption('Reversible Dynamics')

# Define a clock for consistent timestepping at a fixed framerate
clock = pygame.time.Clock()
framerate = 60

# Initialize the active simulation "event loop" flag
simulation_running = True

# Local function to animate the current deformed state of the mesh
def animate_state():
    # Loop through any/all simulation "events":
    for event in pygame.event.get():
        # Check if the user has closed the window to "QUIT" the simulation:
        if event.type == pygame.QUIT:
            # If the user quit, stop the simulation from running:
            # (this will terminate the event loop and end the program)
            simulation_running = False
            quit()
        
    # "Erase" the content of the display window
    window.fill("WhiteSmoke") # reset window to display a full white screen

    # Draw the mesh in its currently deformed configuration:
    max_pressure = 1.0
    xy_deformed, pressure, system_state, eqps, is_dead = REMAT.deform_geometry(coordinates)
    for e in range(0,Nelems):
        points = 100*xy_deformed[connectivity[e,:],:]
        points[:,1] = 600 - points[:,1]
        p1 = max(-pressure[e],0)/max_pressure
        p2 = max(+pressure[e],0)/max_pressure
        color = (255*(1.0-p2),255*(1.0-p1-p2),255*(1.0-p1))
        pygame.draw.polygon(window, color, points)
        pygame.draw.polygon(window, "Black", points, 1)

    # Draw contact interactions:
    for contact in contacts:
        points = 100*xy_deformed[contact[0],:]
        points[:,1] = 600 - points[:,1]
        for point in points:
            pygame.draw.circle(window, "Red", point, 5, 2)
        for segment in contact[1]:
            points = 100*xy_deformed[segment,:]
            points[:,1] = 600 - points[:,1]
            pygame.draw.line(window, "Green", points[0,:], points[1,:], 2)
            
    # Update the window to display the current state of the simulation
    pygame.display.update()

    # Wait enough time for the simulation to update at a fixed framerate
    clock.tick(framerate)

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

# Generate the problem data from the Model object
coordinates, velocities, fixity, connectivity, contacts, truss_connectivity = model.generate_problem()
Nnodes = coordinates.shape[0]
Nelems = connectivity.shape[0]
        
# Define the problem geometry
REMAT.create_geometry(coordinates,velocities,fixity,connectivity,contacts,truss_connectivity)

# Run analysis -------------------------------------------------------------

# Encapsulate the simulation loop inside of an anync function called "main()"
async def main():

    # initialization
    time = 0.0 # [s] starting time
    REMAT.API.initialize()

    # set analysis time-stepping parameters
    dt = 0.25e-2 # [s] time increment
    step_id = 0
    Nsteps = 300
    Nsub_steps = 100

    # perform the reversible analysis
    while simulation_running:

        # Get the "pressed" status of all keys on the keyboard:
        keys = pygame.key.get_pressed()

        # Update the simulation state backward-in-time
        if (keys[pygame.K_LEFT] and (step_id > 0)):
            step_id = step_id - 1
            time = REMAT.API.update_state(-dt,Nsub_steps)

        # Update the simulation state forward-in-time
        if (keys[pygame.K_RIGHT] and (step_id < Nsteps)):
            step_id = step_id + 1
            time = REMAT.API.update_state(+dt,Nsub_steps)

        # Animate the current state
        animate_state()

        # Call "asyncio.sleep(0)" within the simulation loop
        await asyncio.sleep(0)

# Use asyncio to run your "main()" simulation loop
asyncio.run(main())

# --------------------------------------------------------------------------

