# Append the location of the locally installed REMAT package to sys.path
import sys
sys.path.append("../install/package/")

# Load the REMAT and MeshUtils packages
import REMAT
from GeometryFactory import *
from Model import *

# Python package for reading/writing data in the Exodus mesh database format
# NOTE: PYEXODUS V0.1.5 NEEDS TO BE MODIFIED TO WORK CORRECTLY WITH PYTHON 3.12
# https://pypi.org/project/pyexodus/
import pyexodus

# Other needed Python packages
from math import *
import numpy as np
import matplotlib.pyplot as plt

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
    max_pressure = 0.3
    xy_deformed, pressure, system_state = REMAT.deform_geometry(coordinates)
    for e in range(0,Nelems):
        points = 100*xy_deformed[connectivity[e,:],:]
        points[:,1] = 600 - points[:,1]
        p1 = min(max(-pressure[e],0)/max_pressure,1.0)
        p2 = min(max(+pressure[e],0)/max_pressure,1.0)
        color = (255*(1.0-p2),255*(1.0-p1-p2),255*(1.0-p1))
        pygame.draw.polygon(window, color, points)
        pygame.draw.polygon(window, "Black", points, 1)

    # Draw contact interactions (primarily for debugging purposes):
    for contact in contacts:
        points = 100*xy_deformed[contact[0],:]
        points[:,1] = 600 - points[:,1]
        for point in points:
            pygame.draw.circle(window, "Red", point, 5, 2)
        for segment in contact[1]:
            points = 100*xy_deformed[segment,:]
            points[:,1] = 600 - points[:,1]
            pygame.draw.line(window, "Green", points[0,:], points[1,:], 2)

    # Draw Play, Pause, and Rewind buttons
    bx = 20
    by = 20
    bd = 30
    if (keys[pygame.K_LEFT] and (step_id > 0)):
        pygame.draw.polygon(window, "Grey", ((bx+0*bd,by+1*bd),(bx+2*bd,by+0*bd),(bx+2*bd,by+2*bd)), 0)
    elif (keys[pygame.K_RIGHT] and (step_id < Nsteps)):
        pygame.draw.polygon(window, "Grey", ((bx+0*bd,by+0*bd),(bx+2*bd,by+1*bd),(bx+0*bd,by+2*bd)), 0)
    else:
        pygame.draw.polygon(window, "Grey", ((bx+0*bd,by+0*bd),(bx+0.7*bd,by+0*bd),(bx+0.7*bd,by+2*bd),(bx+0*bd,by+2*bd)), 0)
        pygame.draw.polygon(window, "Grey", ((bx+1.3*bd,by+0*bd),(bx+2*bd,by+0*bd),(bx+2*bd,by+2*bd),(bx+1.3*bd,by+2*bd)), 0)
            
    # Update the window to display the current state of the simulation
    pygame.display.update()

    # Wait enough time for the simulation to update at a fixed framerate
    clock.tick(framerate)

# --------------------------------------------------------------------------

# Call C/C++ library API functions from Python:

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

# Pre-process mesh/geometry ------------------------------------------------

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
grid = geom_factory.cartesian_grid([0.0,0.0],[10.0,3.0],[100,30])
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

# Generate the problem data from the Model object
coordinates, velocities, fixity, connectivity, contacts = model.generate_problem()
Nnodes = coordinates.shape[0]
Nelems = connectivity.shape[0]
        
# Define the problem geometry
REMAT.create_geometry(coordinates,velocities,fixity,connectivity,contacts)

# Run analysis -------------------------------------------------------------

# initialization
time = 0.0 # [s] starting time
REMAT.API.initialize()

# set analysis time-stepping parameters
dt = 1.0e-2 # [s] time increment
step_id = 0
Nsteps = 1000
Nsub_steps = 10

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

# --------------------------------------------------------------------------

