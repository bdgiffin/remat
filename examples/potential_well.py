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

# Initialize a new font for drawing text
pygame.font.init()
font_type  = 'PT Mono'
font_size  = 16 # (pixels)
display_font = pygame.font.SysFont(font_type, font_size)
anti_aliasing = False

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
        pygame.draw.polygon(window, "LightGrey", points)
        pygame.draw.polygon(window, "Black", points, 1)

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
        
    # Display system state
    display_pos_x = 700
    window.blit(display_font.render("    Strain energy: {:.3e}".format(system_state[0]),anti_aliasing,"Black"),(display_pos_x,1*font_size))
    window.blit(display_font.render("   Kinetic energy: {:.3e}".format(system_state[1]),anti_aliasing,"Black"),(display_pos_x,2*font_size))
    window.blit(display_font.render(" Potential energy: {:.3e}".format(system_state[2]),anti_aliasing,"Black"),(display_pos_x,3*font_size))
    window.blit(display_font.render("     Total energy: {:.3e}".format(system_state[0]+system_state[1]+system_state[2]),anti_aliasing,"Black"),(display_pos_x,4*font_size))
            
    # Update the window to display the current state of the simulation
    pygame.display.update()

    # Wait enough time for the simulation to update at a fixed framerate
    clock.tick(framerate)

# --------------------------------------------------------------------------

# Call C/C++ library API functions from Python:

# Define global parameters
REMAT.API.define_parameter(b"body_force_y",       -1.0e-1)
REMAT.API.define_parameter(b"initial_velocity_x", +0.0e-1)
REMAT.API.define_parameter(b"initial_velocity_y", -0.0e-1)
REMAT.API.define_parameter(b"mass_damping_factor", 2.0e-2)
REMAT.API.define_parameter(b"contact_stiffness",   5.0e+0)
REMAT.API.define_parameter(b"search_radius",       1.0e+0)

# Define material parameters
REMAT.API.define_parameter(b"density",        1.0)
REMAT.API.define_parameter(b"youngs_modulus", 50.0)
REMAT.API.define_parameter(b"poissons_ratio", 0.28)

# Set the integrator type: "float" (default), or "fixed"
REMAT.API.set_integrator_type(b"float")

# Pre-process mesh/geometry ------------------------------------------------

# Create geometry factory
geom_factory = GeometryFactory()

# create block
block = geom_factory.cartesian_grid([1.0,4.5],[2.0,5.5],[3,3],0.1)
boundary_nodes = block.select_boundary_nodes()

with pygmsh.geo.Geometry() as geom:
    lcar = 1.0
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
Nsteps = 200
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

# --------------------------------------------------------------------------

