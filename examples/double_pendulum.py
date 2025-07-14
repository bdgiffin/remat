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
    xy_deformed, pressure = REMAT.deform_geometry(coordinates)
    for e in range(0,Nelems):
        points = 100*xy_deformed[connectivity[e,:],:]
        points[:,1] = 600 - points[:,1]
        p1 = min(max(-pressure[e],0)/max_pressure,1.0)
        p2 = min(max(+pressure[e],0)/max_pressure,1.0)
        color = (255*(1.0-p2),255*(1.0-p1-p2),255*(1.0-p1))
        pygame.draw.polygon(window, color, points)
        pygame.draw.polygon(window, "Black", points, 1)

    # Draw truss elements
    for e in range(0,Ntruss):
        points = 100*xy_deformed[truss_connectivity[e,:],:]
        points[:,1] = 600 - points[:,1]
        pygame.draw.line(window, "Grey", points[0,:], points[1,:], 4)

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

    # Draw point masses
    for e in range(0,Npoints):
        point = 100*xy_deformed[point_ids[e],:]
        point[1] = 600 - point[1]
        pygame.draw.circle(window, "Black", point, 4*point_mass[e], 0)
        # create tracer
        if (keys[pygame.K_LEFT] and (step_id > 0)):
            tracer.pop()
        elif (keys[pygame.K_RIGHT] and (step_id < Nsteps)):
            tracer.append(point)

    # Draw tracer
    for e in range(0,len(tracer)-1):
        print(tracer[e])
        print(tracer[e+1])
        pygame.draw.line(window, "Blue", tracer[e], tracer[e+1], 2)

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
REMAT.API.define_parameter(b"body_force_y",       -1.0e-1)
REMAT.API.define_parameter(b"initial_velocity_x", +0.0e-1)
REMAT.API.define_parameter(b"initial_velocity_y", -0.0e-1)
REMAT.API.define_parameter(b"mass_damping_factor", 1.0e-2)
REMAT.API.define_parameter(b"contact_stiffness",   0.0e+0)

# Define material parameters
REMAT.API.define_parameter(b"density",         1.0)
REMAT.API.define_parameter(b"youngs_modulus", 1000.0)
REMAT.API.define_parameter(b"area",            1.0)
REMAT.API.define_parameter(b"poissons_ratio", 0.28)

# Pre-process mesh/geometry ------------------------------------------------

# Manually define nodal coordinates and element geometry
Nnodes = 3
Nelems = 0
Ntruss = 2
coordinates  = np.array([[5.0, 3.0],[5.00001,4.25],[5.0,5.5]])
velocities   = np.zeros((3,2))
fixity       = np.array([[True,True],[False,False],[False,False]],dtype=np.bool_)
connectivity = np.empty((0,0),dtype=np.int32)
truss_connectivity = np.array([[0,1],[1,2]],dtype=np.int32)
contacts = []
        
# Define the problem geometry
REMAT.create_geometry(coordinates,velocities,fixity,connectivity,contacts)

# Define truss elements
REMAT.API.define_truss_elements(truss_connectivity,Ntruss)

# Define point masses
Npoints = 1
point_ids = np.array([2],dtype=np.int32)
point_mass = np.array([2.0])
REMAT.API.define_point_mass(point_ids,point_mass,Npoints)
tracer = []

# Run analysis -------------------------------------------------------------

# initialization
time = 0.0 # [s] starting time
REMAT.API.initialize()

# set analysis time-stepping parameters
dt = 1.0e-2 # [s] time increment
step_id = 0
Nsteps = 1000
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

