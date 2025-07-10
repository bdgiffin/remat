# Append the location of the locally installed REMAT package to sys.path
import sys
sys.path.append("../install/package/")

# Load the REMAT package
import REMAT

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

# --------------------------------------------------------------------------

# Initialize pygame window for interactive visualization
pygame.init()
window = pygame.display.set_mode((1000,600))

# Set the window title and icon
pygame.display.set_caption('Reversible Dynamics')


clock = pygame.time.Clock()
framerate = 60

# Initialize the active simulation "event loop" flag
simulation_running = True

# (consider using this to show reversibility in a step-by-step fashion)
# (include a "progress bar" or a clock to show advancement of time forward or backward)

# --------------------------------------------------------------------------

# Call C/C++ library API functions from Python:

# Define global parameters
REMAT.API.define_parameter(b"body_force_y",       -1.0e-2)
REMAT.API.define_parameter(b"initial_velocity_y", -1.0e-1)
REMAT.API.define_parameter(b"mass_damping_factor", 1.0e-1)
REMAT.API.define_parameter(b"contact_stiffness",   2.0e+0)

# Define material parameters
REMAT.API.define_parameter(b"density",        1.0)
REMAT.API.define_parameter(b"youngs_modulus", 1.0)
REMAT.API.define_parameter(b"poissons_ratio", 0.28)

# Pre-process mesh/geometry ------------------------------------------------

Nnodes = 9
Nelems = 4
coordinates = np.array([ 1.0, 1.1,
                         2.0, 1.0,
                         3.0, 0.9,
		         1.0, 2.0,
                         2.0, 2.0,
                         3.0, 2.0,
		 	 1.0, 3.0,
                         2.0, 3.0,
                         3.0, 3.0 ])
connectivity = np.array([[ 0, 1, 4, 3 ],
                         [ 1, 2, 5, 4 ],
                         [ 3, 4, 7, 6 ],
                         [ 4, 5, 8, 7 ]], dtype=np.int32);

# Define the problem geometry
filename = "output.exo"
REMAT.create_geometry(coordinates,connectivity,Nnodes,Nelems,filename)

# Run analysis -------------------------------------------------------------

# initialization
time = 0.0 # [s] starting time
REMAT.API.initialize()
REMAT.output_state()

# set analysis time-stepping parameters
dt = 2.0e-1 # [s] time increment
step_id = 0
Nsteps = 44
Nsub_steps = 10

# perform the forward-in-time analysis
while (simulation_running and (step_id < Nsteps)):
    # Update the simulation state
    step_id = step_id + 1
    time = REMAT.API.update_state(dt,Nsub_steps)
    REMAT.output_state()
    
    # Loop through any/all simulation "events":
    for event in pygame.event.get():
        # Check if the user has closed the window to "QUIT" the simulation:
        if event.type == pygame.QUIT:
            # If the user quit, stop the simulation from running:
            # (this will terminate the event loop and end the program)
            simulation_running = False
        
    # "Erase" the content of the display window
    window.fill("White") # reset window to display a full white screen

    # Draw the mesh in its currently deformed configuration:
    xy_deformed = REMAT.deform_geometry(coordinates)
    for e in range(0,Nelems):
        points = 100*xy_deformed[connectivity[e,:],:]
        points[:,1] = 600 - points[:,1]
        color  = "Blue"
        #points = ((325,75),(376,125),(275,200)) # ((x1,y1), (x2,y2), (x3,y3), ...)
        pygame.draw.polygon(window, color, points, 1)
            
    # Update the window to display the current state of the simulation
    pygame.display.update()

    # Wait enough time for the simulation to update at a fixed framerate
    clock.tick(framerate)

# preformed the backward-in-time analysis
while (simulation_running and (step_id > 0)):
    # Update the simulation state
    step_id = step_id - 1
    time = REMAT.API.update_state(-dt,Nsub_steps)
    REMAT.output_state()
    
    # Loop through any/all simulation "events":
    for event in pygame.event.get():
        # Check if the user has closed the window to "QUIT" the simulation:
        if event.type == pygame.QUIT:
            # If the user quit, stop the simulation from running:
            # (this will terminate the event loop and end the program)
            simulation_running = False
        
    # "Erase" the content of the display window
    window.fill("White") # reset window to display a full white screen

    # Draw the mesh in its currently deformed configuration:
    xy_deformed = REMAT.deform_geometry(coordinates)
    for e in range(0,Nelems):
        points = 100*xy_deformed[connectivity[e,:],:]
        points[:,1] = 600 - points[:,1]
        color  = "Blue"
        pygame.draw.polygon(window, color, points, 1)
            
    # Update the window to display the current state of the simulation
    pygame.display.update()

    # Wait enough time for the simulation to update at a fixed framerate
    clock.tick(framerate)

# finalize the REMAT module (close the Exodus files)
REMAT.finalize()

# --------------------------------------------------------------------------
