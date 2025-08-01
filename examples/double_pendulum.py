# Append the location of the locally installed REMAT package to sys.path
import sys
#sys.path.append("../../install/package/")

# Import the "asyncio" Python package
import asyncio

# Load the REMAT and MeshUtils packages
import REMAT

# Other needed Python packages
from math import *
import numpy as np
import random

# (Optional) use PyGame to visualize the problem as it is running
import pygame
import pygame_widgets
from pygame_widgets.progressbar import ProgressBar

# --------------------------------------------------------------------------

# Encapsulate the simulation loop inside of an anync function called "main()"
async def main():

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
        events = pygame.event.get()
        for event in events:
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
        xy_deformed, pressure, system_state, eqps, is_dead = REMAT.deform_geometry(coordinates)
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
            pygame.draw.line(window, "Grey", points[0,:], points[1,:], 5)
            pygame.draw.circle(window, "Dark Grey", points[0,:], 6, 0)
            pygame.draw.circle(window, "Dark Grey", points[1,:], 6, 0)

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
            point    = 100*xy_deformed[point_ids[e],:]
            point[1] = 600 - point[1]
            point_x = int(point[0])
            point_y = int(point[1])
            radius  = int(6*point_mass[e])
            color   = (0,0,0)
            pygame.draw.circle(window, color, point, radius, 0)
            # create tracers
            if (keys[pygame.K_LEFT] and (step_id > 0)):
                #tracer.pop()
                reverse_tracer.append(point)
            elif (keys[pygame.K_RIGHT] and (step_id < Nsteps)):
                if (step_id == 1):
                    tracer.clear()
                    reverse_tracer.clear()
                tracer.append(point)

        # Draw tracers
        if (len(tracer) > 2):
            pygame.draw.aalines(window, "Blue", False, tracer)
        if (len(reverse_tracer) > 2):
            pygame.draw.aalines(window, "Magenta", False, reverse_tracer)

        # Draw Play, Pause, and Rewind buttons
        bx = 20
        by = 520
        bd = 30
        if (keys[pygame.K_LEFT] and (step_id > 0)):
            pygame.draw.polygon(window, "Grey", ((bx+0*bd,by+1*bd),(bx+2*bd,by+0*bd),(bx+2*bd,by+2*bd)), 0)
        elif (keys[pygame.K_RIGHT] and (step_id < Nsteps)):
            pygame.draw.polygon(window, "Grey", ((bx+0*bd,by+0*bd),(bx+2*bd,by+1*bd),(bx+0*bd,by+2*bd)), 0)
        else:
            pygame.draw.polygon(window, "Grey", ((bx+0*bd,by+0*bd),(bx+0.7*bd,by+0*bd),(bx+0.7*bd,by+2*bd),(bx+0*bd,by+2*bd)), 0)
            pygame.draw.polygon(window, "Grey", ((bx+1.3*bd,by+0*bd),(bx+2*bd,by+0*bd),(bx+2*bd,by+2*bd),(bx+1.3*bd,by+2*bd)), 0)

        # Draw different widgets within the game window:
        pygame_widgets.update(events)

        # Update the window to display the current state of the simulation
        pygame.display.update()

        # Wait enough time for the simulation to update at a fixed framerate
        clock.tick(framerate)

    # --------------------------------------------------------------------------

    # Call C/C++ library API functions from Python:

    # Define global parameters
    REMAT.API.define_parameter(b"body_force_y",       -2.0e-2)
    REMAT.API.define_parameter(b"initial_velocity_x", +0.5e-1)
    REMAT.API.define_parameter(b"initial_velocity_y", -0.0e-1)
    REMAT.API.define_parameter(b"mass_damping_factor", 0.0e-2)
    REMAT.API.define_parameter(b"contact_stiffness",   0.0e+0)

    # Define material parameters
    REMAT.API.define_parameter(b"density",         1.0)
    REMAT.API.define_parameter(b"youngs_modulus", 1000.0)
    REMAT.API.define_parameter(b"yield_stress",   1.0e+9)
    REMAT.API.define_parameter(b"viscosity",       1.0)
    REMAT.API.define_parameter(b"area",            1.0)
    REMAT.API.define_parameter(b"poissons_ratio", 0.28)

    # Set the integrator type: "float" (default), or "fixed"
    REMAT.API.set_integrator_type(b"fixed")

    # Pre-process mesh/geometry ------------------------------------------------

    # Manually define nodal coordinates and element geometry
    length = 3.5
    Nnodes = 4
    Nelems = 0
    Ntruss = Nnodes-1
    coordinates  = np.zeros((Nnodes,2))
    velocities   = np.zeros((Nnodes,2))
    fixity       = np.zeros((Nnodes,2),dtype=np.bool_)
    truss_connectivity = np.zeros((Ntruss,2),dtype=np.int32)
    random.seed(8)
    for i in range(0,Nnodes):
        coordinates[i,:] = [5.0+0.001*random.random(), 3.0+(length/Nnodes)*i]
        if (i == 0):
            fixity[i,:] = [True,True]
        else:
            truss_connectivity[i-1,:] = [i-1,i]

    # Define the problem geometry
    connectivity = np.empty((0,0),dtype=np.int32)
    contacts = []
    REMAT.create_geometry(coordinates,velocities,fixity,connectivity,contacts,truss_connectivity)

    # Define point masses
    Npoints = 1
    point_ids = np.array([Nnodes-1],dtype=np.int32)
    point_mass = np.array([3.0])
    REMAT.API.define_point_mass(point_ids,point_mass,Npoints)
    tracer = []
    reverse_tracer = []

    # Run analysis -------------------------------------------------------------

    # initialization
    time = 0.0 # [s] starting time
    REMAT.API.initialize()

    # set analysis time-stepping parameters
    dt = 1.0e-2 # [s] time increment
    step_id = 0
    Nsteps = 400
    Nsub_steps = 50

    # Create a progress bar to show advancement of time forward or backward
    progressBarColour1=(0,0,200)
    progressBarColour2=(158,128,220)
    progressBar = ProgressBar(window, 100, 570, 800, 10,
                              lambda: step_id / Nsteps, curved=True,
                              completedColour=progressBarColour1,incompletedColour=progressBarColour2)

    # Ensure we have somewhere for the frames
    #try:
    #    os.makedirs("Snaps_fixed")
    #except OSError:
    #    pass
    #file_num = 0

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

        ## Update file number
        #file_num = file_num + 1
        #
        ## Save every frame
        #screenshot_filename = "Snaps_fixed/%04d.png" % file_num
        #pygame.image.save(window, screenshot_filename)

        # Call "asyncio.sleep(0)" within the simulation loop
        await asyncio.sleep(0)

# --------------------------------------------------------------------------

# Use asyncio to run your "main()" simulation loop
asyncio.run(main())


