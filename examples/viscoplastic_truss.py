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
import pygame_widgets
from pygame_widgets.progressbar import ProgressBar

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
font_size  = 32 # (pixels)
display_font = pygame.font.SysFont(font_type, font_size)
anti_aliasing = True

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

    # Draw the original shape for comparison
    if (step_id == 0):
        # Draw truss elements
        for e in range(0,Ntruss):
            points = 100*coordinates[truss_connectivity[e,:],:]
            points[:,1] = 600 - points[:,1]
            pygame.draw.line(window, "LightBlue", points[0,:], points[1,:], 4)

    # Draw the mesh in its currently deformed configuration:
    max_pressure = 0.3
    xy_deformed, pressure, system_state, eqps, is_dead = REMAT.deform_geometry(coordinates)
    for e in range(0,Nelems):
        points = 100*xy_deformed[connectivity[e,:],:]
        points[:,1] = 600 - points[:,1]
        p1 = min(max(-pressure[e],0)/max_pressure,1.0)
        p2 = min(max(+pressure[e],0)/max_pressure,1.0)
        #color = (255*(1.0-p2),255*(1.0-p1-p2),255*(1.0-p1))
        color = "Black"
        pygame.draw.polygon(window, color, points)
        #pygame.draw.polygon(window, "Black", points, 1)

    # Draw truss elements
    for e in range(0,Ntruss):
        if (not is_dead[e]):
            points = 100*xy_deformed[truss_connectivity[e,:],:]
            points[:,1] = 600 - points[:,1]
            p1 = min(max(eqps[e]/eps_fail,0.0),1.0)
            color = (220,220*(1.0-p1),220*(1.0-p1))
            pygame.draw.line(window, "Dark Grey", points[0,:], points[1,:], 8)
            pygame.draw.line(window, color, points[0,:], points[1,:], 4)
            pygame.draw.circle(window, "Dark Grey", points[0,:], 5, 0)
            pygame.draw.circle(window, "Dark Grey", points[1,:], 5, 0)

    # Draw Play, Pause, and Rewind buttons
    bx = 20
    by = 20
    bd = 30
    color = "Dark Blue"
    if (keys[pygame.K_LEFT] and (step_id > 0)):
        pygame.draw.polygon(window, color, ((bx+0*bd,by+1*bd),(bx+2*bd,by+0*bd),(bx+2*bd,by+2*bd)), 0)
    elif (keys[pygame.K_RIGHT] and (step_id < Nsteps)):
        pygame.draw.polygon(window, color, ((bx+0*bd,by+0*bd),(bx+2*bd,by+1*bd),(bx+0*bd,by+2*bd)), 0)
    else:
        pygame.draw.polygon(window, color, ((bx+0*bd,by+0*bd),(bx+0.7*bd,by+0*bd),(bx+0.7*bd,by+2*bd),(bx+0*bd,by+2*bd)), 0)
        pygame.draw.polygon(window, color, ((bx+1.3*bd,by+0*bd),(bx+2*bd,by+0*bd),(bx+2*bd,by+2*bd),(bx+1.3*bd,by+2*bd)), 0)
        
    # Display system state
    #display_pos_x = 300
    #window.blit(display_font.render("    Strain energy: {:.4f}".format(system_state[0]),anti_aliasing,"Black"),(display_pos_x,1*font_size))
    #window.blit(display_font.render("   Kinetic energy: {:.4f}".format(system_state[1]),anti_aliasing,"Black"),(display_pos_x,2*font_size))
    #window.blit(display_font.render(" Potential energy: {:.4f}".format(system_state[2]),anti_aliasing,"Black"),(display_pos_x,3*font_size))
    #window.blit(display_font.render("     Total energy: {:.4f}".format(system_state[0]+system_state[1]+system_state[2]),anti_aliasing,"Blue"),(display_pos_x,4*font_size))
    
    # Draw different widgets within the game window:
    pygame_widgets.update(events)
            
    # Update the window to display the current state of the simulation
    pygame.display.update()

    # Wait enough time for the simulation to update at a fixed framerate
    clock.tick(framerate)

# --------------------------------------------------------------------------

# Call C/C++ library API functions from Python:

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

# Set the integrator type: "float" (default), or "fixed"
REMAT.API.set_integrator_type(b"float")

# Pre-process mesh/geometry ------------------------------------------------

# Manually define nodal coordinates and element geometry
#Nnodes = 3
#Nelems = 0
#Ntruss = 3
#coordinates  = np.array([[5.0, 1.0],[6.0,2.25],[5.0,3.5]])
#velocities   = np.zeros((3,2))
#fixity       = np.zeros((3,2),dtype=np.bool_)
#connectivity = np.empty((0,0),dtype=np.int32)
#truss_connectivity = np.array([[0,1],[1,2],[2,0]],dtype=np.int32)

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

# Generate the problem data from the Model object
coordinates, velocities, fixity, connectivity, contacts, truss_connectivity = model.generate_problem()
Nnodes = coordinates.shape[0]
Nelems = connectivity.shape[0]
Ntruss = truss_connectivity.shape[0]
        
# Define the problem geometry
REMAT.create_geometry(coordinates,velocities,fixity,connectivity,contacts,truss_connectivity)

# Run analysis -------------------------------------------------------------

# initialization
time = 0.0 # [s] starting time
REMAT.API.initialize()

# set analysis time-stepping parameters
dt = 1.0e-3 # [s] time increment
step_id = 0
Nsteps = 150
Nsub_steps = 10

# Create a progress bar to show advancement of time forward or backward
progressBarColour1=(0,0,200)
progressBarColour2=(158,128,220)
progressBar = ProgressBar(window, 100, 20, 800, 10,
                          lambda: step_id / Nsteps, curved=True,
                          completedColour=progressBarColour1,incompletedColour=progressBarColour2)

# Ensure we have somewhere for the frames
try:
    os.makedirs("Snaps_truss_float")
except OSError:
    pass
file_num = 0

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

    # Update file number
    file_num = file_num + 1

    # Save every frame
    screenshot_filename = "Snaps_truss_float/%04d.png" % file_num
    pygame.image.save(window, screenshot_filename)

# --------------------------------------------------------------------------

