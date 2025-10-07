# Module for enabling interactive animation of REMAT simulations

# Import the "asyncio" Python package (for compatibility with pygbag)
import asyncio

# Use REMAT Python API to retrieve data to be displayed
import REMAT

# Use PyGame to visualize the problem as it is running
import pygame
import pygame_widgets
from pygame_widgets.progressbar import ProgressBar

# Other needed Python packages
import os
import numpy as np

class Animation:
    # The Animation object facilitates interactive REMAT simulation and visualization

    def __init__(self,Nsteps,Nsub_steps,dt,
                 display_energy=False,
                 element_color="Grey",   element_field_max=1.0, edge_color="",
                 truss_color="Dark Grey",truss_field_max=1.0,
                 point_color="Black",    point_field_max=1.0):

        # set analysis time-stepping parameters
        self.dt = dt # 1.0e-2 # [s] time increment
        self.step_id = 0
        self.Nsteps = Nsteps # 220
        self.Nsub_steps = Nsub_steps # 80

        # set display parameters
        self.display_energy = display_energy
        self.element_color = element_color
        self.element_field_max = element_field_max
        self.edge_color = edge_color
        self.truss_color = truss_color
        self.truss_field_max = truss_field_max
        self.point_color = point_color
        self.point_field_max = point_field_max
        
        # Initialize pygame window for interactive visualization
        pygame.init()
        self.window = pygame.display.set_mode((1000,600))
        
        # Set the window title and icon
        pygame.display.set_caption('REMAT')
        
        # Initialize a new font for drawing text
        pygame.font.init()
        font_type  = 'PT Mono'
        self.font_size  = 32 # (pixels)
        self.display_font = pygame.font.SysFont(font_type, self.font_size)
        self.anti_aliasing = True

        # Create a progress bar to show advancement of time forward or backward
        progressBarColour1=(0,0,200)
        progressBarColour2=(158,128,220)
        self.progressBar = ProgressBar(self.window, 100, 570, 800, 10,
                                       lambda: self.step_id / self.Nsteps, curved=True,
                                       completedColour=progressBarColour1,incompletedColour=progressBarColour2)
        
        # Define a clock for consistent timestepping at a fixed framerate
        self.clock = pygame.time.Clock()
        self.framerate = 60
        
        # Initialize the active simulation "event loop" flag
        self.simulation_running = True

        # Initialize discrete point tracers
        self.tracer = []
        self.reverse_tracer = []

    # --------------------------------------------------------------------------

    # Function to display the current deformed state of the simulation
    async def display_state(self):
        # Get the "pressed" status of all keys on the keyboard:
        keys = pygame.key.get_pressed()
        
        # Loop through any/all simulation "events":
        events = pygame.event.get()
        for event in events:
            # Check if the user has closed the window to "QUIT" the simulation:
            if event.type == pygame.QUIT:
                # If the user quit, stop the simulation from running:
                # (this will terminate the event loop and end the program)
                self.simulation_running = False
        
        # "Erase" the content of the display window
        self.window.fill("WhiteSmoke") # reset window to display a full white screen

        # get mesh totals
        num_dim   = REMAT.API.get_num_dim()
        num_nodes = REMAT.API.get_num_entities(b"node")
        num_elems = REMAT.API.get_num_entities(b"element")
        num_truss = REMAT.API.get_num_entities(b"truss")
        num_point = REMAT.API.get_num_entities(b"point")

        # Get the currently deformed coordinates of all nodes
        coords = np.zeros((num_nodes,num_dim),dtype=np.double)
        REMAT.API.get_node_coords(coords,True)

        # Draw solid elements
        connectivity = np.zeros((num_elems,4),dtype=np.int32)
        REMAT.API.get_connectivity(b"element",connectivity)
        element_field = REMAT.get_field(b"element",self.element_color)
        for e in range(0,num_elems):
            points = 100*coords[connectivity[e,:],:]
            points[:,1] = 600 - points[:,1]
            if (element_field is not None):
                p1 = min(max(-element_field[e],0)/self.element_field_max,1.0)
                p2 = min(max(+element_field[e],0)/self.element_field_max,1.0)
                color = (255*(1.0-p2),255*(1.0-p1-p2),255*(1.0-p1))
            else:
                color = self.element_color
            pygame.draw.polygon(self.window, color, points)
            if (self.edge_color):
                pygame.draw.polygon(self.window, self.edge_color, points, 1)

        # Draw truss elements
        connectivity = np.zeros((num_truss,2),dtype=np.int32)
        REMAT.API.get_connectivity(b"truss",connectivity)
        truss_field = REMAT.get_field(b"truss",self.truss_color)
        is_dead = REMAT.get_field(b"truss","steps_since_element_death")
        for e in range(0,num_truss):
            if (not is_dead[e]):
                points = 100*coords[connectivity[e,:],:]
                points[:,1] = 600 - points[:,1]
                if (truss_field is not None):
                    p1 = min(max(truss_field[e]/self.truss_field_max,0.0),1.0)
                    color = (220,220*(1.0-p1),220*(1.0-p1))
                else:
                    color = self.truss_color
                pygame.draw.line(self.window, "Dark Grey", points[0,:], points[1,:], 8)
                pygame.draw.line(self.window, color, points[0,:], points[1,:], 4)
                pygame.draw.circle(self.window, "Dark Grey", points[0,:], 5, 0)
                pygame.draw.circle(self.window, "Dark Grey", points[1,:], 5, 0)

        # Draw discrete points
        connectivity = np.zeros((num_point,1),dtype=np.int32)
        REMAT.API.get_connectivity(b"point",connectivity)
        point_mass = REMAT.get_field(b"point","mass")
        for e in range(0,num_point):
            point    = 100*coords[connectivity[e,0],:]
            point[1] = 600 - point[1]
            radius  = int(6*(point_mass[e]/self.point_field_max))
            pygame.draw.circle(self.window, self.point_color, point, radius, 0)
            # create tracers
            if (keys[pygame.K_LEFT] and (self.step_id > 0)):
                self.reverse_tracer.append(point)
            elif (keys[pygame.K_RIGHT] and (self.step_id < self.Nsteps)):
                if (self.step_id == 1):
                    self.tracer.clear()
                    self.reverse_tracer.clear()
                self.tracer.append(point)

        # Draw discrete point tracers
        if (len(self.tracer) > 2):
            pygame.draw.aalines(self.window, "Blue", False, self.tracer)
        if (len(self.reverse_tracer) > 2):
            pygame.draw.aalines(self.window, "Magenta", False, self.reverse_tracer)

        # Draw Play, Pause, and Rewind buttons
        bx = 20
        by = 520
        bd = 30
        color = "Dark Blue"
        if (keys[pygame.K_LEFT] and (self.step_id > 0)):
            pygame.draw.polygon(self.window, color, ((bx+0*bd,by+1*bd),(bx+2*bd,by+0*bd),(bx+2*bd,by+2*bd)), 0)
        elif (keys[pygame.K_RIGHT] and (self.step_id < self.Nsteps)):
            pygame.draw.polygon(self.window, color, ((bx+0*bd,by+0*bd),(bx+2*bd,by+1*bd),(bx+0*bd,by+2*bd)), 0)
        else:
            pygame.draw.polygon(self.window, color, ((bx+0*bd,by+0*bd),(bx+0.7*bd,by+0*bd),(bx+0.7*bd,by+2*bd),(bx+0*bd,by+2*bd)), 0)
            pygame.draw.polygon(self.window, color, ((bx+1.3*bd,by+0*bd),(bx+2*bd,by+0*bd),(bx+2*bd,by+2*bd),(bx+1.3*bd,by+2*bd)), 0)
        
        # Display system state
        if (self.display_energy):
            display_pos_x = 300
            self.window.blit(self.display_font.render("    Strain energy: {:.4f}".format(REMAT.get_field(b"global","elastic_strain_energy")[0]),
                                                      self.anti_aliasing,"Black"),(display_pos_x,1*self.font_size))
            self.window.blit(self.display_font.render("   Kinetic energy: {:.4f}".format(REMAT.get_field(b"global","kinetic_energy")[0]),
                                                      self.anti_aliasing,"Black"),(display_pos_x,2*self.font_size))
            self.window.blit(self.display_font.render(" Potential energy: {:.4f}".format(REMAT.get_field(b"global","potential_energy")[0]),
                                                      self.anti_aliasing,"Black"),(display_pos_x,3*self.font_size))
            self.window.blit(self.display_font.render("     Total energy: {:.4f}".format(REMAT.get_field(b"global","total_energy")[0]),
                                                      self.anti_aliasing,"Black"),(display_pos_x,4*self.font_size))
        
        # Draw different widgets within the game window:
        pygame_widgets.update(events)
    
        # Update the window to display the current state of the simulation
        pygame.display.update()

        # Wait enough time for the simulation to update at a fixed framerate
        self.clock.tick(self.framerate)
    
    # --------------------------------------------------------------------------

    async def start(self,recording=""):
        # Optionally record the animation and output all frames to the "recording" directory
        if recording:
            # Create output directory for recording (if it doesn't already exist)
            try:
                os.makedirs(recording)
            except OSError:
                pass

        # Initialize the animation frame number and step_id
        self.step_id = 0
        frame_num = 0
        
        # Initialize REMAT
        time = 0.0 # [s] starting time
        REMAT.API.initialize()

        # perform the reversible analysis
        while self.simulation_running:

            # Get the "pressed" status of all keys on the keyboard:
            keys = pygame.key.get_pressed()
    
            # Update the simulation state backward-in-time
            if (keys[pygame.K_LEFT] and (self.step_id > 0)):
                self.step_id = self.step_id - 1
                time = REMAT.API.update_state(-self.dt,self.Nsub_steps)
    
            # Update the simulation state forward-in-time
            if (keys[pygame.K_RIGHT] and (self.step_id < self.Nsteps)):
                self.step_id = self.step_id + 1
                time = REMAT.API.update_state(+self.dt,self.Nsub_steps)
        
            # Display the current state of the simulation
            await self.display_state()

            # Update animation frame number
            frame_num = frame_num + 1

            if recording:
                # Save current animation frame
                screenshot_filename = recording + "/%04d.png" % frame_num
                pygame.image.save(self.window, screenshot_filename)
                
            # Call "asyncio.sleep(0)" within the simulation loop (for compatibility with pygbag)
            await asyncio.sleep(0)

        # Consolidate animated frames into a movie (requires ffmpeg)
        if recording:
            os.system("ffmpeg -framerate " + str(self.framerate) + " -i " + recording + "/%04d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p " + recording + ".mp4")
            os.system("rm -rf " + recording)

        # Quit the program
        quit()

    # --------------------------------------------------------------------------
