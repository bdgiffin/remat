# Module for calling the C/C++ REMAT API functions from Python

# Python package for calling C/C++ functions from Python
from ctypes import CDLL, POINTER
from ctypes import c_size_t, c_double, c_int, c_char_p

# Package for reading/writing mesh files
import pyexodus

# Other needed Python packages
import sys
import os
import math
import time as timer
from datetime import timedelta
import numpy as np
from argparse import ArgumentParser

# ---------------------------------------------------------------------------- #

# Module initialization:

# Load the pre-compiled external C/C++ "shared object" libraries
library_name = os.path.join(os.path.dirname(__file__), 'libREMAT.so')
API = CDLL(library_name)

# Define types to convert Numpy arrays into C arrays:

# C-types corresponding to pointers to int/double values
c_int_p    = POINTER(c_int)
c_double_p = POINTER(c_double)

# C-type corresponding to 1D numpy array
ND_POINTER_1 = np.ctypeslib.ndpointer(dtype=np.float64, 
                                      ndim=1,
                                      flags="C")
NI_POINTER_2 = np.ctypeslib.ndpointer(dtype=np.int32, 
                                      ndim=2,
                                      flags="C")

# Define all C/C++ library API function signatures
API.define_parameter.argtypes = [c_char_p, c_double]
API.define_parameter.restype  = None
API.define_geometry.argtypes = [ND_POINTER_1, NI_POINTER_2, c_size_t, c_size_t]
API.define_geometry.restype  = None
API.initialize.argtypes = [c_double]
API.initialize.restype  = None
API.update_state.argtypes = [c_double, c_double]
API.update_state.restype  = None
API.get_field_data.argtypes = [c_int_p, c_double_p, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1]
API.get_field_data.restype  = None

# ---------------------------------------------------------------------------- #

# Generate mesh geometry and initialize the REMAT object prior to initialization
def create_geometry(x,connectivity,Nnodes,Nelems,filename):

    global num_nodes
    global num_elems
    num_nodes = Nnodes
    num_elems = Nelems

    # Call REMAT initialization API function
    API.define_geometry(x,connectivity,Nnodes,Nelems)

# ---------------------------------------------------------------------------- #

    # Create the Exodus file containing the problem info:

    # create a new Exodus file
    try:
        os.remove(filename)
    except OSError:
        pass
    global exo
    exo = pyexodus.exodus(file=filename, mode='w', array_type='numpy', title='Produced by REMAT module', numDims=2, numNodes=Nnodes, numElems=Nelems, numBlocks=1, numNodeSets=0, numSideSets=0, io_size=0, compression=None)
    
    # put node coordinates (copied into 1D arrays)
    position_x = np.zeros(num_nodes)
    position_y = np.zeros(num_nodes)
    position_z = np.zeros(num_nodes)
    for i in range(0,num_nodes):
        position_x[i] = x[2*i+0]
        position_y[i] = x[2*i+1]
    exo.put_coords(xCoords=position_x,yCoords=position_y,zCoords=position_z)
    
    # put element block info for all elements
    exo.put_elem_blk_info(id=1, elemType='QUAD', numElems=Nelems, numNodesPerElem=4, numAttrsPerElem=0)
    exo.put_elem_connectivity(id=1, connectivity=connectivity, shift_indices=1, chunk_size_in_mb=128)

    # set the number of output node variables and their names
    num_node_variables = 9
    exo.set_node_variable_number(num_node_variables)
    exo.put_node_variable_name("displacement_x", 1)
    exo.put_node_variable_name("displacement_y", 2)
    exo.put_node_variable_name("displacement_z", 3)
    exo.put_node_variable_name("velocity_x",     4)
    exo.put_node_variable_name("velocity_y",     5)
    exo.put_node_variable_name("velocity_z",     6)
    exo.put_node_variable_name("force_x",        7)
    exo.put_node_variable_name("force_y",        8)
    exo.put_node_variable_name("force_z",        9)
    
    # set the number of output element variables and their names
    num_elem_variables = 4
    exo.set_element_variable_number(num_elem_variables)
    exo.put_element_variable_name("stress_xx", 1)
    exo.put_element_variable_name("stress_yy", 2)
    exo.put_element_variable_name("stress_xy", 3)
    exo.put_element_variable_name("pressure",  4)

# ---------------------------------------------------------------------------- #

# Write data to the Exodus file containing particle info for the current time state
def output_state():
    
    # Pre-allocate node data arrays for use during output
    ux = np.zeros(num_nodes)
    uy = np.zeros(num_nodes)
    uz = np.zeros(num_nodes)
    vx = np.zeros(num_nodes)
    vy = np.zeros(num_nodes)
    vz = np.zeros(num_nodes)
    fx = np.zeros(num_nodes)
    fy = np.zeros(num_nodes)
    fz = np.zeros(num_nodes)
    
    # Pre-allocate element data arrays for use during output
    sxx      = np.zeros(num_elems)
    syy      = np.zeros(num_elems)
    sxy      = np.zeros(num_elems)
    pressure = np.zeros(num_elems)
    
    # retrieve the simulation state info at the current time
    step_id = c_int(0)
    time    = c_double(0.0)
    API.get_field_data(step_id,time,ux,uy,vx,vy,fx,fy,sxx,syy,sxy,pressure)
    
    # shift step ID by 1 for consistency with the Exodus file format conventions
    step_id = step_id.value + 1;
    
    # create a new output time state
    exo.put_time(step_id, time)
    
    # write nodal variable values at the current time state
    exo.put_node_variable_values("displacement_x", step_id, ux)
    exo.put_node_variable_values("displacement_y", step_id, uy)
    exo.put_node_variable_values("displacement_z", step_id, uz)
    exo.put_node_variable_values("velocity_x",     step_id, vx)
    exo.put_node_variable_values("velocity_y",     step_id, vy)
    exo.put_node_variable_values("velocity_z",     step_id, vy)
    exo.put_node_variable_values("force_x",        step_id, fx)
    exo.put_node_variable_values("force_y",        step_id, fy)
    exo.put_node_variable_values("force_z",        step_id, fz)
    
    # write element variable values at the current time state
    exo.put_element_variable_values(1, "stress_xx", step_id, sxx)
    exo.put_element_variable_values(1, "stress_yy", step_id, syy)
    exo.put_element_variable_values(1, "stress_xy", step_id, sxy)
    exo.put_element_variable_values(1, "pressure",  step_id, pressure)

# ---------------------------------------------------------------------------- #

# Close the connection to the REMAT library and any Exodus files
def finalize():

    exo.close()

# ---------------------------------------------------------------------------- #
