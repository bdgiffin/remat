# Module for calling the C/C++ REMAT API functions from Python

# Python package for calling C/C++ functions from Python
from ctypes import CDLL, POINTER, CFUNCTYPE
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

# C-types corresponding to a function that accepts 2 double arguments and returns 1 double
c_function_2d = CFUNCTYPE(c_double, c_double, c_double)

# C-type corresponding to 1D numpy array
ND_POINTER_1 = np.ctypeslib.ndpointer(dtype=np.float64, 
                                      ndim=1,
                                      flags="C")
ND_POINTER_2 = np.ctypeslib.ndpointer(dtype=np.float64, 
                                      ndim=2,
                                      flags="C")
NB_POINTER_2 = np.ctypeslib.ndpointer(dtype=np.bool_, 
                                      ndim=2,
                                      flags="C")
NI_POINTER_1 = np.ctypeslib.ndpointer(dtype=np.int32, 
                                      ndim=1,
                                      flags="C")
NI_POINTER_2 = np.ctypeslib.ndpointer(dtype=np.int32, 
                                      ndim=2,
                                      flags="C")

# Define all C/C++ library API function signatures
API.set_integrator_type.argtypes = [c_char_p]
API.set_integrator_type.restype  = None
API.define_parameter.argtypes = [c_char_p, c_double]
API.define_parameter.restype  = None
API.define_geometry.argtypes = [ND_POINTER_2, ND_POINTER_2, NB_POINTER_2, NI_POINTER_2, c_size_t, c_size_t]
API.define_geometry.restype  = None
API.define_truss_elements.argtypes = [NI_POINTER_2, c_size_t]
API.define_truss_elements.restype  = None
API.define_contact_interaction.argtypes = [NI_POINTER_1, NI_POINTER_2, c_size_t, c_size_t]
API.define_contact_interaction.restype  = None
API.define_point_mass.argtypes = [NI_POINTER_1, ND_POINTER_1, c_size_t]
API.define_point_mass.restype  = None
API.initialize.argtypes = None
API.initialize.restype  = None
API.initialize_variable_properties.argtypes = [c_function_2d]
API.initialize_variable_properties.restype  = None
API.update_state.argtypes = [c_double, c_int]
API.update_state.restype  = c_double
API.get_field_data.argtypes = [ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1, ND_POINTER_1]
API.get_field_data.restype  = c_double

# ---------------------------------------------------------------------------- #

# Generate mesh geometry and initialize the REMAT object prior to initialization
def create_geometry(x,v,fixity,connectivity,contacts,filename=""):

    global num_nodes
    global num_elems
    global exo_file
    num_nodes = x.shape[0]
    num_elems = connectivity.shape[0]
    exo_file  = filename

    # Call REMAT initialization API function
    API.define_geometry(x,v,fixity,connectivity,num_nodes,num_elems)

    # Define contacts
    for contact in contacts:
        num_contact_nodes = contact[0].shape[0]
        num_contact_segs  = contact[1].shape[0]
        API.define_contact_interaction(contact[0],contact[1],num_contact_nodes,num_contact_segs)

# ............................................................................ #

    # If requested, create the Exodus file containing the problem info:
    if (exo_file != ""):

        # create a new Exodus file
        try:
            os.remove(filename)
        except OSError:
            pass
        global exo
        exo = pyexodus.exodus(file=filename, mode='w', array_type='numpy', title='Produced by REMAT module', numDims=2, numNodes=num_nodes, numElems=num_elems, numBlocks=1, numNodeSets=0, numSideSets=0, io_size=0, compression=None)
    
        # put node coordinates (copied into 1D arrays)
        position_x = x[:,0]
        position_y = x[:,1]
        position_z = np.zeros(num_nodes)
        exo.put_coords(xCoords=position_x,yCoords=position_y,zCoords=position_z)
    
        # put element block info for all elements
        exo.put_elem_blk_info(id=1, elemType='QUAD', numElems=num_elems, numNodesPerElem=4, numAttrsPerElem=0)
        exo.put_elem_connectivity(id=1, connectivity=connectivity, shift_indices=1, chunk_size_in_mb=128)

        # set the number of output global variables and their names
        num_global_variables = 4
        exo.set_global_variable_number(num_global_variables)
        exo.put_global_variable_name("elastic_strain_energy", 1)
        exo.put_global_variable_name("kinetic_energy",        2)
        exo.put_global_variable_name("potential_energy",      3)
        exo.put_global_variable_name("total_energy",          4)

        # set the number of output node variables and their names
        num_node_variables = 15
        exo.set_node_variable_number(num_node_variables)
        exo.put_node_variable_name("displacement_X",       1)
        exo.put_node_variable_name("displacement_Y",       2)
        exo.put_node_variable_name("displacement_Z",       3)
        exo.put_node_variable_name("velocity_X",           4)
        exo.put_node_variable_name("velocity_Y",           5)
        exo.put_node_variable_name("velocity_Z",           6)
        exo.put_node_variable_name("force_X",              7)
        exo.put_node_variable_name("force_Y",              8)
        exo.put_node_variable_name("force_Z",              9)
        exo.put_node_variable_name("dual_displacement_X", 10)
        exo.put_node_variable_name("dual_displacement_Y", 11)
        exo.put_node_variable_name("dual_displacement_Z", 12)
        exo.put_node_variable_name("dual_velocity_X",     13)
        exo.put_node_variable_name("dual_velocity_Y",     14)
        exo.put_node_variable_name("dual_velocity_Z",     15)
        
        # set the number of output element variables and their names
        num_elem_variables = 8
        exo.set_element_variable_number(num_elem_variables)
        exo.put_element_variable_name("stress_XX", 1)
        exo.put_element_variable_name("stress_YY", 2)
        exo.put_element_variable_name("stress_ZZ", 3)
        exo.put_element_variable_name("stress_XY", 4)
        exo.put_element_variable_name("stress_XZ", 5)
        exo.put_element_variable_name("stress_YZ", 6)
        exo.put_element_variable_name("pressure",  7)
        exo.put_element_variable_name("stiffness_scaling_factor", 8)
        
        # Define the output step counter
        global output_step
        output_step = 0

# ---------------------------------------------------------------------------- #

# Write data to the Exodus file containing particle info for the current time state
def output_state():

    if (exo_file != ""):

        # Update the output step counter
        global output_step
        output_step = output_step + 1

        # Pre-allocate node data arrays for use during output
        ux      = np.zeros(num_nodes)
        uy      = np.zeros(num_nodes)
        vx      = np.zeros(num_nodes)
        vy      = np.zeros(num_nodes)
        fx      = np.zeros(num_nodes)
        fy      = np.zeros(num_nodes)
        dual_ux = np.zeros(num_nodes)
        dual_uy = np.zeros(num_nodes)
        dual_vx = np.zeros(num_nodes)
        dual_vy = np.zeros(num_nodes)
        z       = np.zeros(num_nodes)

        # Pre-allocate element data arrays for use during output
        sxx      = np.zeros(num_elems)
        syy      = np.zeros(num_elems)
        sxy      = np.zeros(num_elems)
        pressure = np.zeros(num_elems)
        stiffness_scaling_factor = np.zeros(num_elems)
        zz       = np.zeros(num_elems)

        # Pre-allocate system state data
        system_state = np.zeros(4)

        # retrieve the simulation state info at the current time
        time = API.get_field_data(ux,uy,vx,vy,fx,fy,dual_ux,dual_uy,dual_vx,dual_vy,sxx,syy,sxy,pressure,stiffness_scaling_factor,system_state)

        # create a new output time state
        exo.put_time(output_step, output_step)

        # write global variable values at the current time state
        exo.put_global_variable_value("elastic_strain_energy", output_step, system_state[0])
        exo.put_global_variable_value("kinetic_energy",        output_step, system_state[1])
        exo.put_global_variable_value("potential_energy",      output_step, system_state[2])
        exo.put_global_variable_value("total_energy",          output_step, system_state[3])

        # write nodal variable values at the current time state
        exo.put_node_variable_values("displacement_X",      output_step,      ux)
        exo.put_node_variable_values("displacement_Y",      output_step,      uy)
        exo.put_node_variable_values("displacement_Z",      output_step,       z)
        exo.put_node_variable_values("velocity_X",          output_step,      vx)
        exo.put_node_variable_values("velocity_Y",          output_step,      vy)
        exo.put_node_variable_values("velocity_Z",          output_step,       z)
        exo.put_node_variable_values("force_X",             output_step,      fx)
        exo.put_node_variable_values("force_Y",             output_step,      fy)
        exo.put_node_variable_values("force_Z",             output_step,       z)
        exo.put_node_variable_values("dual_displacement_X", output_step, dual_ux)
        exo.put_node_variable_values("dual_displacement_Y", output_step, dual_uy)
        exo.put_node_variable_values("dual_displacement_Z", output_step,       z)
        exo.put_node_variable_values("dual_velocity_X",     output_step, dual_vx)
        exo.put_node_variable_values("dual_velocity_Y",     output_step, dual_vy)
        exo.put_node_variable_values("dual_velocity_Z",     output_step,       z)

        # write element variable values at the current time state
        exo.put_element_variable_values(1, "stress_XX", output_step,      sxx)
        exo.put_element_variable_values(1, "stress_YY", output_step,      syy)
        exo.put_element_variable_values(1, "stress_ZZ", output_step,       zz)
        exo.put_element_variable_values(1, "stress_XY", output_step,      sxy)
        exo.put_element_variable_values(1, "stress_XZ", output_step,       zz)
        exo.put_element_variable_values(1, "stress_YZ", output_step,       zz)
        exo.put_element_variable_values(1, "pressure",  output_step, pressure)
        exo.put_element_variable_values(1, "stiffness_scaling_factor",  output_step, stiffness_scaling_factor)

# ---------------------------------------------------------------------------- #

# Define spatially varying properties
def define_variable_properties(py_function_xy):

    # Convert python function to C callback
    c_function_xy = c_function_2d(py_function_xy)

    # Call REMAT API function to adjust variable properties
    API.initialize_variable_properties(c_function_xy)
    
# ---------------------------------------------------------------------------- #

# Request the deformed coordinates of all nodes of all elements
def deform_geometry(x):
    
    # Pre-allocate node data arrays for use during output
    ux      = np.zeros(num_nodes)
    uy      = np.zeros(num_nodes)
    vx      = np.zeros(num_nodes)
    vy      = np.zeros(num_nodes)
    fx      = np.zeros(num_nodes)
    fy      = np.zeros(num_nodes)
    dual_ux = np.zeros(num_nodes)
    dual_uy = np.zeros(num_nodes)
    dual_vx = np.zeros(num_nodes)
    dual_vy = np.zeros(num_nodes)
    
    # Pre-allocate element data arrays for use during output
    sxx      = np.zeros(num_elems)
    syy      = np.zeros(num_elems)
    sxy      = np.zeros(num_elems)
    pressure = np.zeros(num_elems)
    stiffness_scaling_factor = np.zeros(num_elems)

    # Pre-allocate system state data
    system_state = np.zeros(4)
    
    # retrieve the simulation state info at the current time
    time = API.get_field_data(ux,uy,vx,vy,fx,fy,dual_ux,dual_uy,dual_vx,dual_vy,sxx,syy,sxy,pressure,stiffness_scaling_factor,system_state)
    
    # sum displacements to nodal coordinates
    coordinates = np.zeros((num_nodes,2))
    for i in range(0,num_nodes):
        coordinates[i,0] = x[i,0] + ux[i]
        coordinates[i,1] = x[i,1] + uy[i]

    return coordinates, pressure, system_state

# ---------------------------------------------------------------------------- #

# Close the connection to the REMAT library and any Exodus files
def finalize():

    exo.close()

# ---------------------------------------------------------------------------- #
