# Module for calling the C/C++ REMAT API functions from Python

# Python package for calling C/C++ functions from Python
from ctypes import CDLL, POINTER, CFUNCTYPE
from ctypes import c_size_t, c_double, c_int, c_char_p, c_bool

# Other needed Python packages
import sys, platform
import os
import numpy as np

# ---------------------------------------------------------------------------- #

# REMAT dynamic library linking:

# Load the pre-compiled external C/C++ "shared object" libraries
if sys.platform == "emscripten":
    library_name = './libREMAT.wasm'
elif sys.platform == "darwin":
    library_name = os.path.join(os.path.dirname(__file__), 'libREMAT.dylib')
elif sys.platform == "win32":
    library_name = os.path.join(os.path.dirname(__file__), 'libREMAT.dll')
else:
    library_name = os.path.join(os.path.dirname(__file__), 'libREMAT.so')
API = CDLL(library_name)

# ---------------------------------------------------------------------------- #

# Define types to convert Numpy arrays into C arrays:

# C-types corresponding to pointers to int/double values
c_int_p    = POINTER(c_int)
c_double_p = POINTER(c_double)

# C-types corresponding to a function that accepts 2 double arguments and returns 1 double
c_function_2d = CFUNCTYPE(c_double, c_double, c_double)

# C-type corresponding to 1D numpy array
ND_POINTER_1 = np.ctypeslib.ndpointer(dtype=np.double,
                                      ndim=1,
                                      flags="C")
ND_POINTER_2 = np.ctypeslib.ndpointer(dtype=np.double,
                                      ndim=2,
                                      flags="C")
NB_POINTER_1 = np.ctypeslib.ndpointer(dtype=np.bool_,
                                      ndim=1,
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

# ---------------------------------------------------------------------------- #

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
API.get_num_dim.argtypes      = None
API.get_num_dim.restype       = c_int
API.get_num_entities.argtypes = [c_char_p]
API.get_num_entities.restype  = c_int
API.get_node_coords.argtypes  = [ND_POINTER_2,c_bool]
API.get_node_coords.restype   = None
API.get_connectivity.argtypes = [c_char_p,NI_POINTER_2]
API.get_connectivity.restype  = None
API.get_num_fields.argtypes   = [c_char_p]
API.get_num_fields.restype    = c_int
API.get_field_name.argtypes   = [c_char_p, c_int]
API.get_field_name.restype    = c_char_p
API.get_fields.argtypes       = [c_char_p, ND_POINTER_2]
API.get_fields.restype        = None
API.get_time.argtypes         = None
API.get_time.restype          = c_double

# ---------------------------------------------------------------------------- #

# Generate mesh geometry and initialize the REMAT object prior to initialization
def create_geometry(x,v,fixity,connectivity,contacts,truss_connectivity):
    num_nodes = x.shape[0]
    num_elems = connectivity.shape[0]
    num_truss = truss_connectivity.shape[0]

    # Call REMAT initialization API function
    API.define_geometry(x,v,fixity,connectivity,num_nodes,num_elems)
    
    # Conditionally define truss elements
    if (num_truss > 0):
        API.define_truss_elements(truss_connectivity,num_truss)

    # Define contacts
    for contact in contacts:
        num_contact_nodes = contact[0].shape[0]
        num_contact_segs  = contact[1].shape[0]
        API.define_contact_interaction(contact[0],contact[1],num_contact_nodes,num_contact_segs)

# ---------------------------------------------------------------------------- #

# Define spatially varying properties
def define_variable_properties(py_function_xy):

    # Convert python function to C callback
    c_function_xy = c_function_2d(py_function_xy)

    # Call REMAT API function to adjust variable properties
    API.initialize_variable_properties(c_function_xy)
    
# ---------------------------------------------------------------------------- #

# Request the specified field by entity type and field name
def get_field(entity_type,field_name):
    # Attempt to find the indicated field
    num_entities = API.get_num_entities(entity_type)
    num_fields = API.get_num_fields(entity_type)
    field_data = np.zeros((num_entities,num_fields),dtype=np.double)
    API.get_fields(entity_type,field_data)
    for i in range(0,num_fields):
        if (field_name == API.get_field_name(entity_type,i).decode('utf-8')):
            return field_data[:,i]

    # Return nothing if the indicated field does not exist
    return None

# ---------------------------------------------------------------------------- #
