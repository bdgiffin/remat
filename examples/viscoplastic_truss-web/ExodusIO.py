import numpy as np
import sys, platform
import os

# Use REMAT Python API to communicate data to/from Exodus file
import REMAT

# Packages for reading/writing mesh files
if (not sys.platform == "emscripten"):
    import pyexodus

class ExodusIO:
    # The ExodusIO object facilitates reading/writing mesh and field data to/from Exodus files

    def __init__(self):
        self.exo = None
        self.output_step = 0
        
    def __del__(self):
        self.finalize() # ensure that the Exodus file is properly closed
   
    # ---------------------------------------------------------------------------- #

    def create(self,filename):
        # get mesh totals
        num_dim   = REMAT.API.get_num_dim()
        num_nodes = REMAT.API.get_num_entities(b"node")
        num_elems = REMAT.API.get_num_entities(b"element")
        num_truss = REMAT.API.get_num_entities(b"truss")
        num_point = REMAT.API.get_num_entities(b"point")

        # determine the number of element blocks
        num_blocks = 0
        if (num_elems > 0):
            num_blocks = num_blocks + 1
        if (num_truss > 0):
            num_blocks = num_blocks + 1
        if (num_point > 0):
            num_blocks = num_blocks + 1
        
        # attempt to create a new Exodus file
        try:
            os.remove(filename)
        except OSError:
            pass
        self.exo = pyexodus.exodus(file=filename,
                                   mode='w',
                                   array_type='numpy',
                                   title='Produced by REMAT ExodusIO module',
                                   numDims=num_dim,
                                   numNodes=num_nodes,
                                   numElems=(num_elems+num_truss),
                                   numBlocks=num_blocks,
                                   numNodeSets=0,
                                   numSideSets=0,
                                   io_size=0,
                                   compression=None)
    
        # put initial (undeformed) node coordinates (copied into 1D arrays)
        coords = np.zeros((num_nodes,num_dim),dtype=np.double)
        REMAT.API.get_node_coords(coords,False)
        position_x = coords[:,0]
        position_y = coords[:,1]
        position_z = np.zeros(num_nodes)
        self.exo.put_coords(xCoords=position_x,
                            yCoords=position_y,
                            zCoords=position_z)

        # reset the element block counter
        num_blocks = 0
    
        # put element block info for all solid elements
        if (num_elems > 0):
            element_type = 'QUAD'
            num_blocks = num_blocks + 1
            element_block_id = num_blocks
            num_nodes_per_elem = 4
            connectivity = np.zeros((num_elems,num_nodes_per_elem),dtype=np.int32)
            REMAT.API.get_connectivity(b"element",connectivity)
            self.exo.put_elem_blk_info(id=element_block_id,
                                       elemType=element_type,
                                       numElems=num_elems,
                                       numNodesPerElem=num_nodes_per_elem,
                                       numAttrsPerElem=0)
            self.exo.put_elem_connectivity(id=element_block_id,
                                           connectivity=connectivity,
                                           shift_indices=1,
                                           chunk_size_in_mb=128)

        # put element block info for all truss elements
        if (num_truss > 0):
            element_type = 'BAR'
            num_blocks = num_blocks + 1
            element_block_id = num_blocks
            num_nodes_per_elem = 2
            connectivity = np.zeros((num_truss,num_nodes_per_elem),dtype=np.int32)
            REMAT.API.get_connectivity(b"truss",connectivity)
            self.exo.put_elem_blk_info(id=element_block_id,
                                       elemType=element_type,
                                       numElems=num_truss,
                                       numNodesPerElem=num_nodes_per_elem,
                                       numAttrsPerElem=0)
            self.exo.put_elem_connectivity(id=element_block_id,
                                           connectivity=connectivity,
                                           shift_indices=1,
                                           chunk_size_in_mb=128)

        # put element block info for all discrete points
        if (num_point > 0):
            element_type = 'CIRCLE'
            num_blocks = num_blocks + 1
            element_block_id = num_blocks
            num_nodes_per_elem = 1
            connectivity = np.zeros((num_point,num_nodes_per_elem),dtype=np.int32)
            REMAT.API.get_connectivity(b"point",connectivity)
            self.exo.put_elem_blk_info(id=element_block_id,
                                       elemType=element_type,
                                       numElems=num_point,
                                       numNodesPerElem=num_nodes_per_elem,
                                       numAttrsPerElem=0)
            self.exo.put_elem_connectivity(id=element_block_id,
                                           connectivity=connectivity,
                                           shift_indices=1,
                                           chunk_size_in_mb=128)

        # set the number of output global variables and their names
        num_global_variables = REMAT.API.get_num_fields(b"global")
        self.exo.set_global_variable_number(num_global_variables)
        for i in range(0,num_global_variables):
            self.exo.put_global_variable_name(REMAT.API.get_field_name(b"global",i).decode('utf-8'), i+1)

        # set the number of output node variables and their names
        num_node_variables = REMAT.API.get_num_fields(b"node")
        self.exo.set_node_variable_number(num_node_variables)
        for i in range(0,num_node_variables):
            self.exo.put_node_variable_name(REMAT.API.get_field_name(b"node",i).decode('utf-8'), i+1)

        # set the number of output element, truss, and discrete point variables and their names
        num_elem_variables  = REMAT.API.get_num_fields(b"element")
        num_truss_variables = REMAT.API.get_num_fields(b"truss")
        num_point_variables = REMAT.API.get_num_fields(b"point")
        self.exo.set_element_variable_number(num_elem_variables+num_truss_variables+num_point_variables)
        for i in range(0,num_elem_variables):
            self.exo.put_element_variable_name(REMAT.API.get_field_name(b"element",i).decode('utf-8'), i+1)
        for i in range(num_elem_variables,num_elem_variables+num_truss_variables):
            self.exo.put_element_variable_name(REMAT.API.get_field_name(b"truss",i).decode('utf-8'), i+1)
        for i in range(num_elem_variables+num_truss_variables,num_elem_variables+num_truss_variables+num_point_variables):
            self.exo.put_element_variable_name(REMAT.API.get_field_name(b"point",i).decode('utf-8'), i+1)
        
        # Reset the output step counter
        self.output_step = 0

    # ---------------------------------------------------------------------------- #

    # Write data to the Exodus file containing info for the current time state
    def output_state(self):

        # Update the output step counter
        self.output_step = self.output_step + 1
        
        # Create a new output time state
        output_time = REMAT.API.get_time()
        self.exo.put_time(self.output_step, self.output_step)

        # write global variable values at the current time state
        entity_type = b"global"
        num_entities = REMAT.API.get_num_entities(entity_type)
        num_fields = REMAT.API.get_num_fields(entity_type)
        field_data = np.zeros((num_entities,num_fields),dtype=np.double)
        REMAT.API.get_fields(entity_type,field_data)
        for i in range(0,num_fields):
            field_variable = field_data[0,i]
            self.exo.put_global_variable_value(REMAT.API.get_field_name(entity_type,i).decode('utf-8'),
                                               self.output_step,
                                               field_variable)

        # write nodal variable values at the current time state
        entity_type = b"node"
        num_entities = REMAT.API.get_num_entities(entity_type)
        num_fields = REMAT.API.get_num_fields(entity_type)
        field_data = np.zeros((num_entities,num_fields),dtype=np.double)
        REMAT.API.get_fields(entity_type,field_data)
        for i in range(0,num_fields):
            field_variable = field_data[:,i]
            self.exo.put_node_variable_values(REMAT.API.get_field_name(entity_type,i).decode('utf-8'),
                                              self.output_step,
                                              field_variable)

        # reset the element block counter
        num_blocks = 0
            
        # write solid element variable values at the current time state
        entity_type = b"element"
        num_entities = REMAT.API.get_num_entities(entity_type)
        if (num_entities > 0):
            num_blocks = num_blocks + 1
            element_block_id = num_blocks
            num_fields = REMAT.API.get_num_fields(entity_type)
            field_data = np.zeros((num_entities,num_fields),dtype=np.double)
            REMAT.API.get_fields(entity_type,field_data)
            for i in range(0,num_fields):
                field_variable = field_data[:,i]
                self.exo.put_element_variable_values(element_block_id,
                                                     REMAT.API.get_field_name(entity_type,i).decode('utf-8'),
                                                     self.output_step,
                                                     field_variable)

        # write truss element variable values at the current time state
        entity_type = b"truss"
        num_entities = REMAT.API.get_num_entities(entity_type)
        if (num_entities > 0):
            num_blocks = num_blocks + 1
            element_block_id = num_blocks
            num_fields = REMAT.API.get_num_fields(entity_type)
            field_data = np.zeros((num_entities,num_fields),dtype=np.double)
            REMAT.API.get_fields(entity_type,field_data)
            for i in range(0,num_fields):
                field_variable = field_data[:,i]
                self.exo.put_element_variable_values(element_block_id,
                                                     REMAT.API.get_field_name(entity_type,i).decode('utf-8'),
                                                     self.output_step,
                                                     field_variable)

        # write discrete point element variable values at the current time state
        entity_type = b"point"
        num_entities = REMAT.API.get_num_entities(entity_type)
        if (num_entities > 0):
            num_blocks = num_blocks + 1
            element_block_id = num_blocks
            num_fields = REMAT.API.get_num_fields(entity_type)
            field_data = np.zeros((num_entities,num_fields),dtype=np.double)
            REMAT.API.get_fields(entity_type,field_data)
            for i in range(0,num_fields):
                field_variable = field_data[:,i]
                self.exo.put_element_variable_values(element_block_id,
                                                     REMAT.API.get_field_name(entity_type,i).decode('utf-8'),
                                                     self.output_step,
                                                     field_variable)

    # ---------------------------------------------------------------------------- #

    # Close the connection to the Exodus file being written
    def finalize(self):
        self.exo.close()
    
    # ---------------------------------------------------------------------------- #
