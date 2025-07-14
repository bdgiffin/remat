from Part import *
import sys
import os
import pyexodus
import numpy as np

class Model:

    def __init__(self):
        self.parts                = []
        self.initial_conditions   = []
        self.boundary_conditions  = []
        self.contact_interactions = []

    def add_part(self, part):
        self.parts.append(part)

    def add_initial_condition(self, node_set, velocity):
        self.initial_conditions.append((node_set,velocity))

    def add_boundary_condition(self, node_set, fixity):
        self.boundary_conditions.append((node_set,fixity))

    def add_contact_interaction(self, node_set, face_set):
        self.contact_interactions.append((node_set,face_set))
    
    def write_exodus(self, filename):
        # get part totals
        num_parts = len(self.parts)
        num_dims = 3
        num_nodes = 0
        num_elements = 0
        for part in self.parts:
            num_dims     = min(num_dims,part.geometry.vertices.shape[1])
            num_nodes    = num_nodes    + part.geometry.vertices.shape[0]
            num_elements = num_elements + part.geometry.connectivity.shape[0]
        coordinates = np.zeros((num_nodes,num_dims))
        index = 0
        for part in self.parts:
            for i in range(part.geometry.vertices.shape[0]):
                coordinates[index,:] = part.geometry.vertices[i,:]
                index = index + 1
            
        # create a new Exodus file
        try:
            os.remove(filename)
        except OSError:
            pass
        exo = pyexodus.exodus(file=filename, mode='w', array_type='numpy', title=None, numDims=num_dims, numNodes=num_nodes, numElems=num_elements, numBlocks=num_parts, numNodeSets=0, numSideSets=0, io_size=8, compression=None)

        # put node coordinates
        exo.put_coords(xCoords=coordinates[:,0],yCoords=coordinates[:,1],zCoords=coordinates[:,2])

        # put element block data
        num_nodes = 0
        for i,part in enumerate(self.parts):
            exo.put_elem_blk_info(id=i+1, elemType='HEX8', numElems=part.geometry.connectivity.shape[0], numNodesPerElem=part.geometry.connectivity.shape[1], numAttrsPerElem=0)

            # put element block 1 connectivity (shift indicies by 1 to ensure compatibility with 1-based indexing in Exodus files, whereas Python uses 0-based indexing)
            exo.put_elem_connectivity(id=i+1, connectivity=part.geometry.connectivity, shift_indices=num_nodes+1, chunk_size_in_mb=128)

            # update node offset for current part
            num_nodes = num_nodes + part.geometry.vertices.shape[0]

        # close the Exodus object (avoids error message)
        exo.close()
    
    def generate_problem(self):
        # get part totals
        num_parts = len(self.parts)
        num_dims = 3
        num_nodes = 0
        num_elements = 0
        num_nodes_per_element = 0
        for part in self.parts:
            num_dims     = min(num_dims,part.geometry.vertices.shape[1])
            num_nodes    = num_nodes    + part.geometry.vertices.shape[0]
            num_elements = num_elements + part.geometry.connectivity.shape[0]
            num_nodes_per_element = max(num_nodes_per_element,part.geometry.connectivity.shape[1])

        # concatenate all nodal coordinates into a single array, and determine part ID offsets
        coordinates = np.zeros((num_nodes,num_dims))
        num_nodes = 0
        for part in self.parts:
            part.geometry.offset_id = num_nodes
            for i in range(part.geometry.vertices.shape[0]):
                coordinates[num_nodes+i,:] = part.geometry.vertices[i,:]
            num_nodes = num_nodes + part.geometry.vertices.shape[0]

        # define initial conditions
        velocities = np.zeros((num_nodes,num_dims))
        for ic in self.initial_conditions:
            for i in ic[0].global_node_ids():
                velocities[i,:] = ic[1]

        # define boundary conditions
        fixity = np.zeros((num_nodes,num_dims),dtype=np.bool_)
        for bc in self.boundary_conditions:
            for i in bc[0].global_node_ids():
                fixity[i,:] = bc[1]

        # concatenate all element connectivities
        connectivity = np.zeros((num_elements,num_nodes_per_element),dtype=np.int32)
        num_elements = 0
        for part in self.parts:
            part_connectivity = part.geometry.global_connectivity()
            for i, elem_connectivity in enumerate(part_connectivity):
                connectivity[num_elements+i,:] = elem_connectivity
            num_elements = num_elements + part.geometry.connectivity.shape[0]

        # define contact interactions
        contacts = []
        for contact in self.contact_interactions:
            contacts.append((contact[0].global_node_ids().astype(np.int32),contact[1].global_node_ids().astype(np.int32)))
        
        # return the data necessary to instantiate the problem
        return coordinates, velocities, fixity, connectivity, contacts
