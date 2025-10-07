from Part import *
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
    
    def generate_problem(self):
        # get part totals
        num_parts = len(self.parts)
        num_dims = 3
        num_nodes = 0
        num_elements = 0
        num_truss = 0
        num_nodes_per_element = 0
        for part in self.parts:
            num_dims     = min(num_dims,part.geometry.vertices.shape[1])
            num_nodes    = num_nodes    + part.geometry.vertices.shape[0]
            if (part.geometry.connectivity.shape[1] == 2):
                num_truss = num_truss + part.geometry.connectivity.shape[0]
            else:
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
        truss_connectivity = np.zeros((num_truss,2),dtype=np.int32)
        num_truss = 0
        num_elements = 0
        for part in self.parts:
            part_connectivity = part.geometry.global_connectivity()
            if (part_connectivity.shape[1] == 2):
                for i, elem_connectivity in enumerate(part_connectivity):
                    truss_connectivity[num_truss+i,:] = elem_connectivity
                num_truss = num_truss + part.geometry.connectivity.shape[0]
            else:
                for i, elem_connectivity in enumerate(part_connectivity):
                    connectivity[num_elements+i,:] = elem_connectivity
                num_elements = num_elements + part.geometry.connectivity.shape[0]

        # define contact interactions
        contacts = []
        for contact in self.contact_interactions:
            contacts.append((contact[0].global_node_ids().astype(np.int32),contact[1].global_node_ids().astype(np.int32)))
        
        # return the data necessary to instantiate the problem
        return coordinates, velocities, fixity, connectivity, contacts, truss_connectivity
