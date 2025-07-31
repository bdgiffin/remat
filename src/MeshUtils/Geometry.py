from NodeSet import *
from FaceSet import *
from Selector import *
from scipy.spatial.transform import Rotation as Rot
import numpy as np
from abc import ABC, abstractmethod

# The total number of created Geometry objects
geometry_count = 0

# Geometry class declaration
class Geometry(ABC):
    # The geometry object consists of an independent collection of elements

    def __init__(self,vertices,connectivity):
        self.vertices     = vertices
        self.connectivity = connectivity
        global geometry_count
        self.geometry_id = geometry_count
        self.offset_id = 0
        geometry_count = geometry_count + 1
    
    def union(self,geometry2,tol=1.0e-6):
        num_unique = 0
        vertex_id_map = np.full(geometry2.vertices.shape[0],-1,dtype=int)
        for j in range(geometry2.vertices.shape[0]):
            for i in range(self.vertices.shape[0]):
                if (np.linalg.norm(self.vertices[i,:] - geometry2.vertices[j,:]) < tol):
                    vertex_id_map[j] = i
            if (vertex_id_map[j] == -1):
                vertex_id_map[j] = self.vertices.shape[0] + num_unique
                num_unique = num_unique + 1
        shape = self.vertices.shape
        self.vertices = np.resize(self.vertices,(shape[0]+num_unique,shape[1]))
        num_unique = 0
        for j in range(geometry2.vertices.shape[0]):
            if (vertex_id_map[j] >= shape[0]):
                self.vertices[shape[0]+num_unique,:] = geometry2.vertices[j,:]
                num_unique = num_unique + 1

        shape = self.connectivity.shape
        self.connectivity = np.resize(self.connectivity,(shape[0]+geometry2.connectivity.shape[0],shape[1]))
        for j in range(geometry2.connectivity.shape[0]):
            self.connectivity[shape[0]+j,:] = vertex_id_map[geometry2.connectivity[j,:]]
    
    def select_nodes(self,selector):
        vertex_ids = []
        for i in range(self.vertices.shape[0]):
            if (selector(self.vertices[i,:])):
                vertex_ids.append(i)
        return NodeSet(self,np.array(vertex_ids,dtype=int))
    
    def select_faces(self,selector):
        selected = np.full(self.vertices.shape[0], False)
        faces = self.face_topology()
        for i in range(self.vertices.shape[0]):
            if (selector(self.vertices[i,:])):
                selected[i] = True
        selected_faces = []
        for e in range(self.connectivity.shape[0]):
            for f in range(faces.shape[0]):
                # only select non-degenerate faces
                if (len(np.unique(self.connectivity[e,faces[f,:]])) > (self.vertices.shape[1]-1)):
                    if (all(selected[self.connectivity[e,faces[f,:]]])):
                        selected_faces.append(self.connectivity[e,faces[f,:]].tolist())
        return FaceSet(self,np.array(selected_faces,dtype=int))

    def select_boundary_nodes(self):
        edge_counts = {}
        for element in self.connectivity:
            # Extract edges of the element
            if (element[2] == element[3]):
                edges = [(element[0], element[1]), (element[1], element[2]), (element[2], element[0])]
            else:
                edges = [(element[0], element[1]), (element[1], element[2]), (element[2], element[3]), (element[3], element[0])]
            for edge in edges:
                # Ensure consistent order for edge (e.g., (min_idx, max_idx))
                orientation = +1
                if (edge[0] > edge[1]):
                    orientation = -1
                sorted_edge = tuple(sorted(edge))
                edge_counts[sorted_edge] = (edge_counts.get(sorted_edge, (0,orientation))[0] + 1, orientation)

        boundary_edges = []
        for edge, count_and_orientation in edge_counts.items():
            if count_and_orientation[0] == 1:
                if count_and_orientation[1] == +1:
                    boundary_edges.append(edge)
                elif count_and_orientation[1] == -1:
                    boundary_edges.append((edge[1], edge[0]))
        boundary_edges = np.array(boundary_edges,dtype=int)
        boundary_nodes = np.unique(boundary_edges.flatten())
        return NodeSet(self,boundary_nodes)

    def select_boundary_faces(self):
        edge_counts = {}
        for element in self.connectivity:
            # Extract edges of the element
            if (element[2] == element[3]):
                edges = [(element[0], element[1]), (element[1], element[2]), (element[2], element[0])]
            else:
                edges = [(element[0], element[1]), (element[1], element[2]), (element[2], element[3]), (element[3], element[0])]
            for edge in edges:
                # Ensure consistent order for edge (e.g., (min_idx, max_idx))
                orientation = +1
                if (edge[0] > edge[1]):
                    orientation = -1
                sorted_edge = tuple(sorted(edge))
                edge_counts[sorted_edge] = (edge_counts.get(sorted_edge, (0,orientation))[0] + 1, orientation)

        boundary_edges = []
        for edge, count_and_orientation in edge_counts.items():
            if count_and_orientation[0] == 1:
                if count_and_orientation[1] == +1:
                    boundary_edges.append(edge)
                elif count_and_orientation[1] == -1:
                    boundary_edges.append((edge[1], edge[0]))
        boundary_edges = np.array(boundary_edges,dtype=int)
        return FaceSet(self,boundary_edges)
    
    def translate(self,distance):
        for i in range(self.vertices.shape[0]):
            self.vertices[i,:] = self.vertices[i,:] + distance
    
    def rotate(self,base_point,axis,angle):
        axis = axis/np.linalg.norm(axis)
        rotation = Rot.from_rotvec(axis*angle)
        for i in range(self.vertices.shape[0]):
            radius = self.vertices[i,:] - base_point
            self.vertices[i,:] = base_point + rotation.apply(radius)
    
    def mirror(self,base_point,plane_normal,copy=False):
        if (copy):
            mirrored_geometry = self.copy()
            mirrored_geometry.mirror(base_point,plane_normal,copy=False)
            self.union(mirrored_geometry)
        else:
            for i in range(self.vertices.shape[0]):
                shifted_node = self.vertices[i,:] - base_point
                dist = np.dot(shifted_node,plane_normal)
                self.vertices[i,:] = base_point + shifted_node - 2*dist*plane_normal
            self.connectivity = self.connectivity[:,self.inverted_topology()]

    def global_connectivity(self):
        return self.connectivity + self.offset_id

    @abstractmethod               
    def copy(self):
        pass

    @abstractmethod
    def face_topology(self):
        pass

    @abstractmethod
    def inverted_topology(self):
        pass

class HexGeometry(Geometry):

    def __init__(self,vertices,connectivity):
        self.vertices     = vertices
        self.connectivity = connectivity
        global geometry_count
        self.geometry_id = geometry_count
        self.offset_id = 0
        geometry_count = geometry_count + 1
                    
    def copy(self):
        return HexGeometry(self.vertices.copy(),self.connectivity.copy())

    def face_topology(self):
        return np.array([[3,2,1,0],[4,5,6,7],[0,1,5,4],[1,2,6,5],[2,3,7,6],[3,0,4,7]],dtype=int)

    def inverted_topology(self):
        return np.array([4,5,6,7,0,1,2,3],dtype=int)

class QuadGeometry(Geometry):

    def __init__(self,vertices,connectivity):
        self.vertices     = vertices
        self.connectivity = connectivity
        global geometry_count
        self.geometry_id = geometry_count
        self.offset_id = 0
        geometry_count = geometry_count + 1
                    
    def copy(self):
        return QuadGeometry(self.vertices.copy(),self.connectivity.copy())

    def face_topology(self):
        return np.array([[0,1],[1,2],[2,3],[3,0]],dtype=int)

    def inverted_topology(self):
        return np.array([3,2,1,0],dtype=int)

class TrussGeometry(Geometry):

    def __init__(self,vertices,connectivity):
        self.vertices     = vertices
        self.connectivity = connectivity
        global geometry_count
        self.geometry_id = geometry_count
        self.offset_id = 0
        geometry_count = geometry_count + 1
                    
    def copy(self):
        return TrussGeometry(self.vertices.copy(),self.connectivity.copy())

    def face_topology(self):
        return np.array([[0,1]],dtype=int)

    def inverted_topology(self):
        return np.array([1,0],dtype=int)
