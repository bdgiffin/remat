import numpy as np

# NodeSet class declaration
class NodeSet:

    def __init__(self,geometry,vertex_ids):
        self.geometry   = geometry
        self.vertex_ids = vertex_ids

    def global_node_ids(self):
        return self.vertex_ids + self.geometry.offset_id
