import numpy as np

# FaceSet class declaration
class FaceSet:

    def __init__(self,geometry,connectivity):
        self.geometry     = geometry
        self.connectivity = connectivity

    def global_node_ids(self):
        return self.connectivity + self.geometry.offset_id
