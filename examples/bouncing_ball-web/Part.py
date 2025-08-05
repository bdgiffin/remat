from Geometry import *
from Material import *

class Part:
    # The part object consists of the geometry, element type, and material

    def __init__(self,geometry,material):
        self.geometry = geometry
        self.material = material
    
