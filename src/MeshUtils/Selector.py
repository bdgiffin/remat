from abc import ABC, abstractmethod
import math

class Selector(ABC):
    @abstractmethod
    def __call__(self,vertex):
        pass

class Select_all(Selector):
    def __call__(self,vertex):
        return True

class Select_X_eq(Selector):
    def __init__(self, x, tol=1.0e-6):
        self.x   = x
        self.tol = tol
    
    def __call__(self,vertex):
        return (abs(self.x - vertex[0]) < self.tol)

class Select_Y_eq(Selector):
    def __init__(self, y, tol=1.0e-6):
        self.y   = y
        self.tol = tol
    
    def __call__(self,vertex):
        return (abs(self.y - vertex[1]) < self.tol)

class Select_Z_eq(Selector):
    def __init__(self, z, tol=1.0e-6):
        self.z   = z
        self.tol = tol
    
    def __call__(self,vertex):
        return (abs(self.z - vertex[2]) < self.tol)

class Select_XY_eq(Selector):
    def __init__(self, x, y, tol=1.0e-6):
        self.x   = x
        self.y   = y
        self.tol = tol
    
    def __call__(self,vertex):
        return (math.sqrt((self.x - vertex[0])*(self.x - vertex[0])+(self.y - vertex[1])*(self.y - vertex[1])) < self.tol)

class Select_XY_window(Selector):
    def __init__(self, lower_left, upper_right, tol=1.0e-6):
        self.x_min = lower_left[0]  - tol
        self.x_max = upper_right[0] - tol
        self.y_min = lower_left[1]  + tol
        self.y_max = upper_right[1] + tol
    
    def __call__(self,vertex):
        return ((vertex[0] > self.x_min) and (vertex[0] < self.x_max) and (vertex[1] > self.y_min) and (vertex[1] < self.y_max))
