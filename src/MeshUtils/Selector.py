from abc import ABC, abstractmethod

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
