from abc import ABC
import numpy as np
from enum import Enum

def cotan(ang):
    return 1/np.tan(ang)

class Orientation(Enum):
    HORIZONTAL = 1
    VERTICAL = 2
    
class RecTriangle:
    def __init__(self, x, y, ang_crt, orientation):
        self.x = x
        self.y = y
        if orientation == Orientation.HORIZONTAL:
            self.ang_low = 0
            self.ang_high = ang_crt
        else:
            self.ang_low = ang_crt
            self.ang_high = np.pi/2
        self.orientation = orientation
        self.ang_crt = ang_crt
    def change_ang(self, new_ang):
        if self.orientation == Orientation.HORIZONTAL:
            self.ang_low = new_ang
        else:
            self.ang_high = new_ang
    def reset_ang(self):
        self.change_ang(self.ang_crt)
    def get_area(self):
        return self.x*self.y/2
    def __str__(self):
        return f'Triangle: {self.x}, {self.y}, {self.ang_crt}'

class ArbitraryTriangle:
    def __init__(self, x, y, ang_low, ang_high, orientation):
        self.x = x
        self.y = y 
        self.ang_low = ang_low 
        self.ang_high = ang_high
        self.orientation = orientation
        self._max_radius()
    def _max_radius(self):
        avg_ang = (self.ang_low+self.ang_high)/2
        
        if self.orientation == Orientation.HORIZONTAL:
            self.max_r = abs(self.x/np.cos(avg_ang))
        else:
            self.max_r = abs(self.y/np.sin(avg_ang))
    def get_area(self):
        if self.orientation == Orientation.HORIZONTAL:
            low_tr = abs(self.x**2*np.tan(self.ang_low))
            high_tr = abs(self.x**2*np.tan(self.ang_high))
            return abs(high_tr-low_tr)
        else:
            low_tr = abs(self.y**2*cotan(self.ang_low))
            high_tr = abs(self.y**2*cotan(self.ang_high))
            return abs(low_tr-high_tr)
    def __str__(self):
        return f'Triangle: {self.x}, {self.y}, {self.max_r}, {self.get_area()}, {self.ang_low}, {self.ang_high}, {self.orientation}'

class Rectangle(ABC):
    """
       A rectangle is just a space with X and Y 
    """
    def __init__(self, X, Y):
    
    
        self.X = X
        self.Y = Y
    def __iter__(self):
        for trig in self.triangles:
            yield trig
    def __str__(self):
        str_out = ""
        for trig in self.triangles:
            str_out += str(trig)+"\n"
        return str_out

class EightRectangle(Rectangle):
    """
    An eight rectangle is a rectangle divided in eight
    sections. 
    """
    def __init__(self, X, Y, x_center, y_center):
        super().__init__(X, Y)
        self._gen_triangles(x_center, y_center)
    
    def _gen_triangles(self, xc, yc):
        triangles = []
        x_series = [self.X-xc, xc, xc, self.X-xc]
        y_series = [self.Y-yc, self.Y-yc, yc, yc]
        orientation_series = [Orientation.HORIZONTAL, Orientation.VERTICAL]
        for i, x in enumerate(x_series):
            for j in range(2):
                y = y_series[i]
                orientation = orientation_series[j]
                ang_crt = np.arctan(y/x)
                triangles.append(RecTriangle(x,y,ang_crt,orientation))
        self.triangles = triangles

class UniformRectangle(Rectangle):
    """
    A rectangle divided by over N triangles with same angle. 
    """
    def __init__(self, X, Y, x_center, y_center, N):
        super().__init__(X, Y)
        self.N = N
        self._gen_triangles(x_center, y_center, N)
    
    def _gen_triangles(self, xc, yc, N):
        triangles = []
        x_series = [self.X-xc, xc, xc, self.X-xc]
        y_series = [self.Y-yc, self.Y-yc, yc, yc]
        orientation_series = [Orientation.HORIZONTAL, Orientation.VERTICAL]
        delta = 2*np.pi/N
        curr = 0
        rev_ori = False
        for i, x in enumerate(x_series):
            y = y_series[i]
            if not rev_ori:
                ang_crt = np.arctan(y/x)+curr
            else:
                ang_crt = curr+np.pi/2-np.arctan(y/x)
            for j in range(2):
                if not j:
                    ang_seq = np.arange(curr, ang_crt+delta, delta)
                    if ang_seq[-1]>ang_crt:
                        ang_seq[-1] = ang_crt
                else:
                    ang_seq = np.arange(ang_crt, curr+np.pi/2+delta, delta)
                    if ang_seq[-1]>curr+np.pi/2:
                        ang_seq[-1] = curr+np.pi/2
                if not rev_ori:
                    orientation = orientation_series[j]
                else:
                    orientation = orientation_series[1-j]
                try:
                    for i, ang in enumerate(ang_seq[:-1]):
                        triangles.append(ArbitraryTriangle(x,y,ang_seq[i],ang_seq[i+1],orientation))
                except IndexError:
                    continue
            curr += np.pi/2
            rev_ori = not rev_ori
        self.triangles = triangles

class IntegrationLimit:
    def __init__(self, low, high, const) -> None:
        self.low = low
        self.high = high
        self.const = const
    def set_high(self, high):
        self.high = high
    def set_low(self, low):
        self.low = low
    def sort_radius(self):
        if self.low is None:
            return self.high
        return self.low
    def __str__(self):
        return f'Low: {self.low}, High: {self.high}, type: {self.const}'

class ProbabilityCalculator(ABC):
    def __init__(self,fov, beta, h, r):
        self.r = r
        self.sinbeta = np.sin(beta)
        self.h = h
        self.hcos = h*np.cos(beta)
        self.b = np.sqrt(r**2+h**2+2*h*r*np.cos(beta))
        self.a = self.hcos + self.r 
        self.cosfov = np.cos(fov)


if __name__=="__main__":
    trg = UniformRectangle(5,3,1,1,360)

