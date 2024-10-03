from abc import ABC
import numpy as np
from enum import Enum

def cotan(ang):
    return 1/np.tan(ang)

class Orientation(Enum):
    HORIZONTAL = 1
    VERTICAL = 2
class Bound(Enum):
    LOWER = 1
    UPPER = 2



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
    def get_Dn(self):
        Dn = self.x if self.orientation == Orientation.HORIZONTAL else self.y
        return Dn 

class ArbitraryTriangle:
    def __init__(self, x, y, ang_low, ang_high, orientation):
        self.x = x
        self.y = y 
        self.ang_low = ang_low 
        self.ang_high = ang_high
        self.avg_ang = (ang_low+ang_high)/2
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
        x_series = [xc, self.X-xc, self.X-xc, xc]
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

class Interval:
    def __init__(self, offset_lb, offset_ub, lb, ub):
        """
        This is the amount of info needed for the Interval solution

        If both are offset or not offset then we use ub < lb
        If ub is offset but lb is not offset or viceversa then we use [pi, lb]U[ub,lb]
        If ub is offset then it is actually 2pi-ub
        If lb is offset then it is actually 2pi-lb
        """
        self.offset_ub = offset_ub
        self.offset_lb = offset_lb
        self.lb = lb
        self.ub = ub
        self.decreasing = False
        if self.ub < self.lb:
            self.lb = self.ub
    def divide_interval(self, divider):
        divided_ub = self.ub
        self.ub = divider
        return Interval(self.offset_lb, self.offset_ub, divider, divided_ub)
class OffsetInterval(Interval):
    def __init__(self, offset_lb, offset_ub, lb, ub, pivoted=False):
        """
        This is the amount of info needed for the Interval solution

        If both are offset or not offset then we use ub < lb
        If ub is offset but lb is not offset or viceversa then we use [pi, lb]U[ub,lb]
        If ub is offset then it is actually 2pi-ub
        If lb is offset then it is actually 2pi-lb
        """
        super().__init__(offset_lb, offset_ub, lb, ub)
        self.pivoted = pivoted
class ProbabilityCalculator(ABC):
    def __init__(self,fov, beta, h, r):
        self.r = r
        self.sinbeta = np.sin(beta)
        self.h = h
        self.hcos = h*np.cos(beta)
        self.b = np.sqrt(r**2+h**2+2*h*r*np.cos(beta))
        self.a = self.hcos + self.r 
        self.cosfov = np.cos(fov)

class 

class IntervalSolver:
    def __init__(self, Llow, Lhigh, interval: Interval, offset: OffsetInterval):
        self.llow = Llow
        self.lhigh = Lhigh
        self.interval = interval
        self.offset = offset

    def _build_solution_root(self):
        if self.interval.offset_lb:
            if self.interval.offset_ub: 
                self._lb_ub_offset_node()
            else:
                self._lb_offset_node()
        else:
            self._no_offset_node()

    def _lb_ub_offset_node(self):
        if self.offset.offset_lb:
            if self.offset.offset_ub:
                """
                  [max(lb,lbd),min(ub,ubd)]
                """
                max_intervs = self._solve_max_equation(...)
                min_intervs = self._solve_min_equation(...)
                res_intervs = self._check_max_min(min_intervs, max_intervs)
            else:
                """
                  [lb*, min(ub*, ubd)] + [max(lb, lbd), ub*]
                """
                max_intervs = self._solve_max_equation(...)
                min_intervs = self._solve_min_equation(...)
                upper_res_intervs = self._check_max_min(...)
                lower_res_intervs = self._check_max_min(...)
                center_res_intervs = self._check_center(...)

        else:
            """
                [max(lb*,lbd),min(ub*,ubd)]
            """
            max_intervs = self._solve_max_equation(...)
            min_intervs = self._solve_min_equation(...)
            res_intervs = self._check_max_min(min_intervs, max_intervs)

    def _lb_offset_node(self):
        if self.offset.offset_lb:
            if self.offset.offset_ub:
                """
                  ([lb,min(ub,ubd)] + [max(lb,lbd),ub])
                """
                max_intervs = self._solve_max_equation(...)
                min_intervs = self._solve_min_equation(...)
                upper_res_intervs = self._check_max_min(...)
                lower_res_intervs = self._check_max_min(...)
                center_res_intervs = self._check_center(...)


            else:
                """
                    [-2pi, min(ub, ubd)] + [max(lb*, lbd*), 2pi]
                """
                max_intervs = self._solve_max_equation(...)
                min_intervs = self._solve_min_equation(...)
                upper_res_intervs = self._check_max_min(...)
                lower_res_intervs = self._check_max_min(...)
                center_res_intervs = self._check_center(...)

               
        else:
            """
                     [lbd, min(ub, ubd)] + [max(lb*, lbd), ubd]
            """
            max_intervs = self._solve_max_equation(...)
            min_intervs = self._solve_min_equation(...)
            upper_res_intervs = self._check_max_min(...)
            lower_res_intervs = self._check_max_min(...)
            center_res_intervs = self._check_center(...)

    def _no_offset_node(self):
        if self.offset.offset_lb:
            if self.offset.offset_ub:
                """
                  [max(lb,lbd*),min(ub,ubd*)]
                """
                max_intervs = self._solve_max_equation(...)
                min_intervs = self._solve_min_equation(...)
                res_intervs = self._check_max_min(min_intervs, max_intervs)
            else:
                """
                   [lb, min(ub, ubd)] + [max(lb, lbd*), ub]
                """
                max_intervs = self._solve_max_equation(...)
                min_intervs = self._solve_min_equation(...)
                upper_res_intervs = self._check_max_min(...)
                lower_res_intervs = self._check_max_min(...)
                center_res_intervs = self._check_center(...)

        else:
            """
                [max(lb,lbd),min(ub,ubd)]
            """
            max_intervs = self._solve_max_equation(...)
            min_intervs = self._solve_min_equation(...)
            res_intervs = self._check_max_min(min_intervs, max_intervs)



class SecantSolution:
    def __init__(self, Llow, Lhigh):
        self.llow = Llow
        self.lhigh = Lhigh
    def solve_equation(self, f):
        solved = False
        an = self.llow
        bn = self.lhigh
        
        epsilon = 1e-3
        if f(an)*f(bn)>0:
            return self._parabolic_solver(f)
        else:
            return self._solve_equation(f,an,bn)
        

    def _parabolic_solver(self, f):
        an = self.llow
        bn = self.lhigh
        candidate = (an+bn)/2
        rn = candidate 
        if  f(an)*f(candidate) > 0 and f(bn)*f(candidate) > 0:
            """
            Here we create two different potential sides, something like a tree
            One interval which is [an, rn] and other that is [rn, bn].
            Exploiting the properties of this, we will select the middle point and either the lower or higher side depending
            on the function value. 
            Lets say mn = (an+rn)/2
            Then if f(an)<0, f(mn)<0, f(rn)<0, and f(an)<f(mn)<f(rn) then the solution is never in [an, mn]
            Same case if f(bn)<f(mn)<f(rn) holds for [mn, bn]
            If f(bn)<f(rn)<f(mn) then the solution is for sure in [rn, bn] and we need to solely look in this interval using an = rn, rn = mn, ...
            If f(an)<f(rn)<f(mn) then the solution is for sure in [an, rn]
            If f(mn)>0 then you can do the parabolic solver on the case
            """
            solved = False
            while not solved:
                mna = (an+rn)/2
                mnb = (bn+rn)/2
                if f(mna)*f(an) < 0:
                    rn = mna
                    bn = rn
                    break
                elif f(mnb)*f(bn) < 0:
                    an = rn 
                    rn = mnb
                    break
                else:
                    if f(rn) < f(mnb):
                        an = rn 
                        rn = mnb 
                    elif f(rn) < f(mna):
                        an = mna 
                    else:
                        an = mna 
                        bn = mnb
        return self._solve_equation(f, an, rn), self._solve_equation(f, bn, rn)
    def _solve_equation(self, f, lb, ub):
        an = lb
        bn = ub
        epsilon = 1e-3
        solved = False
        while not solved:
            rn = (an+bn)/2
            if f(an)*f(rn)>0:
                an = rn
            if f(bn)*f(rn)>0:
                bn = rn
            if f(rn)<epsilon:
                return rn

def my_equation(x):
    term_one = np.arccos((np.cos(20*np.pi/180)*np.sqrt(x**2+1.551)-1.137)/(x*0.422)) - 1.47
    u = (x+0.1)
    term_two = -np.arccos((np.cos(20*np.pi/180)*np.sqrt(u**2+2.5411)-1.137)/(np.sqrt(u**2+0.994**2)*0.422)) - np.arctan((x*0.994)/(0.1*x+1))
    return term_one-term_two
if __name__=="__main__":
    secant = SecantSolution(0.1, 0.695)
    print(secant.solve_equation(my_equation))


