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
    def __init__(self, offset_lb, offset_ub, lb, ub, upper_func = None, lower_func = None, pi_interval = False):
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
        self.upper_func = upper_func
        self.lower_func = lower_func
        self.pi_interval = pi_interval
        self.decreasing = False
        if self.ub < self.lb:
            self.lb = self.ub
    def divide_interval(self, divider):
        divided_ub = self.ub
        self.ub = divider
        return Interval(self.offset_lb, self.offset_ub, divider, divided_ub)
    def inverse_divided_interval(self, divider):
        divided_lb = self.lb
        self.lb = divider
        return Interval(self.offset_lb, self.offset_ub, divided_lb, divider)
    def get_upper_bound_interval(self):
        ret_int = FunctionalInterval(self.offset_lb, self.offset_ub, self.lb, self.ub)
        ret_int.set_function(self.upper_func)
        return ret_int 
    def get_lower_bound_interval(self):
        ret_int = FunctionalInterval(self.offset_lb, self.offset_ub, self.lb, self.ub)
        ret_int.set_function(self.lower_func)
        return ret_int 
class FunctionalInterval(Interval):
    def __init__(self, offset_lb, offset_ub, lb, ub):
        super().__init__(offset_lb, offset_ub, lb, ub)
        self.function = None
    def set_function(self, func):
        self.function = func

class OffsetInterval(Interval):
    def __init__(self, offset_lb, offset_ub, lb, ub, pivoted=False, upper_func = None, lower_func = None):
        """
        This is the amount of info needed for the Interval solution

        If both are offset or not offset then we use ub < lb
        If ub is offset but lb is not offset or viceversa then we use [pi, lb]U[ub,lb]
        If ub is offset then it is actually 2pi-ub
        If lb is offset then it is actually 2pi-lb
        """
        super().__init__(offset_lb, offset_ub, lb, ub, upper_func, lower_func)
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


class IntervalSolver:
    def __init__(self, Llow, Lhigh, interval: Interval, offset: OffsetInterval):
        self.llow = Llow
        self.lhigh = Lhigh
        self.interval = interval
        self.offset = offset
        self.constant_interval = Interval(False, False, Llow, Lhigh, upper_bound = lambda x: np.pi, lower_bound = lambda x: -np.pi, pi_interval = True)
    def _build_solution_root(self):
        if self.interval.offset_lb:
            if self.interval.offset_ub: 
                self._lb_ub_offset_node()
            else:
                self._lb_offset_node()
        else:
            self._no_offset_node()
    def _interval_sorter_generation(self, first_set, second_set):
        output_sets = []
        j = 0
        for i, interv in enumerate(first_set):
            if j < len(second_set):
                available_sc = True
            while available_sc:   
                sc_interv = second_set[j]
                if interv.ub < sc_interv.ub:
                    new_interval = sc_interv.inverse_divided_interval(interv.ub)
                    output_sets.append([new_interval, interv])
                    break
                elif interv.ub == sc_interv.ub:
                    output_sets.append([sc_interv, interv])
                    j+=1
                    break
                else:
                    new_interval = interv.inverse_divided_interval(sc_interv.ub)
                    output_sets.append([sc_interv, new_interval])
                    j+=1
                available_sc = j < len(second_set)
    def _min_max_finder(self):
        max_interval_one, max_interval_two, breaking_point = self._solve_max_equation(self.llow, self.lhigh, self.interval.get_lower_bound_interval(), 
        self.offset.get_lower_bound_interval())
        lower_set = [max_interval_one]
        if breaking_point is not None:
            max_interval_one.ub = breaking_point
            max_interval_two.lb = breaking_point
            lower_set.append(max_interval_two)
            
        min_intervs_one, min_intervs_two, breaking_point = self._solve_min_equation(self.llow, self.lhigh, self.interval.get_upper_bound_interval(), 
        self.offset.get_upper_bound_interval())
        higher_set = [min_intervs_one]
        if breaking_point is not None:
            min_intervs_one.ub = breaking_point
            min_intervs_two.lb = breaking_point
            higher_set.append(min_intervs_two) 
        output_sets = self._interval_sorter_generation(higher_set, lower_set)
        return output_sets  
    def _lb_min_finder(self, lower_bound):
        min_intervs_one, min_intervs_two, breaking_point = self._solve_min_equation(self.llow, self.lhigh, self.interval.get_upper_bound_interval(), 
        self.offset.get_upper_bound_interval())
        higher_set = [min_intervs_one]
        if breaking_point is not None:
            min_intervs_one.ub = breaking_point
            min_intervs_two.lb = breaking_point
            higher_set.append(min_intervs_two) 
        res_intervs = [self._check_max_min(maxi.lb, maxi.ub, lower_bound.get_lower_bound_interval(), maxi) for maxi in higher_set]
        return res_intervs
    def _ub_max_finder(self, upper_bound):
        max_interval_one, max_interval_two, breaking_point = self._solve_max_equation(self.llow, self.lhigh, self.interval.get_lower_bound_interval(), 
        self.offset.get_lower_bound_interval())
        lower_set = [max_interval_one]
        if breaking_point is not None:
            max_interval_one.ub = breaking_point
            max_interval_two.lb = breaking_point
            lower_set.append(max_interval_two)
        res_intervs = [self._check_max_min(mini.lb, mini.ub, mini, upper_bound.get_upper_bound_interval()) for mini in higher_set]
        return res_intervs
    def _interval_finder_merger(self, res_ub, res_lb):
        pivot = None
        j = 0
        output_array = []
        for i, interv in enumerate(res_ub):
            cap = j < len(res_lb)
            if j >= len(res_lb):
                output_array.append([None, None, interv[0], interv[1]])
            while not_capped:
                if interv[0].lb > res_lb[j][0].lb:
                    if res_lb[j][0].ub < interv[0].lb:
                        lower_lb_solo = res_lb[j][0].inverse_divide_interval(interv[0].lb)
                        upper_lb_solo = res_lb[j][1].inverse_divide_interval(interv[0].lb)
                        output_array.append([lower_lb_solo,upper_lb_solo, None, None])
                    else:
                        output_array.append([res_lb[j][0],res_lb[j][1], None, None])
                        j+=1
                        cap = j < len(res_lb)
                    if interv[0].ub < res_lb[j][0].ub:
                        lower_lb_solo = res_lb[j][0].inverse_divide_interval(interv[0].ub)
                        upper_lb_solo = res_lb[j][1].inverse_divide_interval(interv[0].ub)
                        output_array.append([lower_lb_solo,upper_lb_solo, interv[0], interv[1]])
                        break

                    else:
                        lower_lb_paired = interv[j][0].inverse_divide_interval(res_lb[j][0].ub)
                        upper_lb_paired = interv[j][1].inverse_divide_interval(res_lb[j][0].ub)
                        output_array.append([res_lb[j][0],res_lb[j][1], lower_lb_paired, upper_lb_paired])
                        j+=1 
                        cap = j < len(res_lb)
                        
                elif interv[0].lb < res_lb[j][0].lb:
                    if res_lb[j][0].ub > interv[0].lb:
                        lower_lb_solo = interv[0].inverse_divide_interval(res_lb[j][0].lb)
                        upper_lb_solo = interv[1].inverse_divide_interval(res_lb[j][0].lb)
                        output_array.append([None,None, lower_lb_solo, upper_lb_solo])
                    else:
                        output_array.append([None, None, interv[0], interv[1]])
                        break
                    if interv[0].ub > res_lb[j][0].ub:
                        lower_lb_solo = interv[0].inverse_divide_interval(res_lb[j][0].ub)
                        upper_lb_solo = interv[1].inverse_divide_interval(res_lb[j][0].ub)
                        output_array.append([res_lb[j][0],res_lb[j][1], lower_lb_solo, upper_lb_solo])
                        j+=1
                        cap = j < len(res_lb)
                    else:
                        lower_lb_paired = res_lb[j][0].inverse_divide_interval(interv[0].ub)
                        upper_lb_paired = res_lb[j][1].inverse_divide_interval(interv[0].ub)
                        output_array.append([lower_lb_paired,upper_lb_paired, interv[0], interv[1]])
                        break
                else:
                    if interv[0].ub > res_lb[j][0].ub:
                        lower_lb_solo = interv[0].inverse_divide_interval(res_lb[j][0].ub)
                        upper_lb_solo = interv[1].inverse_divide_interval(res_lb[j][0].ub)
                        output_array.append([res_lb[j][0],res_lb[j][1], lower_lb_solo, upper_lb_solo])
                        j+=1
                        cap = j < len(res_lb)
                    else:
                        lower_lb_paired = res_lb[j][0].inverse_divide_interval(interv[0].ub)
                        upper_lb_paired = res_lb[j][1].inverse_divide_interval(interv[0].ub)
                        output_array.append([lower_lb_paired,upper_lb_paired, interv[0], interv[1]])
                        break
        if j < len(res_lb):
            for k in range(j, len(res_lb)):
                output_array.append([res_lb[j][0],res_lb[j][1],None,None])
        return output_array
    def _check_centers_list(self, interval_list):
        output_array = []
        for interv in interval_list:
            if interv[1] is None or interv[2] is None:
                if interv[1] is None:
                    output_array.append([interv[2], interv[3]])
                    continue
                output_array.append([interv[0], interv[1]])
            else:
                center_checked_interval = self._check_center(interv[1].lb, interv[1].ub, interv[1], interv[2])
                interval_gen = self._center_interval_generator(interv[0].lb, interv[0].ub, interv[0], interv[3], center_checked_interval)
                output_array.append(interval_gen)
                

    def _lb_ub_offset_node(self):
        if self.offset.offset_lb:
            if self.offset.offset_ub:
                """
                  [max(lb,lbd),min(ub,ubd)]
                """
                output_sets = self._min_max_finder()
                res_intervs = [self._check_max_min(mini.lb, mini.ub, mini, maxi) for mini, maxi in output_sets]
            else:
                """
                  [lb*, min(ub*, ubd)] + [max(lb, lbd), ub*]
                """
                res_ub = self._ub_max_finder(self.offset)
                res_lb = self._lb_min_finder(self.offset)
                intervals_merged = self._interval_finder_merger(res_ub, res_lb)
                res_intervs = self._check_centers_list(intervals_merged)

        else:
            """
                [max(lb*,lbd),min(ub*,ubd)]
            """
            output_sets = self._min_max_finder()
            res_intervs = [self._check_max_min(mini.lb, mini.ub, mini, maxi) for mini, maxi in output_sets]

    def _lb_offset_node(self):
        if self.offset.offset_lb:
            if self.offset.offset_ub:
                """
                  ([lb,min(ub,ubd)] + [max(lb,lbd),ub])
                """
                res_ub = self._ub_max_finder(self.interval)
                res_lb = self._lb_min_finder(self.interval)
                intervals_merged = self._interval_finder_merger(res_ub, res_lb)
                res_intervs = self._check_centers_list(intervals_merged)


            else:
                """
                    [-pi, min(ub, ubd)] + [max(lb*, lbd*), pi]
                """
                res_ub = self._ub_max_finder(self.constant_interval)
                res_lb = self._lb_min_finder(self.constant_interval)
                intervals_merged = self._interval_finder_merger(res_ub, res_lb)
                res_intervs = self._check_centers_list(intervals_merged)

               
        else:
            """
                     [lbd, min(ub, ubd)] + [max(lb*, lbd), ubd]
            """
            res_ub = self._ub_max_finder(self.interval)
            res_lb = self._lb_min_finder(self.interval)
            intervals_merged = self._interval_finder_merger(res_ub, res_lb)
            res_intervs = self._check_centers_list(intervals_merged)

    def _no_offset_node(self):
        if self.offset.offset_lb:
            if self.offset.offset_ub:
                """
                  [max(lb,lbd*),min(ub,ubd*)]
                """
                output_sets = self._min_max_finder()
                res_intervs = [self._check_max_min(mini.lb, mini.ub, mini, maxi) for mini, maxi in output_sets]
            else:
                """
                   [lb, min(ub, ubd)] + [max(lb, lbd*), ub]
                """
                res_ub = self._ub_max_finder(self.interval)
                res_lb = self._lb_min_finder(self.interval)
                intervals_merged = self._interval_finder_merger(res_ub, res_lb)
                res_intervs = self._check_centers_list(intervals_merged)

        else:
            """
                [max(lb,lbd),min(ub,ubd)]
            """
            output_sets = self._min_max_finder()
            res_intervs = [self._check_max_min(mini.lb, mini.ub, mini, maxi) for mini, maxi in output_sets]
    def _solve_max_equation(self, low, high, interv_one, interv_two):
        """
        Output would be:

        One, None, None: No breaking points, the one is the maximum
        Two, None, None: No breaking points, the two is the maximum
        One, Two, breaking point: Re-do the intervals with One from low to br
        Two, One, breaking point: Re-do the intervals with One from br to high
        """
        f_one = interv_one.function
        f_two = interv_two.function
        secant = SecantSolution(low, high)

        
        if f_one(low) >= f_two(high) and f_one(high) >= f_two(high):
            return interv_one, None, None
        elif f_one(low) <= f_two(high) and f_one(high) <= f_two(high):
            return interv_two, None, None
        else:
            diff_fn = lambda x: f_one(x)-f_two(x)
            br_point = secant.solve_equation(diff_fn)
            if f_one(low) > f_two(low):
                return interv_one, interv_two, br_point
            else:
                return interv_two, interv_one, br_point
    def _solve_min_equation(self, low, high, interv_one, interv_two):
        """
        Output would be:

        One, None, None: No breaking points, the one is the minimum
        Two, None, None: No breaking points, the two is the minimum
        One, Two, breaking point: Re-do the intervals with One from low to br
        Two, One, breaking point: Re-do the intervals with One from br to high
        """
        f_one = interv_one.function
        f_two = interv_two.function
        secant = SecantSolution(low, high)

        
        if f_one(low) <= f_two(high) and f_one(high) <= f_two(high):
            return interv_one, None, None
        elif f_one(low) >= f_two(high) and f_one(high) >= f_two(high):
            return interv_two, None, None
        else:
            diff_fn = lambda x: f_one(x)-f_two(x)
            br_point = secant.solve_equation(diff_fn)
            if f_one(low) <= f_two(low):
                return interv_one, interv_two, br_point
            else:
                return interv_two, interv_one, br_point

    def _check_max_min(self, low, high, interv_bottom, interv_top):
        """
        Output would be:

        True True None: No breaking points, bottom is always bottom and top is always top
        False False None: No breaking points, bottom is always top and top is always bottom
        False, Interval, breaking point: To create the new intervals from br to high
        Interval, False, breaking point: To create the new intervals from low to br
        """
        f_bottom = interv_bottom.function
        f_top = interv_top.function
        if f_top(low) >= f_bottom(low) and f_top(high) >= f_bottom(high):
            return [interv_bottom, interv_top]
        elif f_top(low) < f_bottom(low) and f_top(high) < f_bottom(high):
            return []
        else:
            diff_fn = lambda x: f_top(x)-f_bottom(x)
            br_point = secant.solve_equation(diff_fn)
            if f_bottom(low) > f_top(low):
                interv_bottom.lb = breaking_point
                interv_top.lb = breaking_point
                return [interv_bottom, interv_top]
            else:
                interv_bottom.ub = breaking_point
                interv_top.ub = breaking_point

                return [interv_bottom, interv_top]

    def _check_center(self, low, high, interv_bottom, interv_top):
        """
            The function is the same, but the output should be understood differently
            In the case of both true, then both stay as the lower and upper limits
            In the case of both false, then the lower and upper limits are the ones from the outer side (The outer side is always true)
            In any other case, we look for the breaking point and then obtain three intervals, one before (after) the breakpoint where
            the false case holds and two after the breakpoint
        """
        
        return self._check_max_min(low, high, interv_bottom, interv_top)

    def _center_interval_generator(self, low, high, interv_bottom, interv_top, interv_middle):
        if interv_middle == []:
            return [interv_bottom, interv_top]
        else:
            if interv_middle[0].lb != interv_bottom.lb:
                centerless_bot = interv_bottom.divide_interval(interv_middle[0].ub)
                centerless_top = interv_top.divide_interval(interv_middle[0].ub)
                return [[centerless_bot, centerless_top], [interv_bottom, interv_middle[0]], [interv_middle[1], interv_top]]

            elif interv_middle[0].ub != interv_bottom.ub:
                centerless_bot = interv_bottom.inverse_divide_interval(interv_middle[0].ub)
                centerless_top = interv_top.inverse_divide_interval(interv_middle[0].ub)
                return [[centerless_bot, centerless_top], [interv_bottom, interv_middle[0]], [interv_middle[1], interv_top]]
            else:  
                return [[interv_bottom, interv_middle[0]], [interv_middle[1], interv_top]]
    







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


