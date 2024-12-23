from aux import cotan, UniformRectangle, EightRectangle, ProbabilityCalculator, IntegrationLimit, Orientation
from interval import Interval, Bound, OffsetInterval, OutOfUnitaryBound
import numpy as np
from .interval_generator_solver import IntervalOffsetSolver
class ArccosEquationSolver:
    def __init__(self, thresh):
        self.thresh = thresh

    def solve_offset_equations(self, parameters, theta, Lmin):
        interval_solver = IntervalOffsetSolver()
        output = []
        d = parameters.d
        dsin = parameters.d**2*np.sin(theta)**2
        b1 = np.sqrt(parameters.b**2+dsin)
        costh = np.cos(theta)
        if theta >= np.pi/2 and theta <= 3/2*np.pi:
            pivot = -d/np.cos(theta)
        else:
            pivot = None
        a = parameters.cosfov**2-parameters.sinbeta**2*costh**2
        b = (2*d*parameters.cosfov**2+2*parameters.sinbeta*parameters.a-2*d*parameters.sinbeta**2)*costh
        c = parameters.cosfov**2*d**2*costh**2+parameters.cosfov**2*b1**2-d**2*parameters.sinbeta**2-parameters.a**2+2*parameters.a*d*parameters.sinbeta
        L1ap, L2ap, flagu = None, None, None
        L1bp, L2bp, flagl = self._solve_quadratic_offset(a, b, c, theta,parameters, pivot_point = -d/np.cos(theta), lmin = Lmin)
        
        if pivot is not None:
            b = (2*d*parameters.cosfov**2-2*parameters.sinbeta*parameters.a-2*d*parameters.sinbeta**2)*costh
            c = parameters.cosfov**2*d**2*costh**2+parameters.cosfov**2*b1**2-d**2*parameters.sinbeta**2-parameters.a**2-2*parameters.a*d*parameters.sinbeta
            L1ap, L2ap, flagu = self._solve_quadratic_offset(a, b, c, theta, parameters, True, pivot_point= -d/np.cos(theta), lmin = -d/np.cos(theta))
            if flagu is None:
                pivot = None
        output = interval_solver.offset_intervals_solver(L1bp, L2bp, L1ap, L2ap, pivot, theta, parameters, flagl, flagu)

        return output

    def solve_base_equations(self, parameters, theta, Lmin):
        interval_solver = IntervalOffsetSolver()
        output = []
        costh = np.cos(theta)
        print(costh)
        b = 2*parameters.sinbeta*costh*parameters.a
        a = parameters.cosfov**2-parameters.sinbeta**2*costh**2
        c = parameters.b**2*parameters.cosfov**2-parameters.a**2
        sol1, sol2 = self._solve_quadratic_base(a,b,c,theta,parameters, Lmin)
        output = interval_solver.base_intervals_solver(sol1, sol2, theta, parameters)
        if parameters.cosfov*np.sqrt(parameters.b**2)-parameters.a > 0:
            x_switch = np.sqrt(parameters.cosfov**2*parameters.b**4/parameters.a**2-parameters.b**2)
            for interv in output:
                if interv.lb < x_switch:
                    if interv.ub > x_switch:
                        new_interval = interv.divide_interval(x_switch)
                        new_interval.decreasing = True
                        output.append(new_interval)
                else:
                    interv.decreasing = True
        return output
    def solve_equations_triangle_wrapper(func):
        def triangle_wrapper(*args, **kwargs):
            sol_equations = {}
            for triangle in kwargs.get("triangles"):
                lmin = kwargs.get("lims")[triangle][0].low
                theta = triangle.avg_ang
                sol_equations[triangle] = func(args[0],kwargs.get("parameters"), theta, lmin)
            return sol_equations 
        return triangle_wrapper
    @solve_equations_triangle_wrapper
    def solve_base_equations_triangles(self, parameters, theta, lmin):
        
        return self.solve_base_equations(parameters, theta, lmin)
    @solve_equations_triangle_wrapper
    def solve_offset_equations_triangles(self, parameters, theta, lmin):
            
        return self.solve_offset_equations(parameters, theta, lmin)
    
    def _solve_quadratic_base(self, a, b, c, theta, parameters, lmin):
        """
        Quick analysis
        if a > 0 => L1 > L2
        
        """
        epsilon = 0.001

        if b**2-4*a*c < 0:
            if parameters.cosfov*np.sqrt(parameters.d**2+parameters.b**2)-parameters.a < 0:
                ub = parameters.eq_base(lmin, theta)
                lb = parameters.eq_base(lmin, theta,neg=-1)
                return lb < -np.pi, ub < -np.pi

                
            
        sqrt = np.sqrt(b**2-4*a*c)
        if a>0:
            L1 = (-b+sqrt)/(2*a)
            L2 = (-b-sqrt)/(2*a)
        
        elif a<0:
            L2 = (-b+sqrt)/(2*a)
            L1 = (-b-sqrt)/(2*a)
        
        else:
            L1 = -a/b
            L2 = None
        sols = [L1, L2]
        i = 0
        while i < 2:
            if sols[i] is not None and (np.abs(parameters.arg_acos(sols[i]+0.001)) >= 1 or np.abs(parameters.arg_acos(sols[i]-0.001)) >= 1):   
                if len(sols) > i+1 and sols[i+1] is not None:
                    sols[i] = sols[i+1]
                    sols[i+1] = None
                else:
                    sols[i] = None
                    i += 2
            else:
                i += 1
        L1, L2 = sols[0], sols[1]
        if L1 is None or L1 < 0:
            ub = parameters.eq_base(lmin, theta)
            lb = parameters.eq_base(lmin, theta,neg=-1)
            """
                As we saw before, it can only be below np.pi!
            """
            return lb < -np.pi, ub < -np.pi
        else:
            L2_is_ub = False
            L2_is_lb = False
            ub_L1 = parameters.eq_base(L1, theta)
            lb_L1 = parameters.eq_base(L1, theta,neg=-1)
            L1_is_ub = abs(ub_L1+np.pi) <= epsilon
            L1_is_lb = abs(lb_L1+np.pi) <= epsilon
            if L2 is not None and L2 > 0:
                ub_L2 = parameters.eq_base(L2, theta)
                lb_L2 = parameters.eq_base(L2, theta,neg=-1)
                L2_is_ub = abs(ub_L2+np.pi) <= epsilon
                L2_is_lb = abs(lb_L2+np.pi) <= epsilon
                return SolWrapper(L1, L1_is_lb, L1_is_ub), SolWrapper(L2, L2_is_lb, L2_is_ub)
            return SolWrapper(L1, L1_is_lb, L1_is_ub), None
        


    def _solve_quadratic_offset(self, a, b, c, theta, parameters, pivot = False, pivot_point = None, lmin = None):
        """
            The equation we are solving is of the sort
                    ax^2+bx+c = 0
            However, the solution assumes (x+dcos(theta)) >= 0
            The solutions of the quadratic equation are symmetric around the pivot point (x+dcos(theta))
            However, not all solutions are valid. 
            L1, L2, L3, L4 are the raw solutions, with L1-L2 being below the pivot and L3-L4 being above the pivot
            

        """
        epsilon = 0.001
        if b**2-4*a*c < 0:
            """
                In this case, the logic is the same! Except for a single thing,
                we have no theta because it is a variable depending on L :/
                Solution: Just evaluate this up at zero
                (Guaranteed to be under pivot), and then
                return the truth value

            """
            ub = parameters.eq_offset(lmin if not pivot else pivot_point+0.01, theta, pivot = pivot)
            lb = parameters.eq_offset(lmin if not pivot else pivot_point+0.01, theta,neg=-1, pivot = pivot)
            """
                As we saw before, it can only be below np.pi!
            """
            return lb < -np.pi, ub < -np.pi, 0

        sqrt = np.sqrt(b**2-4*a*c)
        if a>0:
            L1 = (-b+sqrt)/(2*a)
            L2 = (-b-sqrt)/(2*a)
        
        elif a<0:
            L2 = (-b+sqrt)/(2*a)
            L1 = (-b-sqrt)/(2*a)
        else:
            L1 = -a/b
        """
        Checking for fake solutions.
        """
        sol_is_lb = False
        L2_is_ub = False
        L2_is_lb = False
        L1_is_lb = False
        L1_is_ub = False
        print(theta, L1, L2, pivot_point)
        if L1 < lmin if not pivot else pivot_point:
            """
                As we saw before, it can only be below np.pi!
            """
            try:
                ub = parameters.eq_offset(lmin if not pivot else pivot_point+0.01, theta, pivot = pivot)
                lb = parameters.eq_offset(lmin if not pivot else pivot_point+0.01, theta,neg=-1, pivot = pivot)
            except OutOfUnitaryBound:
                return None, None, None
            return lb < -np.pi, ub < -np.pi, 0
        else:
            if np.abs(parameters.eq_offset(L1, theta, pivot = pivot)+np.pi)>epsilon and np.abs(parameters.eq_offset(L1, theta, pivot = pivot)-np.pi)>epsilon:
                if np.abs(parameters.eq_offset(L1, theta, neg=-1, pivot = pivot)+np.pi)>epsilon and np.abs(parameters.eq_offset(L1, theta, neg=-1, pivot = pivot)-np.pi)>epsilon:
                    L1 = None
                else:
                    L1_is_lb = True
            else:
                L1_is_ub = True
        if L2 < lmin if not pivot else pivot_point:
            L2 = None
        else:
            if np.abs(parameters.eq_offset(L2, theta, pivot = pivot)+np.pi)>epsilon and np.abs(parameters.eq_offset(L2, theta, pivot = pivot)-np.pi)>epsilon:
                if np.abs(parameters.eq_offset(L2, theta, neg=-1, pivot = pivot)+np.pi)>epsilon and np.abs(parameters.eq_offset(L2, theta, neg=-1, pivot = pivot)-np.pi)>epsilon:
                    L2 = None 
                else:
                    L2_is_lb = True
            else:
                L2_is_ub = True
            
        
        sol_is_lb = L1_is_lb or L2_is_lb
        sol = L1 if L1 is not None else L2
        if L1 is not None and L2 is not None:
            return SolWrapper(L1, L1_is_lb, L1_is_ub), SolWrapper(L2, L2_is_lb, L2_is_ub), 3
        if sol is not None and sol>=0:

            if sol_is_lb:
                """
                If sol is lb then it is impossible for ub to be outside of -pi,pi
                """

                return sol, False, 1
            else:
                """
                If sol is ub then it can be either ub being above pi or under pi, then lb can be inside or outside
                If lb is outside then ub sol is under, if lb is inside then ub sol is over!
                """
            
                lb = parameters.eq_offset(sol, theta, neg = -1, pivot = pivot)
                return lb < -np.pi, sol, 2
        else:
            if not pivot:
                ub = parameters.eq_offset(lmin, theta, pivot = pivot)
                lb = parameters.eq_offset(lmin, theta,neg=-1, pivot = pivot)
                """
                    As we saw before, it can only be below np.pi!
                """
                return lb < -np.pi, ub < -np.pi, 0
            else:
                ub = parameters.eq_offset(pivot_point+0.01, theta, pivot = pivot)
                lb = parameters.eq_offset(pivot_point+0.01, theta,neg=-1, pivot = pivot)
                """
                    As we saw before, it can only be below np.pi!
                """
                return lb < -np.pi, ub < -np.pi, 0

class SolWrapper:
    def __init__(self, sol, lb, ub):
        self.sol = sol
        self.lb = lb
        self.ub = ub
        self.is_wrapper = True
