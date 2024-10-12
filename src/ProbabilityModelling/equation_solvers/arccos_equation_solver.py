from aux import cotan, UniformRectangle, EightRectangle, ProbabilityCalculator, IntegrationLimit, Orientation, Interval, Bound, OffsetInterval
import numpy as np
from .interval_generator_solver import IntervalOffsetSolver
class ArccosEquationSolver:
    def __init__(self, thresh):
        self.thresh = thresh

    def solve_offset_equations(self, parameters, theta):
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
        L1bp, L2bp, flagl = self._solve_quadratic_offset(a,b,c,theta,parameters)
        if pivot is not None:
            L1ap, L2ap, flagu = self._solve_quadratic_offset(a,b,c,theta,parameters,True)
        output = interval_solver.offset_intervals_solver(L1bp, L2bp, L1ap, L2ap, pivot, theta, parameters, flagl, flagu)

        return output

    def solve_base_equations(self, parameters, theta):
        interval_solver = IntervalOffsetSolver()
        output = []
        costh = np.cos(theta)
        b = 2*parameters.sinbeta*costh*parameters.a
        a = parameters.cosfov**2-parameters.sinbeta**2*costh**2
        c = parameters.b**2*parameters.cosfov**2-parameters.a**2
        sol1, sol2 = self._solve_quadratic_base(a,b,c,theta,parameters)
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

    def solve_base_equations_triangles(self, triangles, parameters):
        sol_base_equations = {}
        for triangle in triangles:
            theta = triangle.avg_ang
            sol_base_equations[triangle] = self.solve_base_equations(parameters, theta)
        return sol_base_equations 

    def solve_offset_equations_triangles(self, triangles, parameters):
        sol_offset_equations = {}
        for triangle in triangles:
            theta = triangle.avg_ang
            sol_offset_equations[triangle] = self.solve_offset_equations(parameters, theta)
            
        return sol_offset_equations
    
    def _solve_quadratic_base(self, a, b, c, theta, parameters):
        """
        Quick analysis
        if a > 0 => L1 > L2
        
        """
        if b**2-4*a*c < 0:
            if parameters.cosfov*np.sqrt(parameters.d**2+parameters.b**2)-parameters.a < 0:
                if theta > np.pi:
                    return False, True
                else:
                    return False, False
            else:
                if theta > np.pi:
                    return True, True
                else:
                    return False, False
                
            
        sqrt = np.sqrt(b**2-4*a*c)
        if a>0:
            L1 = (-b+sqrt)/(2*a)
            L2 = (-b-sqrt)/(2*a)
        
        elif a<0:
            L2 = (-b+sqrt)/(2*a)
            L1 = (-b-sqrt)/(2*a)
        
        else:
            return -a/b, None

        if L1 < 0:
            if parameters.cosfov*np.sqrt(parameters.b**2)-parameters.a < 0:
                if theta > np.pi:
                    return False, True
                else:
                    return False, False
            else:
                if theta > np.pi:
                    return True, True
                else:
                    return False, False
        else:
            if L2 < 0:
                
                return L1, None
            return L1, L2
    def _solve_quadratic_offset(self, a, b, c, theta, parameters, pivot = False):
        """
        In this case, elevating to a square can create *fake* solutions, as the intersection
        between a quadratic and an absolute function have at most two solutions,
        but both reside in opposite halfs 
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
            if parameters.cosfov*np.sqrt(parameters.d**2+parameters.b**2)-parameters.a < 0:
                ub = parameters.eq_offset(0, theta, pivot = pivot)
                lb = parameters.eq_offset(0, theta,neg=-1, pivot = pivot)
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
        if np.abs(parameters.eq_offset(L2, theta, pivot = pivot)+np.pi)>epsilon and np.abs(parameters.eq_offset(L2, theta, pivot = pivot)-np.pi)>epsilon:
            if np.abs(parameters.eq_offset(L2, theta, neg=-1, pivot = pivot)+np.pi)>epsilon and np.abs(parameters.eq_offset(L2, theta, neg=-1, pivot = pivot)-np.pi)>epsilon:
                L2 = None 
            else:
                sol_is_lb = True
            
        if np.abs(parameters.eq_offset(L1, theta, pivot = pivot)+np.pi)>epsilon and np.abs(parameters.eq_offset(L1, theta, pivot = pivot)-np.pi)>epsilon:
            if np.abs(parameters.eq_offset(L1, theta, neg=-1, pivot = pivot)+np.pi)>epsilon and np.abs(parameters.eq_offset(L1, theta, neg=-1, pivot = pivot)-np.pi)>epsilon:
                L1 = None
            else:
                sol_is_lb = True
        sol = L1 if L1 is not None else L2
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
                ub = parameters.eq_offset(0, theta, pivot = pivot)
                lb = parameters.eq_offset(0, theta,neg=-1, pivot = pivot)
                """
                    As we saw before, it can only be below np.pi!
                """
                return lb < -np.pi, ub < -np.pi, 0
            else:
                ub = parameters.eq_offset(np.sqrt(parameters.X**2+parameters.Y**2), theta, pivot = pivot)
                lb = parameters.eq_offset(np.sqrt(parameters.X**2+parameters.Y**2), theta,neg=-1, pivot = pivot)
                """
                    As we saw before, it can only be below np.pi!
                """
                return lb < -np.pi, ub < -np.pi, 0