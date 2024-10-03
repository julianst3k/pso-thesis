from aux import cotan, UniformRectangle, EightRectangle, ProbabilityCalculator, IntegrationLimit, Orientation, Interval, Bound, OffsetInterval
import numpy as np
from analytical_prob_siso import AnalyticalProbability

class AnalyticalMISO(AnalyticalProbability):
    def __init__(self, X, Y, x_center, y_center, fov, beta, h, r, threshs, d):
        super().__init__(X, Y, x_center, y_center, fov, beta, h, r, threshs)
        self.lims_base = self.lims
        self.d = d
        self.rect = UniformRectangle(X, Y, x_center, y_center, 360)

    def _solve_base_equations(self):
        self.sol_base_equations = {}
        for triangle in self.rect:
            costh = np.cos(triangle.avg_ang)
            b = 2*self.sinbeta*costh*self.a
            a = self.cosfov**2-self.sinbeta**2*costh**2
            c = self.b**2*self.cosfov**2-self.a**2
            sol1, sol2 = self._solve_quadratic_base(a,b,c,triangle.avg_ang)
            self.sol_base_equations[triangle] = self._base_intervals_solver(sol1, sol2, triangle.avg_ang)
            
            if self.cosfov*np.sqrt(self.b**2)-self.a > 0:
                x_switch = np.sqrt(self.cosfov**2*self.b**4/self.a**2-self.b**2)
                for interv in self.sol_base_equations[triangle]:
                    if interv.lb < x_switch:
                        if interv.ub > x_switch:
                            new_interval = interv.divide_interval(x_switch)
                            new_interval.decreasing = True
                            self.sol_base_equations[triangle].append(new_interval)
                        else:
                            interv.decreasing = True
            print(sol1, sol2, triangle.avg_ang)
    def _solve_offset_equations(self):
        self.sol_offset_equations = {}
        for triangle in self.rect:
            theta = triangle.avg_ang
            d = self.d
            dsin = self.d**2*np.sin(theta)**2
            b1 = np.sqrt(self.b**2+dsin)
            costh = np.cos(theta)
            if theta >= np.pi/2 and theta <= 3/2*np.pi:
                pivot = -d/np.cos(theta)
            else:
                pivot = None
            a = self.cosfov**2-self.sinbeta**2*costh**2
            b = (2*d*self.cosfov**2+2*self.sinbeta*self.a-2*d*self.sinbeta**2)*costh
            c = self.cosfov**2*d**2*costh**2+self.cosfov**2*b1**2-d**2*self.sinbeta**2-self.a**2+2*self.a*d*self.sinbeta
            L1ap, L2ap, flagu = None, None, None
            L1bp, L2bp, flagl = self._solve_quadratic_offset(a,b,c,theta)
            if pivot is not None:
                L1ap, L2ap, flagu = self._solve_quadratic_offset(a,b,c,theta,True)
            self.sol_offset_equations[triangle] = self._offset_intervals_solver(L1bp, L2bp, L1ap, L2ap, pivot, theta, flagl, flagu)
            if self.cosfov*np.sqrt(self.b**2)-self.a > 0:
                try:
                    x_switch = np.sqrt(self.cosfov**2*self.b**4/self.a**2-self.b1**2)-self.d*costh
                    if x_switch > pivot:
                        x_switch = -np.sqrt(self.cosfov**2*self.b**4/self.a**2-self.b1**2)+self.d*costh
                except:
                    pass
                for interv in self.sol_base_equations[triangle]:
                    if interv.lb < x_switch:
                        if interv.ub > x_switch:
                            new_interval = interv.divide_interval(x_switch)
                            new_interval.decreasing = True
                            self.sol_base_equations[triangle].append(new_interval)
                        else:
                            interv.decreasing = True
                

            

    def _solve_quadratic_offset(self, a, b, c, theta, pivot = False):
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
            if self.cosfov*np.sqrt(self.d**2+self.b**2)-self.a < 0:
                ub = self._eq_offset(0, theta, pivot = pivot)
                lb = self._eq_offset(0, theta,neg=-1, pivot = pivot)
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
        if np.abs(self._eq_offset(L2, theta, pivot = pivot)+np.pi)>epsilon and np.abs(self._eq_offset(L2, theta, pivot = pivot)-np.pi)>epsilon:
            if np.abs(self._eq_offset(L2, theta, neg=-1, pivot = pivot)+np.pi)>epsilon and np.abs(self._eq_offset(L2, theta, neg=-1, pivot = pivot)-np.pi)>epsilon:
                L2 = None 
            else:
                sol_is_lb = True
            
        if np.abs(self._eq_offset(L1, theta, pivot = pivot)+np.pi)>epsilon and np.abs(self._eq_offset(L1, theta, pivot = pivot)-np.pi)>epsilon:
            if np.abs(self._eq_offset(L1, theta, neg=-1, pivot = pivot)+np.pi)>epsilon and np.abs(self._eq_offset(L1, theta, neg=-1, pivot = pivot)-np.pi)>epsilon:
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
            
                lb = self._eq_offset(sol, theta, neg = -1, pivot = pivot)
                return lb < -np.pi, sol, 2
        else:
            if not pivot:
                ub = self._eq_offset(0, theta, pivot = pivot)
                lb = self._eq_offset(0, theta,neg=-1, pivot = pivot)
                """
                    As we saw before, it can only be below np.pi!
                """
                return lb < -np.pi, ub < -np.pi, 0
            else:
                ub = self._eq_offset(np.sqrt(self.X**2+self.Y**2), theta, pivot = pivot)
                lb = self._eq_offset(np.sqrt(self.X**2+self.Y**2), theta,neg=-1, pivot = pivot)
                """
                    As we saw before, it can only be below np.pi!
                """
                return lb < -np.pi, ub < -np.pi, 0
            
        
    def _eq_offset(self, L, theta, neg=1, pivot = False):
        off = np.arctan(L*np.sin(theta)/(L*np.cos(theta)+d))+pivot*np.pi

        return neg*np.arccos((self.cosfov*np.sqrt(L**2+2*d*L*np.cos(theta)+d**2+self.b**2)-self.a)/(np.sqrt(L**2+d**2+2*d*L*np.cos(theta))*self.sinbeta))-off
    def _solve_quadratic_base(self, a, b, c, theta):
        """
        Quick analysis
        if a > 0 => L1 > L2
        
        """
        if b**2-4*a*c < 0:
            if self.cosfov*np.sqrt(self.d**2+self.b**2)-self.a < 0:
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
            if self.cosfov*np.sqrt(self.b**2)-self.a < 0:
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
    def _base_intervals_solver(self, L1, L2, theta):
        epsilon = 0.001
        if not L1:
            interv = [Interval(L1, L2, 0, np.sqrt(self.X**2+self.Y**2))]
        if L1 and L2:
            interv = [Interval(True, True, 0, np.sqrt(self.X**2+self.Y**2))]
        if L2 == None:
            if self.cosfov*np.sqrt(self.b**2)-self.a < 0:
                """
                In this case it only crosses once
                Positive + -np.pi => False to True, and lb would be True True
                Negative + -np.pi => True to False, and ub would be False False
                It should never cross np.pi because the maximum is np.pi, 
                """
                if np.abs(self._eq_base(L1, theta)+np.pi) < epsilon:
                    interv = [Interval(True, False, 0, L1), Interval(True, True, L1, np.sqrt(self.X**2+self.Y**2))]
                elif np.abs(self._eq_base(L1, theta, -1)+np.pi) < epsilon:
                    interv = [Interval(True, False, 0, L1), Interval(False, False, L1, np.sqrt(self.X**2+self.Y**2))]
            else:
                """
                In this case if it cross only once in the positives it would be 
                Positive + -np.pi => True to False, and lb would be True True 
                Negative + -np.pi => False to True, and ub would be False False
                It should never cross np.pi because the maximum is np.pi, 
                """
                if np.abs(self._eq_base(L1, theta)+np.pi) < epsilon:
                    interv = [Interval(True, True, 0, L1), Interval(True, False, L1, np.sqrt(self.X**2+self.Y**2))]
                elif np.abs(self._eq_base(L1, theta, -1)+np.pi) < epsilon:
                    interv = [Interval(False, False, 0, L1), Interval(True, False, L1, np.sqrt(self.X**2+self.Y**2))]

        elif L1<np.sqrt(self.X**2+self.Y**2):
            """
                The only case it crosses twice is in the case of self.cosfov*np.sqrt(self.b**2)-self.a > 0
                In this case it would be
                Positive + -np.pi => True to False, False to True, and True True True
                Negative + -np.pi => False to True, True to False, and False, False, False
            """
            if np.abs(self._eq_base(L2, theta)+np.pi) < epsilon:
                interv = [Interval(True, True, 0, L2), Interval(True, False, L2, L1), Interval(True, True, L1, np.sqrt(self.X**2+self.Y**2))]
            elif np.abs(self._eq_base(L2, theta, -1)+np.pi) < epsilon:
                interv = [Interval(False, False, 0, L2), Interval(True, False, L2, L1), Interval(False, False, L1, np.sqrt(self.X**2+self.Y**2))]
        elif L2<np.sqrt(self.X**2+self.Y**2):
            if np.abs(self._eq_base(L2, theta)+np.pi) < epsilon:
                interv = [Interval(True, True, 0, L2), Interval(True, False, L2, L1)]
            elif np.abs(self._eq_base(L2, theta, -1)+np.pi) < epsilon:
                interv = [Interval(False, False, 0, L2), Interval(True, False, L2, L1)]



    def _offset_intervals_solver(self, L1l, L2l, L1u, L2u, pivot, theta, flagl = None, flagu = None):
        if pivot is None:
            """
            In this case there is no transition so we just need to generate the intervals
            This case is equivalent to the non offset solver! But since we coded the results differently, we need to
            do another function
            """
            interv = self._offset_intervals_no_pivot_solver(L1, L2, theta)
        else:
            """
            In this case, there is a transition between pivot and not pivot that we need to take into account
            To do this, we use the following scheme, but 
            Pre-pivot -> The function acts in reverse (Increasing over distance).
            """
            """
            There is a transition so we need to get the interval AND the transition
            """
            interv_sol_one = self._offset_intervals_pivot_solver_lower(L1l, L2l, pivot, theta, flagl)
            interv_sol_two = self._offset_intervals_pivot_solver_upper(L1u, L2u, pivot, theta, flagu)
            interv_sol_one.extend(interv_sol_two)
            interv = interv_sol_one
        return interv
    def _offset_intervals_pivot_solver_lower(self, L1, L2, pivot, theta, flag):
        """
                Since True and False evaluate to numerals, we need a flag to know anything
                If flag == 0 => Both are boolean
                If flag == 1 => The second is boolean
                If flag == 2 => The first is boolean
        """
        epsilon = 0.001
        if flag == 0:
            
            interv = [OffsetInterval(L1, L2, 0, pivot)]
        
        if flag == 1:

            interv = [OffsetInterval(False, L2, 0, L1), OffsetInterval(True, L2, L1, pivot)]

        if flag == 2:

            interv = [OffsetInterval(L1, False, 0, L2), OffsetInterval(L1, True, L2, pivot)]
        
        return interv
    def _offset_intervals_pivot_solver_upper(self, L1, L2, pivot, theta, flag):
        """
                Since True and False evaluate to numerals, we need a flag to know anything
                If flag == 0 => Both are boolean
                If flag == 1 => The second is boolean
                If flag == 2 => The first is boolean
        """
        maxr = np.sqrt(self.X**2+self.Y**2)
        epsilon = 0.001
        if flag == 0:
            
            interv = [OffsetInterval(L1, L2, 0, pivot)]
        
        if flag == 1:

            interv = [OffsetInterval(True, L2, pivot, L1), OffsetInterval(False, L2, L1, maxr)]

        if flag == 2:

            interv = [OffsetInterval(L1, True, pivot, L2), OffsetInterval(L1, False, L2, maxr)]
        
        return interv
    
    
    def _offset_intervals_no_pivot_solver(self, L1, L2, theta):
        maxr = np.sqrt(self.X**2+self.Y**2)
        epsilon = 0.001
        if flag == 0:
            
            interv = [OffsetInterval(L1, L2, 0, maxr)]
        
        if flag == 1:

            interv = [OffsetInterval(True, L2, 0, L1), OffsetInterval(False, L2, L1, maxr)]

        if flag == 2:

            interv = [OffsetInterval(L1, True, 0, L2), OffsetInterval(L1, False, L2, maxr)]
        
        return interv
        
    def _eq_base(self, L, theta, neg=1):
        return neg*np.arccos((self.cosfov*np.sqrt(L**2+self.b**2)-self.a)/(L*self.sinbeta))-theta


        
    
if __name__ == "__main__":
    beta = np.pi/180*45
    fov = np.pi/180*60
    r = 0.05
    h = 1.2
    x_c = 1
    y_c = 1
    X = 5
    Y = 3
    d = 1
    threshs = [{"thr": -1, "consts": 1},
               {"thr": -0.9, "consts": {"a":-3.2, "b": -0.2}},
               {"thr": -0.6, "consts": {"a":-1.51, "b": 1.3}},
               {"thr": 0.6, "consts": {"a":-1, "b":np.pi/2}},
               {"thr": 0.9, "consts": {"a":-1.51, "b": 1.85}},
               {"thr": 1, "consts": {"a":-3.2, "b": -3.3}}]
    an_prob = AnalyticalMISO(X, Y, x_c, y_c, fov, beta, h, r, threshs, d)
    
    an_prob._solve_offset_equations()
    
