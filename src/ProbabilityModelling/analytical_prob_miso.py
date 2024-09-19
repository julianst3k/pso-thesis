from aux import cotan, UniformRectangle, EightRectangle, ProbabilityCalculator, IntegrationLimit, Orientation, Interval, Bound
import numpy as np
from analytical_prob_siso import AnalyticalProbability

class AnalyticalMISO(AnalyticalProbability):
    def __init__(self, X, Y, x_center, y_center, fov, beta, h, r, threshs):
        super().__init__(X, Y, x_center, y_center, fov, beta, h, r, threshs)
        self.lims_base = self.lims
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
            print(sol1, sol2, triangle.avg_ang)
    def _solve_offset_equations(self):
        self.sol_offset_equations = {}
        for triangle in self.rect:
            theta = triangle.avg_ang
            dsin = self.d**2*np.sin(theta)**2
            b1 = self.b**2+dsin
            costh = np.cos(theta)
            if theta >= np.pi and theta <= 3/2*np.pi:
                pivot = -d/np.cos(theta)
                a = self.cosfov**2-self.sinbeta**2*costh**2
                b = (2*d*self.cosfov**2+2*self.sinbeta*self.a-2*d*self.sinbeta)*costh
                c = self.cosfov**2*d**2*costh**2+self.cosfov**2*b1**2-d**2*self.sinbeta-a**2+2*a*d*self.sinbeta
                L1up, L2up = self._solve_quadratic_offset_pivot(a,b,c,False)
                L1op, L2op = self._solve_quadratic_offset_pivot(a,b,c,True)

                

            
    def _solve_quadratic_offset(self, a, b, c, phi)        

    def _solve_quadratic_offset_pivot(self, a, b, c, over_pivot)

    def _solve_quadratic_base(self, a, b, c, theta):
        """
        Quick analysis
        if a > 0 => L1 > L2
        """
        if b**2-4*a*c < 0:
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
        if L1 == False and L2 == False:
            interv = [Interval(False, False, 0, np.sqrt(self.X**2+self.Y**2))]
        if L1 == True and L2 == True:
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
            if np.abs(self._eq_base(L1, theta)+np.pi) < epsilon:
                interv = [Interval(True, True, 0, L1), Interval(True, False, L1, L2), Interval(True, True, L2, np.sqrt(self.X**2+self.Y**2))]
            elif np.abs(self._eq_base(L1, theta, -1)+np.pi) < epsilon:
                interv = [Interval(False, False, 0, L1), Interval(True, False, L1, L2), Interval(False, False, L2, np.sqrt(self.X**2+self.Y**2))]

    def _offset_intervals_solver(self, L1, L2, theta):
        ...
        


    return interv
    def _eq_base(self, L, theta, neg=1):
        return neg*np.arccos((self.cosfov*np.sqrt(L**2+self.b**2)-self.a)/(L*self.sinbeta))-theta


        
    
if __name__ == "__main__":
    beta = np.pi/180*45
    fov = np.pi/180*25
    r = 0.8
    h = 1.2
    x_c = 1
    y_c = 1
    X = 5
    Y = 3
    threshs = [{"thr": -1, "consts": 1},
               {"thr": -0.9, "consts": {"a":-3.2, "b": -0.2}},
               {"thr": -0.6, "consts": {"a":-1.51, "b": 1.3}},
               {"thr": 0.6, "consts": {"a":-1, "b":np.pi/2}},
               {"thr": 0.9, "consts": {"a":-1.51, "b": 1.85}},
               {"thr": 1, "consts": {"a":-3.2, "b": -3.3}}]
    an_prob = AnalyticalMISO(X, Y, x_c, y_c, fov, beta, h, r, threshs)
    an_prob._solve_base_equations()
    