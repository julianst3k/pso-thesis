from aux import cotan, UniformRectangle, EightRectangle, ProbabilityCalculator, IntegrationLimit, Orientation
import numpy as np
import equation_solvers as eq

class AnalyticalProbability(ProbabilityCalculator):
    def __init__(self, X, Y, x_center, y_center, fov, beta, h, r, threshs):
        super().__init__(fov, beta, h, r)
        self.rect = EightRectangle(X, Y, x_center, y_center)
        self.X = X
        self.Y = Y
        self.threshs = threshs
        self.from_one = False
        self.solve_thresholds()


    def solve_thresholds(self):
        thresh_solver = eq.ThresholdSolver(self.threshs)
        self.lims = thresh_solver.solve_equations(self)
        print(self.lims)
    def print_lims(self, lims = None):
        if lims is None:
            lims = self.lims
        print(lims)
        for lim in lims:
            print(lim)

    
    def calculate_probability(self):
        ...
    def calculate_probability_unitary(self, L_max):
        tot_sum = 0
        triangles = self.rect.triangles
        for triangle in triangles:
            Dn = triangle.get_Dn()
            area = triangle.get_area()
            if triangle.orientation == Orientation.HORIZONTAL:
                if L_max < Dn:
                    tot_sum += area*(triangle.ang_crt-triangle.ang_low)*L_max**2/(Dn**2*np.tan(triangle.ang_crt))
                else:
                    new_low = np.arccos(Dn/L_max)
                    if new_low > triangle.ang_crt:
                        triangle.change_ang(triangle.ang_crt)
                        new_low = triangle.ang_crt
                    tot_sum += area*(np.tan(new_low)/np.tan(triangle.ang_crt)+
                        (triangle.ang_crt-new_low)*(L_max**2/(Dn**2*np.tan(triangle.ang_crt))))
            else:
                if L_max < Dn:
                    tot_sum += area*(triangle.ang_high-triangle.ang_crt)*L_max**2/(Dn**2*cotan(triangle.ang_crt))
                else:
                    new_high = np.arcsin(Dn/L_max)
                    if new_high < triangle.ang_crt:
                        triangle.change_ang(triangle.ang_crt)
                        new_high = triangle.ang_crt
                    tot_sum += area*(cotan(new_high)/cotan(triangle.ang_crt)+
                        (new_high-triangle.ang_crt)*(L_max**2/(Dn**2*cotan(triangle.ang_crt))))
        return tot_sum

    def calculate_probability_ring(self):
        ...
if __name__ == "__main__":
    beta = np.pi/180*55
    fov = np.pi/180*45
    r = 0.05
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
    an_prob = AnalyticalProbability(X, Y, x_c, y_c, fov, beta, h, r, threshs)
    an_prob.print_lims()
