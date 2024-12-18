from aux import cotan, UniformRectangle, EightRectangle, ProbabilityCalculator, IntegrationLimit, Orientation
from integral_funcs import RectangleIntegrator
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
    def print_lims(self, lims = None):
        if lims is None:
            lims = self.lims
        for lim in lims:
            print(lim)

    
    def calculate_probability(self):
        integrator = RectangleIntegrator(self.rect)
        summ = 0
        for lim in self.lims:
            if lim.const != 1:
                summ += integrator.non_origin_integrator(lim.low, lim.high, lim.const, self)
            else:
                summ += integrator.origin_integrator(lim.high)
        return summ
            

if __name__ == "__main__":
    beta = np.pi/180*45
    fov = np.pi/180*30
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
    an_prob = AnalyticalProbability(X, Y, x_c, y_c, fov, beta, h, r, threshs)
    print(an_prob.calculate_probability())
