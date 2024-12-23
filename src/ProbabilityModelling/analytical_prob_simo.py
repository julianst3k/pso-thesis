from aux import cotan, UniformRectangle, EightRectangle, ProbabilityCalculator, IntegrationLimit, Orientation
from integral_funcs import RectangleIntegrator
from interval import Interval
import numpy as np
import equation_solvers as eq
from analytical_prob_siso import AnalyticalProbability

class AnalyticalSIMO(AnalyticalProbability):
    def __init__(self, X, Y, x_center, y_center, fov, beta, h, r, threshs, alpha):
        super().__init__(X, Y, x_center, y_center, fov, beta, h, r, threshs)
        self.alpha = alpha
        self.solve_thresholds()
        self._solve_offset()
    def _solve_offset(self):
        base_solver = eq.ArccosEquationSolver(self.threshs)
        self.simo_intervals = base_solver.solve_base_equations(self, self.alpha, 0)
    def do_pairings(self):
        base_interval = Interval(False, False, 0, np.sqrt(self.X**2+self.Y**2))
        self.sol_base_equations = self.divide_by_lims([base_interval], self.lims)
        self.sol_simo_equations = self.divide_by_lims(self.simo_intervals, self.lims)
        self.base_simo_pairs = self.interval_pairing(self.sol_base_equations, self.sol_simo_equations)
        self.base_simo_pairs = [[pair[0], pair[1]] for pair in self.base_simo_pairs if np.abs(pair[0].lb-pair[0].ub)>1e-4]
        
        for pair in self.base_simo_pairs:
            base, rotated = pair[0], pair[1]
            pair_generator = eq.SIMOPairGenerator(self.alpha, base, rotated, self)
            extra_pair_generated = pair_generator.solve()
            if extra_pair_generated is not None:
                self.base_simo_pairs.append(extra_pair_generated)
        self.base_simo_pairs.sort(key =  lambda x: x[0].lb)
        for pair in self.base_simo_pairs:
            print(pair[0])
            print(pair[1]) 
    
if __name__ == "__main__":
    beta = np.pi/180*45
    fov = np.pi/180*45
    r = 0.05
    h = 1.2
    x_c = 1
    y_c = 1
    X = 5
    Y = 3
    d = 1
    alpha = np.pi/4
    threshs = [{"thr": -1, "consts": 1},
               {"thr": -0.85, "consts": {"a":-3.2, "b": -0.2}},
               {"thr": -0.6, "consts": {"a":-1.51, "b": 1.3}},
               {"thr": 0.6, "consts": {"a":-1, "b":np.pi/2}},
               {"thr": 0.85, "consts": {"a":-1.51, "b": 1.85}},
               {"thr": 1, "consts": {"a":-3.2, "b": 3.3}}]
    an_prob = AnalyticalSIMO(X, Y, x_c, y_c, fov, beta, h, r, threshs, alpha)
    #print([triang.max_r for triang in an_prob.rect])
    an_prob.do_pairings()
