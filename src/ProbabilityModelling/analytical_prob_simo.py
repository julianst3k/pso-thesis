from aux import cotan, UniformRectangle, EightRectangle, ProbabilityCalculator, IntegrationLimit, Orientation
from integral_funcs import RectangleIntegrator
from interval import Interval
import numpy as np
import equation_solvers as eq
from analytical_prob_siso import AnalyticalProbability
from montecarlo_prob import MonteCarloIntegrator
import copy

class AnalyticalSIMO(AnalyticalProbability):
    def __init__(self, X, Y, x_center, y_center, fov, beta, h, r, threshs, alpha):
        super().__init__(X, Y, x_center, y_center, fov, beta, h, r, threshs)
        self.alpha = alpha
        self.solve_thresholds()
        self._solve_offset()
    def _solve_offset(self):
        base_solver = eq.ArccosEquationSolver(self.threshs)
        lmin = self.lims[0].low if self.lims[0].const != 1 else self.lims[1].low
        self.simo_intervals = base_solver.solve_base_equations(self, self.alpha, lmin)
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
    def integrate(self):
        self.do_pairings()
        integral = 0
        integrator = RectangleIntegrator(self.rect)
        for pair in self.base_simo_pairs:
            if pair[1].offset_lb and pair[1].consts != 1:
                copy_rect = copy.deepcopy(integrator.rect)
                val_up = integrator.pair_integrator(pair[0].lb, pair[1].ub, pair[1].consts, self, False)
                integral += val_up
                integrator.rect = copy_rect
                val_down = integrator.pair_integrator(pair[0].lb, pair[1].ub, pair[1].consts, self, True)
                integral += val_down
            else:
                integral += integrator.pair_integrator(pair[0].lb, pair[1].ub, pair[1].consts, self, False)
            #print(integral, pair[0])
        return integral

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
    beta_arr = np.linspace(20,80,60)
    fov_arr = np.linspace(20,80,60)
    miso_arr = np.zeros((60,60,2))

    for u, beta in enumerate(beta_arr):
        for v, fov in enumerate(fov_arr):
            an_prob = AnalyticalSIMO(X, Y, x_c, y_c, fov*np.pi/180, beta*np.pi/180, h, r, threshs, alpha)
            integrator = MonteCarloIntegrator()
            mont = integrator.simo_integrator(10000, beta = beta, fov = fov, alpha=alpha)
            miso_arr[u, v, :] = np.array([an_prob.integrate(), mont])
            an_prob = AnalyticalSIMO(X, Y, x_c, y_c, fov*np.pi/180, beta*np.pi/180, h, r, threshs, alpha)
            print(an_prob.integrate(), mont, beta, fov)
        
    miso_arr.tofile("simo_arr.csv", sep=",")

