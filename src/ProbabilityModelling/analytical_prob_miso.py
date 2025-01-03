from aux import cotan, UniformRectangle, EightRectangle, ProbabilityCalculator, IntegrationLimit, Orientation
from interval import Interval, Bound, OffsetInterval, OutOfUnitaryBound
from integral_funcs import TriangleIntegrator
import numpy as np
from analytical_prob_siso import AnalyticalProbability
import equation_solvers as eq

class AnalyticalMISO(AnalyticalProbability):
    def __init__(self, X, Y, x_center, y_center, fov, beta, h, r, threshs, d):
        super().__init__(X, Y, x_center, y_center, fov, beta, h, r, threshs)
        self.lims_base = self.lims
        self.d = d
        self.rect = UniformRectangle(X, Y, x_center, y_center, 360)
        self._solve_lims_offset()

    def _solve_lims_offset(self, theta=None):
        thresh_solver = eq.ThresholdSolver(self.threshs)
        self.offset_lims = thresh_solver.solve_lims_offset(self, theta)
        
    def _modulus_solver(func):
        def _modulus_wrapper(self, *args, **kwargs):
            output, divider = func(self, *args, **kwargs)
            if divider is not None:
                for triangle in self.rect:
                    self._sort_intervals_by_lb(output[triangle])
                    output[triangle] = self.divide_by_lims(output[triangle], divider[triangle] if triangle in divider else divider)
                    for sol in output[triangle]:
                        if sol.lb == sol.ub:
                            output[triangle].remove(sol)
                    if triangle.avg_ang > 3.13 and triangle.avg_ang < 3.14:
                        for sol in output[triangle]:
                            print(sol)

            return output
        return _modulus_wrapper    
    def _generate_sol_offset_equations(self, theta = None):
        self.sol_offset_equations = self._solve_offset_wrapper(theta)
        if theta is not None:
            self._sort_intervals_by_lb(self.sol_offset_equations)
        else:
            for triangle in self.rect:
                self._sort_intervals_by_lb(self.sol_offset_equations[triangle])

    def _generate_sol_base_equations(self, theta = None):
        self.sol_base_equations = self._solve_base_wrapper(theta)
        if theta is not None:
            self._sort_intervals_by_lb(self.sol_base_equations)
        else:
            for triangle in self.rect:
                self._sort_intervals_by_lb(self.sol_base_equations[triangle])


    
    def _interval_limits_generator(self, theta = None):
        self._solve_lims_offset(theta)
        self._generate_sol_base_equations(theta)
        self._generate_sol_offset_equations(theta)
        if theta is not None:
            self.sol_base_equations = self.divide_by_lims(self.sol_base_equations, self.lims)
            self.sol_offset_equations = self.divide_by_lims(self.sol_offset_equations, self.offset_lims)

    def _offset_base_pair_generator(self, theta = None):
        self._interval_limits_generator(theta)
        debug = False
        if theta is not None:
            self.base_offset_pairs = self.interval_pairing(self.sol_base_equations, self.sol_offset_equations)
        else:
            self.base_offset_pairs = {}
            for triangle in self.rect:
                if triangle.avg_ang > 2.87 and triangle.avg_ang < 2.89:
                    for sol in self.sol_offset_equations[triangle]:
                        print(sol)
                self.base_offset_pairs[triangle] = self.interval_pairing(self.sol_base_equations[triangle], self.sol_offset_equations[triangle], debug = debug)
    @_modulus_solver
    def _solve_base_wrapper(self, theta = None):
        base_solver = eq.ArccosEquationSolver(self.threshs)
        
        if theta is not None:
            lmin = self.lims[0].low
            return base_solver.solve_base_equations(self, theta, lmin), None
        else:
            self._remap_lims_to_triangles()
            return base_solver.solve_base_equations_triangles(triangles=self.rect, parameters=self, lims=self.lims), self.lims
    @_modulus_solver
    def _solve_offset_wrapper(self, theta = None):
        base_solver = eq.ArccosEquationSolver(self.threshs)
        if theta is not None:
            lmin = self.offset_lims[0].low
            return base_solver.solve_offset_equations(self, theta, lmin), None
        else:
            return base_solver.solve_offset_equations_triangles(triangles=self.rect, parameters=self, lims=self.offset_lims), self.offset_lims
    def _remap_lims_to_triangles(self):
        aux_array = [lim for lim in self.lims]
        self.lims = {}
        for triangle in self.rect:
            self.lims[triangle] = aux_array
        
    def eq_offset(self, L, theta, neg=1, pivot = False, is_offset = False, negative_mod = False):
        epsilon = 0.01 # Relevant to avoid floating problems
        acos_arg = (self.cosfov*np.sqrt(L**2+2*d*L*np.cos(theta)+d**2+self.b**2)-self.a)/(np.sqrt(L**2+d**2+2*d*L*np.cos(theta))*self.sinbeta)
        if acos_arg > 1+epsilon:
            raise OutOfUnitaryBound
        off = np.arctan(L*np.sin(theta)/(L*np.cos(theta)+d))+pivot*np.pi
        acos_arg = max(min(1, acos_arg),-1)
        return 2*np.pi*is_offset*(-1)**negative_mod + neg*np.arccos(acos_arg)-off
    def eq_offset_int(self, L, interval, theta, offset):
        sign = (-1)**(interval.is_neg)
        return self.eq_offset(L, theta, sign, interval.pivoted, offset, interval.ub_over_pi)



    def print_intervals(self, intervals):
        for interv in intervals:
            print(interv)
    
    

    
    def _lower_upper_pairs_generator(self, interval_pairs, theta):
        output = []
        for interval_pair in interval_pairs:
            base_pair = interval_pair[0]
            offset_pair = interval_pair[1]
            solver = eq.PairGenerator(theta, base_pair, offset_pair, self.interval_diff, self.interval_diff_d, self)
            solutions = solver.solve()
            output.extend(solutions)
        return output 

    def interval_diff(self, u, theta, base_interval, offset_interval, is_lb = False):
        offset_u = np.sqrt(u**2+self.d**2+2*self.d*u*np.cos(theta))
        arg_acos_base = self.arg_acos(u)
        arg_acos_offset = self.arg_acos(offset_u)
        
        try:
            atan = np.arctan((u*np.sin(theta))/(u*np.cos(theta)+self.d))+offset_interval.pivoted*np.pi
        except RuntimeWarning:
            if np.sin(theta) > 0:
                atan = np.pi/2 
            else:
                atan = 3*np.pi/2
        if is_lb:
            offset = 2*np.pi*offset_interval.offset_lb-np.arccos(arg_acos_offset)-atan 
            base = 2*np.pi*base_interval.offset_lb-np.arccos(arg_acos_base)-theta
        else:
            offset = 2*np.pi*offset_interval.offset_ub+np.arccos(arg_acos_offset)-atan 
            base = 2*np.pi*base_interval.offset_ub+np.arccos(arg_acos_base)-theta
        #if not np.isnan(base-offset):
        #    print(base, offset, base-offset, offset_interval.pivoted, base_interval.pivoted, u)
        return base-offset
    def base_derivative(self, u, theta = None):
        arg_acos_base = self.arg_acos(u)
        return -1/np.sqrt(1-arg_acos_base**2)*(u**2/np.sqrt(u**2+self.b**2)*self.cosfov*self.sinbeta-
        ((self.cosfov*np.sqrt(u**2+self.b**2)-self.a)*self.sinbeta))/(u**2*self.sinbeta**2)
    def base_derivative_int(self, u, interval, theta = None):
        return (-1)**(interval.is_neg)*self.base_derivative(u, theta)
    def offset_derivative_int(self, u, interval, theta):
        offset_u = np.sqrt(u**2+self.d**2+2*self.d*u*np.cos(theta))
        offset_derivative_tan = -1/(offset_u**2)*(self.d*np.sin(theta))
        return (-1)**(interval.is_neg)*self.offset_derivative(u, theta)+2*offset_derivative_tan*(interval.is_neg)
    def offset_derivative(self, u, theta):
        
        offset_u = np.sqrt(u**2+self.d**2+2*self.d*u*np.cos(theta))
        offset_derivative_arccos = self.base_derivative(offset_u)*(u+self.d*np.cos(theta))/offset_u
        offset_derivative_tan = -1/(offset_u**2)*(self.d*np.sin(theta))
        return offset_derivative_arccos - offset_derivative_tan

    def interval_diff_d(self, u, theta, base_interval, offset_interval, is_lb = False):
        offset_u = np.sqrt(u**2+self.d**2+2*self.d*u*np.cos(theta))
        arg_acos_base = self.arg_acos(u)
        arg_acos_offset = self.arg_acos(offset_u)
        base_derivative = self.base_derivative(u)
        offset_derivative_arccos = self.base_derivative(offset_u)*(u+self.d*np.cos(theta))/offset_u
        offset_derivative_tan = 1/(offset_u**2)*(self.d*np.sin(theta))
        if is_lb:
            return -base_derivative+offset_derivative_arccos+offset_derivative_tan
        else:
            return base_derivative-offset_derivative_arccos+offset_derivative_tan

    def integrate(self):
        self._offset_base_pair_generator()
        integrator = TriangleIntegrator(self.rect)
        integral = 0
        pairs_dict = {}
        for triangle in self.rect:
            #for pair in self.base_offset_pairs[triangle]:
            #    print(pair[0])
            #    print(pair[1])
            print(triangle.avg_ang)
            if triangle.avg_ang > 3.13 and triangle.avg_ang < 3.14:
                for pair in self.base_offset_pairs[triangle]:
                    print(pair[0])
                    print(pair[1], pair[1].pivoted)
            pairs = self._lower_upper_pairs_generator(self.base_offset_pairs[triangle], triangle.avg_ang)
            pairs_dict[triangle] = pairs
            #for pair in pairs_dict[triangle]:
            #    print(pair[0])
            #    print(pair[1])
        integral = integrator(pairs_dict, self)
        return integral




    def integral_debug(self):
        self._interval_limits_generator()
        for triang in self.rect:
            print(f"Average angle: {triang.avg_ang}, Top Angle: {triang.ang_high}, Low Angle: {triang.ang_low}, Max_r: {triang.max_r}")
            for interv in self.sol_offset_equations[triang]:
                if interv.lb <= triang.max_r:
                    print(str(interv)+f'Integrate {interv._integrate_debug(triang, self)}')
    
    
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
               {"thr": -0.85, "consts": {"a":-3.2, "b": -0.2}},
               {"thr": -0.6, "consts": {"a":-1.51, "b": 1.3}},
               {"thr": 0.6, "consts": {"a":-1, "b":np.pi/2}},
               {"thr": 0.85, "consts": {"a":-1.51, "b": 1.85}},
               {"thr": 1, "consts": {"a":-3.2, "b": 3.3}}]
    an_prob = AnalyticalMISO(X, Y, x_c, y_c, fov, beta, h, r, threshs, d)
    #print([triang.max_r for triang in an_prob.rect])
    print(an_prob.integrate())


