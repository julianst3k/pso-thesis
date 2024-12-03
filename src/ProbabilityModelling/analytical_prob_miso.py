from aux import cotan, UniformRectangle, EightRectangle, ProbabilityCalculator, IntegrationLimit, Orientation
from interval import Interval, Bound, OffsetInterval, OutOfUnitaryBound
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
                    self._sort_intervals_by_lb(output)
                    output[triangle] = self.divide_by_lims(output[triangle], divider[triangle] if triangle in divider else divider)
            return output
        return _modulus_wrapper    
    def _generate_sol_offset_equations(self, theta = None):
        self.sol_offset_equations = self._solve_offset_wrapper(theta)
        self._sort_intervals_by_lb(self.sol_offset_equations)

    def _generate_sol_base_equations(self, theta = None):
        self.sol_base_equations = self._solve_base_wrapper(theta)
        self._sort_intervals_by_lb(self.sol_base_equations)
    
    def _interval_divide(self, theta = None):
        self._solve_lims_offset(theta)
        self._generate_sol_base_equations(theta)
        self._generate_sol_offset_equations(theta)
        if theta is not None:
            self.sol_base_equations = self.divide_by_lims(self.sol_base_equations, self.lims)
            self.sol_offset_equations = self.divide_by_lims(self.sol_offset_equations, self.offset_lims)

    def _interval_fit(self, theta = None):
        self._interval_divide(theta)
        if theta is not None:
            self.base_offset_pairs = self.interval_fitting(self.sol_base_equations, self.sol_offset_equations)
        else:
            self.base_offset_pairs = {}
            for triangle in self.rect:
                self.base_offset_pairs[triangle] = self.interval_fitting(self.sol_base_equations[triangle], self.sol_offset_equations[triangle])

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
        if np.abs((self.cosfov*np.sqrt(L**2+2*d*L*np.cos(theta)+d**2+self.b**2)-self.a)/(np.sqrt(L**2+d**2+2*d*L*np.cos(theta))*self.sinbeta)) > 1:
            raise OutOfUnitaryBound
        off = np.arctan(L*np.sin(theta)/(L*np.cos(theta)+d))+pivot*np.pi

        return 2*np.pi*is_offset*(-1)**negative_mod + neg*np.arccos((self.cosfov*np.sqrt(L**2+2*d*L*np.cos(theta)+d**2+self.b**2)-self.a)/(np.sqrt(L**2+d**2+2*d*L*np.cos(theta))*self.sinbeta))-off
    def eq_offset_int(self, L, interval, theta, offset):
        sign = (-1)**(interval.is_neg)
        return self.eq_offset(L, theta, sign, interval.pivoted, offset, interval.ub_over_pi)
    def eq_base(self, L, theta, neg=1, pivot = False, is_offset = False):
        return  2*np.pi*is_offset + neg*np.arccos((self.cosfov*np.sqrt(L**2+self.b**2)-self.a)/(L*self.sinbeta))-theta
    def eq_base_int(self, L, interval, theta, offset):
        sign = (-1)**(interval.is_neg)
        return self.eq_base(L, theta, sign, interval.pivoted, offset)

    def _sort_intervals_by_lb(self, list_of_intervals):
        list_of_intervals.sort(key = lambda interval: interval.lb)


    def print_intervals(self, intervals):
        for interv in intervals:
            print(interv)
    
    def divide_by_lims(self, list_of_intervals, list_of_thresholds):
        """
        Input:
        list_of_intervals: The list of intervals generated by Arccos_Equation_Solver
        list_of_thresholds: The list of thresholds generated by Treshold_Equation_Solver
        Result: A new list of intervals
        """
        current_interval = 0
        return_interval = []
        for i, lim in enumerate(list_of_thresholds):
            while current_interval < len(list_of_intervals):
                curr = list_of_intervals[current_interval]
                if curr.lb < lim.low:
                    _ = curr.inverse_divided_interval(lim.low)
                if lim.low <= curr.ub:
                    if lim.high < curr.lb:
                        break
                    curr.set_consts(lim.const)
                    if lim.high < curr.ub and lim.high > curr.lb:
                        added_interval = curr.inverse_divided_interval(lim.high)
                        return_interval.append(added_interval)
                        break                                
                    else:
                        return_interval.append(curr)
                        current_interval += 1
                        if lim.high == curr.ub:
                            break
                else:
                    current_interval += 1
                    continue
        self._sort_intervals_by_lb(return_interval)
        return return_interval

    def interval_fitting(self, first_set, second_set):
        """
        It is assumed that both first and second sets are ordered
        """
        current_interval = 0
        return_interval = []
        for i, inter in enumerate(first_set):
            not_completed = current_interval < len(second_set)
            while not_completed:
                curr = second_set[current_interval]
                if inter.ub < curr.ub:
                    if inter.ub < curr.lb:
                        """
                        In this case, the interval is dismissed because no interval is intersecting
                        """
                        break
                    else:
                        if inter.lb <= curr.lb:
                            _ = inter.inverse_divided_interval(curr.lb) # We discard the lower half
                        new_curr = curr.inverse_divided_interval(inter.ub)
                        return_interval.append([inter, new_curr])

                        break
                else:
                    if inter.lb > curr.ub:
                        current_interval += 1 
                        not_completed = current_interval < len(second_set)
                    else:
                        new_inter = inter.inverse_divided_interval(curr.ub) # We discard the lower half
                        if inter.lb >= curr.lb:
                            _ = curr.inverse_divided_interval(new_inter.lb)
                        return_interval.append([new_inter, curr])
                        current_interval += 1
                        not_completed = current_interval < len(second_set)
        return return_interval
    def _interval_solver(self, interval_pairs, theta):
        output = []
        for interval_pair in interval_pairs:
            print(interval_pair[0], interval_pair[1])
        for interval_pair in interval_pairs:

            base_pair = interval_pair[0]
            offset_pair = interval_pair[1]
            solver = eq.PairGenerator(theta, base_pair, offset_pair, self.interval_diff, self.interval_diff_d, self)
            solutions = solver.solve()
            for inter in solutions:
                print("Lower", inter[0])
                print("Upper", inter[1])
            output.extend(solutions)
    def arg_acos(self, u):
        return (self.cosfov*np.sqrt(u**2+self.b**2)-self.a)/(u*self.sinbeta)

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



    def integral_debug(self):
        self._interval_divide()
        for triang in self.rect:
            print([interv for interv in self.sol_offset_equations[triangle]])

    
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
    an_prob = AnalyticalMISO(X, Y, x_c, y_c, fov, beta, h, r, threshs, d)
    an_prob.integral_debug()


