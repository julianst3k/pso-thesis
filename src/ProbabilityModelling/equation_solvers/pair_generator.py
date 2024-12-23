from aux import NewtonRaphson
from interval import Interval, OffsetInterval, OutOfUnitaryBound
import numpy as np



class PairGenerator:
    def __init__(self, theta, interval: Interval, offset: OffsetInterval, func, dfunc, parameters):
        self.llow = max(interval.lb, offset.lb)
        self.lhigh = min(interval.ub, offset.ub)
        self.theta = theta
        self.func = func
        self.dfunc = dfunc
        self.parameters = parameters
        self.pair_tree = PairTreeGenerator(self, interval, offset)
    def solve(self):
        return self.pair_tree._build_solution_root()

class SIMOPairGenerator:
    """
    El PG no es lo suficiente general para el caso SIMO (LO hice pensando para el MISO), por lo que tendria que
    escribirlo nuevamente, so prefiero escribir una extension
    """
    def __init__(self, alpha, base: Interval, rotated: Interval, parameters):
        self.llow = max(base.lb, rotated.lb)
        self.lhigh = min(base.ub, rotated.ub)
        self.alpha = alpha
        self.rotated = rotated
        self.base = base
        self.parameters = parameters
    def solve(self):
        if self.rotated.offset_lb:
            self._solve(reverse=True) # Lower
            self._solve() # Higher
        else:
            self._solve()
    def _solve(self, reverse = False):
        parameters = self.parameters
        upper_angle = self.alpha if not reverse else 0
        lower_angle = self.alpha if reverse else 0
        low_int = self.base if not reverse else self.rotated
        high_int = self.base if reverse else self.rotated
        offset = True if reverse else False
        costh = np.cos(self.alpha/2)
        b = 2*parameters.sinbeta*costh*parameters.a if reverse else -2*parameters.sinbeta*costh*parameters.a
        a = parameters.cosfov**2-parameters.sinbeta**2*costh**2
        c = parameters.b**2*parameters.cosfov**2-parameters.a**2
        sol1, sol2 = self._solve_quadratic_base(a,b,c,self.alpha,parameters, self.llow)
        tol = 1e-4
        low_ph = None
        up_ph = None
        sol2 = sol2 if sol2 >= self.llow and sol2 <= self.lhigh else None
        sol1 = sol1 if sol1 >= self.llow and sol1 <= self.lhigh else None
        print(sol1, sol2)
        if sol2 is not None:
            try:
                lower = parameters.eq_base(sol2, lower_angle, -1, is_offset = offset)
                upper = parameters.eq_base(sol2, upper_angle, 1)
                if np.abs(lower-upper) < tol:
                    # Aceptar solucion, corta arriba o abajo?
                    try:
                        lower = parameters.eq_base(sol2+0.001, lower_angle, -1, is_offset = offset)
                        upper = parameters.eq_base(sol2+0.001, upper_angle, 1)
                        if lower < upper:
                            low_int.inverse_divide_interval(sol2)
                            high_int.inverse_divide_interval(sol2)
                        else:
                            # Cortar la parte de arriba, pero mantenemos el intervalo para el caso de sol1
                            low_ph = low_int.divide_interval(sol1)
                            up_ph = high_int.divide_interval(sol1)                            
                    except OutOfUnitaryBound:
                        # sol2+0.01 fuera de la solucion, hay que ver si hay que cortar abajo
                        lower = parameters.eq_base(sol2-0.001, lower_angle, -1, is_offset = offset)
                        upper = parameters.eq_base(sol2-0.001, upper_angle, 1)
                        if lower > upper:
                            low_int.inverse_divide_interval(sol1)
                            high_int.inverse_divide_interval(sol1)
            except OutOfUnitaryBound:
                pass
            
        if sol1 is not None:
            try:
                lower = parameters.eq_base(sol1, lower_angle, -1, is_offset = offset)
                upper = parameters.eq_base(sol1, upper_angle, 1)
                if np.abs(lower-upper) < tol:
                    # Aceptar solucion, corta arriba o abajo?
                    try:
                        lower = parameters.eq_base(sol1+0.001, lower_angle, -1, is_offset = offset)
                        upper = parameters.eq_base(sol1+0.001, upper_angle, 1)
                        if lower < upper:
                            # Cortar la parte de abajo
                            if up_ph is not None:
                                up_ph.inverse_divide_interval(sol1)
                                low_ph.inverse_divide_interval(sol1)
                            else:
                                low_int.inverse_divide_interval(sol1)
                                high_int.inverse_divide_interval(sol1)
                        else:
                            low_int.divide_interval(sol1)
                            high_int.divide_interval(sol1)                            
                    except OutOfUnitaryBound:
                        # sol2+0.01 fuera de la solucion, hay que ver si hay que cortar abajo
                        lower = parameters.eq_base(sol1-0.001, lower_angle, -1, is_offset = offset)
                        upper = parameters.eq_base(sol1-0.001, upper_angle, 1)
                        if lower > upper:
                            low_int.inverse_divide_interval(sol1)
                            high_int.inverse_divide_interval(sol1)
            except OutOfUnitaryBound:
                pass
        if up_ph is None:
            return None
        else:
            return [low_ph,up_ph]
    def _solve_quadratic_base(self, a, b, c, theta, parameters, lmin):
            """
            Quick analysis
            if a > 0 => L1 > L2
            
            """
            epsilon = 0.001

            if b**2-4*a*c < 0:
                return None, None
                    
                
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
            return L1, L2

class PairWrapper:
    def __init__(self, interval, offset):
        self.interval = interval
        self.offset = offset
    def get_offset(self):
        return self.offset
    def get_interval(self):
        return self.interval
    def __getitem__(self, index):
        if index == 0:
            return self.interval
        return self.offset




class PairTreeGenerator:
    def __init__(self, pair_generator, interval, offset):
        self.pg = pair_generator
        self.constant_interval = Interval(False, False, self.pg.llow, self.pg.lhigh, upper_func = lambda x: np.pi, lower_func = lambda x: -np.pi, pi_interval = True)
        self.pair_solver = PairSolver(self.pg)
        self.int_gen = PairIntervalGenerator()
        self.pairs = PairWrapper(interval, offset)
    def _build_solution_root(self):
        """
        We can understand the solution as a tree.
        The root is the conditions of the Lower Bound of the Base Interval
        The Tree would look like this:

                                Base
                    LB Offset           LB No Offset
                UB Offset UB No Offset  UB No Offset -> Since LB < UB, UB will never have offset without LB   
        """
        if self.pairs[0].offset_lb:
            if self.pairs[0].offset_ub: 
                return self._lb_ub_offset_node()
            else:
                return self._lb_offset_node()
        else:
            return self._no_offset_node()

 

    def _lb_ub_offset_node(self):
        if self.pairs[1].offset_lb:
            if self.pairs[1].offset_ub:
                """
                  [max(lb,lbd),min(ub,ubd)]
                """
                output_sets = self.pair_solver._min_max_finder(self.pairs)
                res_intervs = [self.pair_solver._check_max_min(max(mini.lb, maxi.lb), min(mini.ub, maxi.ub), mini, maxi) for mini, maxi in output_sets]
            else:
                """
                  [lb*, min(ub*, ubd)] + [max(lb, lbd), ub*]
                        _lb_min              _ub_max
                """
                res_ub = self.pair_solver._ub_max_finder(self.pairs.get_interval(), self.pairs)
                res_lb = self.pair_solver._lb_min_finder(self.pairs.get_interval(), self.pairs)
                intervals_merged = self.int_gen.interval_finder_merger(res_ub, res_lb)
                res_intervs = self.pair_solver._check_centers_list(intervals_merged)

        else:
            """
                [max(lb*,lbd),min(ub*,ubd)]
            """
            output_sets = self.pair_solver._min_max_finder(self.pairs)
            
            res_intervs = [self.pair_solver._check_max_min(max(mini.lb, maxi.lb), min(mini.ub, maxi.ub), mini, maxi) for mini, maxi in output_sets]
        filtered_output = self.int_gen.pair_filtering(res_intervs)
        return filtered_output
    def _lb_offset_node(self):
        if self.pairs[1].offset_lb:
            if self.pairs[1].offset_ub:
                """
                  ([lbd,min(ub,ubd)] + [max(lb,lbd),ubd])
                """
                res_ub = self.pair_solver._ub_max_finder(self.pairs.get_offset(), self.pairs)
                res_lb = self.pair_solver._lb_min_finder(self.pairs.get_offset(), self.pairs)
                intervals_merged = self.int_gen.interval_finder_merger(res_ub, res_lb)
                res_intervs = self.pair_solver._check_centers_list(intervals_merged)


            else:
                """
                    [-pi, min(ub, ubd)] + [max(lb*, lbd*), pi]
                """
                res_ub = self.pair_solver._ub_max_finder(self.constant_interval, self.pairs)
                res_lb = self.pair_solver._lb_min_finder(self.constant_interval, self.pairs)
                intervals_merged = self.int_gen.interval_finder_merger(res_ub, res_lb)
                res_intervs = self.pair_solver._check_centers_list(intervals_merged)

               
        else:
            """
                     [lbd, min(ub, ubd)] + [max(lb*, lbd), ubd]
            """
            res_ub = self.pair_solver._ub_max_finder(self.pairs.get_offset(), self.pairs)
            res_lb = self.pair_solver._lb_min_finder(self.pairs.get_offset(), self.pairs)
            intervals_merged = self.int_gen.interval_finder_merger(res_ub, res_lb)
            res_intervs = self.pair_solver._check_centers_list(intervals_merged)
        filtered_output = self.int_gen.pair_filtering(res_intervs)
        return filtered_output

    def _no_offset_node(self):
        if self.pairs[1].offset_lb:
            if self.pairs[1].offset_ub:
                """
                  [max(lb,lbd*),min(ub,ubd*)]
                """
                output_sets = self.pair_solver._min_max_finder(self.pairs)
                res_intervs = [self.pair_solver._check_max_min(max(mini.lb, maxi.lb), min(mini.ub, maxi.ub), mini, maxi) for mini, maxi in output_sets]
                filtered_output = self.int_gen.pair_filtering(res_intervs)
            else:
                """
                   [lb, min(ub, ubd)] + [max(lb, lbd*), ub]
                """
                res_ub = self.pair_solver._ub_max_finder(self.pairs.get_interval(), self.pairs)
                res_lb = self.pair_solver._lb_min_finder(self.pairs.get_interval(), self.pairs)
                intervals_merged = self.int_gen.interval_finder_merger(res_ub, res_lb)
                res_intervs = self.pair_solver._check_centers_list(intervals_merged)

        else:
            """
                [max(lb,lbd),min(ub,ubd)]
            """
            output_sets = self.pair_solver._min_max_finder(self.pairs)
            res_intervs = [self.pair_solver._check_max_min(max(mini.lb, maxi.lb), min(mini.ub, maxi.ub), mini, maxi) for mini, maxi in output_sets]
        filtered_output = self.int_gen.pair_filtering(res_intervs)
        return filtered_output
    
class PairSolver:
    def __init__(self, pair_generator):
        self.theta = pair_generator.theta
        self.func = pair_generator.func
        self.dfunc = pair_generator.dfunc
        self.parameters = pair_generator.parameters
        self.llow = pair_generator.llow
        self.lhigh = pair_generator.lhigh
        self.int_gen = PairIntervalGenerator()
    def _solve_max_equation(self, low, high, interv_one, interv_two):
        """
        Output would be:

        One, None, None: No breaking points, the one is the maximum
        Two, None, None: No breaking points, the two is the maximum
        One, Two, breaking point: Re-do the intervals with One from low to br
        Two, One, breaking point: Re-do the intervals with One from br to high
        """
        newton = NewtonRaphson(low, high, self.func, self.dfunc, self.parameters)
        is_lb = True
        s1, s2 = newton.solve(self.theta, interv_one, interv_two, is_lb)
        func_val = self.func(low+0.001, self.theta, interv_one, interv_two, is_lb)
        if s1 and s2:
            """
            Two solutions
            """
            if func_val >= 0:
                return interv_one, s1, s2
            else:
                return interv_two, s1, s2
        else:
            sol = s1 if s1 else s2 
            if not sol:
                if func_val >= 0:
                    return interv_one, None, None
                else:
                    return interv_two, None, None
            else:
                if func_val >= 0:
                    return interv_one, interv_two, sol
                else:
                    return interv_two, interv_one, sol

    def _solve_min_equation(self, low, high, interv_one, interv_two):
        """
        Output would be:

        One, None, None: No breaking points, the one is the minimum
        Two, None, None: No breaking points, the two is the minimum
        One, Two, breaking point: Re-do the intervals with One from low to br
        Two, One, breaking point: Re-do the intervals with One from br to high
        """
        newton = NewtonRaphson(low, high, self.func, self.dfunc, self.parameters)
        is_lb = False
        print("Solving Min Equation")
        s1, s2 = newton.solve(self.theta, interv_one, interv_two, is_lb)
        func_val = self.func(low+0.001, self.theta, interv_one, interv_two, is_lb)
        if s1 and s2:
            """
            Two solutions
            """
            if func_val <= 0:
                return interv_one, s1, s2
            else:
                return interv_two, s1, s2
        else:
            sol = s1 if s1 else s2 
            if not sol:
                if func_val <= 0:
                    return interv_one, None, None
                else:
                    return interv_two, None, None
            else:
                if func_val <= 0:
                    return interv_one, interv_two, sol
                else:
                    return interv_two, interv_one, sol

    def _check_max_min(self, low, high, interv_bottom, interv_top, reverse = False):
        """
        Output would be:

        True True None: No breaking points, bottom is always bottom and top is always top
        False False None: No breaking points, bottom is always top and top is always bottom
        False, Interval, breaking point: To create the new intervals from br to high
        Interval, False, breaking point: To create the new intervals from low to br
        """
        func = lambda x, *args: interv_top.gen_func(self.parameters)(x, interv_top, self.theta, interv_top.offset_lb if reverse else interv_top.offset_ub) - interv_bottom.gen_func(self.parameters)(x, interv_bottom, self.theta, interv_bottom.offset_ub if reverse else interv_bottom.offset_lb)
        dfunc = lambda x, *args: interv_top.gen_dfunc(self.parameters)(x, interv_top, self.theta) - interv_bottom.gen_dfunc(self.parameters)(x, interv_bottom, self.theta)
        newton = NewtonRaphson(low, high, func, dfunc, self.parameters)
        s1, s2 = newton.solve(self.theta, interv_bottom, interv_top)
        func_val = func(low+0.001)
        if s1 and s2:
            """
            Two solutions
            """
            if max_interval_one == base_min:
                new_one = base_min.divide_interval(s1)
                offset_min.ub = s2
                new_one.lb = s2
                offset_min.lb = s1
                higher_set.append(offset_min, new_one)
            if max_interval_one == offset_min:
                new_one = max_interval_one.divide_interval(s1)
                base_min.ub = s2
                new_one.lb = s2
                base_min.lb = s1
                higher_set.append(base_min, new_one)
            if func_val >= 0:
                last_interv_top = interv_top.divide_interval(s1)
                last_interv_bot = interv_top.divide_interval(s1)
                last_interv_top.lb = s2
                last_interv_bot.lb = s2 
                return [interv_bottom, interv_top, last_interv_top, last_interv_bot]
            else:
                interv_top.lb = s1
                interv_bot.lb = s1
                interv_top.ub = s2
                interv_top.ub = s2
                return [interv_bottom, interv_top]
        else:
            sol = s1 if s1 else s2 
            if not sol:
                if func_val >= 0:
                    return [interv_bottom, interv_top]
                else:
                    return []
            else:
                if func_val <= 0:
                    interv_bottom.lb = sol
                    interv_top.lb = sol
                    return [interv_bottom, interv_top]
                else:
                    interv_bottom.ub = sol
                    interv_top.ub = sol
                    return [interv_bottom, interv_top]
    def _min_max_finder(self, pair):
        """
        We need to get the maximum and minimum for the intersection. Min Max finder
        takes both intervals (Self.interval, Self.offset), and push a function (So we don't
        modify the intervals but their children) and we solve the equation
        f(x) = g(x)
        The equation can have from 0 (No intersection, one is over the other all the time) to 2
        (It alternates) 
        """
        base_min = pair[0].push_functional_interval(True)
        offset_min = pair[1].push_functional_interval(True)
        max_interval_one, max_interval_two, breaking_point = self._solve_max_equation(self.llow, self.lhigh, base_min, 
        offset_min)
        lower_set = self._breaking_point_insertion(max_interval_one, max_interval_two, breaking_point, base_min, offset_min)
        base_max = pair[0].push_functional_interval(False)
        offset_max = pair[1].push_functional_interval(False)
            
        min_intervs_one, min_intervs_two, breaking_point = self._solve_min_equation(self.llow, self.lhigh, base_max, 
        offset_max)
        print("Breaking point", breaking_point)
        higher_set = self._breaking_point_insertion(min_intervs_one, min_intervs_two, breaking_point, base_max, offset_max)
        output_sets = self.int_gen.interval_sorter_generation(higher_set, lower_set)
        
        return output_sets  
    def _breaking_point_insertion(self, first_interval, second_interval, breaking_point, base, offset):
        output_set = [first_interval]
        if second_interval is not None:
            try:
                if second_interval.is_interval():
                    first_interval.ub = breaking_point
                    second_interval.lb = breaking_point
                    output_set.append(second_interval) 
            except AttributeError:
                s1 = second_interval
                s2 = breaking_point 
                if first_interval == base:
                    new_one = base.divide_interval(s1)
                    offset.ub = s2
                    new_one.lb = s2
                    offset.lb = s1
                    output_set.append(offset)
                    output_set.append(new_one)
                if first_interval == offset:
                    new_one = first_interval.divide_interval(s1)
                    base.ub = s2
                    new_one.lb = s2
                    base.lb = s1
                    output_set.append(base)
                    output_set.append(new_one) 
        return output_set
    def _lb_min_finder(self, lower_bound, pair):
        """
        Similar to the previous case but we check the output to see if it is above Lower Bound
        """
        base_max = pair[0].push_functional_interval(False)
        offset_max = pair[1].push_functional_interval(False)
        min_intervs_one, min_intervs_two, breaking_point = self._solve_min_equation(self.llow, self.lhigh, base_max, 
        offset_max)
        higher_set = self._breaking_point_insertion(min_intervs_one, min_intervs_two, breaking_point, base_max, offset_max)
        print("Higher Set", higher_set[0])
        unf_res_intervs = [self._check_max_min(max(maxi.lb, lower_bound.lb), min(maxi.ub, lower_bound.ub), lower_bound.push_functional_interval(True), maxi) for maxi in higher_set]
        res_intervs = self.int_gen.pair_filtering(unf_res_intervs)
        return res_intervs
    def _ub_max_finder(self, upper_bound, pair):
        base_min = pair[0].push_functional_interval(True)
        offset_min = pair[1].push_functional_interval(True)
        max_interval_one, max_interval_two, breaking_point = self._solve_max_equation(self.llow, self.lhigh, base_min, offset_min)
        lower_set = self._breaking_point_insertion(max_interval_one, max_interval_two, breaking_point, base_min, offset_min)
        print("Lower Set", lower_set[0])
        unf_res_intervs = [self._check_max_min(max(mini.lb, upper_bound.lb), min(mini.ub, upper_bound.ub), mini, upper_bound.push_functional_interval(False)) for mini in lower_set]
        res_intervs = self.int_gen.pair_filtering(unf_res_intervs)
        return res_intervs
    def _check_centers_list(self, interval_list):
        output_array = []
        for interv in interval_list:
            if interv[1] is None or interv[2] is None:
                if interv[1] is None:
                    output_array.append([interv[2], interv[3]])
                    continue
                output_array.append([interv[0], interv[1]])
            else:
                center_checked_interval = self._check_center(interv[1].lb, interv[1].ub, interv[1], interv[2])
                interval_gen = self.center_interval_generator(interv[0].lb, interv[0].ub, interv[0], interv[3], center_checked_interval)
                output_array.append(interval_gen)
        return output_array
    
    def _check_center(self, low, high, interv_bottom, interv_top):
        """
            The function is the same, but the output should be understood differently
            In the case of both true, then both stay as the lower and upper limits
            In the case of both false, then the lower and upper limits are the ones from the outer side (The outer side is always true)
            In any other case, we look for the breaking point and then obtain three intervals, one before (after) the breakpoint where
            the false case holds and two after the breakpoint
        """
        
        return self._check_max_min(low, high, interv_bottom, interv_top, reverse = True)
class PairIntervalGenerator:
    def interval_sorter_generation(self, first_set, second_set):
        """
        The purpose of this class is to make the pairs that are going to be integrated
        This method sorts the pairs to make them fit with each other to have proper 
        upper and lower radius!

        """
        output_sets = []
        j = 0
        for i, interv in enumerate(first_set):
            if j < len(second_set):
                available_sc = True
            while available_sc:   
                sc_interv = second_set[j]
                if interv.ub < sc_interv.ub:
                    new_interval = sc_interv.inverse_divided_interval(interv.ub)
                    output_sets.append([new_interval, interv])
                    break
                elif interv.ub == sc_interv.ub:
                    output_sets.append([sc_interv, interv])
                    j+=1
                    break
                else:
                    new_interval = interv.inverse_divided_interval(sc_interv.ub)
                    output_sets.append([sc_interv, new_interval])
                    j+=1
                available_sc = j < len(second_set)
        return output_sets
    def center_interval_generator(self, low, high, interv_bottom, interv_top, interv_middle):
        if interv_middle == []:
            return [interv_bottom, interv_top]
        else:
            if interv_middle[0].lb != interv_bottom.lb:
                centerless_bot = interv_bottom.divide_interval(interv_middle[0].ub)
                centerless_top = interv_top.divide_interval(interv_middle[0].ub)
                return [[centerless_bot, centerless_top], [interv_bottom, interv_middle[0]], [interv_middle[1], interv_top]]

            elif interv_middle[0].ub != interv_bottom.ub:
                centerless_bot = interv_bottom.inverse_divide_interval(interv_middle[0].ub)
                centerless_top = interv_top.inverse_divide_interval(interv_middle[0].ub)
                return [[centerless_bot, centerless_top], [interv_bottom, interv_middle[0]], [interv_middle[1], interv_top]]
            else:  
                return [[interv_bottom, interv_middle[0]], [interv_middle[1], interv_top]]
    def pair_filtering(self, intervals):
        out = []
        for interval in intervals:
            if len(interval) == 0:
                continue
            elif len(interval) == 4:
                out.append([interval[0], interval[1]])
                out.append([interval[2], interval[3]])
            else:
                out.append(interval)
        return out  
    def interval_finder_merger(self, res_ub, res_lb):
        starting_index = 0
        output_array = []
        for i, ub_interv in enumerate(res_ub):
            if starting_index >= len(res_lb):
                output_array.append([None, None, ub_interv[0], ub_interv[1]])
            starting_index = self._interval_finder_merger_inner_iterator(res_lb, ub_interv, starting_index, output_array)
        if starting_index < len(res_lb):
            for k in range(starting_index, len(res_lb)):
                output_array.append([res_lb[k][0],res_lb[k][1],None,None])
        return output_array
    def _interval_finder_merger_inner_iterator(self, res_lb, ub_interv, starting_index, output_array):
        j = starting_index
        not_capped = j < len(res_lb)
        while not_capped:
            lb_interv = res_lb[j] 
            if ub_interv[0].lb > lb_interv[0].lb:
                if lb_interv[0].ub > ub_interv[0].lb:
                    lower_lb_solo, upper_lb_solo = self._inverse_divide_pair_wrapper(ub_interv[0].lb, lb_interv)
                    output_array.append([lower_lb_solo,upper_lb_solo, None, None])
                else:
                    output_array.append([lb_interv[0],lb_interv[1], None, None])
                    j+=1
                    cap = j < len(res_lb)
                if ub_interv[0].ub < lb_interv[0].ub:
                    lower_lb_paired, upper_lb_paired = self._inverse_divide_pair_wrapper(ub_interv[0].ub, lb_interv)
                    output_array.append([lower_lb_paired,upper_lb_paired, ub_interv[0], ub_interv[1]])
                    break

                else:
                    lower_lb_paired, upper_lb_paired = self._inverse_divide_pair_wrapper(lb_interv[0].ub, ub_interv)
                    output_array.append([lb_interv[0],lb_interv[1], lower_lb_paired, upper_lb_paired])
                    j+=1 
                    cap = j < len(res_lb)
                    
            elif ub_interv[0].lb < lb_interv[0].lb:
                if lb_interv[0].lb < ub_interv[0].ub:
                    lower_ub_solo, upper_ub_solo = self._inverse_divide_pair_wrapper(lb_interv[0].lb, ub_interv)
                    output_array.append([None,None, lower_ub_solo, upper_ub_solo])
                else:
                    output_array.append([None, None, ub_interv[0], ub_interv[1]])
                    break
                if ub_interv[0].ub > lb_interv[0].ub:
                    lower_lb_paired, upper_lb_paired = self._inverse_divide_pair_wrapper(lb_interv[0].ub, ub_interv)
                    output_array.append([lb_interv[0],lb_interv[1], lower_lb_paired, upper_lb_paired])
                    j+=1
                    cap = j < len(res_lb)
                else:
                    lower_lb_paired, upper_lb_paired = self._inverse_divide_pair_wrapper(ub_interv[0].ub, lb_interv)
                    output_array.append([lower_lb_paired, upper_lb_paired, ub_interv[0], ub_interv[1]])
                    break
            else:
                if ub_interv[0].ub > lb_interv[0].ub:
                    lower_ub_paired, upper_ub_paired = self._inverse_divide_pair_wrapper(lb_interv[0].ub, ub_interv)
                    output_array.append([lb_interv[0],lb_interv[1], lower_ub_paired, upper_ub_paired])
                    j+=1
                    cap = j < len(res_lb)
                else:
                    lower_lb_paired, upper_lb_paired = self._inverse_divide_pair_wrapper(ub_interv[0].ub, lb_interv)
                    output_array.append([lower_lb_paired,upper_lb_paired, ub_interv[0], ub_interv[1]])
                    break
        return j
    def inverse_divide_pair_wrapper(self, pivot, sets):
        lower = sets[0].inverse_divided_interval(pivot)
        upper = sets[1].inverse_divided_interval(pivot)
        return lower, upper
    