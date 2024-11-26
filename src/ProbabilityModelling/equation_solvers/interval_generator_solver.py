from aux import cotan, UniformRectangle, EightRectangle, ProbabilityCalculator, IntegrationLimit, Orientation
from interval import Interval, Bound, OffsetInterval, OutOfUnitaryBound
import numpy as np

class IntervalOffsetSolver:
    """
    The responsibility of this class is managing intervals from the offset solving.

    """
    def __init__(self):
        ...
    def base_intervals_solver(self, L1, L2, theta, parameters):
        maxr = np.sqrt(parameters.X**2+parameters.Y**2)
        try:
            if L1.is_wrapper:
                if L2 is not None and (L2.ub or L2.lb):
                    if (L1.ub or L1.lb):
                        """
                        Remember, L2 < L1
                        """

                        if L2.ub and L1.ub:
                            interv = [Interval(True, True, 0, L2.sol), Interval(True, False, L2.sol, L1.sol), Interval(True, True, L1.sol, maxr)]

                        elif L2.lb and L1.lb:
                            interv = [Interval(False, False, 0, L2.sol), Interval(True, False, L2.sol, L1.sol), Interval(False, False, L1.sol, maxr)]
                        elif L2.lb and L1.ub:
                            bool_off_lb, val = self._get_base_bool(L2.sol, theta, parameters, -1)
                            bool_off_ub, val = self._get_base_bool(L1.sol, theta, parameters, 1)
                            interv = [Interval(not bool_off_lb, not bool_off_ub, 0, L2.sol), Interval(bool_off_lb, not bool_off_ub, L2.sol, L1.sol), Interval(bool_off_lb, bool_off_ub, L1.sol, maxr)]

                        else:
                            bool_off_lb, val = self._get_base_bool(L1.sol, theta, parameters, -1)
                            bool_off_ub, val = self._get_base_bool(L2.sol, theta, parameters, 1)
                            interv = [Interval(not bool_off_lb, not bool_off_ub, 0, L2.sol), Interval(not bool_off_lb, bool_off_ub, L2.sol, L1.sol), Interval(bool_off_lb, bool_off_ub, L1.sol, maxr)]


                    else:

                        if L2.ub:
                            interv = [Interval(True, False, 0, L2.sol), Interval(True, True, L2.sol, maxr)]

                        else:
                            interv = [Interval(False, False, 0, L2.sol), Interval(True, False, L2.sol, maxr)]
                else:
                    if (L1.ub or L1.lb):
                        if L1.ub:
                            interv = [Interval(True, False, 0, L1.sol), Interval(True, True, L1.sol, maxr)]

                        else:
                            interv = [Interval(False, False, 0, L1.sol), Interval(True, False, L1.sol, maxr)]

                    else:
                        bool_off_lb, val = self._get_base_bool(L1.sol, theta, parameters, -1)
                        bool_off_ub, val = self._get_base_bool(L1.sol, theta, parameters, 1)

                        interv = [Interval(bool_off_lb, bool_off_ub, 0, maxr)]

        except AttributeError:
            interv = [Interval(L1, L2, 0, maxr)]
        return interv

    
    def offset_intervals_solver(self, L1l, L2l, L1u, L2u, pivot, theta, parameters, flagl = None, flagu = None):
        if pivot is None:
            """
            In this case there is no transition so we just need to generate the intervals
            This case is equivalent to the non offset solver! But since we coded the results differently, we need to
            do another function
            """
            interv = self._offset_generic_solver(L1l, L2l, theta, parameters, flagl)
        else:
            """
            In this case, there is a transition between pivot and not pivot that we need to take into account
            To do this, we use the following scheme, but 
            Pre-pivot -> The function acts in reverse (Increasing over distance).
            """
            """
            There is a transition so we need to get the interval AND the transition
            """
            interv_sol_one = self._offset_generic_solver(L1l, L2l, theta, parameters, flagl, pivot = pivot)
            interv_sol_two = self._offset_generic_solver(L1u, L2u, theta, parameters, flagu, True, pivot = pivot)
            interv_sol_one.extend(interv_sol_two)
            interv = interv_sol_one
        return interv

    
    def _offset_generic_solver(self, L1, L2, theta, parameters, flag, pivoted = False, pivot = None):
        """
                Since True and False evaluate to numerals, we need a flag to know anything
                If flag == 0 => Both are boolean
                If flag == 1 => The second is boolean
                If flag == 2 => The first is boolean
                If flag == 3 => Neither are boolean
                In the case of flag 3, both L1 and L2 will be dictionaries that will tell us the direction
                and the starting point
        """
        maxr = np.sqrt(parameters.X**2+parameters.Y**2) if pivot is None or pivoted else pivot
        floor = pivot if pivoted and pivot > 0 else 0
        if flag == 0:
            
            interv = [OffsetInterval(L1, L2, floor, maxr, pivoted=pivoted)]
        
        if flag == 1:

            bool_off, val = self._get_offset_bool(L1, theta, parameters, -1, pivoted = pivoted)
            interv = [OffsetInterval(not bool_off, L2, floor, L1, pivoted=pivoted, over_pi = val), OffsetInterval(bool_off, L2, L1, maxr, pivoted=pivoted, over_pi = val)]

        if flag == 2:
                
            bool_off, val = self._get_offset_bool(L2, theta, parameters, 1, pivoted = pivoted)

            interv = [OffsetInterval(L1, not bool_off, floor, L2, pivoted=pivoted, over_pi = val), OffsetInterval(L1, bool_off, L2, maxr, pivoted=pivoted, over_pi = val)]
            
        if flag == 3:
            if L1.lb and L2.lb:
                interv = [OffsetInterval(False, False, floor, L2.sol, pivoted=pivoted, over_pi = val), OffsetInterval(True, False, L2.sol, L1.sol, pivoted=pivoted), OffsetInterval(False, False, L1.sol, maxr, pivoted=pivoted, over_pi = val)]

            elif L1.ub and L2.ub:
                interv = [OffsetInterval(True, True, floor, L2.sol, pivoted=pivoted, over_pi = val), OffsetInterval(True, False, L2.sol, L1.sol, pivoted=pivoted, over_pi = val), OffsetInterval(True, True, L1.sol, maxr, pivoted=pivoted, over_pi = val)]
        
            elif L1.ub:
                bool_off_lb, val = self._get_offset_bool(L2.sol, theta, parameters, -1, pivoted = pivoted)
                bool_off_ub, val = self._get_offset_bool(L1.sol, theta, parameters, 1, pivoted = pivoted)
                interv = [OffsetInterval(not bool_off_lb, not bool_off_ub, floor, L2.sol, pivoted=pivoted, over_pi = val), OffsetInterval(bool_off_lb, not bool_off_ub, L2.sol, L1.sol, pivoted=pivoted, over_pi = val), 
                OffsetInterval(bool_off_lb, bool_off_ub, L1.sol, maxr, pivoted=pivoted, over_pi = val)]
            else:
                bool_off_lb, val = self._get_offset_bool(L1.sol, theta, parameters, -1, pivoted = pivoted)
                bool_off_ub, val = self._get_offset_bool(L2.sol, theta, parameters, 1, pivoted = pivoted)
                interv = [OffsetInterval(not bool_off_lb, not bool_off_ub, floor, L2.sol, pivoted=pivoted, over_pi = val), OffsetInterval(not bool_off_lb, bool_off_ub, L2.sol, L1.sol, pivoted=pivoted, over_pi = val), 
                OffsetInterval(bool_off_lb, bool_off_ub, L1.sol, maxr, pivoted=pivoted, over_pi = val)]

        return interv
    
    def _get_offset_bool(self, radius, theta, parameters, neg, pivoted):
        try: 
            val = parameters.eq_offset(radius+0.01, theta, neg, pivoted)
        except OutOfUnitaryBound:
            val = parameters.eq_offset(radius-0.01, theta, neg, pivoted)
        return abs(val) > np.pi, val < 0

    def _get_base_bool(self, radius, theta, parameters, neg, pivoted = None):
        try: 
            val = parameters.eq_base(radius+0.01, theta, neg)
        except OutOfUnitaryBound:
            val = parameters.eq_base(radius-0.01, theta, neg)
        return abs(val) > np.pi, val < 0

