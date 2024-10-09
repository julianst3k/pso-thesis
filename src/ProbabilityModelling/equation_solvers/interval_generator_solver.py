from aux import cotan, UniformRectangle, EightRectangle, ProbabilityCalculator, IntegrationLimit, Orientation, Interval, Bound, OffsetInterval
import numpy as np

class IntervalOffsetSolver:
    """
    The responsibility of this class is managing intervals from the offset solving.

    """
    def __init__(self):
        ...

    def base_intervals_solver(self, L1, L2, theta, parameters):
        """
        The method takes two parameters, L1 and L2, which code information given certain requirements
        L1 can be:
        - False: In this case the offset never happens and it stays permanently 
        - True: In this case the offset always happens and it stays permanently
        - A number: In this case the offset happens after a value
        L2 can be:
        - False, True: According to the previous case
        - A number: In this case, two offset happen after a value. As we know the behavior of the equation,
        we know that these two offsets exist in an specific case
        - None: In this case, the L1 offset happens after a value and we need to know which side is crossing (Upper or lower)
        """
        epsilon = 0.001
        if not L1:
            interv = [Interval(L1, L2, 0, np.sqrt(parameters.X**2+parameters.Y**2))]
        if L1 and L2:
            interv = [Interval(True, True, 0, np.sqrt(parameters.X**2+parameters.Y**2))]
        if L2 == None:
            if parameters.cosfov*np.sqrt(parameters.b**2)-parameters.a < 0:
                """
                In this case it only crosses once
                Positive + -np.pi => False to True, and lb would be True True
                Negative + -np.pi => True to False, and ub would be False False
                It should never cross np.pi because the maximum is np.pi, 
                """
                if np.abs(parameters.eq_base(L1, theta)+np.pi) < epsilon:
                    interv = [Interval(True, False, 0, L1), Interval(True, True, L1, np.sqrt(parameters.X**2+parameters.Y**2))]
                elif np.abs(parameters.eq_base(L1, theta, -1)+np.pi) < epsilon:
                    interv = [Interval(True, False, 0, L1), Interval(False, False, L1, np.sqrt(parameters.X**2+parameters.Y**2))]
            else:
                """
                In this case if it cross only once in the positives it would be 
                Positive + -np.pi => True to False, and lb would be True True 
                Negative + -np.pi => False to True, and ub would be False False
                It should never cross np.pi because the maximum is np.pi, 
                """
                if np.abs(parameters.eq_base(L1, theta)+np.pi) < epsilon:
                    interv = [Interval(True, True, 0, L1), Interval(True, False, L1, np.sqrt(parameters.X**2+parameters.Y**2))]
                elif np.abs(parameters.eq_base(L1, theta, -1)+np.pi) < epsilon:
                    interv = [Interval(False, False, 0, L1), Interval(True, False, L1, np.sqrt(parameters.X**2+parameters.Y**2))]

        elif L1<np.sqrt(parameters.X**2+parameters.Y**2):
            """
                The only case it crosses twice is in the case of self.cosfov*np.sqrt(self.b**2)-self.a > 0
                In this case it would be
                Positive + -np.pi => True to False, False to True, and True True True
                Negative + -np.pi => False to True, True to False, and False, False, False
            """
            if np.abs(parameters.eq_base(L2, theta)+np.pi) < epsilon:
                interv = [Interval(True, True, 0, L2), Interval(True, False, L2, L1), Interval(True, True, L1, np.sqrt(parameters.X**2+parameters.Y**2))]
            elif np.abs(parameters.eq_base(L2, theta, -1)+np.pi) < epsilon:
                interv = [Interval(False, False, 0, L2), Interval(True, False, L2, L1), Interval(False, False, L1, np.sqrt(parameters.X**2+parameters.Y**2))]
        elif L2<np.sqrt(parameters.X**2+parameters.Y**2):
            if np.abs(parameters.eq_base(L2, theta)+np.pi) < epsilon:
                interv = [Interval(True, True, 0, L2), Interval(True, False, L2, L1)]
            elif np.abs(parameters.eq_base(L2, theta, -1)+np.pi) < epsilon:
                interv = [Interval(False, False, 0, L2), Interval(True, False, L2, L1)]
    def offset_intervals_solver(self, L1l, L2l, L1u, L2u, pivot, theta, parameters, flagl = None, flagu = None):
        if pivot is None:
            """
            In this case there is no transition so we just need to generate the intervals
            This case is equivalent to the non offset solver! But since we coded the results differently, we need to
            do another function
            """
            interv = self._offset_intervals_no_pivot_solver(L1, L2, theta, parameters)
        else:
            """
            In this case, there is a transition between pivot and not pivot that we need to take into account
            To do this, we use the following scheme, but 
            Pre-pivot -> The function acts in reverse (Increasing over distance).
            """
            """
            There is a transition so we need to get the interval AND the transition
            """
            interv_sol_one = self._offset_intervals_pivot_solver_lower(L1l, L2l, pivot, theta, parameters, flagl)
            interv_sol_two = self._offset_intervals_pivot_solver_upper(L1u, L2u, pivot, theta, parameters, flagu)
            interv_sol_one.extend(interv_sol_two)
            interv = interv_sol_one
        return interv
    def _offset_intervals_pivot_solver_lower(self, L1, L2, pivot, theta, parameters, flag):
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
    def _offset_intervals_pivot_solver_upper(self, L1, L2, pivot, theta, parameters, flag):
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
    
    
    def _offset_intervals_no_pivot_solver(self, L1, L2, theta, parameters):
        maxr = np.sqrt(self.X**2+self.Y**2)
        epsilon = 0.001
        if flag == 0:
            
            interv = [OffsetInterval(L1, L2, 0, maxr)]
        
        if flag == 1:

            interv = [OffsetInterval(True, L2, 0, L1), OffsetInterval(False, L2, L1, maxr)]

        if flag == 2:

            interv = [OffsetInterval(L1, True, 0, L2), OffsetInterval(L1, False, L2, maxr)]
        
        return interv
