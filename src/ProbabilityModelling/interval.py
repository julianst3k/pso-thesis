from integral_funcs import MISOOffsetIntegrator, MISOBaseIntegrator
from enum import Enum

class Bound(Enum):
    LOWER = 1
    UPPER = 2
class Interval:
    def __init__(self, offset_lb, offset_ub, lb, ub, upper_func = None, lower_func = None, consts = {}, pi_interval = False):
        """
        This is the amount of info needed for the Interval solution

        If both are offset or not offset then we use ub < lb
        If ub is offset but lb is not offset or viceversa then we use [pi, lb]U[ub,lb]
        If ub is offset then it is actually 2pi-ub
        If lb is offset then it is actually 2pi-lb
        """
        self.offset_ub = offset_ub
        self.offset_lb = offset_lb
        self.lb = self.low = lb
        self.ub = self.high = ub
        self.upper_func = upper_func
        self.lower_func = lower_func
        self.pi_interval = pi_interval
        self.decreasing = False
        self.consts = consts 
        self.is_offset = False
        if self.ub < self.lb:
            self.lb = self.ub
    def divide_interval(self, divider):
        divided_ub = self.ub
        self.ub = divider
        return Interval(self.offset_lb, self.offset_ub, divider, divided_ub, self.upper_func, self.lower_func, self.consts, self.pi_interval)
    def inverse_divided_interval(self, divider):
        divided_lb = self.lb
        self.lb = divider
        return Interval(self.offset_lb, self.offset_ub, divided_lb, divider, self.upper_func, self.lower_func, self.consts, self.pi_interval)
    def push_functional_interval(self, is_neg):
        return FunctionalInterval(self.offset_lb, self.offset_ub, self.lb, self.ub, consts = self.consts, pi_interval = self.pi_interval, is_offset = self.is_offset, is_neg = is_neg)
    def set_consts(self, consts):
        self.consts = consts
    def __str__(self):
        return f'Offset Ub: {self.offset_ub}, Offset Lb: {self.offset_lb}, Lb: {self.lb}, Ub: {self.ub}, Consts: {self.consts}, Is_Offset: {self.is_offset}'
    def is_interval(self):
        return True
    def integrate_lb(self, triang, parameters):
        return self._integrate(triang, parameters, offset_lb, True)
    def integrate_ub(self, triang, parameters):
        return self._integrate(triang, parameters, offset_ub, False)
    def _integrate(self, triang, parameters, offset, is_lb):
        integrator = MISOBaseIntegrator(self.lb, min(self.ub, triang.max_r), self.consts, parameters)
        acos_integral = (-1)**(is_lb)*integrator.acos_integrator(triang)
        angle_integral = -integrator.angle_integrator(triang)
        if offset:
            pi_const_integral = integrator.pi_const_integrator(triang)
            return (-1)**(is_lb)*acos_integral+2*pi_const_integral+angle_integral
        return (-1)**(is_lb)*acos_integral+angle_integral

class OutOfUnitaryBound(Exception):
    ...
class FunctionalInterval(Interval):
    def __init__(self, offset_lb, offset_ub, lb, ub, consts, pi_interval, is_offset, is_neg, is_pivoted = False, ub_over_pi = False):
        super().__init__(offset_lb, offset_ub, lb, ub, consts = consts, pi_interval = pi_interval)
        self.pivoted = is_pivoted
        self.is_offset = is_offset
        self.is_neg = is_neg
        self.ub_over_pi = ub_over_pi
        self.function = None
    def gen_func(self, params):
        if self.is_offset:
            return params.eq_offset_int
        elif self.pi_interval:
            return lambda x: (-1)**(not self.is_neg)*np.pi
        else:
            return params.eq_base_int
    def gen_dfunc(self, params):
        if self.is_offset:
            return params.offset_derivative_int
        elif self.pi_interval:
            return 0 
        else:
            return params.base_derivative_int
    def divide_interval(self, divider):
        divided_ub = self.ub
        self.ub = divider
        return FunctionalInterval(self.offset_lb, self.offset_ub, divider, divided_ub, self.consts, self.pi_interval, self.is_offset, self.is_neg, self.pivoted, self.ub_over_pi)
    def inverse_divided_interval(self, divider):
        divided_lb = self.lb
        self.lb = divider
        return FunctionalInterval(self.offset_lb, self.offset_ub, divided_lb, divider, self.consts, self.pi_interval, self.is_offset, self.is_neg, self.pivoted, self.ub_over_pi)

class OffsetInterval(Interval):
    def __init__(self, offset_lb, offset_ub, lb, ub, pivoted=False, upper_func = None, lower_func = None, consts = {}, over_pi = False):
        """
        This is the amount of info needed for the Interval solution

        If both are offset or not offset then we use ub < lb
        If ub is offset but lb is not offset or viceversa then we use [pi, lb]U[ub,lb]
        If ub is offset then it is actually 2pi-ub
        If lb is offset then it is actually 2pi-lb
        """
        super().__init__(offset_lb, offset_ub, lb, ub, upper_func, lower_func, consts)
        self.pivoted = pivoted
        self.is_offset = True 
        self.ub_over_pi = over_pi 
    def push_functional_interval(self, is_neg):
        return FunctionalInterval(self.offset_lb, self.offset_ub, self.lb, self.ub, consts = self.consts, pi_interval = self.pi_interval, is_offset = self.is_offset, is_neg = is_neg, is_pivoted = self.pivoted, ub_over_pi = self.ub_over_pi)
    def divide_interval(self, divider):
        divided_ub = self.ub
        self.ub = divider
        return OffsetInterval(self.offset_lb, self.offset_ub, divider, divided_ub, self.pivoted, self.upper_func, self.lower_func, self.consts, self.ub_over_pi)
    def inverse_divided_interval(self, divider):
        divided_lb = self.lb
        self.lb = divider
        return OffsetInterval(self.offset_lb, self.offset_ub, divided_lb, divider, self.pivoted, self.upper_func, self.lower_func, self.consts, self.ub_over_pi)
    def integrate_lb(self, triang, parameters):
        return self._integrate(triang, parameters, offset_lb, False, True)
    def integrate_ub(self, triang, parameters):
        return self._integrate(triang, parameters, offset_ub, over_pi, False)
    def _integrate(self, triang, parameters, is_lb):
        offset = self.is_offset
        over_pi = self.ub_over_pi 
        integrator = MISOOffsetIntegrator(self.lb, min(self.ub, triang.max_r), self.consts, parameters)
        acos_integral = integrator.acos_integrator(triang)
        atan_integral = integrator.atan_integrator(triang)
        if offset:
            pi_const_integral = integrator.pi_const_integrator(triang)
            if over_pi:
                return (-1)**(is_lb)*acos_integral-atan_integral-(self.pivoted)*pi_const_integral-2*pi_const_integral
            else:
                return (-1)**(is_lb)*acos_integral-atan_integral-(self.pivoted)*pi_const_integral+2*pi_const_integral
        return (-1)**(is_lb)*acos_integral-atan_integral-(self.pivoted)*pi_const_integral
    def _integrate_debug(self, triang, parameters):
        integrator = MISOOffsetIntegrator(self.lb, min(self.ub, triang.max_r), self.consts, parameters)
        acos_integral = integrator.acos_integrator(triang)
        return acos_integral