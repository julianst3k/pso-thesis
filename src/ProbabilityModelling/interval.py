from integral_funcs import MISOOffsetIntegrator, MISOBaseIntegrator
from enum import Enum
import numpy as np

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
            return lambda x, *args: (-1)**(self.is_neg)*np.pi
        else:
            return params.eq_base_int
    def set_lb_ub(self, lb, ub):
        self.lb = lb 
        self.ub = ub
        return self
    def gen_dfunc(self, params):
        if self.is_offset:
            return params.offset_derivative_int
        elif self.pi_interval:
            return lambda x, *args: 0
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
    def integrate_lb(self, triang, parameters):
        return self._integrate(triang, parameters, self.offset_lb, False, True)
    def integrate_ub(self, triang, parameters):
        return self._integrate(triang, parameters, self.offset_ub, self.ub_over_pi, False)
    def riemman_integral(self, triang, params):
        angs = np.linspace(triang.ang_low, triang.ang_high, 100)
        radius = np.linspace(self.lb, min(self.ub, triang.max_r), 100)
        diff_ang = angs[1]-angs[0]
        diff_rad = radius[1]-radius[0]
        summ_tot = 0
        summ_arc = 0
        summ_atan = 0
        for r in radius:
            for ang in angs:
                func = params.eq_offset_int_debug
                func_ret = func(r, self, ang, self.offset_lb if self.is_neg else self.offset_ub)
                summ_tot += diff_ang*diff_rad*func_ret[0]*r
                summ_arc += diff_ang*diff_rad*func_ret[1]*r
                summ_atan += diff_ang*diff_rad*func_ret[2]*r
        return summ_tot, summ_arc, summ_atan
    def _integrate(self, triang, parameters, offset, over_pi, is_lb):
        if self.is_offset:
            integrator = MISOOffsetIntegrator(self.lb, min(self.ub, triang.max_r), self.consts, parameters)
            if self.consts == 1:
                pi_const_integral = integrator.pi_const_integrator(triang)
                return (-1)**(is_lb)*pi_const_integral
            acos_integral = integrator.acos_integrator(triang)
            atan_integral = integrator.atan_integral(triang)
            if np.abs(triang.avg_ang-1.13) < 0.01:
                print(f"Analytical Integral {acos_integral}, {atan_integral}")
            #print(acos_integral, atan_integral, self.riemman_integral(triang, parameters), self.lb, self.ub, triang.ang_high, triang.ang_low)
            pi_const_integral = integrator.pi_const_integrator(triang)
            if offset:
                if over_pi:
                    return (-1)**(is_lb)*acos_integral-atan_integral-2*pi_const_integral-(self.pivoted)*pi_const_integral
                else:
                    return (-1)**(is_lb)*acos_integral-atan_integral+2*pi_const_integral-(self.pivoted)*pi_const_integral
            return (-1)**(is_lb)*acos_integral-atan_integral-(self.pivoted)*pi_const_integral
        if self.pi_interval:
            integrator = MISOBaseIntegrator(self.lb, min(self.ub, triang.max_r), self.consts, parameters)
            pi_const_integral = integrator.pi_const_integrator(triang)
            return (-1)**(is_lb)*pi_const_integral            
        else:
            integrator = MISOBaseIntegrator(self.lb, min(self.ub, triang.max_r), self.consts, parameters)
            if self.consts == 1:
                pi_const_integral = integrator.pi_const_integrator(triang)
                return (-1)**(is_lb)*pi_const_integral
            acos_integral = integrator.acos_integrator(triang)
            angle_integral = -integrator.angle_integrator(triang)
            if offset:
                pi_const_integral = integrator.pi_const_integrator(triang)
                return (-1)**(is_lb)*acos_integral+2*pi_const_integral+angle_integral
            return (-1)**(is_lb)*acos_integral+angle_integral

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
    def _integrate_debug(self, triang, parameters):
        integrator = MISOOffsetIntegrator(self.lb, min(self.ub, triang.max_r), self.consts, parameters)
        acos_integral = integrator.acos_integrator(triang)
        return acos_integral