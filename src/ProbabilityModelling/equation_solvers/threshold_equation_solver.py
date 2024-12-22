from aux import cotan, UniformRectangle, EightRectangle, ProbabilityCalculator, IntegrationLimit, Orientation
from interval import Interval, Bound, OffsetInterval
import numpy as np

class ThresholdSolver:
    def __init__(self, thresh):
        self.threshs = thresh
    def solve_equations(self, parameters):
        lims = []
        if parameters.cosfov*parameters.b-parameters.a > 0:
            parameters.from_one = True
            self.threshs = self.threshs[::-1]
            parameters.threshs = self.threshs
        print(parameters.from_one)
        filled = False
        for thresh in self.threshs:
            u = thresh["thr"]
            a = parameters.cosfov**2-u**2*parameters.sinbeta**2
            if np.abs(a) < 1e-8:
                L1 = (parameters.cosfov**2*parameters.b**2 - parameters.a**2)/(2*u*parameters.sinbeta)
                L2 = None
            else:
                b = -2*u*parameters.sinbeta*parameters.a
                c = parameters.cosfov**2*parameters.b**2 - parameters.a**2
                L1, L2 = self._solve_quadratic(a,b,c,u,parameters)
            if L1 != None and L1 >= 0:
                if parameters.from_one:
                    lims.append(IntegrationLimit(L1, None, thresh["consts"]))
                    xL1 = L1
                    if L2 != None:
                        lims.append(IntegrationLimit(None, L2, thresh["consts"]))
                        xL2 = L2
                else:
                    lims.append(IntegrationLimit(None, L1, thresh["consts"]))
                    if L2 != None:
                        lims.append(IntegrationLimit(L2, None, thresh["consts"]))
            if L1 == None and not filled and not parameters.from_one:
                lims.append(IntegrationLimit(None, np.sqrt(parameters.X**2+parameters.Y**2), thresh["consts"]))
                filled = True
            print(lims[-1], L1, L2)
        lims.sort(key= lambda x: x.sort_radius())


        
        remove_index = None
        for i, integr in enumerate(lims):
            if i == 0 and not parameters.from_one:
                integr.set_low(0)
            else:
                if integr.high == None:
                    try:
                        if lims[i+1].low is None:
                            integr.set_high(lims[i+1].high)
                            remove_index = i+1
                        else:
                            integr.set_high(lims[i+1].low)
                    except IndexError:
                        integr.set_high(np.sqrt(parameters.X**2+parameters.Y**2))
                else:
                    integr.set_low(lims[i-1].high)
        if lims[-1].high < np.sqrt(parameters.X**2+parameters.Y**2):
            """
            This occurs if there is no solutions after self.lims[-1].high. In this case, the 
            solution will converge to cos(FoV)/sin(beta), so we just need to match with the constants

            """
            converger = parameters.cosfov/parameters.sinbeta
            constant = self._find_limit_constants(converger, parameters)
            if constant is not None:
                lims.append(IntegrationLimit(lims[-1].high, np.sqrt(parameters.X**2+parameters.Y**2), constant))
        if remove_index is not None:
            del lims[remove_index]
        return lims
    def _find_limit_constants(self, converger, parameters):
        for i, constant in enumerate(self.threshs):
            if parameters.from_one:
                """
                 This case starts by 1 so first question would be, is bigger than one? 
                """
                if converger >= constant["thr"] and i==0:
                    return None
                else:
                    """
                    It will never get into negatives in infinity so we dont worry about i+1 being bigger than the
                    length
                    """
                    if converger <= constant["thr"] and converger >= self.threshs[i+1]["thr"]:
                        return constant["consts"]
            else:
                if converger >= constant["thr"] and converger <= self.threshs[i+1]["thr"]:
                    return constant["consts"]



    def _eq_offset_lims(self, L, theta, parameters):         
        u = np.sqrt(L**2+parameters.d**2+2*parameters.d*L*np.cos(theta))
        return (parameters.cosfov*np.sqrt(u**2+parameters.b**2)-parameters.a)/(u*parameters.sinbeta)

    def _gen_threshold_solutions(self, parameters):
        thresh_res = []
        for i, thresh in enumerate(self.threshs):
            u = thresh["thr"]
            a = parameters.cosfov**2-u**2*parameters.sinbeta**2
            if np.abs(a) < 1e-8:
                L1 = (parameters.cosfov**2*parameters.b**2 - parameters.a**2)/(2*u*parameters.sinbeta)
                L2 = None
            else:
                b = -2*u*parameters.sinbeta*parameters.a
                c = parameters.cosfov**2*parameters.b**2 - parameters.a**2
                L1, L2 = self._solve_quadratic(a,b,c,u,parameters)
            thresh_res.append([L1, L2])
        return thresh_res
    def solve_lims_offset(self, parameters, theta = None):
        offset_lims = {}
        if parameters.cosfov*parameters.b-parameters.a > 0:
            parameters.from_one = True
        thresh_res = self._gen_threshold_solutions(parameters)
        if theta is not None:
            output_array = []
            self._solve_lims_offset_theta(theta, thresh_res, parameters, output_array)
            return output_array
        else:
            for triangle in parameters.rect:
                offset_lims[triangle] = []
                theta = triangle.avg_ang
                self._solve_lims_offset_theta(theta, thresh_res, parameters, offset_lims[triangle])
        return offset_lims
    def _offset_array_generator(self, sols, thresh, output_array):
        if len(sols) == 0:
            return
        if len(sols) == 2:
            output_array.append(IntegrationLimit(None, sols[1], thresh["consts"]))
            output_array.append(IntegrationLimit(sols[0], None, thresh["consts"]))

        if len(sols) == 4:
            left_sols = sols[0:2]
            right_sols = sols[2:4]
            output_array.append(IntegrationLimit(left_sols[0], None, thresh["consts"]))
            output_array.append(IntegrationLimit(None, left_sols[1], thresh["consts"]))
            output_array.append(IntegrationLimit(right_sols[0], None, thresh["consts"]))
            output_array.append(IntegrationLimit(None, right_sols[1], thresh["consts"]))

        output_array.sort(key= lambda x: x.sort_radius())

    def _solve_lims_offset_theta(self, theta, thresh_res, parameters, output_array):
        filled = False
        costh = np.cos(theta)
        center = -parameters.d*costh # Cute center
        probe_value = self._eq_offset_lims(center, theta, parameters) if center**2+parameters.d**2+2*parameters.d*center != 0 else None # Lets test the value in center
        if probe_value is not None and probe_value >= -1 and probe_value <=1:
            consts_probe = self._find_limit_constants(probe_value, parameters)
            output_array.append(IntegrationLimit(center, None, consts_probe))
            output_array.append(IntegrationLimit(None, center, consts_probe))
        elif not parameters.from_one:
            output_array.append(IntegrationLimit(center, None, self.threshs[0]["consts"]))
            output_array.append(IntegrationLimit(None, center, self.threshs[0]["consts"]))
        for i, thr in enumerate(thresh_res):
            thresh = self.threshs[i]
            sols = []
            if thr[0] == None and not filled and not parameters.from_one:
                output_array.append(IntegrationLimit(None, np.sqrt(parameters.X**2+parameters.Y**2), thresh["consts"]))
            elif thr[0] is not None:
                u = thr[0]
                a = 1
                b = 2*parameters.d*np.cos(theta)
                c = parameters.d**2-u**2
                sol_one, sol_two = self._solve_quadratic(a,b,c,thresh["thr"],parameters,func=self._eq_offset_lims,theta=theta,ignore_negatives=False)
                sols.append(sol_one) if sol_one is not None else ...
                sols.append(sol_two)  if sol_two is not None else ...
                if thr[1] is not None:
                    u = thr[1]
                    c = parameters.d**2-u**2

                    sol_one, sol_two = self._solve_quadratic(a,b,c,thresh["thr"],parameters,func=self._eq_offset_lims,theta=theta,ignore_negatives=False)
                    sols.append(sol_one) if sol_one is not None else ...
                    sols.append(sol_two) if sol_two is not None else ...
            if len(sols) > 0:
                """
                We can have Zero, Two or Four solutions, depending on the case.
                If we have Two solutions, one of each side of center, then we are going to simply append at each side as usual
                If we have Four solutions, two of each side of center, then we need to append both solutions using one as the starting point and 
                one as the ending point
                """
                # If probe_value is None then we need to look for its direction, Downtrend or Uptrend. self.from_one gives us that information
                # If we have a prob_value, then our probe helps us to see their direction
                sols.sort()
                self._offset_array_generator(sols, thresh, output_array)
        remove_index = None
        for i, integr in enumerate(output_array):
            if i == 0 and not parameters.from_one and integr.low is None:
                integr.set_low(0)
            else:
                if integr.high == None:
                    try:
                        if output_array[i+1].low is None:
                            integr.set_high(output_array[i+1].high)
                            remove_index = i+1
                        else:
                            integr.set_high(output_array[i+1].low)
                    except IndexError:
                        integr.set_high(np.sqrt(parameters.X**2+parameters.Y**2))

                else:
                    integr.set_low(output_array[i-1].high)
        if output_array[-1].high < np.sqrt(parameters.X**2+parameters.Y**2):
            """
            This occurs if there is no solutions after self.lims[-1].high. In this case, the 
            solution will converge to cos(FoV)/sin(beta), so we just need to match with the constants

            """
            converger = parameters.cosfov/parameters.sinbeta
            constant = self._find_limit_constants(converger, parameters)
            if constant is not None:
                output_array.append(IntegrationLimit(output_array[-1].high, np.sqrt(parameters.X**2+parameters.Y**2), constant))
        if remove_index is not None:
            del output_array[remove_index]
        output_array = filter(lambda x: np.abs(x.high - x.low) > 0.001 and x.high > x.low, output_array)

    def _solve_quadratic(self, a, b, c, u, parameters, func=None, theta=None, ignore_negatives = True):
        args = {}
        args["parameters"] = parameters
        if func is None:
            func = self._eq
        if b**2-4*a*c<0:
            return None, None
        sqrt = np.sqrt(b**2-4*a*c)
        if theta is not None:
            args["theta"] = theta
        L1 = (-b+sqrt)/(2*a)
        args["L"] = L1
        if abs(func(**args)-u) > 0.01:
            L1 = None
        L2 = (-b-sqrt)/(2*a)
        args["L"] = L2

        if abs(func(**args)-u) > 0.01:
            L2 = None

        if ignore_negatives:
            if L1 is not None and L1 > 0:
                if L2 is None or L2 < 0:
                    return L1, None
                elif L1 < L2:
                    return L1, L2
                else:
                    return L2, L1
            else:
                if L2 is None or L2 < 0:
                    return None, None
                else:
                    return L2, None
        else:
            if L1 is not None:
                if L2 is None:
                    return L1, None
                elif L1 < L2:
                    return L1, L2
                else:
                    return L2, L1
            else:
                if L2 is None:
                    return None, None
                else:
                    return L2, None

    def _eq(self, L, parameters):
        return (parameters.cosfov*np.sqrt(L**2+parameters.b**2)-parameters.a)/(L*parameters.sinbeta)