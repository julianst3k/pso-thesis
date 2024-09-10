from aux import cotan, UniformRectangle, EightRectangle, ProbabilityCalculator, IntegrationLimit
import numpy as np


class AnalyticalProbability(ProbabilityCalculator):
    def __init__(self, X, Y, x_center, y_center, fov, beta, h, r, threshs):
        super().__init__(fov, beta, h, r)
        self.rect = EightRectangle(X, Y, x_center, y_center)
        self.X = X
        self.Y = Y
        self.threshs = threshs
        self.from_one = False
        self._solve_equations()
    def _solve_equations(self):
        self.lims = []
        if self.cosfov*self.b-self.a > 0:
            self.from_one = True
            self.threshs = threshs[::-1]
        filled = False
        for thresh in self.threshs:
            u = thresh["thr"]
            a = self.cosfov**2-u**2*self.sinbeta**2
            if a == 0:
                L1 = (self.cosfov**2*self.b**2 - self.a**2)/(2*u*self.sinbeta)
                L2 = None
            else:
                b = -2*u*self.sinbeta*self.a
                c = self.cosfov**2*self.b**2 - self.a**2
                L1, L2 = self._solve_quadratic(a,b,c,u)
            if L1 != None:
                if self.from_one:
                    self.lims.append(IntegrationLimit(L1, None, thresh["consts"]))
                    xL1 = L1
                    if L2 != None:
                        self.lims.append(IntegrationLimit(None, L2, thresh["consts"]))
                        xL2 = L2
                else:
                    self.lims.append(IntegrationLimit(None, L1, thresh["consts"]))
                    if L2 != None:
                        self.lims.append(IntegrationLimit(L2, None, thresh["consts"]))
            if L1 == None and not filled and not self.from_one:
                self.lims.append(IntegrationLimit(None, np.sqrt(self.X**2+self.Y**2), thresh["consts"]))
                filled = True
        self.lims.sort(key= lambda x: x.sort_radius())


        
        self.print_lims()
        remove_index = None
        for i, integr in enumerate(self.lims):
            if i == 0 and not self.from_one:
                integr.set_low(0)
            else:
                if integr.high == None:
                    try:
                        if self.lims[i+1].low is None:
                            integr.set_high(self.lims[i+1].high)
                            remove_index = i+1
                        else:
                            integr.set_high(self.lims[i+1].low)
                    except IndexError:
                        integr.set_high(np.sqrt(self.X**2+self.Y**2))

                else:
                    integr.set_low(self.lims[i-1].high)
        if remove_index is not None:
            del self.lims[remove_index]
            
    def print_lims(self):
        for lim in self.lims:
            print(lim)

    def _solve_quadratic(self, a, b, c, u):
        if b**2-4*a*c<0:
            return None, None
        sqrt = np.sqrt(b**2-4*a*c)
        L1 = (-b+sqrt)/(2*a)
        if abs(self._eq(L1)-u) > 0.01:
            L1 = None
        L2 = (-b-sqrt)/(2*a)
        if abs(self._eq(L2)-u) > 0.01:
            L2 = None

        print(L1, L2)
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
    def _eq(self, L):
        return (self.cosfov*np.sqrt(L**2+self.b**2)-self.a)/(L*self.sinbeta)
    def calculate_probability(self):
        ...
    def calculate_probability_unitary(self):
        ...
    def calculate_probability_ring(self):
        ...
if __name__ == "__main__":
    beta = np.pi/180*45
    fov = np.pi/180*25
    r = 0.8
    h = 1.2
    x_c = 1
    y_c = 1
    X = 5
    Y = 3
    threshs = [{"thr": -1, "consts": 1},
               {"thr": -0.9, "consts": {"a":-3.2, "b": -0.2}},
               {"thr": -0.6, "consts": {"a":-1.51, "b": 1.3}},
               {"thr": 0.6, "consts": {"a":-1, "b":np.pi/2}},
               {"thr": 0.9, "consts": {"a":-1.51, "b": 1.85}},
               {"thr": 1, "consts": {"a":-3.2, "b": -3.3}}]
    an_prob = AnalyticalProbability(X, Y, x_c, y_c, fov, beta, h, r, threshs)
    an_prob.print_lims()