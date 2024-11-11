import numpy as np 
import scipy as sp 
from aux import EightRectangle
def cotan(x):
    return np.cos(x)/np.sin(x)

def sin_power_primitive(n: int):
    if n == 1:
        return -np.cos(x)
    elif n == 0:
        return x
    else:
        return -np.cos(x)*np.sin(n-1)/n - (n-1)/n * sin_power_primitive(n-2, x)
def cos_power_primitive(n: int, x):
    if n == 1:
        return np.sin(x)
    elif n == 0:
        return x
    else:
        return np.sin(x)*np.cos(x)**(n-1)/n - (n-1)/n * cos_power_primitive(n-2, x)
class NonOriginIntegrator:
    def __init__(self, lb, ub, consts, parameters):
        self.lb = lb
        self.ub = ub
        self.consts = consts
        self.params = parameters
    def integral_even(self, Dn, theta_bot, theta_top, n):
        theta_cr = theta_top
        max_Dn = Dn/np.cos(theta_cr)
        if self.ub > max_Dn:
            return 1/(Dn**2*np.tan(theta_cr))*self.integral_cosine(Dn, theta_bot, theta_top, n), theta_top
        elif self.ub < Dn:
            return 1/(Dn**2*np.tan(theta_cr))*self.integral_linear(Dn, theta_bot, theta_top, n), None
        else:
            theta_mid = np.arccos(Dn/self.ub)
            cosine_int = 1/(Dn**2*np.tan(theta_cr))*self.integral_cosine(Dn, theta_mid, theta_top, n)
            linear_int = 1/(Dn**2*np.tan(theta_cr))*self.integral_linear(Dn, theta_bot, theta_mid, n)
            return cosine_int+linear_int, theta_mid
    def integral_cosine(self, Dn, theta_bot, theta_top, n):
        sum_int_a, sum_int_b, sum_int_c = self.integral_lambdas()
        ...
    def integral_linear(self, Dn, theta_bot, theta_top, n):
        sum_int_a, sum_int_b, sum_int_c = self.integral_lambdas()
        sum_int = sum_int_a(self.ub)-sum_int_a(self.lb)+sum_int_b(self.ub)-sum_int_b(self.lb)+sum_int_c(self.ub)-sum_int_c(self.lb)
        return (theta_top-theta_bot)*sum_int
    def integral_lambdas_cosine(self, Dn):
        a = self.consts["a"]
        b = self.consts["b"]
        cosfov = self.params.cosfov
        sinbeta = self.params.sinbeta
        sigma = Dn**2/self.params.b**2
        sum_int_one = lambda x, t: cosfov/sinbeta*a*Dn*b*(np.tan(t)*np.sqrt(np.cos(2*t)+1+2*sigma**2)/np.sqrt(2))+
        np.sqrt(1+sigma**2)*sp.special.ellipkinc(t, 1/(1+sigma**2))-np.sqrt(1+sigma**2)*sp.special.ellipeinc(t, 1/(1+sigma**2))
        sum_int_two = lambda x, t: self.sinh_int_cos(x, t)
    def sinh_int_cos(self, x, t):
        ...
    def integral_lambdas_linear(self):
        a = self.consts["a"]
        b = self.consts["b"]
        cosfov = self.params.cosfov
        sinbeta = self.params.sinbeta
        sum_int_a = lambda x: a*(0.5*cosfov/sinbeta*(x*np.sqrt(x**2+self.params.b**2)+self.params.b**2*np.ln(np.sqrt(x**2+self.params.b**2)+x)))
        sum_int_b = lambda x: self.params.a/sinbeta*X
        sum_int_c = lambda x: b*x**2
        return sum_int_a, sum_int_b, sum_int_c
    def integral_odd(self, Dn, theta_bot, theta_top, n):

class OriginIntegrator:
    def __init__(self, radius):
        self.radius = radius
    def unitary_integral_even(self, Dn, theta_bot, theta_top, n):
        L1 = self.radius
        theta_cr = theta_top
        max_Dn = Dn/np.cos(theta_cr)
        if L1 > max_Dn:
            return 1, None
        elif L1 < Dn:
            return (theta_top-theta_bot)*L1**2/((Dn**2)*np.tan(theta_cr)), None
        else:
            theta_mid = np.arccos(Dn/L1)
            return 1/np.tan(theta_cr)*np.tan(theta_mid)-1/np.tan(theta_cr)*np.tan(theta_bot)+(theta_cr-theta_mid)*L1**2/((Dn**2)*np.tan(theta_cr)), theta_mid
    def unitary_integral_odd(self, Dn, theta_bot, theta_top, n):
        L1 = self.radius
        theta_cr = theta_bot
        max_Dn = Dn/np.sin(theta_cr)
        if L1 > max_Dn:
            return 1, None
        elif L1 < Dn:
            return (theta_top-theta_bot)*L1**2/((Dn**2)*cotan(theta_cr)), None
        else:
            theta_mid = np.arcsin(Dn/L1)
            return 1/cotan(theta_cr)*cotan(theta_mid)-1/cotan(theta_cr)*cotan(theta_top)+(theta_mid-theta_cr)*L1**2/((Dn**2)*cotan(theta_cr)), theta_mid

class RectangleIntegrator:
    def __init__(self, rec):
        self.rect = rec
    
    def _integrate(func):
        def _integrator_wrapper(self, *args):
            integrator = func(self, *args)
            X, Y = self.rect.X, self.rect.Y
            int_sum = 0 
            for i, tri in enumerate(self.rect):
                if i % 2 == 0:
                    Dn = tri.x if i%4 == 0 else tri.y
                    integr, crt = integrator.unitary_integral_even(Dn, tri.ang_crt, tri.ang_high, i)
                else:
                    Dn = tri.y if i%4 == 1 else tri.x
                    integr, crt = integrator.unitary_integral_odd(Dn, tri.ang_crt, tri.ang_high, i)
                int_sum += integr*tri.get_area()/(X*Y)
                if crt is not None:
                    tri.change_ang(crt)

            return int_sum
        return _integrator_wrapper   
    @_integrate
    def unitary_integrator(self, L1):
        integrator = OriginIntegrator(L1)
        return integrator

    


if __name__=="__main__":
    X, Y, xc, yc = 5, 3, 1, 1
    rec = EightRectangle(X, Y, xc, yc)
    intrec = RectangleIntegrator(rec)
    L1 = 1
    int_sum = intrec.unitary_integrator(L1)
    print(int_sum)
    print(np.pi/(15))
