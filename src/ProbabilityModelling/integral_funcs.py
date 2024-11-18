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
        sum_int_one, sum_int_two, sum_int_three, sum_int_four = self.integral_lambdas_cosine(Dn)
        inner_sum = (-sum_int_a(self.lb)-sum_int_b(self.lb)-sum_int_c(self.lb))*(theta_top-theta_bot)
        outer_sum = sum_int_one(theta_top)-sum_int_one(theta_bot)+sum_int_two(theta_bot, theta_top)+sum_int_three(theta_top)-sum_int_three(theta_bot)+sum_int_four(theta_top)-sum_int_four(theta_bot)
        return outer_sum+inner_sum
    def integral_sine(self, Dn, theta_bot, theta_top, n):
        sum_int_a, sum_int_b, sum_int_c = self.integral_lambdas()
        sum_int_one, sum_int_two, sum_int_three, sum_int_four = self.integral_lambdas_sine(Dn)
        inner_sum = (-sum_int_a(self.lb)-sum_int_b(self.lb)-sum_int_c(self.lb))*(theta_top-theta_bot)
        outer_sum = sum_int_one(theta_top)-sum_int_one(theta_bot)+sum_int_two(theta_bot, theta_top)+sum_int_three(theta_top)-sum_int_three(theta_bot)+sum_int_four(theta_top)-sum_int_four(theta_bot)
        return outer_sum+inner_sum
    
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
        sum_int_one = lambda t: (cosfov/sinbeta*a*Dn*b*(np.tan(t)*np.sqrt(np.cos(2*t)+1+2*sigma**2)/np.sqrt(2))+ 
        np.sqrt(1+sigma**2)*sp.special.ellipkinc(t, 1/(1+sigma**2))-np.sqrt(1+sigma**2)*sp.special.ellipeinc(t, 1/(1+sigma**2)))
        sum_int_two = lambda tb, tt: cosfov/sinbeta*a*self.sinh_int_cos(x, tb, tt, Dn)
        sum_int_three = lambda t: self.params.a/sinbeta*Dn*self.log_cos(t)
        sum_int_four = lambda t: b*Dn/2*np.tan(t)
        return sum_int_one, sum_int_two, sum_int_three, sum_int_four
    def integral_lambdas_sine(self, Dn):
        a = self.consts["a"]
        b = self.consts["b"]
        cosfov = self.params.cosfov
        sinbeta = self.params.sinbeta
        sigma = Dn**2/self.params.b**2
        sum_int_one = lambda t: (cosfov/sinbeta*a*Dn*b*(-cotan(t)*np.sqrt(-np.cos(2*t)+1+2*sigma**2)/np.sqrt(2))+
        (sigma+1/sigma)*sp.special.ellipkinc(t, -1/sigma**4)-sigma*sp.special.ellipeinc(t, -1/sigma**4))
        sum_int_two = lambda tb, tt: cosfov/sinbeta*a*self.sinh_int_cos(x, tb, tt, Dn)
        sum_int_three = lambda t: self.params.a/sinbeta*Dn*np.log(np.sin(t)/np.cos(t))
        sum_int_four = lambda t: b*Dn/2*cotan(t)
        return sum_int_one, sum_int_two, sum_int_three, sum_int_four
    def sinh_int_cos(self, tb, tt, Dn, max_order = 10):
        b = self.params.b 
        summation = 0
        multiplier = b**2
        if Dn < b:
            tm = np.arccos(Dn/b)
            if tt < tm:
                tm = tt
            for order in range(max_order):
                factor = 1/sp.special.factorial(order)
                for j in range(order):
                    factor *= (-1/2-j)
                summation += 1/(2*order+1)*factor*(self.f_sec(tm, order)-self.f_sec(tb, order))*(Dn/b)**(2*order+1)
            for order in range(1,max_order):
                factor = 1/sp.special.factorial(order)
                for j in range(order):
                    factor *= (-1/2-j)
                summation -= 1/(2*order)*factor*(self.f_cos(tt, order)-self.f_cos(tm, order))

            summation += np.log(2*Dn/b)*(tt-tm)-(self.f_logcos(tt)-self.f_logcos(tm))
        else:
            summation += np.log(2*Dn/b)*(tt-tb)-(self.f_logcos(tt)-self.f_logcos(tb))
            for order in range(1,max_order):
                factor = 1/sp.special.factorial(order)
                for j in range(order):
                    factor *= (-1/2-j)
                summation -= 1/(2*order)*factor*(self.f_cos(tt, order)-self.f_cos(tb, order))
        summation += np.log(b)*(tt-tb)
        return summation*multiplier
    def sinh_int_sin(self, tb, tt, Dn, max_order = 10, debug = False):
        """
        TODO: Test against Wolfram values
        integral of arcsinh(Dn/(bsin(x))) 
        """
        b = self.params.b if not debug else 1
        summation = 0
        multiplier = b**2
        if Dn < b:
            tm = np.arccos(Dn/b)
            if tt < tm:
                tm = tt
            for order in range(max_order):
                factor = 1/sp.special.factorial(order)
                for j in range(order):
                    factor *= (-1/2-j)
                summation += 1/(2*order+1)*factor*(self.f_cosec(tt, order)-self.f_cosec(tm, order))*(Dn/b)**(2*order+1)
            for order in range(1,max_order):
                factor = 1/sp.special.factorial(order)
                for j in range(order):
                    factor *= (-1/2-j)
                summation -= 1/(2*order)*factor*(self.f_sin(tm, order)-self.f_sin(tb, order))*(Dn/b)**(-2*order)

            summation += np.log(2*Dn/b)*(tt-tm)-(self.f_logsin(tt)-self.f_logsin(tm))
        else:
            summation += np.log(2*Dn/b)*(tt-tb)-(self.f_logsin(tt)-self.f_logsin(tb))
            for order in range(1,max_order):
                factor = 1/sp.special.factorial(order)
                for j in range(order):
                    factor *= (-1/2-j)
                summation -= 1/(2*order)*factor*(self.f_sin(tt, order)-self.f_sin(tb, order))
        summation += np.log(b)*(tt-tb)
        return summation*multiplier
    def f_sec(self, t, order):
        """
            Tested: Yes
            Integral of 1/cos^(2*order+1)(x)
        """
        sum_base = np.tan(t)
        if order > 0:
            sum_base *= 1/np.cos(t)**(2*order-1)
        else:
            print(np.log((np.sin(t/2)+np.cos(t/2))/(np.cos(t/2)-np.sin(t/2))))
            return np.log((np.sin(t/2)+np.cos(t/2))/(np.cos(t/2)-np.sin(t/2)))
        summation = 0
        multiplier = 1
        for i in range(order):
            if i==0:
                multiplier *= 1/(2*order-2*i)
            else:
                multiplier *= (order-(2*i-1)/2)/(order-i)
            summation += multiplier*np.cos(t)**(2*i)
        summation *= sum_base
        summation += multiplier*np.log((np.sin(t/2)+np.cos(t/2))/(np.cos(t/2)-np.sin(t/2)))
        return summation
    def f_cosec(self, t, order):
        """
            Tested: Yes
            Integral of 1/sin^(2*order+1)(x)
        """
        sum_base = -cotan(t)
        if order > 0:
            sum_base *= 1/np.sin(t)**(2*order-1)
        else:
            return np.log(np.sin(t/2)/np.cos(t/2))
        summation = 0
        multiplier = 1
        for i in range(order):
            if i==0:
                multiplier *= 1/(2*order-2*i)
            else:
                multiplier *= (order-(2*i-1)/2)/(order-i)
            summation += multiplier*np.sin(t)**(2*i)
        summation *= sum_base
        summation += multiplier*np.log(np.sin(t/2)/np.cos(t/2))
        return summation
    def f_cos(self, t, order):
        """
            Tested: Yes
            Integral of cos^(2*order)(x)
        """
        sum_base = np.sin(t)*np.cos(t)**(2*order-1)
        summation = 0
        multiplier = 1
        if order == 0:
            return t
        for i in range(order):
            if i==0:
                multiplier *= 1/(2*order-2*i)
            else:
                multiplier *= (order-(2*i-1)/2)/(order-i)
            summation += multiplier*np.cos(t)**(-2*i)
        summation *= sum_base
        summation += multiplier*t
        return summation
    def f_sin(self, t, order):
        """
        Tested: Yes
        Integral of sin^(2*order)(x)
        """
        sum_base = -np.cos(t)*np.sin(t)**(2*order-1)
        summation = 0
        multiplier = 1
        if order == 0:
            return t
        for i in range(order):
            if i==0:
                multiplier *= 1/(2*order-2*i)
            else:
                multiplier *= (order-(2*i-1)/2)/(order-i)
            summation += multiplier*np.sin(t)**(-2*i)
        summation *= sum_base
        summation += multiplier*t
        return summation
    def f_logcos(self, t):
        """
        Tested: Yes
        Integral of log(cos(x))
        """
        return -(np.pi/2-t)*np.log(np.pi/2-t)+(np.pi/2-t)+1/(3*6)*(np.pi/2-t)**3
    def f_logsin(self, t):
        """
        Tested: Yes
        Integral of log(sin(x))
        """
        return t*np.log(t)-t-1/(3*6)*t**3
    def log_cos(self, t):
        return np.log((np.sin(t/2)+np.cos(t/2))/(np.cos(t/2)-np.sin(t/2)))
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
        ...

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
    intrec = NonOriginIntegrator(1,2, None, None)
    for Dn in np.arange(0.8,1.4, 0.2):
        print([(x, intrec.sinh_int_sin(x, np.pi/2, Dn, debug = True)) for x in np.linspace(np.pi/10, np.pi/2, 10)])
