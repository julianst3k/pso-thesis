import numpy as np 
import scipy as sp 
from aux import UniformRectangle
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
class ParamMocker:
    def __init__(self, d, b):
        self.d = d
        self.b = b
class MISOOffsetIntegrator:
    def __init__(self, lb, ub, consts, parameters):
        self.lb = lb
        self.ub = ub
        self.consts = consts
        self.params = parameters
    def pi_const_integrator(self, triang):
        theta_bot, theta_top = triang.ang_low, triang.ang_high
        return (theta_top-theta_bot)*(self.ub**2-self.lb**2)*np.pi/2
    def acos_integrator(self, triang):
        # Los lambdas tienen distintas entradas por como se generan las soluciones :(
        lambda_over_one_wrap, lambda_over_two_wrap, lambda_under_one_wrap, lambda_under_two_wrap = self._acos_lambda_wrapper() 
        xb = self.lb
        xt = self.ub
        d = self.params.d
        b = self.params.b
        a = self.params.a
        avg_t = triang.avg_ang
        tb, tt = triang.ang_low, triang.ang_high
        x_thresh_u, x_thresh_d = self.solve_b_u_threshold(avg_t, d, b)
        """
        Se genera un intervalo entre x_u y x_d, el cual se intersecta con [xb, xt]
        """
        integr = 0
        if x_u is None or xb >= x_u or xt <= x_d: # No hay interseccion, como es creciente entonces estamos por sobre b
            integr += lambda_under_one_wrap(xt, xb, tt, tb)
            integr += lambda_under_two_wrap(xt, xb, tt, tb)
        elif xt < x_u and xb >= x_d: # Hay full interseccion, como es creciente entonces estamos por bajo b
            integr += lambda_over_one_wrap(xt, xb, tt, tb)
            integr += lambda_over_two_wrap(xt, xb, tt, tb)
        elif xt > x_u and xb < x_u: # [xb, x_u] bajo b, [x_u, xt] sobre b
            integr += lambda_over_one_wrap(x_u, xb, tt, tb)
            integr += lambda_over_two_wrap(x_u, xb, tt, tb)
            integr += lambda_under_one_wrap(xt, x_u, tt, tb)
            integr += lambda_under_two_wrap(xt, x_u, tt, tb)
        elif xt >= x_d and xb < x_d: # [xb, x_d] sobre b, [x_d, xt] bajo b
            integr += lambda_over_one_wrap(xt, x_d, tt, tb)
            integr += lambda_over_two_wrap(xt, x_d, tt, tb)
            integr += lambda_under_one_wrap(x_d, xb, tt, tb)
            integr += lambda_under_two_wrap(x_d, xb, tt, tb)
        lambda_one_a_constant, _ = self.acos_lambda_under_b(use_discerner=False)
        cosfov = self.params.cosfov
        sinbeta = self.params.sinbeta
        integr *= cosfov/sinbeta
        aux_integr = lambda_one_a_constant(xt, tt)-lambda_one_a_constant(xt,tb)-lambda_one_a_constant(xb, tt)+lambda_one_a_constant(xb, tb)
        integr += a*(aux_integr)/(b*sinbeta)
        integr *= self.consts["a"]
        integr += self.consts["b"]/2*(xt**2-xb**2)*(tt-tb)
        return integr
    def _acos_lambda_wrapper(self, N = 10):
        lambda_over_one, lambda_over_two = self.acos_lambda_over_b(N) # f(x,t), f(xt, xb, t)
        lambda_under_one, lambda_under_two = self.acos_lambda_under_b() # f(x,t), f(x, tt, tb) 
        lambda_over_one_wrap = lambda xt, xb, tt, tb: lambda_over_one(xt, tt)-lambda_over_one(xt,tb)-lambda_over_one(xb, tt)+lambda_over_one(xb, tb)
        lambda_over_two_wrap = lambda xt, xb, tt, tb: lambda_over_two(xt, xb, tt)-lambda_over_two(xt, xb, tb)
        lambda_under_one_wrap = lambda xt, xb, tt, tb: lambda_under_one(xt, tt)-lambda_under_one(xt,tb)-lambda_under_one(xb, tt)+lambda_under_one(xb, tb)
        lambda_under_two_wrap = lambda xt, xb, tt, tb: lambda_under_two(xt, tt, tb)-lambda_under_two(xb, tt, tb)
        return lambda_over_one_wrap, lambda_over_two_wrap, lambda_under_one_wrap, lambda_under_two_wrap
    def solve_b_u_threshold(self, avg_t, d, b):
        a = 1
        b = 2*d*np.cos(avg_t)
        c = d**2-b**2
        if b**2-4*a*c >= 0:
            sol_one = (-b+np.sqrt(b**2-4*a*c))/(2*c)
            sol_two = (-b-np.sqrt(b**2-4*a*c))/(2*c)
        else:
            """
            Para que lo de arriba no ocurra, c > 0, lo que implica que siempre es mayor a b**2
            """
            return None, None

    def acos_lambda_over_b(self, N = 10):
        b = self.params.b
        lambda_one = lambda x, t: 1/2*x**2*t
        lambda_two = lambda xt, xb, t: b**2/4*self.arctan_acos_integral(xt,xb,t, N)
        return lambda_one, lambda_two
    def acos_lambda_under_b(self, use_discerner = True):
        """
        Tested: Yes
        """
        def sign_discerner(x,t,d,use_discerner):
            multiplier = d*np.sin(t)
            eta = 4*d*x/(x+d)**2
            if t<np.pi or not use_discerner:
                summ = np.log(x+d*np.cos(t)+np.sqrt(x**2+2*d*x*np.cos(t)+d**2))-1
                summ *= multiplier
                summ += (x+d)*sp.special.ellipeinc(t/2, eta)-(x**2-d**2)/(x+d)*sp.special.ellipkinc(t/2, eta)
                print(sp.special.ellipkinc(t/2, eta), sp.special.ellipeinc(t/2, eta), eta, t/2)
            elif t>=np.pi:
                summ = np.log(-x-d*np.cos(t)+np.sqrt(x**2+2*d*x*np.cos(t)+d**2))-1
                summ *= multiplier
                summ += (x**2-d**2)/(x+d)*sp.special.ellipkinc(t/2, eta)-(x+d)*sp.special.ellipeinc(t/2, eta)
            return summ
        b = self.params.b
        d = self.params.d
        lambda_one = lambda x, t: b*sign_discerner(x,t,d, use_discerner) # <- Tested 
        lambda_two = lambda x, tt, tb: self.x_dcos_ineq_integral(x, tt, tb, d)
        return lambda_one, lambda_two
    def x_dcos_ineq_integral(self, x, tt, tb, d):
        def cos_expansion(x, d, t, N = 10):
            summ = 0
            for n in range(1,N):
                summ += (d/x)**n*(-1)**n/(n*(n+2))*(np.cos(t)**(n+2)*np.sin(t)-self.f_cos(t, n+3))
        lambda_upper = lambda x, t: x**3/3*t + x**2/2*d*np.sin(t)+1/4*d**2*x*t-1/8*x*d**2*np.sin(2*t)-1/2*d**3*(np.sin(t)**3/3*np.log(x)+
        cos_expansion(x,d,t))
        lambda_lower = lambda x, t: -x**2/2*d*np.cos(t)+x**4/(8*d)*np.log(np.tan(t/2))+d*x**3/3*np.log(np.sin(t))+x**2*d**2/(4*d)*(np.cos(t)+np.log(np.tan(t/2)))
        
        if x >= -d*np.cos(tt):
            summ = lambda_upper(x, tt)
            if x >= -d*np.cos(tb):
                summ -= lambda_upper(x, tb)
            else:
                tm = np.arccos(-L/d)
                if tb > tm:
                    tm = 2*np.pi-tm  
                summ -= lambda_upper(x, tm)
                summ -= lambda_lower(x, tm)
                summ += lambda_lower(x, tb)
        else:
            if x >= -d*np.cos(tb):
                summ -= lambda_upper(x, tb)
                if tb > tm:
                    tm = 2*np.pi-tm
                summ += lambda_upper(x, tm)
                summ += lambda_lower(x, tm)
                summ -= lambda_lower(x, tt)

            else:
                summ -= lambda_lower(x, tt)-lambda_lower(x,tb)
        return summ


    def arctan_acos_integral(self, xt, xb, t, N):
        """
        Integral of 
        log(xt^2+d^2+2dxtcos(t))-dcot(t)tan^-1(xt+dcos(t)/dsin(t)) 
        """
        def arctan_acos_sum(self, xt, xb, t, d, N):
            xt_coef = (2*d*xt)/(xt**2+d**2)
            xb_coef = (2*d*xb)/(xb**2+d**2)
            summ = 0
            for n in range(1,N+1):
                summ += xt_coef**n*(-1)**(n+1)*self.f_cos(t, n)/n
                summ -= xb_coef**n*(-1)**(n+1)*self.f_cos(t, n)/n
            return summ
        def arctan_tanh_expr(self, xt, xb, t, d):
            eta = (xt*xb+d**2)/(2*d*(xt+xb))
            multiplier = d/2*(xt-xb)/(xt+xb)
            if eta**2 > 1:
                summ = -2*eta/np.sqrt(eta**2-1)*np.arctanh((eta-1)/np.sqrt(eta**2-1)*np.tan(t/2)**2)    
            else:
                summ = 2*eta/np.sqrt(1-eta**2)*np.arctanh((eta-1)/np.sqrt(eta**2-1)*np.tan(t/2)**2)    
            summ += t
            return summ*multiplier
        d = self.params.d
        return arctan_acos_sum(xt, xb, t, d, N) + t*np.log((xt**2+d**2)/(xb**2+d**2)) - arctan_tanh_expr(xt, xb, t, d)
    def f_cos(self, t, order):
        """
            Tested: Yes
            Integral of cos^(order)(x)
        """
        sum_base = np.sin(t)*np.cos(t)**(order-1)
        summation = 0
        multiplier = 1
        if order == 0:
            return t
        print(order)
        for i in range(order//2):
            print(i)
            if i==0:
                multiplier *= 1/(order)
            else:
                multiplier *= (order-2*i+1)/(order-2*i)
            summation += multiplier*np.cos(t)**(-2*i)
            
        summation *= sum_base
        if order % 2:
            summation += multiplier*np.sin(t)
        else:
            summation += multiplier*t
        return summation    
    def atan_integral(self, triang):
        """
        Tested = Necesito mas ejemplos pero por ahora esta bien
        """
        lambda_list = self.atan_lambdas(triang)
        xt = self.ub
        xb = self.lb
        tt = triang.ang_high
        tb = triang.ang_low            
        summ = 0
        d = self.params.d
        subsum = np.zeros(3)
        for i, lamb in enumerate(lambda_list):
            presum = summ
            if i == 4:
                if np.abs(xt-d) >= 0.1:
                    summ += lamb(xt, tt)-lamb(xt,tb)
                if np.abs(xb-d) >= 0.1:
                    summ -= lamb(xb, tt)-lamb(xb,tb)
                subsum[2] += lamb(xt, tt)-lamb(xt,tb) if np.abs(xt-d) >= 0.1 else 0
            else:
                summ += lamb(xt, tt)-lamb(xt,tb)-lamb(xb, tt)+lamb(xb, tb)
                if i < 2:
                    subsum[0] += lamb(xt, tt)-lamb(xt,tb)
                    
                elif i == 2:
                    subsum[1] += lamb(xt, tt)-lamb(xt,tb)
                else:
                    subsum[2] += lamb(xt, tt)-lamb(xt,tb)
        return summ, subsum
    def atan_lambdas(self, triang):
        """
        TODO: Excepcion si x pasa por d. Si x pasa por d, entonces el logaritmo se indefine,
        pero realmente la expresion general sera 0*infinito.
        Solucion: Si x pasa por d, tomar un epsilon (0.01?) hacia arriba y hacia abajo, y asumir que la integral dentro del epsilon 
        no contiene ese termino (ie termino 4) 
        """
        def approx_val(x, d, tb):
            approx_val = np.arctan(d*np.sin(tb)/(x+d*np.cos(tb)))
            return approx_val
        def approx_der(x, d, tb):
            approx_der = 1/(x**2+2*d*x*np.cos(tb)+d**2)*(x*d*np.cos(tb)+d**2)
            return approx_der
        d = self.params.d
        avg = triang.avg_ang
        lambda_list = []
        lambda_list.append(lambda x, t: (x**2/2)*(approx_val(x,d,avg))*t+(x**2/2)*approx_der(x,d,avg)*(t**2/2-avg*t))
        lambda_list.append(lambda x, t: -(d**2*(np.sin(2*t)/2*(approx_val(x,d,avg)-approx_der(x,d,avg)*avg)+approx_der(x,d,avg)*(t*np.sin(2*t)/2+np.cos(2*t)/4))))
        lambda_list.append(lambda x, t: -1/2*d*np.cos(t)*x)
        lambda_list.append(lambda x, t: -1/4*d**2*(np.cos(2*t)/2*np.log(x**2+d**2+2*d*x*np.cos(t))-d*x*(x**4+d**4)/(4*d**3*x**3)*np.log(x**2+d**2+2*d*x*np.cos(t))))
        lambda_list.append(lambda x, t: 1/4*d**3*x*(-np.cos(2*t)/(4*d*x)+np.cos(t)*(x**2+d**2)/(2*d**2*x**2)))
        return lambda_list
class MISOBaseIntegrator:
    def __init__(self, lb, ub, consts, parameters):
        self.lb = lb
        self.ub = ub
        self.consts = consts
        self.params = parameters
    def acos_integrator(self, triang):
        theta_bot, theta_top = triang.ang_low, triang.ang_high
        sum_int_a, sum_int_b, sum_int_c = self.integral_lambdas()
        sum_int = sum_int_a(self.ub)-sum_int_a(self.lb)+sum_int_b(self.ub)-sum_int_b(self.lb)+sum_int_c(self.ub)-sum_int_c(self.lb)
        return (theta_top-theta_bot)*sum_int
    def integral_lambdas_linear(self):
        a = self.consts["a"]
        b = self.consts["b"]
        cosfov = self.params.cosfov
        sinbeta = self.params.sinbeta
        sum_int_a = lambda x: a*(0.5*cosfov/sinbeta*(x*np.sqrt(x**2+self.params.b**2)+self.params.b**2*np.ln(np.sqrt(x**2+self.params.b**2)+x)))
        sum_int_b = lambda x: self.params.a/sinbeta*X
        sum_int_c = lambda x: b*x**2
        return sum_int_a, sum_int_b, sum_int_c
    def angle_integrator(self, triang):
        theta_bot, theta_top = triang.ang_low, triang.ang_high
        return (theta_top**2-theta_bot**2)/2*(self.ub**2-self.lb**2)/2
    def pi_const_integrator(self, triang):
        theta_bot, theta_top = triang.ang_low, triang.ang_high
        return (theta_top-theta_bot)*(self.ub**2-self.lb**2)*np.pi/2


class NonOriginIntegrator:
    def __init__(self, lb, ub, consts, parameters):
        self.lb = lb
        self.ub = ub
        self.consts = consts
        self.params = parameters
    def integral_even(self, Dn, theta_bot, theta_top, theta_crt, n):
        max_Dn = Dn/np.cos(theta_top)
        if self.ub > max_Dn:
            return 1/(Dn**2*np.tan(theta_cr))*self.integral_cosine(Dn, theta_bot, theta_top, n), theta_top
        elif self.ub < Dn:
            return 1/(Dn**2*np.tan(theta_cr))*self.integral_linear(Dn, theta_bot, theta_top, n), None
        else:
            theta_mid = np.arccos(Dn/self.ub)
            if theta_mid < theta_top:
                theta_mid = theta_top
            cosine_int = 1/(Dn**2*np.tan(theta_cr))*self.integral_cosine(Dn, theta_bot, theta_mid, n)
            linear_int = 1/(Dn**2*np.tan(theta_cr))*self.integral_linear(Dn, theta_mid, theta_top, n)
            return cosine_int+linear_int, theta_mid
    def integral_odd(self, Dn, theta_bot, theta_top, theta_crt, n):
        max_Dn = Dn/np.sin(theta_bot)
        if self.ub > max_Dn:
            return 1/(Dn**2*np.cotan(theta_cr))*self.integral_sine(Dn, theta_bot, theta_top, n), theta_top
        elif self.ub < Dn:
            return 1/(Dn**2*np.cotan(theta_cr))*self.integral_linear(Dn, theta_bot, theta_top, n), None
        else:
            theta_mid = np.arcsin(Dn/self.ub)
            if theta_mid > theta_top:
                theta_mid = theta_top
            sine_int = 1/(Dn**2*np.cotan(theta_cr))*self.integral_sine(Dn, theta_mid, theta_top, n)
            linear_int = 1/(Dn**2*np.cotan(theta_cr))*self.integral_linear(Dn, theta_bot, theta_mid, n)
            return sine_int+linear_int, theta_mid

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
    def integral_lambdas_cosine(self, Dn, debug = False):
        """
            All the functions relation to solving the integral of the function Gn(Dn/cos(x))
            Gn(x) = X(X^2+b^2)^0.5 + ln((X^2+b^2)^0.5+X) - X + X^2
            Tested: Yes
        """
        if not debug:
            a = self.consts["a"]
            b = self.consts["b"]
            cosfov = self.params.cosfov
            sinbeta = self.params.sinbeta
            sigma = Dn/self.params.b
            ap = self.params.a 
        else:
            cosfov = 1
            sinbeta = 1
            a = 1
            b = 1
            sigma = Dn
            ap = 1
        sum_int_one = lambda t: (cosfov/sinbeta*a*Dn*b*(np.tan(t)*np.sqrt(np.cos(2*t)+1+2*sigma**2)/np.sqrt(2)+ 
        np.sqrt(1+sigma**2)*sp.special.ellipkinc(t, (1/(1+sigma**2)))-np.sqrt(1+sigma**2)*sp.special.ellipeinc(t, (1/(1+sigma**2)))))
        sum_int_two = lambda tb, tt: cosfov/sinbeta*a*self.sinh_int_cos(x, tb, tt, Dn)
        sum_int_three = lambda t: ap/sinbeta*Dn*self.log_cos(t)
        sum_int_four = lambda t: b*Dn/2*np.tan(t)
        return sum_int_one, sum_int_two, sum_int_three, sum_int_four
    def integral_lambdas_sine(self, Dn, debug = False):
        """
            All the functions relation to solving the integral of the function Gn(Dn/sin(x))
            Gn(x) = X(X^2+b^2)^0.5 + ln((X^2+b^2)^0.5+X) - X + X^2
            Tested: Yes
        """
        if not debug:
            a = self.consts["a"]
            b = self.consts["b"]
            cosfov = self.params.cosfov
            sinbeta = self.params.sinbeta
            sigma = Dn/self.params.b
            ap = self.params.a 
        else:
            cosfov = 1
            sinbeta = 1
            a = 1
            b = 1
            sigma = Dn
            ap = 1
        sum_int_one = lambda t: (cosfov/sinbeta*a*Dn*b*(-cotan(t)*np.sqrt(-np.cos(2*t)+1+2*sigma**2)/np.sqrt(2)+
        (sigma+1/sigma)*sp.special.ellipkinc(t, -1/sigma**2)-sigma*sp.special.ellipeinc(t, -1/sigma**2)))
        sum_int_two = lambda tb, tt: cosfov/sinbeta*a*self.sinh_int_cos(x, tb, tt, Dn)
        sum_int_three = lambda t: ap/sinbeta*Dn*np.log(np.sin(t)/np.cos(t))
        sum_int_four = lambda t: -b*Dn/2*cotan(t)
        return sum_int_one, sum_int_two, sum_int_three, sum_int_four
    def sinh_int_cos(self, tb, tt, Dn, max_order = 10, debug = True):
        """
        Tested: Yes, small error if Dn/b < 1, but nothing to do ...
        integral of arcsinh(Dn/(bsin(x))) 
        """
        b = self.params.b if not debug else 1
        summation = 0
        multiplier = b**2
        if Dn < b:
            tm = np.arccos(Dn/b)
            if tt < tm:
                tm = tt
            if tm > tb:
                tm = tb
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
        Tested: Yes, small error if Dn/b < 1, but nothing to do ...
        integral of arcsinh(Dn/(bsin(x))) 
        """
        b = self.params.b if not debug else 1
        summation = 0
        multiplier = b**2
        if Dn < b:
            tm = np.arccos(Dn/b)
            if tt < tm:
                tm = tt
            if tm < tb:
                tm = tb
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

            summation += np.log(2*Dn/b)*(tm-tb)-(self.f_logsin(tm)-self.f_logsin(tb))
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

class OriginIntegrator:
    def __init__(self, radius):
        self.radius = radius
    def integral_even(self, Dn, theta_bot, theta_top, theta_cr, n):
        L1 = self.radius
        max_Dn = Dn/np.cos(theta_cr)
        if L1 > max_Dn:
            return 1, None
        elif L1 < Dn:
            return (theta_top-theta_bot)*L1**2/((Dn**2)*np.tan(theta_cr)), None
        else:
            theta_mid = np.arccos(Dn/L1)
            return 1/np.tan(theta_cr)*np.tan(theta_mid)-1/np.tan(theta_cr)*np.tan(theta_bot)+(theta_cr-theta_mid)*L1**2/((Dn**2)*np.tan(theta_cr)), theta_mid
    def integral_odd(self, Dn, theta_bot, theta_top, theta_cr, n):
        L1 = self.radius
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
                    integr, crt = integrator.integral_even(Dn, tri.ang_low, tri.ang_high, tri.ang_crt, i)
                else:
                    Dn = tri.y if i%4 == 1 else tri.x
                    integr, crt = integrator.integral_odd(Dn, tri.ang_low, tri.ang_high, tri.ang_crt, i)
                int_sum += integr*tri.get_area()/(X*Y)
                if crt is not None:
                    tri.change_ang(crt)

            return int_sum
        return _integrator_wrapper   
    @_integrate
    def unitary_integrator(self, L1):
        integrator = OriginIntegrator(L1)
        return integrator
    @_integrate
    def non_origin_integrator(self, L1, L2, consts, parameters):
        integrator = NonOriginIntegrator(L1, L2, consts, parameters)
        return integrator
class TriangleIntegrator:
    def __init__(self, rect):
        self.rect = rect
    
    def _integrate(func):
        def _integrator_wrapper(self, *args):
            integrator = func(self, *args)
            area = self.triang.get_area()
            
    def wrapper_integrator(self, list_of_lims, parameters):
        X, Y = self.rect.X, self.rect.Y
        area = X*Y
        tot_int = 0
        for triang, lims in list_of_lims:
            for lim in lims:
                interv_bot, interv_top = lim
                integral = interv_top.integrate_ub(triang, parameters) - interv_bot.integrate_lb(triang, parameters)
                tot_int += integral * triang.get_area()
        return tot_int

        
            



if __name__=="__main__":
    X, Y, xc, yc = 5, 3, 3, 3
    rec = UniformRectangle(X, Y, xc, yc, 180)
    d = 1
    b = 1.2
    mocker = ParamMocker(d,b)
    lb = 1.1
    ub = 1.5
    intrec = MISOOffsetIntegrator(lb, ub, None, mocker)
    t = 1
    for i in range(10):
        fcos = intrec.f_cos(t, i)-intrec.f_cos(0, i)
        print(fcos)
    #one, two = intrec.acos_lambda_under_b()
    #print(two(1.5, 1.5,1))