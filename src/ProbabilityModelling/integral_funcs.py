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
    def __init__(self, d, b, a, cosfov, sinbeta):
        self.d = d
        self.b = b
        self.a = a
        self.cosfov = cosfov
        self.sinbeta = sinbeta
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
        x_u, x_d = self.solve_b_u_threshold(avg_t, d, b)
        """
        Se genera un intervalo entre x_u y x_d, el cual se intersecta con [xb, xt]
        """
        integr = 0
        if x_u is None or xb >= x_u or xt <= x_d: # No hay interseccion, como es creciente entonces estamos por sobre b
            integr += lambda_over_one_wrap(xt, xb, tt, tb)            
            integr += lambda_over_two_wrap(xt, xb, tt, tb)
        elif xt < x_u and xb >= x_d: # Hay full interseccion, como es creciente entonces estamos por bajo b
            print("2")
            integr += lambda_under_one_wrap(xt, xb, tt, tb)
            integr += lambda_under_two_wrap(xt, xb, tt, tb)
            print(integr, lambda_under_one_wrap(xt, xb, tt, tb), lambda_under_two_wrap(xt, xb, tt, tb))
        elif xt > x_u and xb < x_u: # [xb, x_u] bajo b, [x_u, xt] sobre b
            print("3")
            integr += lambda_under_one_wrap(x_u, xb, tt, tb)
            integr += lambda_under_two_wrap(x_u, xb, tt, tb)
            integr += lambda_over_one_wrap(xt, x_u, tt, tb)
            integr += lambda_over_two_wrap(xt, x_u, tt, tb)
        elif xt >= x_d and xb < x_d: # [xb, x_d] sobre b, [x_d, xt] bajo b
            print("4")
            integr += lambda_over_one_wrap(xt, x_d, tt, tb)
            integr += lambda_over_two_wrap(xt, x_d, tt, tb)
            integr += lambda_under_one_wrap(x_d, xb, tt, tb)
            integr += lambda_under_two_wrap(x_d, xb, tt, tb)
        #print(integr)
        lambda_one_a_constant, _ = self.acos_lambda_under_b(use_discerner=False)
        cosfov = self.params.cosfov
        sinbeta = self.params.sinbeta
        integr *= cosfov/sinbeta
        pre_a = integr*self.consts["a"]
        aux_integr = lambda_one_a_constant(xt, tt)-lambda_one_a_constant(xt,tb)-lambda_one_a_constant(xb, tt)+lambda_one_a_constant(xb, tb)
        integr += -a*(aux_integr)/(b*sinbeta)
        integr *= self.consts["a"]
        integr += self.consts["b"]/2*(xt**2-xb**2)*(tt-tb)
        return integr
    def _acos_lambda_wrapper(self, N = 10):
        lambda_over_one, lambda_over_two = self.acos_lambda_over_b(N) # f(x,t), f(xt, xb, t)
        lambda_under_one, lambda_under_two = self.acos_lambda_under_b() # f(x,t), f(x, tt, tb) 
        lambda_over_one_wrap = lambda xt, xb, tt, tb: lambda_over_one(xt, tt)-lambda_over_one(xt,tb)-lambda_over_one(xb, tt)+lambda_over_one(xb, tb)
        lambda_over_two_wrap = lambda xt, xb, tt, tb: lambda_over_two(xt, xb, tt)-lambda_over_two(xt, xb, tb)
        lambda_under_one_wrap = lambda xt, xb, tt, tb: lambda_under_one(xt, tt)-lambda_under_one(xt,tb)-lambda_under_one(xb, tt)+lambda_under_one(xb, tb)
        lambda_under_two_wrap = lambda xt, xb, tt, tb: lambda_under_two(xt, xb, tt, tb)
        return lambda_over_one_wrap, lambda_over_two_wrap, lambda_under_one_wrap, lambda_under_two_wrap
    def solve_b_u_threshold(self, avg_t, d, b):
        ae = 1
        be = 2*d*np.cos(avg_t)
        ce = d**2-b**2
        if be**2-4*ae*ce >= 0:
            sol_one = (-b+np.sqrt(be**2-4*ae*ce))/(2*ce)
            sol_two = (-b-np.sqrt(be**2-4*ae*ce))/(2*ce)
            return max(sol_one, sol_two), min(sol_two, sol_one)
        else:
            """
            Para que lo de arriba no ocurra, c > 0, lo que implica que siempre es mayor a b**2
            """
            return None, None

    def acos_lambda_over_b(self, N = 10):
        """

        """
        b = self.params.b
        lambda_one = lambda x, t: 1/2*x**2*t
        lambda_two = lambda xt, xb, t: b**2/2*self.arctan_acos_integral(xt,xb,t, N)
        return lambda_one, lambda_two
    def acos_lambda_under_b(self, use_discerner = True):
        """
        Tested: Yes
        """
        def sign_discerner(x,t,d,use_discerner):
            multiplier = d*np.sin(t)
            eta = 4*d*x/(x+d)**2
            base = 2*(x+d)*sp.special.ellipeinc(t/2, eta)
            summ = np.arctanh((x+d*np.cos(t))/np.sqrt(x**2+2*d*x*np.cos(t)+d**2)) if abs(np.sin(t)) > 0.0001 else 0
            summ *= multiplier
            summ += (x+d)*sp.special.ellipeinc(t/2, eta)-(x**2-d**2)/(x+d)*sp.special.ellipkinc(t/2, eta)
            summ *= -1
            summ += base
            return summ
        b = self.params.b
        d = self.params.d
        lambda_one = lambda x, t: b*sign_discerner(x,t,d, use_discerner) # <- Tested 
        lambda_two = lambda xt, xb, tt, tb: self.x_dcos_ineq_integral(xt, xb, tt, tb, d)/(2*b**2) # <- Tested
        return lambda_one, lambda_two
    def x_dcos_ineq_integral_over_theta(self, d):
        """
        Tested: The sum is tested, the rest is assumed to be right
        """
        def lambda_upper_aux(x, tt, tb, d):
            ## Tested, I think..
            def cos_expansion(x, d, t, N = 10):
                ## Integral of d^2sin^2dcos(x)/(L+dcos(x)) // :)
                ## 1/(L+dcos(x)) -> log|L+dcos(x)|
                ## Si L+dcos(x) < 0 -> log(-L-dcos(x)) => No me conviene pq no puedo factorizar para afuera
                ## Solucion: Dividir por -dcos(x) (El cual es negativo)
                if x+d*np.cos(t) < 0:
                    u = (np.sin(t/2)+np.cos(t/2))/(np.cos(t/2)-np.sin(t/2))
                    base = 0
                    base = 1/3*(np.sin(t)**3*np.log(-d*np.cos(t))-(-np.log(np.abs(u))+np.sin(t)+(np.sin(t)**3)/3)) 
                    summ = 0
                    for n in range(1,N):
                        summ += (x/d)**n*(-1)**(n-1)/n*(self.f_cos(t, 1-n)-self.f_cos(t, 3-n))
                    return summ + base
                else:
                    if x == 0: 
                        # En este caso, log(x) = 0, y d/x = inf,
                        # Sin embargo, la solucion real es log(dcos(x)) => usar simplemente el caso base de arriba
                        u = (np.sin(t/2)+np.cos(t/2))/(np.cos(t/2)-np.sin(t/2))
                        base = 1/3*(np.sin(t)**3*np.log(d*np.cos(t))+(np.log(u)-np.sin(t)-np.sin(t)**3/3)) 
                        return base 
                    if x > d*np.cos(t):
                        base = np.sin(t)**3/3*np.log(x)
                        summ = 0
                        for n in range(1,N):
                            summ += (d/x)**n*(-1)**(n)/(n*(n+2))*(np.cos(t)**(n+2)*np.sin(t)-self.f_cos(t, n+3))
                    else:
                        u = (np.sin(t/2)+np.cos(t/2))/(np.cos(t/2)-np.sin(t/2))
                        base = 1/3*(np.sin(t)**3*np.log(d*np.cos(t))-(-np.log(np.abs(u))+np.sin(t)+(np.sin(t)**3)/3)) 
                        summ = 0
                        for n in range(1,N):
                            summ += (x/d)**n*(-1)**(n-1)/n*(self.f_cos(t, 1-n)-self.f_cos(t, 3-n))

                return summ + base

            aux_lambda = lambda x, t: x**3/3*t + x**2/2*d*np.sin(t)+1/4*d**2*x*t-1/8*x*d**2*np.sin(2*t)-1/2*d**3*(cos_expansion(x, d, t))
            if x >= -d*np.cos(tt): # x+dcos(t) > 0
            ## This hasnt been tested but lets assume it is true cuz why not
                summ = aux_lambda(x, tt)
                if x >= -d*np.cos(tb):
                    summ -= aux_lambda(x, tb)
                else:
                    tm = np.arccos(-x/d)
                    if tb > tm:
                        tm = 2*np.pi-tm  
                    summ += aux_lambda(x, tb)
            else:
                summ = 0
                if x >= -d*np.cos(tb):
                    tm = np.arccos(-x/d)
                    summ -= aux_lambda(x, tb)
                    if tb > tm:
                        tm = 2*np.pi-tm
                    summ -= aux_lambda(x, tt)

                else:
                    summ -= aux_lambda(x, tt)-aux_lambda(x,tb)
            return summ
        
        lambda_upper = lambda x, tt, tb: lambda_upper_aux(x, tt, tb, d) 
        lambda_lower_aux = lambda x, t: (-1)**(t<np.pi)*(-x**2/2*d*np.cos(t)+x**4/(8*d)*np.log(np.abs(np.tan(t/2)))+d*x**3/3*np.log(np.abs(np.sin(t)))+x**2*d**2/(4*d)*(np.cos(t)+np.log(np.abs(np.tan(t/2)))))
        lambda_lower = lambda x, tt, tb: lambda_lower_aux(x,tt) - lambda_lower_aux(x,tb)
        return lambda_upper, lambda_lower
    def x_dcos_ineq_integral(self, xt, xb, tt, tb, d):
        """
        Tested: Kinda, lets assume it is all good cuz i cant bother
        """
        lambda_up, lambda_low = self.x_dcos_ineq_integral_over_theta(d)
        ## tt, tb in [0,pi] or [pi,2pi] so...
        avg = (tt + tb)/2
        summ = 0
        sign_of_sin = np.sin(avg) < 0
        if  np.abs(xb + d*np.cos(avg)) < np.abs(d*np.sin(avg)):
            if np.abs(xt + d*np.cos(avg)) < np.abs(d*np.sin(avg)):
                summ += lambda_low(xt, tt, tb)-lambda_low(xb, tt, tb)
            else:
                # [xb, xm] -> Below, [xm, xt] > Above
                summ += lambda_up(xt, tt, tb)-lambda_low(xb, tt,tb)

                if xb+d*np.cos(avg) > 0 or xt + d*np.cos(avg) > 0:
                    xm = d*((-1)**sign_of_sin*np.sin(avg)-np.cos(avg))
                else:                  
                    xm = d*(-(-1)**sign_of_sin*np.sin(avg)-np.cos(avg))
                summ -= lambda_up(xm, tt, tb)-lambda_low(xm, tt,tb)
        else:
            if np.abs(xt + d*np.cos(avg)) > np.abs(np.sin(avg)):
                summ += lambda_up(xt, tt, tb)-lambda_up(xb, tt, tb)
            else:
                summ += lambda_low(xt, tt, tb)-lambda_up(xb, tt,tb)
                # [xb, xm] -> Above, [xm, xt] > Below
                if xb+d*np.cos(avg) > 0:
                    xm = d*((-1)**sign_of_sin*np.sin(avg)-np.cos(avg))
                else:
                    xm = d*(-(-1)**sign_of_sin*np.sin(avg)-np.cos(avg))
                summ -= lambda_low(xm, tt, tb)-lambda_up(xm, tt,tb)
        return summ


    def arctan_acos_integral(self, xt, xb, t, N):
        """
        Integral of 
        log(xt^2+d^2+2dxtcos(t))-dcot(t)tan^-1(xt+dcos(t)/dsin(t))

        Tested: Si  
        """
        def arctan_acos_sum(xt, xb, t, d, N):
            xt_coef = (2*d*xt)/(xt**2+d**2)
            xb_coef = (2*d*xb)/(xb**2+d**2)
            summ = 0
            for n in range(1,N+1):
                summ += xt_coef**n*(-1)**(n+1)*self.f_cos(t, n)/n
                summ -= xb_coef**n*(-1)**(n+1)*self.f_cos(t, n)/n
            return summ
        def arctan_tanh_expr(xt, xb, t, d):
            eta = (xt*xb+d**2)/(d*(xt+xb))
            multiplier = d*(xt-xb)/(xt+xb)
            if np.abs((xt-xb)/(2*(xt+xb))*np.sin(t)/(eta+np.cos(t))) < 1:
                multiplier = d*(xt-xb)/(xt+xb)
                if eta**2 > 1:
                    summ = -2*eta/np.sqrt(eta**2-1)*np.arctan((eta-1)/np.sqrt(eta**2-1)*np.tan(t/2))    
                else:
                    u = (eta-1)/np.sqrt(1-eta**2)*np.tan(t/2)
                    summ = 2*eta/np.sqrt(1-eta**2)*np.arctanh(u) if np.abs(u) <= 1 else 2*eta/np.sqrt(1-eta**2)*np.arctanh(1/u) 
                summ += t
            else:
                multiplier = d
                summ = (np.pi/2-eta/((xt-xb)/(d*(xt+xb))))*np.log(np.abs(np.sin(t)))
                summ += eta*(t+np.cos(t)/np.sin(t))
            return summ*multiplier
        d = self.params.d
        log_int = 1/2*(arctan_acos_sum(xt, xb, t, d, N) + t*np.log((xt**2+d**2)/(xb**2+d**2)))
        #print(f"Arctan: {arctan_tanh_expr(xt, xb, t, d)}, {log_int}, {log_int - arctan_tanh_expr(xt, xb, t, d)}")
        return log_int - arctan_tanh_expr(xt, xb, t, d)
    def arctan_acos_integral_debug(self, xt, xb, tt, tb, N=10):
        """
        Integral of 
        log(xt^2+d^2+2dxtcos(t))-dcot(t)tan^-1(xt+dcos(t)/dsin(t))

        Tested: Si  
        """
        def arctan_acos_sum(xt, xb, t, d, N):
            xt_coef = (2*d*xt)/(xt**2+d**2)
            xb_coef = (2*d*xb)/(xb**2+d**2)
            summ = 0
            for n in range(1,N+1):
                summ += xt_coef**n*(-1)**(n+1)*self.f_cos(t, n)/n
                summ -= xb_coef**n*(-1)**(n+1)*self.f_cos(t, n)/n
            return summ
        def arctan_tanh_expr(xt, xb, t, d):
            eta = (xt*xb+d**2)/(d*(xt+xb))
            multiplier = d*(xt-xb)/(xt+xb)
            if np.abs((xt-xb)/(2*(xt+xb))*np.sin(t)/(eta+np.cos(t))) < 1:
                multiplier = d*(xt-xb)/(xt+xb)
                if eta**2 > 1:
                    summ = -2*eta/np.sqrt(eta**2-1)*np.arctan((eta-1)/np.sqrt(eta**2-1)*np.tan(t/2))    
                else:
                    u = (eta-1)/np.sqrt(1-eta**2)*np.tan(t/2)
                    summ = 2*eta/np.sqrt(1-eta**2)*np.arctanh(u) if np.abs(u) <= 1 else 2*eta/np.sqrt(1-eta**2)*np.arctanh(1/u) 
                print(eta, multiplier)
                summ += t
            else:
                multiplier = d
                summ = (np.pi/2-eta/((xt-xb)/(d*(xt+xb))))*np.log(np.abs(np.sin(t)))
                summ += eta*(t+np.cos(t)/np.sin(t))
            return summ*multiplier
        d = self.params.d
        d = self.params.d
        
        log_int = 1/2*(arctan_acos_sum(xt, xb, tt, d, N))-1/2*(arctan_acos_sum(xt, xb, tb, d, N))
        print(f"Arctan: {arctan_tanh_expr(xt, xb, tt, d)-arctan_tanh_expr(xt, xb, tb, d)}, {log_int}")
  
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
        elif order < 0:
            return self.f_sec(t, order)
        for i in range(order//2 if order % 2 == 0 else order//2+1):
            if i==0:
                multiplier *= 1/(order)
            else:
                multiplier *= (order-2*i+1)/(order-2*i)
            summation += multiplier*np.cos(t)**(-2*i)
        summation *= sum_base
        if order % 2 == 0:
            summation += multiplier*t
        return summation    

    def f_sec(self, t, order):
        sum_base = np.tan(t)
        order = -order
        if order > 1:
            sum_base *= 1/np.cos(t)**(order-2)
        else:
            return np.log(np.abs((np.sin(t/2)+np.cos(t/2))/(np.cos(t/2)-np.sin(t/2))))
        summation = 0
        multiplier = 1
        for i in range((order)//2 if order % 2 == 0 else (order-1)//2):
            if i==0:
                multiplier *= 1/(order-1)
            else:
                multiplier *= (order-2*i)/(order-2*i-1)
            summation += multiplier*np.cos(t)**(2*i)
        pre_sum = summation
        summation *= sum_base
        if order % 2 == 1:
            summation += multiplier*np.log(np.abs((np.sin(t/2)+np.cos(t/2))/(np.cos(t/2)-np.sin(t/2))))
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
        return summ
    def atan_lambdas(self, triang):
        """
        Done!
        Tested: Si
        """
        def approx_val(x, d, tb):
            approx_val = np.arctan(d*np.sin(tb)/(x+d*np.cos(tb)))
            print(approx_val, tb)
            return approx_val
        def approx_der(x, d, tb):
            approx_der = 1/(x**2+2*d*x*np.cos(tb)+d**2)*(x*d*np.cos(tb)+d**2)
            return approx_der
        def fourth_lambda(x, d, t):
            "int d^2sin(x)cos(x)log(x^2+d^2+2dxcos(x)) dt"
            multiplier = d**2/2
            log_x_d = 0
            log_x_d = -(np.log((x**2+d**2))-1+np.log(2))*np.cos(2*t)/4
            expansion = -sum([((2*x*d)/(x**2+d**2))**n*(-1)**(n-1)*np.cos(t)**(n+2)/(n*(n+2)) for n in range(1,15)])
            return (log_x_d+expansion)*multiplier

        d = self.params.d
        avg = triang.avg_ang
        lambda_list = []
        lambda_list.append(lambda x, t: (x**2/2)*(approx_val(x,d,avg))*t+(x**2/2)*approx_der(x,d,avg)*(t**2/2-avg*t))
        lambda_list.append(lambda x, t: -(d**2/2*(np.sin(2*t)/2*(approx_val(x,d,avg)-approx_der(x,d,avg)*avg)+approx_der(x,d,avg)*(t*np.sin(2*t)/2+np.cos(2*t)/4))))
        lambda_list.append(lambda x, t: -1/2*d*np.cos(t)*x)
        lambda_list.append(lambda x, t: -fourth_lambda(x, d, t))
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
        sum_int_a = lambda x: a*(0.5*cosfov/sinbeta*(x*np.sqrt(x**2+self.params.b**2)+self.params.b**2*np.log(np.sqrt(x**2+self.params.b**2)+x)))
        sum_int_b = lambda x: -self.params.a/sinbeta*x
        sum_int_c = lambda x: b*x**2/2
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
            integral = 2/(Dn**2*np.tan(theta_crt))*self.integral_cosine(Dn, theta_bot, theta_top, n)*1/np.pi
            print(f"The cosine integral from {self.lb} to {self.ub}, with max radius {Dn}, angles: [{theta_bot},{theta_top}] amounts to {self.integral_cosine(Dn, theta_bot, theta_top, n)*2/(Dn**2*np.tan(theta_crt))*1/np.pi}", Dn**2*np.tan(theta_crt))
            return integral, theta_top
        elif self.ub < Dn:
            integral = 2/(Dn**2*np.tan(theta_crt))*self.integral_linear(Dn, theta_bot, theta_top, n)*1/np.pi
            print(f"The linear integral from {self.lb} to {self.ub}, with max radius {Dn}, angles: [{theta_bot},{theta_top}] amounts to {self.integral_linear(Dn, theta_bot, theta_top, n)}")
            return integral, None
        else:
            theta_mid = np.arccos(Dn/self.ub)
            if theta_mid < theta_top:
                theta_mid = theta_top
            cosine_int = 2/(Dn**2*np.tan(theta_crt))*self.integral_cosine(Dn, theta_bot, theta_mid, n)*1/np.pi
            linear_int = 2/(Dn**2*np.tan(theta_crt))*self.integral_linear(Dn, theta_mid, theta_top, n)*1/np.pi
            print(f"The linear integral from {self.lb} to {self.ub}, with max radius {Dn}, angles: [{theta_mid},{theta_top}] amounts to {self.integral_linear(Dn, theta_mid, theta_top, n)}")
            print(f"The cosine integral from {self.lb} to {self.ub}, with max radius {Dn}, angles: [{theta_bot},{theta_mid}] amounts to {self.integral_cosine(Dn, theta_bot, theta_mid, n)}")
            return cosine_int+linear_int, theta_mid
    def integral_odd(self, Dn, theta_bot, theta_top, theta_crt, n):
        max_Dn = Dn/np.sin(theta_bot)
        if self.ub > max_Dn:
            integral = 2/(Dn**2*cotan(theta_crt))*self.integral_sine(Dn, theta_bot, theta_top, n)*1/np.pi
            print(f"The sine integral from {self.lb} to {self.ub}, with max radius {Dn}, angles: [{theta_bot},{theta_top}] amounts to {self.integral_sine(Dn, theta_bot, theta_top, n)*2/(Dn**2*cotan(theta_crt))*1/np.pi}", Dn**2*cotan(theta_crt))
            return integral, theta_bot
        elif self.ub < Dn:
            integral = 2/(Dn**2*cotan(theta_crt))*self.integral_linear(Dn, theta_bot, theta_top, n)*1/np.pi
            print(f"The linear integral from {self.lb} to {self.ub}, with max radius {Dn}, angles: [{theta_bot},{theta_top}] amounts to {self.integral_linear(Dn, theta_bot, theta_top, n)}")
            return integral, None
        else:
            theta_mid = np.arcsin(Dn/self.ub)
            if theta_mid > theta_top:
                theta_mid = theta_top
            
            sine_int = 2/(Dn**2*cotan(theta_crt))*self.integral_sine(Dn, theta_mid, theta_top, n)*1/np.pi
            linear_int = 2/(Dn**2*cotan(theta_crt))*self.integral_linear(Dn, theta_bot, theta_mid, n)*1/np.pi
            print(f"The linear integral from {self.lb} to {self.ub}, with max radius {Dn}, angles: [{theta_bot},{theta_mid}] amounts to {self.integral_linear(Dn, theta_bot, theta_mid, n)}")
            print(f"The sine integral from {self.lb} to {self.ub}, with max radius {Dn}, angles: [{theta_mid},{theta_top}] amounts to {self.integral_sine(Dn, theta_mid, theta_top, n)}")

            return sine_int+linear_int, theta_mid

    def integral_cosine(self, Dn, theta_bot, theta_top, n):
        sum_int_a, sum_int_b, sum_int_c = self.integral_lambdas_linear()
        sum_int_one, sum_int_two, sum_int_three, sum_int_four = self.integral_lambdas_cosine(Dn)
        inner_sum = (-sum_int_a(self.lb)-sum_int_b(self.lb)-sum_int_c(self.lb))*(theta_top-theta_bot)
        outer_sum = sum_int_one(theta_top)-sum_int_one(theta_bot)+sum_int_two(theta_bot, theta_top)+sum_int_three(theta_top)-sum_int_three(theta_bot)+sum_int_four(theta_top)-sum_int_four(theta_bot)
        #print(f"The integral from {self.lb} to {self.ub}, with max radius {Dn}, max angle {theta_top} and min angle {theta_bot} amounts to {outer_sum+inner_sum}")
        #print(f"Innter sum {inner_sum}")
        #print(f"Outer sum {outer_sum}, Terms: \n - 1 {(sum_int_one(theta_top)-sum_int_one(theta_bot))+-sum_int_a(self.lb)*(theta_top-theta_bot)+sum_int_two(theta_bot, theta_top)} \n - 2 {sum_int_two(theta_bot, theta_top)} \n - 3 {sum_int_three(theta_top)-sum_int_three(theta_bot)} \n - 4 {sum_int_four(theta_top)-sum_int_four(theta_bot)}")
        return outer_sum+inner_sum
    def integral_sine(self, Dn, theta_bot, theta_top, n):

        sum_int_a, sum_int_b, sum_int_c = self.integral_lambdas_linear()
        sum_int_one, sum_int_two, sum_int_three, sum_int_four = self.integral_lambdas_sine(Dn)
        inner_sum = (-sum_int_a(self.lb)-sum_int_b(self.lb)-sum_int_c(self.lb))*(theta_top-theta_bot)
        outer_sum = sum_int_one(theta_top)-sum_int_one(theta_bot)+sum_int_two(theta_bot, theta_top)+sum_int_three(theta_top)-sum_int_three(theta_bot)+sum_int_four(theta_top)-sum_int_four(theta_bot)
        #print(f"The integral from {self.lb} to {self.ub}, with max radius {Dn}, max angle {theta_top} and min angle {theta_bot} amounts to {outer_sum+inner_sum}")
        #print(f"Innter sum {inner_sum}")
        #print(f"Outer sum {outer_sum}, Terms: \n - 1 {sum_int_one(theta_top)-sum_int_one(theta_bot)} \n - 2 {sum_int_two(theta_bot, theta_top)} \n - 3 {sum_int_three(theta_top)-sum_int_three(theta_bot)} \n - 4 {sum_int_four(theta_top)-sum_int_four(theta_bot)}")
        return outer_sum+inner_sum
    
    def integral_linear(self, Dn, theta_bot, theta_top, n):
        sum_int_a, sum_int_b, sum_int_c = self.integral_lambdas_linear()
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
            bp = self.params.b
        else:
            cosfov = 1
            sinbeta = 1
            a = 1
            b = 1
            sigma = Dn
            ap = 1
            bp = 1
        sum_int_one = lambda t: (cosfov/sinbeta*a*Dn*bp/2*(np.tan(t)*np.sqrt(np.cos(2*t)+1+2*sigma**2)/np.sqrt(2)+ 
        np.sqrt(1+sigma**2)*sp.special.ellipkinc(t, (1/(1+sigma**2)))-np.sqrt(1+sigma**2)*sp.special.ellipeinc(t, (1/(1+sigma**2)))))
        sum_int_two = lambda tb, tt: cosfov/sinbeta*a*self.sinh_int_cos(tb, tt, Dn)/2
        sum_int_three = lambda t: -ap/sinbeta*Dn*self.log_cos(t)*a
        sum_int_four = lambda t: b*Dn**2/2*np.tan(t)
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
            bp = self.params.b
        else:
            cosfov = 1
            sinbeta = 1
            a = 1
            b = 1
            sigma = Dn
            ap = 1
            bp = 1
        sum_int_one = lambda t: (cosfov/sinbeta*a*Dn*bp/2*(-cotan(t)*np.sqrt(-np.cos(2*t)+1+2*sigma**2)/np.sqrt(2)+
        (sigma+1/sigma)*sp.special.ellipkinc(t, -1/sigma**2)-sigma*sp.special.ellipeinc(t, -1/sigma**2)))
        sum_int_two = lambda tb, tt: cosfov/sinbeta*a*self.sinh_int_sin(tb, tt, Dn)/2
        sum_int_three = lambda t: -ap/sinbeta*Dn*np.log(np.sin(t/2)/np.cos(t/2))*a
        sum_int_four = lambda t: -b*Dn**2/2*cotan(t)
        return sum_int_one, sum_int_two, sum_int_three, sum_int_four
    def sinh_int_cos(self, tb, tt, Dn, max_order = 10, debug = False):
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
                summation += 1/(2*order+1)*factor*(self.f_sec(tm, order)-self.f_sec(tb, order))*(Dn/b)**(2*order+1)
            for order in range(1,max_order):
                factor = 1/sp.special.factorial(order)
                for j in range(order):
                    factor *= (-1/2-j)
                summation -= 1/(2*order)*factor*(self.f_cos(tt, order)-self.f_cos(tm, order))*(Dn/b)**(-2*order)

            summation += np.log(2*Dn/b)*(tt-tm)-(self.f_logcos(tt)-self.f_logcos(tm))
        else:
            summation += np.log(2*Dn/b)*(tt-tb)-(self.f_logcos(tt)-self.f_logcos(tb))
            for order in range(1,max_order):
                factor = 1/sp.special.factorial(order)
                for j in range(order):
                    factor *= (-1/2-j)
                summation -= 1/(2*order)*factor*(self.f_cos(tt, order)-self.f_cos(tb, order))*(Dn/b)**(-2*order)
        summation += np.log(b)*(tt-tb)
        #print(f"Integral Arcsinh {summation*multiplier}")

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
            tm = np.arcsin(Dn/b)
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
            partial_summation = 0
            for order in range(1,max_order):
                factor = 1/sp.special.factorial(order)
                for j in range(order):
                    factor *= (-1/2-j)
                summation -= 1/(2*order)*factor*(self.f_sin(tt, order)-self.f_sin(tb, order))*(Dn/b)**(-2*order)
                partial_summation -= 1/(2*order)*factor*(self.f_sin(tt, order)-self.f_sin(tb, order))*(Dn/b)**(-2*order)
            #print(f"Integral Arcsinh {partial_summation*multiplier}, Integral logsin {(self.f_logsin(tt)-self.f_logsin(tb))}")
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
        base = -t*np.log(2)
        summ = 0
        for i in range(1,15):
            summ += np.sin(2*i*t)/(2*i**2)*(-1)**(i-1)
        return base + summ
    def f_logsin(self, t):
        """
        Tested: Yes
        Integral of log(sin(x))
        """
        base = -t*np.log(2)
        summ = 0
        for i in range(1,15):
            summ -= np.sin(2*i*t)/(2*i**2)
        return base + summ 
    def log_cos(self, t):
        return np.log((np.sin(t/2)+np.cos(t/2))/(np.cos(t/2)-np.sin(t/2)))
    def integral_lambdas_linear(self):
        a = self.consts["a"]
        b = self.consts["b"]
        cosfov = self.params.cosfov
        sinbeta = self.params.sinbeta
        sum_int_a = lambda x: a*(0.5*cosfov/sinbeta*(x*np.sqrt(x**2+self.params.b**2)+self.params.b**2*np.log(np.sqrt(x**2+self.params.b**2)+x)))
        sum_int_b = lambda x: -self.params.a/sinbeta*x*a
        sum_int_c = lambda x: b*x**2/2
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
                print(int_sum, tri.ang_crt, i)
            if len(args) > 2:
                print(f"The [{args[0]},{args[1]}] total is {int_sum}")
            return int_sum
        return _integrator_wrapper   
    @_integrate
    def origin_integrator(self, L1):
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
        area = 0
        tot_int = 0
        for triang, lims in list_of_lims:
            for lim in lims:
                interv_bot, interv_top = lim
                integral = interv_top.integrate_ub(triang, parameters)/(2*pi) - interv_bot.integrate_lb(triang, parameters)/(2*pi)
                tot_int += integral * triang.get_area()
            area += triang.get_area()
        return tot_int/area
    def __call__(self, list_of_lims, parameters):
        return self.wrapper_integrator(list_of_lims, parameters)
    


        
            



if __name__=="__main__":
    X, Y, xc, yc = 5, 3, 3, 3
    rec = UniformRectangle(X, Y, xc, yc, 180)
    d = 1
    b = 1.44
    a = 0.65
    cosfov = 0.5
    sinbeta = 0.5
    mocker = ParamMocker(d,b,a,cosfov,sinbeta)
    lb = 1.1
    ub = 1.5
    intrec = MISOOffsetIntegrator(lb, ub, None, mocker)
    t = 1
    N = 10
    for triang in rec:
        if triang.avg_ang > 0.9 and triang.avg_ang < 1.1:
            print(intrec.acos_integrator(triang), triang.ang_low, triang.ang_high, lb, ub)

