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

class OriginIntegrator:
    def __init__(self):
        ...
    def unitary_integral_even(self, L1, Dn, theta_bot, theta_top, n):
        theta_cr = theta_top
        max_Dn = Dn/np.cos(theta_cr)
        if L1 > max_Dn:
            return 1, None
        elif L1 < Dn:
            return (theta_top-theta_bot)*L1**2/((Dn**2)*np.tan(theta_cr)), None
        else:
            theta_mid = np.arccos(Dn/L1)
            return 1/np.tan(theta_cr)*np.tan(theta_mid)-1/np.tan(theta_cr)*np.tan(theta_bot)+(theta_cr-theta_mid)*L1**2/((Dn**2)*np.tan(theta_cr)), theta_mid
    def unitary_integral_odd(self, L1, Dn, theta_bot, theta_top, n):
        theta_cr = theta_bot
        max_Dn = Dn/np.sin(theta_cr)
        if L1 > max_Dn:
            return 1, None
        elif L1 < Dn:
            return (theta_top-theta_bot)*L1**2/((Dn**2)*cotan(theta_cr)), None
        else:
            theta_mid = np.arcsin(Dn/L1)
            return 1/cotan(theta_cr)*cotan(theta_mid)-1/cotan(theta_cr)*cotan(theta_top)+(theta_mid-theta_cr)*L1**2/((Dn**2)*cotan(theta_cr)), theta_mid



if __name__=="__main__":
    X, Y, xc, yc = 5, 3, 1, 1
    rec = EightRectangle(X, Y, xc, yc)
    integrator = OriginIntegrator()
    int_sum = 0
    for i, tri in enumerate(rec):
        if i % 2 == 0:
            Dn = tri.x if i%4 == 0 else tri.y
            int_sum += integrator.unitary_integral_even(1, Dn, 0, tri.ang_crt, i)[0]*tri.get_area()/(X*Y)
            print(integrator.unitary_integral_even(1, Dn, 0, tri.ang_crt, i)[0]*2, tri.ang_crt/tri.get_area(), tri.ang_crt)
        else:
            Dn = tri.y if i%4 == 1 else tri.x
            int_sum += integrator.unitary_integral_odd(1, Dn, tri.ang_crt, np.pi/2, i)[0]*tri.get_area()/(X*Y)
    print(int_sum)
    print(np.pi/(15))
