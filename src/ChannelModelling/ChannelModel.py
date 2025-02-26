import numpy as np
import model as chb
import ChannelSimulation.SimulationInput as sim
from scipy.special import gamma, riccati_jn, riccati_yn
import abc
class WHP:
    def __init__(self, whp_tuple):
        self.whp_tuple = whp_tuple
    def add_tuple(self, w, h, p):
        self.whp_tuple.append({"w":w,"h":h,"p":p})
    def check_valid(self):
        return sum([whp["p"] for whp in self.whp_tuple])==1
    def to_bind(self, ymin, ymax, ymax_real, xmax):
        whp_bind = chb.WH_Probabilities()
        for tup in self.whp_tuple:
            whp_bind.push_probability(tup["w"], tup["h"], tup["p"], ymin, ymax, ymax_real, xmax+tup["w"]/2, -tup["w"]/2)
        return whp_bind
class Shadowing_Parameters:
    def __init__(self, positions, whp):
        self.positions = positions
        self.whp = whp
    def to_bind(self):
        shp = chb.Shadow_Coll()
        for position in self.positions:
            shp.create_collection(position["xt"],position["yt"],position["zi"],position["xr"],position["yr"],position["zj"],self.whp)
        return shp
class WallParameters:
    def __init__(self, reflection_coef):
        self.pw = reflection_coef
    def to_bind(self):
        return chb.WallParameters(self.pw)

class TransmitterParameters:
    def __init__(self, coordinate, beta, alpha, ang_rad):
        self.beta = beta
        self.alpha = alpha
        self.ang_rad = ang_rad
        self.m = -np.log(2) / (np.log(abs(np.cos(self.ang_rad * np.pi / 180))))
        self.coordinate = coordinate

    def to_bind(self):
        chb.TransmitterParameters(np.array(self.coordinate), self.beta, self.alpha, self.ang_rad,
                                  self.m).getCoordinate()

        return chb.TransmitterParameters(np.array(self.coordinate), self.beta, self.alpha, self.ang_rad, self.m)


class TransmitterAggregate:
    def __init__(self):
        self.transmitters = []

    def add_transmitter(self, transmitter: TransmitterParameters):
        self.transmitters.append(transmitter)


class ReceiverParameters:
    def __init__(self, Ap, eta, fov, angle, ele, center):
        self.Ap = Ap
        self.eta = eta
        self.fov = fov
        self.r = 0.05
        self.angle = angle
        self.ele = 90 - ele
        self.center = center
        self.calculate_coord()

    def calculate_coord(self, rot=0):
        self.angle += rot
        angle = self.angle
        self.angle = self.angle%360
        self.coordinate = [self.center[0] + self.r * np.cos(np.deg2rad(angle)) * np.sin(np.deg2rad(self.ele)),
                           self.center[1] + self.r * np.sin(np.deg2rad(angle)) * np.sin(np.deg2rad(self.ele)),
                           self.center[2] + self.r * np.cos(np.deg2rad(self.ele)),
                           angle,
                           90 - self.ele]
        

    def to_bind(self):
        return chb.ReceiverParameters(self.Ap, self.eta, self.fov, self.angle, self.ele, np.array(self.center),
                                      np.array(self.coordinate))


class ReceiverAggregate:
    def __init__(self):
        self.pd = []

    def add_pd(self, pd: ReceiverParameters):
        self.pd.append(pd)

    def rotate(self, rot):
        for pd in self.pd:
            pd.calculate_coord(rot)


class TunnelParameters:
    def __init__(self, x, y, z):
        self.X = x
        self.Y = y
        self.Z = z

    def to_bind(self):
        return chb.TunnelParameters(self.X, self.Y, self.Z)


class SimulationParameters:
    def __init__(self, t, c, sampling_time, t_rise, t_fall):
        self.t = t
        self.c = c
        self.Sampling_Time = sampling_time
        self.time = np.arange(0, 35e-9, self.Sampling_Time)
        self.time = np.round(self.time, decimals=12)
        self.t_rise = t_rise
        self.t_fall = t_fall
        self.h_led = 10 * (np.exp(-self.time / self.t_fall) - np.exp(-self.time / self.t_rise))
    def to_bind(self, scattering):
        sim_chb = chb.SimulationParameters(self.t, self.c, np.array(self.time), np.array(self.h_led), scattering)
        return sim_chb

class ScatteringParameters:
    def __init__(self, p, g, f, gamma, n, kr, km, ka, radius):
        self.p = p
        self.g = g
        self.f = f
        self.gamma = gamma
        self.n = n
        self.kr, self.km, self.ka = kr, km, ka
        self.radius = radius
    def to_bind(self):
        return chb.ScatteringParameters(self.n, self.radius, self.f, self.g, self.gamma, self.kr, self.km, self.p)

class CoefficientCalculator:
    def __init__(self, atmospheric_parameters, **kwargs):
        if "ignore_kr" in kwargs and kwargs["ignore_kr"]:
            self.kr = None
        else:
            self.kr = self.calculate_kr(atmospheric_parameters)
        self.km, self.ka = self.calculate_km_ka(atmospheric_parameters)
    def __call__(self):
        return self.kr, self.km, self.ka
    def calculate_kr(self, parameters):
        wavelength = parameters.wavelength
        wavelength_cm = wavelength*10**3
        density_air = parameters.density_rayleigh
        density_scenario = parameters.density_rayleigh
        depolarization_factor = parameters.depolarization ## 0.031 for Air
        ns = (591817/(238.0185-wavelength**-2)+167909/(57.362-wavelength**-2))/10**-8+1
        return density_scenario*24*np.pi**3/(wavelength_cm**4*density_air**2)*(ns**2-1)**2/(ns**2+2)**2*(6+3*depolarization_factor)/(6-7*depolarization_factor)
    def calculate_km_ka(self, parameters):
        sum_value = 0
        sum_value_ka = 0
        m = parameters.refraction
        x = 2 * np.pi * parameters.mean_radius / parameters.wavelength
        iterations = int(x + 4.05 * x ** 0.33333 + 2.0) + 1
        an, bn = self.calculate_anbn_iterative(m, x, iterations)
        for i in range(1, iterations):
            sum_value += (2 * i + 1) * (abs(an[i - 1]) ** 2 + abs(bn[i - 1]) ** 2)
            sum_value_ka += (2 * i + 1) * (an[i - 1] + bn[i - 1]).real

        km_value = parameters.N_Mie * parameters.wavelength ** 2 / (2*np.pi) * sum_value
        ka_value = parameters.N_Mie * parameters.wavelength ** 2 / (2*np.pi) * sum_value_ka
        return km_value, ka_value

    def calculate_bessel(self, x, iterations):
        bo = np.sin(x)
        b1 = np.sin(x) / x - np.cos(x)
        bessels = np.zeros(iterations, dtype=np.complex128)
        bessels[0] = bo
        bessels[1] = b1
        for i in range(2, iterations):
            bessels[i] = (2 * (i - 1) + 1) / x * bessels[i - 1] - bessels[i - 2]
        return bessels

    def calculate_hankel(self, x, iterations):
        ho = complex(np.sin(x), np.cos(x))
        h1 = complex(np.sin(x) / x - np.cos(x), np.cos(x) / x + np.sin(x))
        hankels = np.zeros(iterations, dtype=np.complex128)
        hankels[0] = ho
        hankels[1] = h1

        for i in range(2, iterations):
            hankels[i] = (2 * (i - 1) + 1) / x * hankels[i - 1] - hankels[i - 2]

        return hankels

    def calculate_Dn(self, x, iterations):
        dn = np.zeros(iterations, dtype=np.complex128)
        n = iterations - 1
        alpha = 1
        beta = 1
        a1 = (-1) ** 2 * (2 * n + 2 * 1 - 1) / x
        a2 = (-1) ** 3 * (2 * n + 2 * 2 - 1) / x
        a = a2
        inverse = 2 / x
        alpha_prev = a1 * (a2 + 1 / a1)
        alpha = alpha_prev * alpha
        beta_prev = a2
        beta = beta_prev * beta
        i = 3
        while np.abs(np.abs(alpha_prev / beta_prev) - 1) > 1e-12:
            beta_prev = a + 1 / beta_prev
            alpha_prev = a + 1 / alpha_prev
            a = -a + (-1) ** (i - 1) * inverse
            alpha = alpha_prev * alpha
            beta = beta_prev * beta
            i += 1
        dn[-1] = -n / x + alpha / beta
        for i in reversed(range(iterations - 1)):
            dn[i] = (i + 1) / x - 1 / (dn[i + 1] + (i + 1) / x)

        return dn

    def calculate_anbn_iterative(self, m, x, iterations):
        dn_mx = self.calculate_Dn(m * x, iterations)
        bessels = self.calculate_bessel(x, iterations)
        hankels = self.calculate_hankel(x, iterations)
        n = np.arange(1, iterations)
        an = ((dn_mx[1:] / m + n / x) * bessels[1:] - bessels[:-1]) / (
                    (dn_mx[1:] / m + n / x) * hankels[1:] - hankels[:-1])
        bn = ((m * dn_mx[1:] + n / x) * bessels[1:] - bessels[:-1]) / (
                    (m * dn_mx[1:] + n / x) * hankels[1:] - hankels[:-1])
        return an, bn
class AtmosphericParameters:
    def __init__(self, parameter_dictionary):
        self.wavelength = parameter_dictionary['wavelength']
        self.depolarization = parameter_dictionary['depolarization']
        self.density_rayleigh = parameter_dictionary['density_rayleigh']
        if "N_Mie" in parameter_dictionary.keys():
            self.N_Mie = parameter_dictionary["N_Mie"]
        else:
            self.N_Mie = self.calculate_N_Mie(parameter_dictionary)
        if "mean_radius" in parameter_dictionary.keys():
            self.mean_radius = parameter_dictionary["mean_radius"]
        else:
            self.mean_radius = self.calculate_mean_radius(parameter_dictionary)
        self.refraction = parameter_dictionary['refraction']
    def calculate_N_Mie(self, parameter_dictionary):
        W = parameter_dictionary['mass_per_unit']
        p = parameter_dictionary['dust_density']
        n = 3
        shape_param = parameter_dictionary['shape_mie']
        mean = parameter_dictionary['diameter_mean']
        Cv = np.pi/6
        integral = mean**n*gamma(n/shape_param+1)
        return W/(p*Cv*integral)
    def calculate_mean_radius(self, parameter_dictionary):
        W = parameter_dictionary['mass_per_unit']
        p = parameter_dictionary['dust_density']
        n = 1
        shape_param = parameter_dictionary['shape_mie']
        mean = parameter_dictionary['diameter_mean']
        Cv = np.pi/6
        integral = mean**n*gamma(n/shape_param+1)
        return integral

def calculate_kr_km_ka(N_Mie, wavelength, mean_radius, refraction):
    parameter_dictionary = {}
    parameter_dictionary['wavelength'] = wavelength
    parameter_dictionary['N_Mie'] = N_Mie
    parameter_dictionary["mean_radius"] = mean_radius
    parameter_dictionary['refraction'] = refraction
    parameter_dictionary['depolarization'] = 0.031
    parameter_dictionary['density_rayleigh'] = 2.547e19
    atmo_params = AtmosphericParameters(parameter_dictionary)
    coefs = CoefficientCalculator(atmo_params, **{"ignore_kr": False})
    kr, km, ka = coefs()
    return kr, km, ka


