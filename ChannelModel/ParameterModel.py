import numpy as np


class WallParameters:
    def __init__(self, reflection_coef):
        self.pw = reflection_coef


class TransmitterParameters:
    def __init__(self, coordinate, beta, alpha, ang_rad, m):
        self.beta = beta
        self.alpha = alpha
        self.ang_rad = ang_rad
        self.m = -np.log(2)/(np.log(abs(np.cos(self.ang_rad*np.pi/180))))
        self.coordinate = coordinate


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
        self.r = 0.1
        self.angle = angle
        self.ele = ele
        self.center = center
        self.coordinate = self.calculate_coord()

    def calculate_coord(self, rot = 0):
        angle = self.angle+rot

        return [self.center[0]+self.r*np.cos(np.deg2rad(angle))*np.sin(self.ele),
                self.center[1]+self.r*np.sin(np.deg2rad(angle))*np.sin(self.ele),
                self.center[2]+self.r*np.cos(self.ele),
                angle,
                90-np.rad2deg(self.ele)]


class ReceiverAggregate:
    def __init__(self):
        self.pd = []

    def add_led(self, pd: ReceiverParameters):
        self.leds.append(pd)

    def rotate(self, rot):
        for pd in self.pd:
            pd.calculate_coord(rot)


class TunnelParameters:
    def __init__(self, x, y, z):
        self.X = x
        self.Y = y
        self.Z = z

class SimulationParameters:
    def __init__(self, t, c, sampling_time, t_rise, t_fall):
        self.t = t
        self.c = c
        self.Sampling_Time = sampling_time
        self.time = np.linspace(0, 35e-9, retstep=self.Sampling_Time)
        self.time = np.round(self.time, decimals=12)
        self.t_rise = t_rise
        self.t_fall = t_fall
        self.h_led = 10*(np.exp(-self.time/self.t_fall)-np.exp(-self.time/self.t_rise))
