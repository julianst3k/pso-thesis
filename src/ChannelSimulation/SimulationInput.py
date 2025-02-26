import numpy as np
import model as chb
import ChannelModelling.ChannelModel as cm
from ChannelModelling.ChannelModel import calculate_kr_km_ka
import ChannelSimulation.ShadowingSimulation as sm
from abc import ABC, abstractmethod
from yaml import load, dump, Loader


class Model(ABC):
    def do_receiver_aggregate(self):
        new_aggregate = chb.ReceiverAggregate()
        for receiver in self.Receiver_Agg.pd:
            receiver_to_bind = receiver.to_bind()
            new_aggregate.pushReceiver(receiver_to_bind)
            self.bind_recv.append(receiver_to_bind)
        return new_aggregate
    def do_transmitter_aggregate(self):
        new_aggregate = chb.TransmitterAggregate()
        for transmitter in self.Transmitter_Agg.transmitters:
            trans_to_bind = transmitter.to_bind()
            new_aggregate.pushTransmitter(trans_to_bind)
            self.bind_trans.append(trans_to_bind)
        return new_aggregate
    def create_pybind_objects(self):
        self.wp_bind = self.wall.to_bind()
        self.tag_bind = self.do_transmitter_aggregate()
        rag = self.do_receiver_aggregate()
        self.scatt_bind = self.scattering.to_bind()
        self.par_bind = self.sim_param.to_bind(self.scatt_bind)
        self.tunn_bind = self.tunnel.to_bind()

        return self.wp_bind, self.tag_bind, rag, self.par_bind, self.tunn_bind

    def create_pybind_rotations(self, num_rots):
        wp, tag, rag, par, tunn = self.create_pybind_objects()
        self.receiver_configurations = chb.ReceiverConfigurations()
        self.rags = []
        self.receiver_configurations.pushReceiver(rag)
        self.rags.append(rag)
        for i in range(1, num_rots):
            rot = 360 / (num_rots + 1)
            self.Receiver_Agg.rotate(rot)
            rag = self.do_receiver_aggregate()
            self.receiver_configurations.pushReceiver(rag)
            self.rags.append(rag)
        self.shadowing_parameters_bind = self.shadowing_parameters.whp_obj.to_bind(0, 3, 3, 5)

        return wp, tag, self.receiver_configurations, par, tunn





class SimulatedModel(Model):
    def __init__(self, config_params):
        self.parameter_factory = ParameterFactory(config_params)
        self.Transmitter_Agg = self.parameter_factory.create_transmitters_agg()
        self.wall = self.parameter_factory.create_walls()
        self.sim_param = self.parameter_factory.create_sim_params()
        self.scattering = self.parameter_factory.create_scattering_params()
        self.tunnel = self.parameter_factory.create_tunnel()
        self.shadowing_parameters = self.parameter_factory.create_shadowing_parameters()
        self.bind_trans = []
        self.bind_recv = []
        self.receiver_configurations = None
        self.rags = []
        self.set_receiver = False



    def generate_receiver(self, coords):
        self.Receiver_Agg = self.parameter_factory.create_receivers_agg(coords)
        self.set_receiver = True

    def create_pybind_rotations(self, num_rots):
        if self.set_receiver:
            super().create_pybind_rotations(num_rots)
    def do_receiver_aggregate(self):
        if self.set_receiver:
            new_aggregate = super().do_receiver_aggregate()
        return new_aggregate

class ParameterFactory:
    def __init__(self, config_params):
        fl = open(config_params, 'r')
        self.config = load(fl, Loader=Loader)
    def create_transmitters_agg(self):
        Transmitter_Agg = cm.TransmitterAggregate()
        alpha_t_outside = self.config["alpha_t_outside"]
        alpha_t_inside = self.config["alpha_t_inside"]
        beta_t = self.config["beta_t"]
        for i in range(4):
            if i == 0 or i == 3:
                alphat = alpha_t_outside * (-1) ** i + 90
            else:
                alphat = alpha_t_inside * (-1) ** i + 90
            Transmitter_Agg.add_transmitter(cm.TransmitterParameters([(i + 1), 0.5, 3], beta_t, alphat, 60))
        return Transmitter_Agg
    def create_receivers_agg(self, coords):
        Receiver_Agg = cm.ReceiverAggregate()
        Ap = self.config["Ap"]
        eta = self.config["eta"]
        fov = self.config["fov"]
        beta_r = self.config["beta_r"]
        for j in range(4):
            Receiver_Agg.add_pd(cm.ReceiverParameters(Ap, eta, fov, j * 90 + 30, beta_r, coords))
        return Receiver_Agg

    def create_walls(self):
        return cm.WallParameters(self.config["reflection_coef"])

    def create_scattering_params(self):
        p = self.config["p"]
        g = self.config["g"]
        f = self.config["f"]
        gamma = self.config["gamma"]
        n = self.config["n"]
        N_Mie = float(self.config["N_Mie"])
        wavelength = float(self.config["wavelength"])
        particle_size = float(self.config["particle_size"])
        kr, km, ka = calculate_kr_km_ka(N_Mie, wavelength, particle_size, complex(1.5, 0.0014))
        radius = self.config["radius"]
        return cm.ScatteringParameters(p, g, f, gamma, n, kr, km, ka, radius)

    def create_sim_params(self):
        t = float(self.config["t"])
        c = float(self.config["c"])
        t_rise = float(self.config["t_rise"])
        t_fall = float(self.config["t_fall"])
        sampling_time = float(self.config["sampling_time"])
        return cm.SimulationParameters(t,c,sampling_time,t_rise,t_fall)

    def create_tunnel(self):
        x = self.config["x"]
        y = self.config["y"]
        z = self.config["z"]
        return cm.TunnelParameters(x,y,z)

    def create_shadowing_parameters(self):
        h_min = self.config["h_min"]
        w_min = self.config["w_min"]
        h_max = self.config["h_max"]
        w_max = self.config["w_max"]
        h_mean = self.config["h_mean"]
        w_mean = self.config["w_mean"]
        std_h = self.config["std_h"]
        std_w = self.config["std_w"]
        wrap = sm.WHP_Wrapper()
        wrap.generate_norm_distribution(h_mean, h_min, h_max, w_mean, w_min, w_max, std_h, std_w)

        return wrap