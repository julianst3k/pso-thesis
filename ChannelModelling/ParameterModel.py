import numpy as np
import ChannelModelling.build.model as chb


class WallParameters:
    def __init__(self, reflection_coef):
        self.pw = reflection_coef


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


class ModelSeed:
    def __init__(self, betat, betar, fov, alphat_outside, alphat_inside, coord_rec):
        self.Transmitter_Agg = TransmitterAggregate()
        self.Receiver_Agg = ReceiverAggregate()
        for i in range(4):
            if i == 0 or i == 3:
                alphat = alphat_outside
            else:
                alphat = alphat_inside
            self.Transmitter_Agg.add_transmitter(TransmitterParameters([(i + 1), 0.5, 3], betat, 90, 60))
        for j in range(4):
            self.Receiver_Agg.add_pd(ReceiverParameters(0.0001, 1.5, fov, j * 90 + 30, betar, coord_rec))
        self.wall = WallParameters(0.6)
        self.sim_param = SimulationParameters(5e-9, 3e8, 0.25e-9, 0.5e-9, 1e-9)
        self.tunnel = TunnelParameters(6, 3, 3.5)
        self.bind_trans = []
        self.bind_recv = []
        self.receiver_configurations = None
        self.rags = []
        self.rotations = False

    def create_pybind_objects(self):
        self.wp_bind = chb.WallParameters(self.wall.pw)
        self.tag_bind = self.do_transmitter_aggregate()
        rag = self.do_receiver_aggregate()
        self.par_bind = chb.SimulationParameters(self.sim_param.t, self.sim_param.c, np.array(self.sim_param.time),
                                                 np.array(self.sim_param.h_led))
        self.tunn_bind = chb.TunnelParameters(self.tunnel.X, self.tunnel.Y, self.tunnel.Z)

        return self.wp_bind, self.tag_bind, rag, self.par_bind, self.tunn_bind

    def do_transmitter_aggregate(self):
        new_aggregate = chb.TransmitterAggregate()
        for transmitter in self.Transmitter_Agg.transmitters:
            trans_to_bind = transmitter.to_bind()
            new_aggregate.pushTransmitter(trans_to_bind)
            self.bind_trans.append(trans_to_bind)
        return new_aggregate

    def do_receiver_aggregate(self):
        new_aggregate = chb.ReceiverAggregate()
        for receiver in self.Receiver_Agg.pd:
            receiver_to_bind = receiver.to_bind()
            new_aggregate.pushReceiver(receiver_to_bind)
            self.bind_recv.append(receiver_to_bind)
        return new_aggregate

    def create_pybind_rotations(self, num_rots):
        if not self.rotations:
            wp, tag, rag, par, tunn = self.create_pybind_objects()
            self.receiver_configurations = chb.ReceiverConfigurations()
            self.rags = []
            self.receiver_configurations.pushReceiver(rag)
            self.rags.append(rag)
            self.rotations = True
            for i in range(1, num_rots):
                rot = 360 / (num_rots + 1)
                self.Receiver_Agg.rotate(rot)
                rag = self.do_receiver_aggregate()
                self.receiver_configurations.pushReceiver(rag)
                self.rags.append(rag)
            return wp, tag, self.receiver_configurations, par, tunn
        else:
            return self.wp_binds, self.tag_bind, self.receiver_configurations, self.par_bind, self.tunn_bind

    def calculate_fitness(self):
        fitness = FitnessCalculator(self)
        fit = fitness.calculate_fitness(8)
        return fit


class FitnessCalculator:
    def __init__(self, model_seed):
        self.model_seed = model_seed

    def calculate_fitness(self, num_rot=8):
        model = self.get_model_seed()
        model.create_pybind_rotations(num_rot)
        resulting_array = chb.channel_calculation(model.wp_bind, model.tag_bind, model.receiver_configurations,
                                                  model.tunn_bind, model.par_bind)
        fitness = self.calculate_non_dominant(resulting_array, 1e-7)
        return fitness

    def get_model_seed(self):
        return self.model_seed

    def calculate_non_dominant(self, input_array, noise):
        solutions = []
        for matrix in input_array:
            h_one = matrix[:, 0] + matrix[:, 3]
            h_two = matrix[:, 1] + matrix[:, 2]
            sinr_one = [(h1 / (h_two[i] + noise), h_two[i] / (h1 + noise)) for i, h1 in enumerate(h_one)]
            queue = SolQueue()
            for i, tp in enumerate(sinr_one):
                queue.add_value(tp[0], tp[1], i)
            solutions.append(queue.pop_solutions())
        final_queue = SolQueue()
        for i, sol in enumerate(solutions):
            final_queue.add_value(-sol[0], -sol[1], i)
        solution = final_queue.pop_solutions_greed()
        print(input_array[0])
        print(solution)
        return [solution[0], solution[1]]


class SolNode:
    def __init__(self, sinr, i):
        self.sinr = sinr
        self.i = i


class SolQueue:
    def __init__(self):
        self.queue_one = [None, None]
        self.queue_two = [None, None]

    def add_value(self, sinr1, sinr2, i):
        self.compare(self.queue_one, sinr1, i)
        self.compare(self.queue_two, sinr2, i)

    def compare(self, queue, sinr, i):
        if queue[0] is None or queue[0].sinr < sinr:
            sol_node = SolNode(sinr, i)
            try:
                if queue[1].sinr < sinr:
                    queue[0] = queue[1]
                    queue[1] = sol_node
                else:
                    queue[0] = sol_node
            except AttributeError:
                queue[1] = sol_node

    def pop_solutions(self):
        if self.queue_one[1].i == self.queue_two[1].i:
            if (self.queue_one[1].sinr - self.queue_one[0].sinr) < (self.queue_two[1].sinr - self.queue_two[0].sinr):
                return [self.queue_one[0].sinr, self.queue_two[1].sinr]
            else:
                return [self.queue_one[1].sinr, self.queue_two[0].sinr]
        else:
            return [self.queue_one[1].sinr, self.queue_two[1].sinr]
    def pop_solutions_greed(self):
        return [self.queue_one[1].sinr, self.queue_two[1].sinr]


