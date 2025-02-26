import model as chb
from abc import ABC, abstractmethod
import ChannelSimulation.SimulationInput as sim
from enum import Enum

class Runner(Enum):
    CHANNEL = 1
    RESPONSE = 2


class GenericRun(ABC):
    @abstractmethod
    def update(self):
        ...
    @abstractmethod
    def __call__(self):
        ...

class ModelRunChannel(GenericRun):
    def __init__(self, config_params):
        self.base_model = sim.SimulatedModel(config_params)
        self.coords_generated = False
    def update(self, coords):
        self.base_model.generate_receiver(coords)
        self.coords_generated = True
    def __call__(self, restype=""):
        if self.coords_generated:
            model = self.base_model
            self.base_model.create_pybind_rotations(8)
            resulting_array = chb.channel_calculation(model.wp_bind, model.tag_bind, model.receiver_configurations,
                                                      model.tunn_bind, model.par_bind, model.shadowing_parameters_bind)

            return resulting_array
        else:
            return ValueError('Coordinates Not Stated')

class ModelRunResponse(GenericRun):
    def __init__(self, config_params):
        self.base_model = sim.SimulatedModel(config_params)
        self.coords_generated = False
    def update(self, coords):
        self.base_model.generate_receiver(coords)
        self.coords_generated = True
    def __call__(self, restype = "", omit_rotations = False):
        if self.coords_generated:
            model = self.base_model
            if not omit_rotations:
                self.base_model.create_pybind_rotations(8)

            if restype == "los":
                resulting_array = chb.response_calculation_los(model.wp_bind, model.tag_bind, model.receiver_configurations,
                                                        model.tunn_bind, model.par_bind, model.shadowing_parameters_bind)
            elif restype == "nlos":
                resulting_array = chb.response_calculation_nlos(model.wp_bind, model.tag_bind, model.receiver_configurations,
                                        model.tunn_bind, model.par_bind, model.shadowing_parameters_bind)
            elif restype == "scatt":
                resulting_array = chb.response_calculation_scattering(model.wp_bind, model.tag_bind, model.receiver_configurations,
                                        model.tunn_bind, model.par_bind, model.shadowing_parameters_bind)
            else:
                resulting_array = chb.response_calculation(model.wp_bind, model.tag_bind, model.receiver_configurations,
                                        model.tunn_bind, model.par_bind, model.shadowing_parameters_bind)

            return resulting_array
        else:
            return ValueError('Coordinates Not Stated')


class ModelRunWrapper(GenericRun):
    def __init__(self, config_params, run_type: Runner):
        if run_type.name == 'CHANNEL':
            self.model_runner = ModelRunChannel(config_params)
        if run_type.name == 'RESPONSE':
            self.model_runner = ModelRunResponse(config_params)
    def update(self, coords):
        self.model_runner.update(coords)
    def __call__(self, restype="", omit_rotations = False):
        return self.model_runner(restype, omit_rotations)

