import numpy as np
from abc import ABC
import logging

FORMAT = '%(asctime)s %(name)-12s %(message)s'
logging.basicConfig(format=FORMAT, level=logging.DEBUG, force=True)
logger = logging.getLogger(__name__)


class IdentityFunction:
    def __init__(self, objective_length):
        self.objective_length = objective_length
    def fitness(self, value):
        return value
    def get_objective_length(self):
        return self.objective_length



class OutputFunction(ABC):

    def __init__(self, lambda_function):
        self.lambda_function = lambda_function

    def generate(self, instance, owner):
        output = np.zeros(owner.get_list().shape)
        for i,value in enumerate(owner.get_list()):
            output[i,:] = self.lambda_function(value)
        return output


class UniformOutput(OutputFunction):

    def __init__(self, lb_limits, ub_limits):
        self.lb_limits = lb_limits
        self.ub_limits = ub_limits

    def generate(self, input_list):
        qty_outputs, qty_params = input_list.shape
        lb = np.tile(self.lb_limits, (qty_outputs, 1))
        ub = np.tile(self.ub_limits, (qty_outputs, 1))
        try:
            output = lb+np.random.rand(qty_outputs, qty_params)*(ub-lb)
            logger.info('Output Example: %s', output[0, :]) # debug purpose
            return output
        except ValueError:
            logger.warning('Shape: %s', lb)


class RandomOutput(OutputFunction):

    def __init__(self, multiplier):
        self.multiplier = multiplier

    def generate(self, input_list):
        return np.random.rand(input_list.shape)






