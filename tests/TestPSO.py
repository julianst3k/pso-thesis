import pytest
from optimization.utils import IdentityFunction, UniformOutput
import optimization.ParticleSwarmOptimization as pso
import numpy as np
import logging


class TestPSO:
    def test_initial_conditions(self):
        initial_function = UniformOutput([0,0,0], [1,1,1])
        initial_array = np.zeros((100,3))
        result = initial_function.generate(initial_array)

        assert result.shape[0] == 100

    def test_optimizer_ic(self):
        numero_de_particulas = 100
        dimension = 2
        parameters = {0: {"lb": 0, "ub": 1}, 1: {"lb":0, "ub":1}}
        turbulence = 0.0001
        hypercubes = 10
        optimizer = pso.ParticleSwarmOptimization(numero_de_particulas, dimension, parameters, turbulence, hypercubes)
        initial_function = UniformOutput([0,0], [1,1])
        optimizer.initial_conditions(initial_function)
        eval = optimizer.evaluation(IdentityFunction(2), first_iter=True)

        assert optimizer.particle[0][0] != 0
        assert optimizer.get_bests()[0][0] != 0



     #   params = {0: {"lb": 0, "ub":1}, 1: {"lb":0, "ub": 1}}
     #   pso_class = pso.ParticleSwarmOptimization(100, 2, params, 0, 10)
     #   pso_class.initial_condition()
