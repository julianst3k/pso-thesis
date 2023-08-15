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
        optimizer.evaluation(IdentityFunction(2), first_iter=True)

        assert optimizer.particle[0][0] != 0
        assert optimizer.get_bests()[0][0] != 0
        assert len(optimizer.repository.return_non_dominant_nodes())>0


    def test_optimizer_ic_iter(self):
        numero_de_particulas = 100
        dimension = 2
        parameters = {0: {"lb": 0, "ub": 1}, 1: {"lb":0, "ub":1}}
        turbulence = 0.0001
        hypercubes = 10
        optimizer = pso.ParticleSwarmOptimization(numero_de_particulas, dimension, parameters, turbulence, hypercubes)
        initial_function = UniformOutput([0,0], [1,1])
        optimizer.initial_conditions(initial_function)
        optimizer.evaluation(IdentityFunction(2), first_iter=True)
        optimizer.update_state()

        assert all([particle[0]>=0 and particle[1]>=0 and particle[0]<=1 and particle[1]<= 1 for particle in optimizer.get_particle()])

    def test_optimizer_ic_check_rep(self):
        numero_de_particulas = 100
        dimension = 2
        parameters = {0: {"lb": 0, "ub": 1}, 1: {"lb": 0, "ub": 1}}
        turbulence = 0.0001
        hypercubes = 10
        optimizer = pso.ParticleSwarmOptimization(numero_de_particulas, dimension, parameters, turbulence, hypercubes)
        initial_function = UniformOutput([0, 0], [1, 1])
        optimizer.initial_conditions(initial_function)
        optimizer.evaluation(IdentityFunction(2), first_iter=True)
        optimizer.update_state()
        nodes = optimizer.repository.return_non_dominant_nodes()
        print([n_two.get_particle() for n_two in nodes])

        optimizer.evaluation(IdentityFunction(2))
        optimizer.update_state()
        nodes_two = optimizer.repository.return_non_dominant_nodes()
        print([n_two.get_particle() for n_two in nodes_two])
        assert all([particle[0] >= 0 and particle[1] >= 0 and particle[0] <= 1 and particle[1] <= 1 for particle in
                    optimizer.get_particle()])
        assert all([n_two.get_particle()[0] <= n_one.get_particle()[0] and n_two.get_particle()[1] <= n_one.get_particle()[1] for n_two in nodes_two for n_one in nodes])

    def test_optimizer_ic_check_rep(self):
        numero_de_particulas = 100
        dimension = 2
        parameters = {0: {"lb": 0, "ub": 1}, 1: {"lb": 0, "ub": 1}}
        turbulence = 0.0001
        hypercubes = 10
        iters = 10
        initial_function = UniformOutput([0, 0], [1, 1])
        optimizer = pso.ParticleOptimizationIterator(numero_de_particulas, dimension, parameters, turbulence, hypercubes, iters, initial_function, IdentityFunction(2))
        optimizer.initial_conditions(initial_function)
        optimizer.evaluation(IdentityFunction(2), first_iter=True)
        optimizer.update_state()
        nodes = optimizer.repository.return_non_dominant_nodes()
        print([n_two.get_particle() for n_two in nodes])

        optimizer.evaluation(IdentityFunction(2))
        optimizer.update_state()
        nodes_two = optimizer.repository.return_non_dominant_nodes()
        print([n_two.get_particle() for n_two in nodes_two])
        assert all([particle[0] >= 0 and particle[1] >= 0 and particle[0] <= 1 and particle[1] <= 1 for particle in
                    optimizer.get_particle()])
        assert all([n_two.get_particle()[0] <= n_one.get_particle()[0] and n_two.get_particle()[1] <= n_one.get_particle()[1] for n_two in nodes_two for n_one in nodes])





     #   params = {0: {"lb": 0, "ub":1}, 1: {"lb":0, "ub": 1}}
     #   pso_class = pso.ParticleSwarmOptimization(100, 2, params, 0, 10)
     #   pso_class.initial_condition()
