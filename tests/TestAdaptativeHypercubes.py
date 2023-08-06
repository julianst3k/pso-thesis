import pytest
import optimization.Archive_Tree as at
import numpy as np
import optimization.AdaptativeHypercubes as ah

class TestAdaptativeHypercubes:

    def test_one(self):
        adap_cubes = ah.AdaptativeHypercubes(2, 10)
        list_values = np.zeros((10,2))
        list_particles = np.zeros((10,2))
        for i in range(10):
            list_values[i][0] = i
            list_values[i][1] = 10-i
            list_particles[i][0] = list_particles[i][1] = i
        tree = at.ArchiveTreeController(list_values, list_particles)
        non_dominated = tree.return_non_dominant_nodes()
        adap_cubes.update_hyper(non_dominated)

        assert adap_cubes.get_max(0) == 9.9

    def test_two(self):
        adap_cubes = ah.AdaptativeHypercubes(2, 10)
        list_values = np.zeros((10,2))
        list_particles = np.zeros((10,2))
        for i in range(10):
            list_values[i][0] = i
            list_values[i][1] = 10-i
            list_particles[i][0] = list_particles[i][1] = i
        tree = at.ArchiveTreeController(list_values, list_particles)
        non_dominated = tree.return_non_dominant_nodes()
        adap_cubes.update_hyper(non_dominated)

        assert adap_cubes.return_index([4.5,4.5])[0] == 5

    def test_two(self):
        adap_cubes = ah.AdaptativeHypercubes(2, 10)
        list_values = np.zeros((10,2))
        list_particles = np.zeros((10,2))
        for i in range(10):
            list_values[i][0] = i
            list_values[i][1] = 10-i
            list_particles[i][0] = list_particles[i][1] = i
        tree = at.ArchiveTreeController(list_values, list_particles)
        non_dominated = tree.return_non_dominant_nodes()
        adap_cubes.hypercube_routine(non_dominated)

        assert adap_cubes.return_index([4.5,4.5])[0] == 5
        assert adap_cubes.get_fitness()[0] == 0.1
