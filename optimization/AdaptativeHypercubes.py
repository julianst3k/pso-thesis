import optimization.Archive_Tree as at
import numpy as np
import logging

FORMAT = '%(asctime)s %(name)-12s %(message)s'
logging.basicConfig(format=FORMAT,level=logging.WARNING,force=True)
logger = logging.getLogger('AdaptativeHypercubes')


class AdaptativeHypercubes:

    def __init__(self, dims, qty):
        self.hypercube_params = {}
        self.qty = qty
        for i in range(dims):
            self.hypercube_params[i] = {"max": 0, "min": 0}
        self.storage = np.empty((10,10), dtype=object)
        self.collection_of_cubes = []
        self.fitness = []
        self.count = []

    def return_index(self, values):
        dimension = len(values)
        coordinates = np.zeros(dimension, dtype=np.int32)
        logger.info("Valor %s", values)
        for i,value in enumerate(values):
            max_val = self.get_max(i)
            min_val = self.get_min(i)
            logger.info("Valor, Valor: %s", value)
            logger.info("Valor, Maximo: %s, Minimo: %s", max_val, min_val)
            try:
                index = round((value-min_val)*self.qty/(max_val-min_val))
            except OverflowError:
                index = 5
            except ValueError:
                index = 5

            coordinates[i] = int(index)
        return coordinates

    def update_hyper(self, nodes):
        logger.info("Update %s", nodes)

        value_array = np.empty((len(nodes),len(self.get_hypercube_params())))
        for i, node in enumerate(nodes):
            value = node.get_value()
            value_array[i,:] = value
        minimal_values = np.min(value_array, axis=0)
        maximal_values = np.max(value_array, axis=0)
        logger.info("Value Array %s, Max: %s, Min: %s", value_array, maximal_values, minimal_values)
        for dim in self.hypercube_params:
            slack_factor = (maximal_values[dim]-minimal_values[dim])/10
            self.hypercube_params[dim]["max"] = maximal_values[dim]+slack_factor
            self.hypercube_params[dim]["min"] = minimal_values[dim]-slack_factor

    def hypercube_store_positions(self, nodes):
        logger.info("Storing %s", nodes)
        self.collection_of_cubes = []
        for node in nodes:
            value = node.get_value()
            coord = self.return_index(value)
            try:
                self.storage[coord[0], coord[1]].add(node)
            except AttributeError:
                self.storage[coord[0], coord[1]] = at.ArchiveTreeNodeCollector()
                self.storage[coord[0], coord[1]].add_node(node)
                self.collection_of_cubes.append(self.storage[coord[0], coord[1]])

    def compute_fitness(self):
        self.fitness = [10/len(cube) if len(cube) > 1 else 10 for cube in self.collection_of_cubes]
        sum_fitness = sum(self.fitness)
        self.fitness = [val/sum_fitness for val in self.fitness]

    def roulette_wheel(self):
        select_value = np.random.choice(len(self.get_fitness()), p=self.get_fitness())
        cube_selected = self.collection_of_cubes[select_value]
        return cube_selected.get_random_node()

    def hypercube_routine(self, nodes):
        self.update_hyper(nodes)
        self.hypercube_store_positions(nodes)
        self.compute_fitness()

    def get_count(self):
        return self.count

    def get_hypercube_params(self):
        return self.hypercube_params

    def get_fitness(self):
        return self.fitness

    def get_max(self, dim):
        logger.info("Maximo %s", self.get_hypercube_params()[dim]["max"])
        return self.get_hypercube_params()[dim]["max"]

    def get_min(self, dim):
        return self.get_hypercube_params()[dim]["min"]

    def get_list_max(self):
        return [self.get_hypercube_params()[key]["max"] for key in self.get_hypercube_params()]

    def get_list_min(self):
        return [self.get_hypercube_params()[key]["min"] for key in self.get_hypercube_params()]
