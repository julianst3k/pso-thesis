import ChannelModelling.ParameterModel as pm
import optimization.ParticleSwarmOptimization as pso
from yaml import load, dump, Loader
from optimization.utils import UniformOutput
import pandas as pd
import numpy as np
import re
import logging
FORMAT = '%(asctime)s %(name)-12s %(message)s'
logging.basicConfig(format=FORMAT,level=logging.INFO,force=True)
logger = logging.getLogger('Optimization')


class FitnessModel:
    def __init__(self, center):
        self.center = center
    def fitness(self, particle):
        mseed = pm.ModelSeed(particle[0], particle[1], particle[2], particle[3], particle[4], self.center)
        return mseed.calculate_fitness()
    def get_objective_length(self):
        return 2

class Optimizator:
    def __init__(self, center, config):
        self.fitness_function = FitnessModel(center)
        self.center = center
        fl = open(config, 'r')
        self.config = load(fl, Loader=Loader)
        self.initial_function = self.create_initial_function()
        self.tracked_results = []
    def create_initial_function(self):

        order = ["betat", "betar", "fov", "alphao", "alphai"]
        lb = [0] * len(order)
        ub = [0] * len(order)
        for param in self.config["parameters"]:
            index_order = order.index(list(param.keys())[0])
            lb[index_order] = param[list(param.keys())[0]]["lb"]
            ub[index_order] = param[list(param.keys())[0]]["ub"]
        ret = UniformOutput(lb,ub)
        return ret
    def do_optimization(self, agg):
        order = ["betat", "betar", "fov", "alphao", "alphai"]
        new_config = {list(param.keys())[0]: {"lb":param[list(param.keys())[0]]["lb"],
                                              "ub":param[list(param.keys())[0]]["ub"]} for param in self.config["parameters"]}
        sorted_dict = {elem: new_config[elem] for elem in order}
        optimizer = pso.ParticleOptimizationIterator(self.config["numero_de_particulas"], self.config["dimension"],
                                                     sorted_dict, self.config["turbulence"],
                                                     self.config["hypercubes"], self.config["iters"],
                                                     self.initial_function, self.fitness_function, track = True, multiprocess=True, agg = agg)
        self.tracked_results, self.good_sols, archive_node = optimizer.get_tracked_results()
        return self.tracked_results, self.good_sols, archive_node
    def save_value(self, name):
        df = pd.DataFrame(self.tracked_results)
        df.to_csv(name)

if __name__=="__main__":
    x_spacing = np.arange(2.8,3,0.5)
    y_spacing = np.arange(0.5,1.3,0.5)
    good_sols_agg = [[45,60,60,0,0],[60,30,60,-20,-20]]
    archive_nodes_agg = []
    logger.info("Starting Optimization")
    for x in x_spacing:
        for y in y_spacing:
            logger.info("Optimizing {}, {}".format(x, y))
            output_file = "pso_final_"+str(x).replace(".","")+str(y).replace(".","")+".csv"
            opt = Optimizator([x,y,1.8], "config.yaml")
            tracked_results, _, _ = opt.do_optimization(good_sols_agg)
            opt.save_value(output_file)
            del opt
            del tracked_results
