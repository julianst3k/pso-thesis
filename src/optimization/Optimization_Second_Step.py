import ChannelModelling.ParameterModel as pm
import optimization.ParticleSwarmOptimization as pso
from yaml import load, dump, Loader
from optimization.utils import UniformOutput
import pandas as pd
import numpy as np
import logging
FORMAT = '%(asctime)s %(name)-12s %(message)s'
logging.basicConfig(format=FORMAT,level=logging.WARNING,force=True)

logger = logging.getLogger("Optimization Second Step")


class FitnessModel:
    def __init__(self, center):
        self.center = center
    def fitness(self, particle):
        mseed = pm.ModelSeed(particle[0], particle[1], particle[2], particle[3], particle[4], self.center)
        return mseed.calculate_fitness()
    def get_objective_length(self):
        return 2


class Evaluator:
    def __init__(self, particle):
        self.particle = particle
    def evaluate(self, x_space, y_space):
        queue = pm.SolQueue()
        i = 0
        for x in x_space:
            for y in y_space:
                center = [x,y,1.8]
                fitness_function = FitnessModel(center)
                ret = fitness_function.fitness(self.particle)
                queue.add_value(ret[0], ret[1], i)
                i+=1

        return queue.pop_solutions_greed()





if __name__=="__main__":
    candidates = pd.read_csv("pso_opt_valfilt.csv")
    sols = []
    dictionary = {}
    for i, c in candidates.iterrows():
        if (c["x"]<=0.8 or c["y"]>=2):
            param_c = [c["betat"], c["betar"], c["fov"], c["alphai"], c["alphao"]]
            eval = Evaluator(param_c)
            eval_res = eval.evaluate(np.arange(0.3,3,0.5), np.arange(0.5,2.6,0.5))
            logger.warning("Ret: %f, Ret: %f", eval_res[0], eval_res[1])
            dictionary["sinr_h1"], dictionary["sinr_h2"] = eval_res[0], eval_res[1]
            dictionary["betat"], dictionary["betar"], dictionary["fov"], dictionary["alphai"], dictionary["alphao"] = c["betat"], c["betar"], c["fov"], c["alphai"], c["alphao"]
            sols.append(dictionary)
    out_df = pd.DataFrame(sols)
    out_df.to_csv("Final_Result.csv")
