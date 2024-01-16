import numpy as np
import optimization.Archive_Tree as at
import optimization.AdaptativeHypercubes as ah
from threading import Thread
from optimization.utils import UniformOutput

class ParticleSwarmOptimization:
    def __init__(self, particle_number, dimension, parameters, turbulence, hypercubes):
        self.parameters = parameters
        self.particle_number = particle_number
        self.turbulence_factor = turbulence
        self.hypercube_number = hypercubes
        self.hypercubes = ah.AdaptativeHypercubes(dimension,self.hypercube_number)
        self.repository = at.ArchiveTreeController([], [], self.hypercubes)
        self.dimension = dimension
        self.initialize_parameters()

    def initialize_parameters(self):
        self.velocidad = np.zeros((self.particle_number,len(self.parameters)))
        self.particle = np.zeros((self.particle_number,len(self.parameters)))
        self.eval_vec = np.zeros(self.particle_number)
        self.bfp = np.zeros((self.particle_number, self.dimension))
        self.bfp_pos = np.zeros((self.particle_number, len(self.parameters)))

        self.rep = np.zeros((self.particle_number, len(self.parameters)))

    def initial_conditions(self, generator, aggs):
        self.aggs = aggs
        self.particle = generator.generate(self.particle)


    def evaluation(self, fit_function, first_iter=False):
        evaluated_values = np.zeros((self.particle_number, fit_function.get_objective_length()))
        for i, particle in enumerate(self.particle):
            evaluated_values[i,:] = fit_function.fitness(particle)
        self.archive_controller(evaluated_values, self.get_hypercubes())
        self.gen_hypercubes()
        self.set_bests(evaluated_values, first_iter=first_iter)
    def multiprocess_evaluation(self, fit_function, first_iter=False):
        agss_values = np.zeros((len(self.aggs), fit_function.get_objective_length()))
        if first_iter:
            processes = []
            for i in range(4):
                processes.append(Thread(target=self.multiprocess_fitness, args=(fit_function, agss_values, i*len(self.aggs)//4, (i+1)*len(self.aggs)//4, True,)))
                processes[i].start()
            for p in processes:
                p.join()
            self.archive_controller(agss_values, self.get_hypercubes(), is_agg=True)

        evaluated_values = np.zeros((self.particle_number, fit_function.get_objective_length()))
        processes = []
        for i in range(4):
            processes.append(Thread(target=self.multiprocess_fitness, args=(fit_function, evaluated_values, i*self.particle_number//4, (i+1)*self.particle_number//4,)))
            processes[i].start()
        for p in processes:
            p.join()

        self.archive_controller(evaluated_values, self.get_hypercubes())
        self.gen_hypercubes()
        self.set_bests(evaluated_values, first_iter=first_iter)

    def multiprocess_fitness(self, fit_function, evaluated_values, start, end, is_agg = False):
        for i in range(start, end):
            if is_agg:
                evaluated_values[i, :] = fit_function.fitness(self.aggs[i])
                print(self.aggs[i], evaluated_values[i,:], "aggs")
            else:
                evaluated_values[i, :] = fit_function.fitness(self.particle[i, :])
                print(self.particle[i,:], evaluated_values[i, :])




    def gen_hypercubes(self):
        nodes = self.repository.return_non_dominant_nodes()
        self.hypercubes.hypercube_routine(nodes)

    def get_hypercubes(self):
        return self.hypercubes

    def get_random_rep(self):
        for i, _ in enumerate(self.rep):
            node = self.hypercubes.roulette_wheel()
            particle = node.get_particle()
            self.rep[i, :] = particle

    def archive_controller(self, evaluated_values, hypercubes, is_agg = False):
        if is_agg:
            archive_tree_controller = at.ArchiveTreeController(evaluated_values, self.aggs, is_agg = is_agg)
            self.agg_archive = archive_tree_controller
        else:
            archive_tree_controller = at.ArchiveTreeController(evaluated_values, self.particle)
        non_dominant_nodes = archive_tree_controller.return_non_dominant_nodes()
        self.repository.update_tree(non_dominant_nodes)

    def set_bests(self, eval, first_iter):
        for i, value in enumerate(eval):
            if first_iter:
                self.bfp[i, :] = value
                self.bfp_pos[i, :] = np.array([val for val in self.particle[i, :]])
            else:
                if all([self.bfp[i, j] >= value[j] for j, _ in enumerate(value)]):
                    self.bfp[i, :] = value
                    self.bfp_pos[i, :] = np.array([val for val in self.particle[i, :]])
                elif any([self.bfp[i,j]>-0.5 for j, _ in enumerate(self.bfp[i,:])]) and all(value[j]<-0.5 for j, _ in enumerate(value)):
                    self.bfp[i,:] = value
                    self.bfp_pos[i, :] = np.array([val for val in self.particle[i, :]])

    def update_state(self):
        self.get_random_rep()
        self.velocidad = 0.4*self.velocidad + np.random.rand(self.particle_number, len(self.parameters))*(self.bfp_pos-self.particle) \
                         + np.random.rand(self.particle_number, len(self.parameters))*(self.rep - self.particle)+(0.5-np.random.rand(self.particle_number, len(self.parameters)))*self.turbulence_factor
        self.particle += self.velocidad
        self.check_boundaries()

    def mutate_state(self, it, max_iter, mutrate = 2.5):
        for particle in self.particle:
            dim = len(particle)
            for i in range(dim):
                dims_affected = np.random.binomial(1, (1-it/max_iter)**(dim/mutrate))
                if dims_affected == 1:
                    param = list(self.parameters.keys())[i]
                    particle[i] = np.random.uniform(self.parameters[param]["lb"], self.parameters[param]["ub"])

    def check_boundaries(self):
        lb_array = []
        ub_array = []
        for param in self.parameters:
            lb_array.append(self.parameters[param]["lb"])
            ub_array.append(self.parameters[param]["ub"])
        lb_comp = np.tile(lb_array, (self.particle_number, 1))
        ub_comp = np.tile(ub_array, (self.particle_number, 1))
        lb_comparison = self.particle < lb_comp
        ub_comparison = self.particle > ub_comp
        self.velocidad += -2*self.velocidad*lb_comparison
        self.velocidad += -2*self.velocidad*ub_comparison
        self.particle += (lb_comp-self.particle)*lb_comparison+(ub_comp-self.particle)*ub_comparison

    def get_bests(self):
        return self.bfp

    def get_particle(self):
        return self.particle

    def get_repository(self):
        return self.repository

class ParticleOptimizationIterator:

    def __init__(self, particle_number, dimension, parameters, turbulence, hypercubes, iters, initial_function, fit_function, track = False, multiprocess = False, agg = []):
        self.parameters = parameters
        self.particle_number = particle_number
        self.turbulence_factor = turbulence
        self.hypercube_number = hypercubes
        self.dimension = dimension
        self.optimizator = ParticleSwarmOptimization(particle_number, dimension, parameters, turbulence, hypercubes)
        self.initial_function = initial_function
        self.function = fit_function
        self.results = []
        self.track, self.sols, self.agg_nodes = self.do_optimization(iters, fit_function, track, multiprocess, agg)

    def reset_optimizator(self):
        self.optimizator = ParticleSwarmOptimization(self.particle_number, self.dimension, self.parameters, self.turbulence, self.hypercubes)

    def do_optimization(self, iters, function, track = False, multiprocess = False, agg = []):
        self.optimizator.initial_conditions(self.initial_function, agg)
        opt_dict = {}
        sol_dict = []
        for i in range(iters):
            if not multiprocess:
                self.get_optimizator().evaluation(function)
            else:
                self.get_optimizator().multiprocess_evaluation(function, first_iter=i==0)
            self.get_optimizator().update_state()
            if track:
                nodes = self.get_best_results().return_non_dominant_nodes()
                if i == 0:
                    opt_dict["iter"] = [i]*len(nodes)
                    opt_dict["particles"] = [node.particle for node in nodes]
                    opt_dict["value"] = [node.value for node in nodes]
                else:
                    opt_dict["iter"].extend([i]*len(nodes))
                    opt_dict["particles"].extend([node.particle for node in nodes])
                    opt_dict["value"].extend([node.value for node in nodes])
                if i == iters-1:
                    sol_dict = [node.particle for node in nodes]
            print(list(zip(opt_dict["particles"], opt_dict["value"])))
            self.get_optimizator().mutate_state(i+1, iters+1)

        return opt_dict, sol_dict, self.get_optimizator().agg_archive.return_non_dominant_nodes()

    def get_optimizator(self):
        return self.optimizator

    def get_best_results(self):
        return self.optimizator.get_repository()

    def get_tracked_results(self):
        return self.track, self.sols, self.agg_nodes

