import numpy as np
import optimization.Archive_Tree as at
import optimization.AdaptativeHypercubes as ah

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
        self.bfp = np.zeros((self.particle_number, len(self.parameters)))
        self.rep = np.zeros((self.particle_number, len(self.parameters)))

    def initial_conditions(self, generator):
        self.particle = generator.generate(self.particle)


    def evaluation(self, fit_function, first_iter=False):
        evaluated_values = np.zeros((self.particle_number, fit_function.get_objective_length()))
        for i, particle in enumerate(self.particle):
            evaluated_values[i,:] = fit_function.fitness(particle)
        self.archive_controller(evaluated_values, self.get_hypercubes())
        self.gen_hypercubes()
        self.set_bests(evaluated_values, first_iter=first_iter)

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

    def archive_controller(self, evaluated_values, hypercubes):
        archive_tree_controller = at.ArchiveTreeController(evaluated_values, self.particle)
        non_dominant_nodes = archive_tree_controller.return_non_dominant_nodes()
        self.repository.update_tree(non_dominant_nodes)

    def set_bests(self, eval, first_iter):
        for i, value in enumerate(eval):
            if first_iter:
                self.bfp[i, :] = value
            else:
                if all([self.bfp[i, j] >= value[j] for j, _ in enumerate(value)]):
                    self.bfp[i, :] = value

    def update_state(self):
        self.velocidad = 0.4*self.velocidad + np.random.rand()*(self.bfp-self.particle) + np.random.rand()*(self.rep - self.particle)+self.turbulence_factor
        self.particle += self.velocidad
        self.check_boundaries()

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
        return self.get_repository()

class ParticleOptimizationIterator:

    def __init__(self, particle_number, dimension, parameters, turbulence, hypercubes, iters, initial_function, fit_function):
        self.parameters = parameters
        self.particle_number = particle_number
        self.turbulence_factor = turbulence
        self.hypercube_number = hypercubes
        self.dimension = dimension
        self.optimizator = ParticleSwarmOptimization(particle_number, dimension, parameters, turbulence, hypercubes)
        self.initial_function = initial_function
        self.function = fit_function
        self.do_optimization(iters, fit_function)

    def reset_optimizator(self):
        self.optimizator = ParticleSwarmOptimization(self.particle_number, self.dimension, self.parameters, self.turbulence, self.hypercubes)

    def do_optimization(self, iters, function):
        self.optimizator.initial_conditions(self.initial_function)
        for i in range(iters):
            self.optimizator.evaluation(function)
            self.update_state()
        return

    def get_optimizator(self):
        return self.optimizator

    def get_best_results(self):
        return self.optimizator.get_repository()


