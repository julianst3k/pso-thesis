import numpy as np
import random
import matplotlib.pyplot as plt
class MetropolisHastings:
    def __init__(self, f, q):
        self.q = q
        self.f = f
    def fit(self, x_o, N):
        x_out = []
        for i in range(N):
            u = np.random.rand()
            x_n = abs(np.random.normal(x_o, 10))
            A = min(1, self.f(x_n, x_o)*self.q(x_o,x_n)/self.q(x_n,x_o))
            x_o = x_n if u<A else x_o
            x_out.append(x_o)
        return x_out
class DataGenerator:
    def __init__(self, u, o, N):
        self.x = [np.random.normal(u, o) for i in range(N)]
    def data(self):
        return self.x

class Function:
    def __init__(self, data):
        self.data = data


class MuProbability(Function):
    def __init__(self, data):
        super().__init__(data)
    def __call__(self, mu, mu_2):
        return np.prod([np.exp(-(d - mu)**2/(2*5))/np.exp(-(d - mu_2)**2/(2*5)) for d in self.data]) * np.exp(-(mu**2)/(2*10))/np.exp(-(mu_2**2)/(2*10))


class QProbability(Function):
    def __init__(self):
        pass
    def __call__(self, x_o, x_n):
        return np.exp(-(x_o-x_n)**2/(2*10)) + np.exp(-(x_o+x_n)**2/(2*10)) 



def simulation(N: int, iters: int):
    mu = abs(np.random.normal(0, 10))
    data = DataGenerator(mu, 5, N).data()
    prob_Mu = MuProbability(data)
    q_prob = QProbability()
    metropolis = MetropolisHastings(prob_Mu, q_prob)
    return metropolis.fit(0, iters), mu


if __name__ == "__main__":
    random.seed(19814697)
    list_met, mu = simulation(100, 10000)
    frq, edges = np.histogram(list_met, 100)
    fig, ax = plt.subplots()
    ax.bar(edges[:-1], frq, width=np.diff(edges))
    print(list_met)
    print(mu)
    plt.show()

