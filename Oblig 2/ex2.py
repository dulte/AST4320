import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
from tqdm import tqdm


class RandomWalk:
    def __init__(self,epsilon,k,max_iter=1e7):
        self.epsilon = epsilon
        self.max_iter = max_iter
        self.k = k


    def reset(self):
        self.S = 2*np.pi/self.k
        self.var = self.get_variance()
        self.path = []
        self.delta = 0


    def get_variance(self):
        return np.pi/self.S**4

    def update_S(self):
        self.S -= self.epsilon



    def update_var(self):
        self.update_S()
        self.var =  self.get_variance()

    def make_step(self):
        self.update_S()
        step_var = self.get_variance() - self.var
        self.delta += np.random.normal(loc=0,scale=sqrt(step_var))
        self.path.append([self.delta,self.S])
        self.var =  self.get_variance()

    def walk(self):
        iter = 0
        while self.S>1 and iter<self.max_iter:
            self.make_step()

            iter += 1

        return np.array(self.path)


class MakeMultipleWalks(RandomWalk):
    def __init__(self,epsilon,k,max_iter=1e7):
        RandomWalk.__init__(self,epsilon,k,max_iter)


    def take_walks(self,nb_walks):
        deltas = np.zeros(int(nb_walks))
        for i in tqdm(range(int(nb_walks))):
            self.reset()
            path = self.walk()

            deltas[i] = path[-1,0]


        return deltas,path[-1,1]


    def get_path_distribution(self,nb_walks,bins=10):
        deltas,S = self.take_walks(nb_walks)


        cont_deltas = np.linspace(np.min(deltas),np.max(deltas),1000)

        plt.hist(deltas,bins="auto",normed=True)
        plt.plot(cont_deltas,self.P(cont_deltas,S))
        plt.show()


    def P(self,delta,S):
        var = np.pi/S**4
        return 1./(sqrt(2*np.pi)*sqrt(var))*np.exp(-(delta**2)/(2*var))


if __name__ == '__main__':
    k = 2*np.pi/((np.pi/1e-4)**(1./4))
    epsilon = 1e-1
    walk = MakeMultipleWalks(epsilon,k)

    walk.get_path_distribution(1e5)
