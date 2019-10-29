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
        while self.S-self.epsilon>=1 and iter<self.max_iter:
            self.make_step()

            iter += 1

        return np.array(self.path)
    
    def restricted_walk(self,delta_crit):
        iter = 0
        while self.S-self.epsilon>=1 and iter<self.max_iter:
            self.make_step()

            iter += 1
            
            if self.delta >= delta_crit:
                return None
        
        if iter < self.max_iter:
            return np.array(self.path)
        else:
            return None


class MakeMultipleWalks(RandomWalk):
    def __init__(self,epsilon,k,max_iter=1e7):
        RandomWalk.__init__(self,epsilon,k,max_iter)
        self.S = None
        self.deltas = None
        self.restriced_deltas = None


    def take_walks(self,nb_walks):
        deltas = np.zeros(int(nb_walks))
        for i in tqdm(range(int(nb_walks))):
            self.reset()
            path = self.walk()

            deltas[i] = path[-1,0]


        self.deltas,self.S = deltas,path[-1,1]


    def take_restricted_walks(self,nb_walks,delta_crit=1):
        deltas = []
        for i in tqdm(range(int(nb_walks))):
            self.reset()

            path = self.walk()
            if not (path[:,0] > delta_crit).any():
                deltas.append(path[-1,0])
                self.S = path[-1,1]
         
        self.restriced_deltas = deltas


    def get_path_distribution(self):
        deltas = self.deltas
        S = self.S
        cont_deltas = np.linspace(np.min(deltas),np.max(deltas),1000)

        plt.hist(deltas,bins="auto",normed=True,label="Distrbution of Random Walk")
        plt.plot(cont_deltas,self.P(cont_deltas,S),"r",label="Normal Distrbution")
        plt.title(r"Distrbution of Density Contrast $\delta$",fontsize=25)
        plt.xlabel(r"$\delta$",fontsize=25)
        plt.ylabel(r"$P(\delta|M)$",fontsize=25)
        plt.legend(loc="best")
        plt.show()


    def get_path_never_cross_distribution(self,delta_crit=1):


        deltas = self.restriced_deltas
        S = self.S
        cont_deltas = np.linspace(np.min(deltas),np.max(deltas),1000)
      
        plt.hist(deltas,bins=100,weights=np.ones_like(100)*1e4,label=r"Distrbution of Random Walk")
        plt.plot(cont_deltas,self.P_nc(cont_deltas,S,delta_crit),"r",label="Normal Distrbution")
        plt.title(r"Distrbution of Density Contrast $\delta$ That Never Crosses $\delta_{crit}$ = %s" %delta_crit,fontsize=25)
        plt.xlabel(r"$\delta$",fontsize=25)
        plt.ylabel(r"$P_{nc}(\delta|M)$",fontsize=25)
        plt.legend(loc="best")
        plt.show()






    def P(self,delta,S):
        var = np.pi/S**4
        return 1./(sqrt(2*np.pi)*sqrt(var))*np.exp(-(delta**2)/(2*var))

    def P_nc(self,delta,S,delta_crit):
        S = 1.0
        var = np.pi/S**4

        

        return 1./(sqrt(2*np.pi)*sqrt(var))*(np.exp(-(delta**2)/(2*var)) - \
                np.exp(-((2*delta_crit-delta)**2)/(2*var))) *2.2 #With 

if __name__ == '__main__':
    k = 2*np.pi/((np.pi/1e-4)**(1./4))
    epsilon = 1.23e-1
    walk = MakeMultipleWalks(epsilon,k)
    walk.take_walks(1e5)
    walk.get_path_distribution()
    walk.take_restricted_walks(1e5,1)
    walk.get_path_never_cross_distribution()

    
