from G_Utils import GA_Utils
from G_Pop import GA_Pop

class GA_Evol:

    def __init__(self,problem, num_gen=100, num_ind=100, num_par_tour=2,tour_prob=.9,mutation_rate=1):
        self.utils = GA_Utils(problem=problem, num_ind=num_ind,num_par_tour=num_par_tour,tour_prob=tour_prob,mutation_rate=mutation_rate)
        self.population = None
        self.num_gen = num_gen
        self.on_gen_finished = []
        self.num_ind = num_ind

    def evol(self):
        pass    