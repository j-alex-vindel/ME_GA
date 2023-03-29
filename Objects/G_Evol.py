from G_Utils import GA_Utils
from G_Pop import GA_Pop
from tqdm import tqdm
from utility_functions import func_name_print

class GA_Evol:

    def __init__(self,problem, num_gen=100, num_ind=30, num_par_tour=2,tour_prob=.9,mutation_rate=1):
        self.utils = GA_Utils(problem=problem, num_ind=num_ind,num_par_tour=num_par_tour,tour_prob=tour_prob,mutation_rate=mutation_rate)
        self.population = None
        self.num_gen = num_gen
        self.on_gen_finished = []
        self.num_ind = num_ind

    def evol(self):
        self.population = self.utils.create_pop()
        self.utils.fast_nondom_sort(self.population)
        for front in self.population.fronts:
            self.utils.crowding_dist(front)

        children = self.utils.create_children(self.population)
        returned_population = None

        for i in tqdm(range(self.num_gen)):
            print(f"\n*** Generation: {i+1}/{self.num_gen}")
            self.population.extend(children)
            self.utils.fast_nondom_sort(self.population)
            new_population = GA_Pop()
            front_num = 0

            while len(new_population) + len(self.population.fronts[front_num]) <= self.num_ind:
                self.utils.crowding_dist(self.population.fronts[front_num])
                new_population.extend(self.population.fronts[front_num])
                front_num += 1
            self.utils.crowding_dist(self.population.fronts[front_num])
            self.population.fronts[front_num].sort(key=lambda individual: individual.crowding_distance)
            new_population.extend(self.population.fronts[front_num][0:self.num_ind - len(new_population)])
            returned_population = self.population
            self.population = new_population
            self.utils.fast_nondom_sort(self.population)
            for front in self.population.fronts:
                self.utils.crowding_dist(front)
            children = self.utils.create_children(self.population)
        
        return returned_population.fronts[0]
