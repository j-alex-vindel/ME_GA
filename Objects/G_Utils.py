from G_Pop import GA_Pop

class GA_utils:

    def __init__(self,problem,num_ind,num_par_tour,):
        self.num_ind = num_ind
        self.problem = problem
        pass

    def create_pop(self):
        population = GA_Pop()
        for _ in range(self.num_ind):
            individual = self.problem.generate_individual()
            self.problem.calc_obj(individual)
            population.append(individual)
        return population
    
    def fast_nondom_sort(self,population):
        population.fronts = [[]]

        for individual in population:
            individual.dominated_cout = 0
            individual.dominated_solutions = []
            for other_individual in population:
                if individual.dominates(other_individual):
                    individual.dominated_solutions.append(other_individual)
                
                elif other_individual.dominates(individual):
                    individual.dominated_count += 1
            
            if individual.dominated_count ==0:
                individual.rank = 0
                population.fronts[0].append(individual)

        i = 0
        while len(population.fronts[i]) > 0: 
            temp = []
            for individual in population.fronts[i]:
                for other_individual in individual.dominated_solutions:
                    other_individual.domination_count -= 1

                    if other_individual.domination_count == 0:
                        other_individual.rank = i+1
                        temp.append(other_individual)
            i = i + 1
            population.fronts.append(temp)

    def crowding_dist(self,front):
        pass


