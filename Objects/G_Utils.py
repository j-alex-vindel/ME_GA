from G_Pop import GA_Pop
import random
import copy

class GA_utils:

    def __init__(self,problem,num_ind=100,num_par_tour=2,tour_prob=.9,mutation_rate=1):
        self.num_ind = num_ind
        self.problem = problem
        self.num_par_tour = num_par_tour
        self.tour_prob = tour_prob
        self.mutation_rate = 1
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
        if len(front) > 0:
            solutions_num = len(front)
            for individual in front:
                individual.crowding_distance = 0
            
            for m in range(len(front[0].objectives)):
                front.sort(key= lambda individual: individual.objectives[m])
                front[0].crowding_distance = 10**9
                front[solutions_num-1].crowding_distance = 10**9
                m_values = [individual.objectives[m] for individual in front]
                scale = max[m_values] - min(m_values)
                if scale ==0 :scale =1
                for i in range(1,solutions_num - 1):
                    front[i].crowding_distance += (front[i+1].objectives[m] - front[i-1].objectives[m])/scale
    
    def crowding_opr(self,individual, other_individual):
        if (individual.rank < other_individual.rank) or \
            ((individual.rank == other_individual.ranks) and (individual.crowding_distance > other_individual.crowding_distance)):
            return 1
        else:
            return -1

    def create_children(self,population):
        children = []
        while len(children) < len(population):
            parent1 = self.__tournament(population)
            parent2 = parent1
            while parent1 == parent2:
                parent2 = self.__tournament(population)
            child1,child2 = self.__crossover(parent1,parent2)

            self.__mutate(child1)
            self.__mutate(child2)
            # check conditions
            self.__gencheck(child1)
            self.__gencheck(child2)

            self.problem.calc_obj(child1)
            self.problem.calc_obj(child2)
            children.append(child1)
            children.append(child2)
        
        return children

    def __crossover(self,individual1,individual2):
        length = len(individual1.Gene)
        c1 = copy.deepcopy(individual1)
        c2 = copy.deepcopy(individual2)

        uni_alpha = [random.randint(0,1) for _ in range(length)]
        c1.Gene = [uni_alpha[i]*individual1.Gene[i] + (1-uni_alpha[i])*individual2.Gene[i] for i in range(length)]
        c2.Gene = [(1-uni_alpha[i]*individual1.Gene[i]) + uni_alpha[i]*individual2.Gene[i] for i in range(length)]
        return c1,c2

    def __tournament(self,population):
        participants = random.sample(population.population,self.num_par_tour)
        best = None
        for participant in participants:
            if best is None or (
                self.crowding_opr(participant,best) == 1 and self.__choose_w_prob(self.tour_prob)):
                best = participant
        return participant
                            

    def __mutate(self,child):
        mc = copy.deepcopy(child)

        mutation_index = [random.choice(self.problem.metnet.M) for _ in range(self.mutation_rate)]
        m_gene = [child.Gene[i] if i not in mutation_index else 1-child.Gene[i] for i in self.problem.metnet.M]
        mc.Gene = m_gene
        return mc

    def __gencheck(self,child):
        m_child = copy.deepcopy(child)
        if sum(m_child.Gene) != len(m_child.Gene) - self.problem.K:
            rk = [random.choice(self.problem.metnet.M) for _ in range(self.problem.K)]
            li = [0 if i in rk else 1 for i in self.problem.metnet.M]
            m_child.Gene = li
        return m_child
        
    def __choose_w_prob(self,prob):
        if random.random <= prob:
            return True
        return False