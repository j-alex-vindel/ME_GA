import copy
import math
import numpy as np
import random
from random import randint
from typing import List, Tuple,Type
from Individual import Individual, Out_Result
from G_Models import mobjective_model,mobjective_solve,true_vij
from Met_Net import Metabolic_Network
from utility_functions import func_timer
from collections import namedtuple

Fit = namedtuple('Fitness',['Cost','Chemical','Biomass','True_chem'])
Genome = List[int]
Population = List[Individual]
Model = Type[object]
MN = Type[Metabolic_Network]
MutationRate = int 
Priority = List[float]
Weights = List[float]

@func_timer
def run_ga(maxit:int,Npop:int,obj:MN,fmodel:Model,smodel:Model,priority:Priority,weights:Weights,k:int,pc:int=1,beta:float=1,mu:float=.1):
    nc = np.round(pc*Npop/2)*2
    # Creating a population
    pop = [Individual() for _ in range(Npop)] 

    # Create a best solution structure
    bestsol = Individual(Cost = -np.Infinity) 

    # Initial Solutions
    for i in range(Npop):
        pop[i].Gene = gene_generator(obj,k)
        # pop[i].Strat =  get_strategy(pop[i].Gene,obj) 
        pop[i].Chem, pop[i].Bio = mobjective_solve(obj,pop[i].Gene,fmodel,obj.M,priority,weights)
        vi = true_vij(obj,pop[i].Gene,smodel,obj.M)
        pop[i].Cost = pop[i].Bio + pop[i].Chem
        pop[i].TrueBio,pop[i].TrueChem = vi[obj.biomass], vi[obj.chemical]
        pop[i].Fitness = Fit(pop[i].Cost,pop[i].Chem,pop[i].Bio,pop[i].TrueChem)
        if pop[i].Cost > bestsol.Cost:
            bestsol = copy.deepcopy(pop[i])

    # Best Cost of Iterations
    bestcost = np.empty(maxit)
    

    # Main Loop     
    for it in range(maxit):

        costs = np.array([x.Cost for x in pop]) # List of costs in the population
        avg_costs = np.mean(costs)
        if avg_costs != 0:
            costs = costs / avg_costs
        probs = np.exp(-beta*costs)



        popc = [] #Children Population
        for _ in range(int(nc//2)):
           
            # Select Parents Randomely (Roulette Wheel)
            # q = np.random.permutation(Npop)
           
            # p1 = pop[roulette_wheel(probs)]
            # p2 = pop[roulette_wheel(probs)]

            p1 = tournament_selection(pop)
            p2 = tournament_selection(pop)
            
            # Crossover
            c1,c2 = two_point_crossover(p1,p2,obj)

            # Mutate
            c1 = mutate(c1,mu)
            c2 = mutate(c2,mu)

            # Define a feasibility function to allow the offsrpings be on the feasible side
            
            # Evaluate Offspring 1
            c1.Chem, c1.Bio = mobjective_solve(obj,c1.Gene,fmodel,obj.M,priority,weights)
            c1.Cost = c1.Bio + c1.Chem
            vs = true_vij(obj,c1.Gene,smodel,obj.M)
            c1.TrueBio,c1.TrueChem = vs[obj.biomass],vs[obj.chemical]
            c1.Fitness = Fit(c1.Cost,c1.Chem,c1.Bio,c1.TrueChem)
            if c1.Cost > bestsol.Cost:
                bestsol = copy.deepcopy(c1)
                
            
            # Evaluate Offspring 2
            c2.Chem, c2.Bio = mobjective_solve(obj,c2.Gene,fmodel,obj.M,priority,weights)
            c2.Cost = c2.Bio + c2.Chem
            vs = true_vij(obj,c2.Gene,smodel,obj.M)
            c2.TrueBio,c2.TrueChem = vs[obj.biomass],vs[obj.chemical]
            c2.Fitness = Fit(c2.Cost,c2.Chem,c2.Bio,c2.TrueChem)
            if c2.Cost > bestsol.Cost:
                bestsol = copy.deepcopy(c2)
         
            popc.append(c1)
            popc.append(c2)

        # Merge Sort and Select
        pop += popc
        pop = sorted(pop,key=lambda x: x.Cost,reverse=True)
        pop = pop[0:Npop]
        

        # Store Best Cost
        bestcost[it] = bestsol.Cost


        # Show Iteration Information
        print("Iteration {}: Best Cost = {}".format(it+1, bestcost[it]))


    result=Out_Result()
    result.Pop = pop
    result.BestCost = bestcost
    result.BestSol = bestsol
 

    return result


def gene_generator(obj:MN,k:int) -> Genome:
    rk = [random.choice(obj.KO) for i in range(k)]
    li = [0 if i in rk else 1 for i in obj.M]
    return li

def uniform_crossover(pa:Type[Individual],pb:Type[Individual]) -> Tuple[Type[Individual],Type[Individual]]:
    if len(pa.Gene) != len(pb.Gene):
        raise ValueError("Genomes a and b must be of same length")
    lenght = len(pa.Gene)
    c1 = copy.deepcopy(pa)
    c2 = copy.deepcopy(pb)
    
    uni_alpha = [randint(0,1) for i in range(lenght)]

    c1.Gene = [uni_alpha[i]*pa.Gene[i] + (1-uni_alpha[i])*pb.Gene[i] for i in range(len(pa.Gene))]

    c2.Gene = [(1-uni_alpha[i]*pa.Gene[i]) + uni_alpha[i]*pb.Gene[i] for i in range(len(pa.Gene))]

    return c1,c2

def mutate(p:Type[Individual],mu:MutationRate) -> Type[Individual]:

    # currently performs under a double mutation rate
    y = copy.deepcopy(p)
    flag = np.random.rand(len(y.Gene))
    ri = [i for i in range(len(y.Gene)) if flag[i] <= mu]
    y.Gene = [y.Gene[i] if i not in ri else 1-y.Gene[i] for i in range(len(y.Gene))]
    return y

def get_strategy(y:Genome,obj:MN) -> List[str]:
    kindex = [index for index,value in enumerate(y) if value == 0]

    return [obj.rxn[i] for i in kindex]

def roulette_wheel(p:List[float]) -> int:
    c = np.cumsum(p)
    r = sum(p)*np.random.rand()
    ind = np.argwhere(r<=c)
    return ind[0][0]

def tournament_selection(pop: Population) -> Genome:
    rk = [random.choice([_ for _ in range(len(pop))]) for i in range(2)] # returns a list with two indices
    f1 = copy.deepcopy(pop[rk[0]])
    f2 = copy.deepcopy(pop[rk[1]])
    if f1.Cost >= f2.Cost:
        return f1
    else:
        return f2

def two_point_crossover(pa:Type[Individual],pb:Type[Individual],obj:MN) -> Tuple[Type[Individual],Type[Individual]]:
    c1 = copy.deepcopy(pa)
    c2 = copy.deepcopy(pb)
    if c1.Cost >= c2.Cost and sum(c1.Gene) <= len(c1.Gene)-2:
        padre = c1.Gene
        a,b = [index for index,value in enumerate(padre) if value !=1][0],[index for index,value in enumerate(padre) if value !=1][1]
    elif sum(c2.Gene) <= len(c1.Gene)-2:
        padre = c2.Gene
        a,b = [index for index,value in enumerate(padre) if value !=1][0],[index for index,value in enumerate(padre) if value !=1][1]
    else:
        a,b = [random.choice(obj.KO) for i in range(2)][0],[random.choice(obj.KO) for i in range(2)][1]


    
    c1.Gene = pa.Gene[0:a] + pb.Gene[a:b] + pa.Gene[b:]
    c2.Gene = pb.Gene[0:a] + pa.Gene[a:b] + pb.Gene[b:]
    return c1,c2

    
def check_feasibility(p:Type[Individual],k:int) -> Genome:
    g = copy.deepcopy(p)
    # if sum(g.Gene) != len(g.Gene) - k:
    #     # g.Gene = 

    pass