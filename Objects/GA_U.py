
import copy
import random
from typing import List,Type,Tuple
from G_Models import MILP_Model, MILP_Solve
from Individual import Individual
from Met_Net import Metabolic_Network

MN = Type[Metabolic_Network]
Genome = List[int]
Model = Type[MILP_Model]


def GA_Run(Network:MN=None,Npop:int=None,MIModel:Model=None,K:int=None):
    # Create population
    pop = [Individual() for _ in range(Npop)]
    
    # Initial Solutions
    for i in range(Npop):
        pop[i].Gene = Gene_Gen(Network=Network,K=K)
        pop[i].Result = MILP_Solve(network=Network,y=pop[i].Gene,model=MIModel)
        pop[i].Cost = pop[i].Result.Cost
        pop[i].Bio = pop[i].Result.Biom
        pop[i].Chem = pop[i].Result.Chem

    # s_pop = sorted(pop,key=lambda x: x.Cost,reverse=True)

    
    
    pass

def Gene_Gen(Network:MN=None,K:int=None) -> Genome:
    rk = [random.choice(Network.KO) for i in range(K)]
    li = [0 if i in rk else 1 for i in Network.M]
    return li

def Uni_Cross(Pa:Type[Individual]=None,Pb:Type[Individual]=None) ->Tuple[Type[Individual],Type[Individual]]:
    lenght = len(Pa.Gene)
    c1 = copy.deepcopy(Pa)
    c2 = copy.deepcopy(Pb)

    uni_alpha = [random.randint(0,1) for i in range(lenght)]
    c1.Gene = [uni_alpha[i]*Pa.Gene[i] + (1-uni_alpha[i])*Pb.Gene[i] for i in range(lenght)]
    c2.Gene = [(1-uni_alpha[i]*Pa.Gene[i]) + uni_alpha[i]*Pb.Gene[i] for i in range(lenght)]
    return c1,c2
