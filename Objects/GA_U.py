
import random
from typing import List,Type
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
        # pop[i].Fitness = MILP_Solve(network=MN,y=pop[i].Gene,model=MIModel)
    print(pop)
    
    pass

def Gene_Gen(Network:MN=None,K:int=None) -> Genome:
    rk = [random.choice(Network.KO) for i in range(K)]
    li = [0 if i in rk else 1 for i in Network.M]
    return li

