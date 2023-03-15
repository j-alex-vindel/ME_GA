from G_Ind import GA_Ind
from G_Models import MILP_Model, MILP_Solve
from utility_functions import gene_generator,func_name_print
from typing import Type
from Met_Net import Met_Net

MetabolicNetwork = Type[Met_Net]
class GA_BLP:

    def __init__(self,MetNet:MetabolicNetwork=None,K:int=None):
        self.metnet = MetNet
        self.K = K
        self.model = MILP_Model(self.metnet)
    

    def generate_individual(self):
        individual = GA_Ind()
        individual.Gene = gene_generator(self.metnet,self.K)
        return individual
    
 
    def calc_obj(self,individual):
        individual.biomass, individual.chemical = MILP_Solve(network=self.metnet,model=self.model,y=individual.Gene)


