from G_Ind import GA_Ind
from G_Models import MILP_Model, MILP_Solve
from utility_functions import gene_generator

class GA_BLP:

    def __init__(self,MetNet,K):
        self.metnet = MetNet
        self.K = K
        self.model = MILP_Model(self.metnet)

    def generate_individual(self):
        individual = GA_Ind()
        individual.Gene = gene_generator(self.metnet,self.K)
        return individual
    
    def calc_obj(self,individual):
        individual.biomass, individual.chemical = MILP_Solve(network=self.metnet,model=self.model,y=individual.Gene)


