from dataclasses import dataclass
from typing import List, Any, Tuple, Type
from collections import namedtuple
from dataclasses_json import dataclass_json


Genome = List[int]
OBJ = namedtuple('OBJ',['chemical','biomass'])

class GA_Ind(object):
    
    def __init__(self):

        self.Gene: Genome = None
        self.rank:int = None
        self.crowding_distance:float = None
        self.dominated_solutions:Type[List] = None
        self.dominated_count:int = None
        self.biomass:float = None
        self.chemical:float = None
        self.features = None
    
    @property
    def objectives(self):
        self._objectives = (self.chemical,self.biomass)
        return self._objectives
    @property
    def cost(self):
        self._cost = self.chemical + self.biomass
        return self._cost  

    @property
    def gene_index(self):
        if self.Gene != None:
            self._index = [i for i,gene in enumerate(self.Gene) if gene == 0]
        else:
            return None
        return self._index

    def __eq__(self, other):
        if isinstance(self,other.__class__):
            return self.features == other.features
        return False

    def dominates(self, other_individual):
        and_condition = True
        or_condition = False
        
        for first, second in zip(self.objectives, other_individual.objectives):
            and_condition = and_condition and first <= second
            or_condition = or_condition or first < second
        
        return (and_condition and or_condition)