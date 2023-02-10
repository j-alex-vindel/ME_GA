from dataclasses import dataclass
from typing import List, Any, Tuple, Type
from collections import namedtuple
from dataclasses_json import dataclass_json


Genome = List[int]


class Individual(object):
    
    def __init__(self):

        self.Gene: Genome = None
        self.Cost: float = None
        self.Biom: float = None
        self.Chem: float = None
        self.rank:int = None
        self.crowding_distance:float = None
        self.dominated_solutions:Type[List] = None
        self.Obj:Type[List] = None

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