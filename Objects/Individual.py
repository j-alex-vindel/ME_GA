from dataclasses import dataclass
from typing import List, Any, Tuple
from collections import namedtuple
from dataclasses_json import dataclass_json


Genome = List[int]

@dataclass
class Individual:
    Gene: Genome = None
    Cost: float = None
    Fitness:namedtuple =  None
    Bio: float = None
    Chem: float = None
    TrueBio: float = None
    TrueChem: float = None
    Strat: List[str] = None

Population = List[Individual]

@dataclass_json
@dataclass
class Out_Result:
    Pop: Population = None
    BestSol: Any = None
    BestCost: float = None 
    