import time
from typing import List,NewType
import gurobipy as gp
from gurobipy import GRB
import copy


FBA = NewType('FBA_vector',List[float])


def set_constructor(list:List[str]) -> List[int]:
    if list != None:
        return [i for i in range(len(list))]


def func_timer(f):
    def wrapper(*args,**kwargs):
        start = time.time()
        v = f(*args,**kwargs)
        total_time = time.time() - start
        print('>>>> Total running time of function:', total_time)
        return v
    return wrapper

def wildtype_FBA(obj,wildtype:bool=True,mutant:bool=False)->FBA:
    LB_wt = copy.deepcopy(obj.LB)
    UB_wt = copy.deepcopy(obj.UB)

    if wildtype and not mutant:
        objective = obj.biomass
        FVA = False
    elif mutant and not wildtype:
        objective = obj.chemical
        FVA = True
    else:
        raise Exception("Both wildtype and mutant cannot be True or False at the same time")

    wt = gp.Model()
    v_wt = wt.addVars(obj.M, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='v_wt')
    
    wt.setObjective(1*v_wt[objective],GRB.MAXIMIZE)
    wt.addMConstr(obj.S,v_wt,'=',obj.b,name='Stoi')
    wt.addConstrs((LB_wt[j] <= v_wt[j] for j in obj.M), name='LBwt')
    wt.addConstrs((UB_wt[j] >= v_wt[j] for j in obj.M), name='UBwt')
    if FVA:
        wt.addConstr((v_wt[obj.biomass] >= obj.minprod), name='minprod')

    wt.Params.OptimalityTol = obj.infeas
    wt.Params.IntFeasTol = obj.infeas
    wt.Params.FeasibilityTol = obj.infeas
    wt.Params.OutputFlag = 0
    wt.optimize()
    if wt.status == GRB.OPTIMAL:
        wt_vs =  [wt.getVarByName('v_wt[%s]'%a).x for a in obj.M] 
    else:
        wt_vs = [-2000 for _ in obj.M]

    return wt_vs