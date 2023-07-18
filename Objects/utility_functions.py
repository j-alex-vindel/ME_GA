import time
from typing import List,NewType
import gurobipy as gp
from gurobipy import GRB
import copy
import random
import pandas as pd

FBA = NewType('FBA_vector',List[float])


def func_name_print(f):
    name = f.__name__
    def wrapper(*args,**kwargs):
        print(f"{' '*2}Running >> {name}")
        v = f(*args,**kwargs)
        return v
    return wrapper

def func_timer(f):
    def wrapper(*args,**kwargs):
        start = time.time()
        v = f(*args,**kwargs)
        total_time = time.time() - start
        print('>>>> Total running time of function:', total_time)
        return v
    return wrapper

def genecheck(f):
    def wrapper(*args,**kwargs):
        v = None
        while v == None:
            v,k = f(*args,**kwargs)
            if sum(v) == len(v) - k:
                v = v
            else:
                v = None
        return v
    return wrapper


def set_constructor(list:List[str]) -> List[int]:
    if list != None:
        return [i for i in range(len(list))]


def gene_generator(network=None,k:int=None):
    rk = random.sample(network.KO,k) 
    li = [0 if i in rk else 1 for i in network.M]
    return li

def gene2name(network=None,gene:List[int]=None) -> List[str]:

    strat = [network.Rxn[i] for i in network.M if gene[i]<.5]

    return strat

def save2df(population=None,network=None,write:bool=False):
    calc = {
        "Bio":[],
        "Che":[],
        "Strat":[],
        "Front":[],
        "Strain":[],
        "K":[]
    }
    for index,front in enumerate(population):
        for ele in front:
            bio = ele.biomass*-1
            che = ele.chemical*-1
            strat = gene2name(network=network,gene=ele.Gene)
            
            calc['Bio'].append(bio)
            calc['Che'].append(che)
            calc['Front'].append(index)
            calc['Strain'].append(network.Name[:3])
            calc['K'].append(len(strat))
            calc['Strat'].append(strat)
   
    df = pd.DataFrame.from_dict(calc)

    if write:
        df.to_csv(f"../Results/GA{network.Name[:3]}_P{len(df)}.csv")
        print(f"File saved!")
    
    print(df.head(5))
    return calc

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
    v_wts = [v_wt[i] for i in obj.M]
    
    wt.setObjective(1*v_wt[objective],GRB.MAXIMIZE)
    wt.addMConstr(obj.S,v_wts,'=',obj.b,name='Stoi')
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