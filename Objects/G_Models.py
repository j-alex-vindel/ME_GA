import gurobipy as gp
from gurobipy import GRB
from typing import List, Tuple,Type
import copy


Vector = List[float]
Model = object
Y = List[int]
MN = object


# Add single level reformulation MILP OptKnock

def flux_balance_analysis(obj:MN) -> Vector:
    LB_WT = obj.LB.copy()
    UB_WT = obj.UB.copy()

    wt = gp.Model()

    v_wt = wt.addVars(obj.M,lb=-GRB.INFINITY, ub= GRB.INFINITY, vtype=GRB.CONTINUOUS, name='v_wt')

    wt.setObjective(1*v_wt[obj.biomass], GRB.MAXIMIZE)

    wt.addConstrs((gp.quicksum(obj.S[i,j]*v_wt[j] for j in obj.M) == 0 for i in obj.N), name='Swt')
    wt.addConstrs((LB_WT[j] <= v_wt[j] for j in obj.M), name='LBwt')
    wt.addConstrs((v_wt[j] <= UB_WT[j] for j in obj.M), name='UBwt')
    wt.Params.OptimalityTol = obj.infeas
    wt.Params.IntFeasTol = obj.infeas
    wt.Params.FeasibilityTol = obj.infeas
    wt.Params.OutputFlag = 0
    wt.optimize()
    if wt.status == GRB.OPTIMAL:
        wt_vs =  [wt.getVarByName('v_wt[%s]'%a).x for a in obj.M] 

    
    return wt_vs

def single_objective_model(obj:MN) -> Model:
    fba = gp.Model('Flux_Balance_Analysis')
    
    vi = fba.addVars(obj.M,lb=-GRB.INFINITY, ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='v')

    fba.setObjective(vi[obj.biomass], GRB.MAXIMIZE)
    fba.addConstrs((gp.quicksum(obj.S[i,j] *vi[j] for j in obj.M) == 0 for i in obj.N), name='S')

    fba.addConstr(vi[obj.biomass] >= obj.minprod, name='minprod')

    fba.update()
    fba._m1  = fba.copy()
    
    return fba._m1

def true_vij(obj:MN,y:Y,model:Model,M:List[int]) -> List[float]:
    model.setAttr('LB',model.getVars(),[obj.LB[j]*y[j] for j in M])
    model.setAttr('UB',model.getVars(),[obj.UB[j]*y[j] for j in M]) 
    
    model.Params.OptimalityTol = obj.infeas
    model.Params.IntFeasTol = obj.infeas
    model.Params.FeasibilityTol = obj.infeas

    model.optimize()

    if model.status == GRB.OPTIMAL:

        vij = [model.getVarByName('v[%s]'%a).x for a in M]
        # vchemical = vij[obj.chemical]

    elif model.status in (GRB.INFEASIBLE, GRB.UNBOUNDED,GRB.INF_OR_UNBD):
        vij = [-1000 for i in M]
    
    model.reset(1)
    return vij
    
    return 


def mobjective_model(obj:MN) -> Model:
    model = gp.Model('Single_Model')

    v = model.addVars(obj.M,lb=-GRB.INFINITY, ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='v')
    model.addConstrs((gp.quicksum(obj.S[i,j] *v[j] for j in obj.M) == 0 for i in obj.N), name='S')

    model.addConstr(v[obj.biomass] >= obj.minprod, name='minprod')

    model.update()
    model._m1  = model.copy()
    
    return model._m1

def mobjective_solve(obj:MN,y:Vector, basic_model:Model,M:List[int],priority:List[float],weight:List[float])->Model:
    
    v = basic_model.getVars()
    basic_model.setAttr('LB',v,[obj.LB[j]*y[j] for j in M])
    basic_model.setAttr('UB',v,[obj.UB[j]*y[j] for j in M])
    basic_model.ModelSense = GRB.MAXIMIZE
    basic_model.setObjectiveN(1*v[obj.chemical],0,priority[0],weight[0])
    basic_model.setObjectiveN(1*v[obj.biomass],1,priority[1],weight[1])

    basic_model.update()
    basic_model.Params.OutputFlag = 0
    basic_model.optimize()

    nObjectives = basic_model.NumObj
    if basic_model.Status == GRB.OPTIMAL:
        results = []
        for o in range(nObjectives):
            # Set which objective we will query
            basic_model.params.ObjNumber = o
            # Query the o-th objective value
            results.append(basic_model.ObjNVal)
    else:
        results = [-1000,-1000]

    basic_model.reset(1)
    return results

def function_optimize(y:Vector,model:Model,obj:MN,M:List[int],index:int) -> float:
    v = model.getVars()
    model.setAttr('LB',v,[obj.LB[j]*y[j] for j in M])
    model.setAttr('UB',v,[obj.UB[j]*y[j] for j in M]) 
    model.setObjective(v[index],GRB.MAXIMIZE)
    model.Params.OptimalityTol = obj.infeas
    model.Params.IntFeasTol = obj.infeas
    model.Params.FeasibilityTol = obj.infeas
    model.optimize()

    if model.status == GRB.OPTIMAL:
        obj = model.getObjective().getValue()

    elif model.status in (GRB.INFEASIBLE, GRB.UNBOUNDED,GRB.INF_OR_UNBD):
  
        obj = -1000 
    
    model.reset(1)
    return obj 


# ================================================== NOT TO BE USED ====================================================
# This functions shpuld not be used!!!!!!!!!! 
# def bi_objective_model(obj:MN,alpha:float) -> Model:

#     model = gp.Model('Multi_Objective')

#     v = model.addVars(obj.M,lb=-GRB.INFINITY, ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='v')

#     model.setObjective(v[obj.biomass] + v[obj.chemical],GRB.MAXIMIZE)

#     model.addConstrs((gp.quicksum(obj.S[i,j] *v[j] for j in obj.M) == 0 for i in obj.N), name='S')

#     model.addConstr(v[obj.biomass] >= obj.minprod, name='minprod')

#     model.update()
#     model._m1  = model.copy()
    
#     return model._m1


def multiobjective_function(y:Vector,obj:MN, model:Model,M:List[int]) -> Tuple[int,Vector]:
    
    model.setAttr('LB',model.getVars(),[obj.LB[j]*y[j] for j in M])
    model.setAttr('UB',model.getVars(),[obj.UB[j]*y[j] for j in M]) 
    
    model.Params.OptimalityTol = obj.infeas
    model.Params.IntFeasTol = obj.infeas
    model.Params.FeasibilityTol = obj.infeas

    model.optimize()

    if model.status == GRB.OPTIMAL:
        obj = model.getObjective().getValue()
        vij = [model.getVarByName('v[%s]'%a).x for a in M]
    elif model.status in (GRB.INFEASIBLE, GRB.UNBOUNDED,GRB.INF_OR_UNBD):
        vij = [-1000 for i in range(len(obj.M))]
        obj = -1000 
    
    model.reset(1)
    return obj, vij

#