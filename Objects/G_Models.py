import gurobipy as gp
from gurobipy import GRB
from typing import List, Tuple,Type
from collections import namedtuple
from Met_Net import Met_Net
import copy


Vector = List[float]
Model = Type[object]
Y = List[int]
MN = Type[Met_Net]
Result = namedtuple('Result',['MetNet','Strategy','Vs','Time','Soltype'])
RMILP = namedtuple('RMILP',['Cost','Chem','Biom'])
K = Type[int]


def MILP_Model(network:MN=None) -> Model:
    '''
    Returns a model object ready to be solved
    ''' 
    m = gp.Model("Single_Level_Reformulation")
    #Variables
    v = m.addVars(network.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='v')
    # Dual Variables
    l = m.addVars(network.N,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='l')
    a1 = m.addVars(network.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='a1')
    b1 = m.addVars(network.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='b1')
    a2 = m.addVars(network.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='a2')
    b2 = m.addVars(network.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='b2')
    a = m.addVars(network.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='a')
    b = m.addVars(network.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='b')
    
    # Objective
    m.setObjective((1*v[network.chemical]),GRB.MAXIMIZE)
    # Constraints
    # Stoichimetric Constrs
    m.addMConstr(network.S,v,'=',network.b,name='S')
    
    # Dual Objective
    m.addConstr((v[network.biomass] >= (sum(a1[j]*network.UB[j] - b1[j]*network.LB[j] for j in network.M)
     + sum(a2[j]*network.UB[j] - b2[j]*network.LB[j] for j in network.M))),name='DO')
    
    # Dual Constraints
    m.addConstrs((gp.quicksum(network.S.transpose()[i,j]*l[j] for j in network.N)
              - b[i]
              + a[i] - b2[i] + a2[i]
               == network.c[i] for i in network.M)
             ,name='SD')

    m.addConstr(v[network.biomass] >= network.minprod, name='min_growth')

    m.update()

    model = m.copy()

    return model

def MILP_Solve(network:MN=None,y:Y=None,model:Model=None) -> Result:
    '''
    Takes a y vector to solve the milp model
    '''
    # Retrieve variables from model
    v = [model.getVarByName('v[%s]'%a) for a in network.M]
    l = [model.getVarByName('l[%s]'%a) for a in network.N]
    a = [model.getVarByName('a[%s]'%a) for a in network.M]
    a1 = [model.getVarByName('a1[%s]'%a) for a in network.M]
    a2 = [model.getVarByName('a2[%s]'%a) for a in network.M]
    b = [model.getVarByName('b[%s]'%a) for a in network.M]
    b1 = [model.getVarByName('b1[%s]'%a) for a in network.M]
    b2 = [model.getVarByName('b2[%s]'%a) for a in network.M]
    
    # Add Linerarization Constraints
    # Linearization
    model.addConstrs((a1[j] <= network.BM*y[j] for j in network.M),name='l1_a1')

    model.addConstrs((a1[j] >= - network.BM*y[j] for j in network.M),name='l2_a1')

    model.addConstrs((a1[j] <= a[j] + network.BM*(1-y[j]) for j in network.M),name='l3_a1')

    model.addConstrs((a1[j] >= a[j] - network.BM*(1-y[j]) for j in network.M),name='l4_a1')

    model.addConstrs((b1[j] <= network.BM*y[j] for j in network.M),name='l1_b1')

    model.addConstrs((b1[j] >= -network.BM*y[j] for j in network.M),name='l2_b1')

    model.addConstrs((b1[j] <= b[j] + network.BM*(1-y[j]) for j in network.M),name='l3_b1')

    model.addConstrs((b1[j] >= b[j] - network.BM*(1-y[j]) for j in network.M),name='l4_b1')

    # Bounds
    model.addConstrs((network.LB[j]*y[j] <= v[j] for j in network.M), name='LB')
    model.addConstrs((v[j] <= network.UB[j]*y[j] for j in network.M), name='UB')

    model.addConstrs((network.LB[j] <= v[j] for j in network.M),name='lb')
    model.addConstrs((v[j] <= network.UB[j] for j in network.M),name='ub')

    # Add Parameters
    model.Params.OptimalityTol = network.infeas
    model.Params.IntFeasTol = network.infeas
    model.Params.FeasibilityTol = network.infeas
    model.Params.NodefileStart = 0.5
    # model.Params.Presolve = 0
    model.Params.OutputFlag = 0
    model.update()
    model.optimize()

    # Update
    # Solve
    if model.status == GRB.OPTIMAL:
        chem = model.getObjective().getValue()
        vs = [model.getVarByName('v[%d]'%j).x for j in network.M]

    elif model.status == GRB.TIME_LIMIT:
        vs = [model.getVarByName('v[%d]'%j).x for j in network.M]
    
    if model.status in (GRB.INFEASIBLE,GRB.INF_OR_UNBD,GRB.UNBOUNDED):
        # print('Model status: *** INFEASIBLE or UNBOUNDED ***')
        vs = [-200 for i in network.M]
    
    cost = vs[network.chemical] + vs[network.biomass]
    chem = vs[network.chemical]
    biom = vs[network.biomass]
    
    return  RMILP(cost,chem,biom)
    
    
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


# def multiobjective_function(y:Vector,obj:MN, model:Model,M:List[int]) -> Tuple[int,Vector]:
    
#     model.setAttr('LB',model.getVars(),[obj.LB[j]*y[j] for j in M])
#     model.setAttr('UB',model.getVars(),[obj.UB[j]*y[j] for j in M]) 
    
#     model.Params.OptimalityTol = obj.infeas
#     model.Params.IntFeasTol = obj.infeas
#     model.Params.FeasibilityTol = obj.infeas

#     model.optimize()

#     if model.status == GRB.OPTIMAL:
#         obj = model.getObjective().getValue()
#         vij = [model.getVarByName('v[%s]'%a).x for a in M]
#     elif model.status in (GRB.INFEASIBLE, GRB.UNBOUNDED,GRB.INF_OR_UNBD):
#         vij = [-1000 for i in range(len(obj.M))]
#         obj = -1000 
    
#     model.reset(1)
#     return obj, vij

#  def MILP_refor(network:Model=None,k:K=None,log:bool=True,speed:bool=False,threads:bool=False) -> Result:
#     '''
#     MILP_solve(network=network,k=k)
#         return Result[MetNet,Strategy,Flows,Time,Soltype]
#     '''
#     print(f'**** Solving ReacKnock k={k} ****')
#     print(f'# Variables (reactions in the network): {len(network.M)}')
#     print('Current Infeasibility:',network.infeas,sep=' -> ')
#     print('KO set: ',len(network.KO), ' reactions')
#     print(f"MN: {network.Name}")

#     lb = copy.deepcopy(network.LB)
#     lb[network.biomass] = network.minprod

#     m = gp.Model()

#     # Variables
#     v = m.addVars(network.M,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='v')
#     y = m.addVars(network.M,vtype=GRB.BINARY,name='y')

#     # Dual Variables
#     l = m.addVars(network.N,lb=-GRB.INFINITY,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='lambda')
#     a1 = m.addVars(network.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='alpha1')
#     b1 = m.addVars(network.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='beta1')
#     a2 = m.addVars(network.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='alpha1')
#     b2 = m.addVars(network.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='beta1')
#     a = m.addVars(network.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='alpha')
#     b = m.addVars(network.M,lb=0,ub=GRB.INFINITY,vtype=GRB.CONTINUOUS,name='beta')
    
#     # Objective
#     m.setObjective((1*v[network.chemical]),GRB.MAXIMIZE)

#     # Knapsack Constrs
#     m.addConstrs((y[j] == 1 for j in network.M if j not in network.KO), name='y_essentials')

#     m.addConstr(sum(1-y[j] for j in network.KO) == k, name='knapsack')

#     # Stoichimetric Constrs
#     m.addMConstr(network.S,v,'=',network.b,name='Stoi')
#     # m.addConstrs((gp.quicksum(network.S[i,j] * v[j] for j in network.M) == network.b[i] for i in network.N),name='Stoichiometry')
    
#     # Dual Objective
#     m.addConstr((v[network.biomass] >= (sum(a1[j]*network.UB[j] - b1[j]*lb[j] for j in network.M)
#      + sum(a2[j]*network.UB[j] - b2[j]*lb[j] for j in network.M))),name='dual-objective')
    
#     # Dual Constraints
#     m.addConstrs((gp.quicksum(network.S.transpose()[i,j]*l[j] for j in network.N)
#               - b[i]
#               + a[i] - b2[i] + a2[i]
#                == network.c[i] for i in network.M)
#              ,name='S_dual')

#     # m.addConstr((gp.quicksum(network.S.transpose()[network.biomas,j]*l[j] for j in network.N)
#     #         - b[network.biomas]
#     #         + a[network.biomas]
#     #         - b2[network.biomas] + a2[network.biomas] == 1), name='Sdual_t')
    
#     # Linearization
#     m.addConstrs((a1[j] <= network.BM*y[j] for j in network.M),name='l1_a1')

#     m.addConstrs((a1[j] >= - network.BM*y[j] for j in network.M),name='l2_a1')

#     m.addConstrs((a1[j] <= a[j] + network.BM*(1-y[j]) for j in network.M),name='l3_a1')

#     m.addConstrs((a1[j] >= a[j] - network.BM*(1-y[j]) for j in network.M),name='l4_a1')

#     m.addConstrs((b1[j] <= network.BM*y[j] for j in network.M),name='l1_b1')

#     m.addConstrs((b1[j] >= -network.BM*y[j] for j in network.M),name='l2_b1')

#     m.addConstrs((b1[j] <= b[j] + network.BM*(1-y[j]) for j in network.M),name='l3_b1')

#     m.addConstrs((b1[j] >= b[j] - network.BM*(1-y[j]) for j in network.M),name='l4_b1')

#     # Bounds
#     m.addConstrs((lb[j]*y[j] <= v[j] for j in network.M), name='LB')
#     m.addConstrs((v[j] <= network.UB[j]*y[j] for j in network.M), name='UB')

#     m.addConstrs((lb[j] <= v[j] for j in network.M),name='lb')
#     m.addConstrs((v[j] <= network.UB[j] for j in network.M),name='ub')

#     m.Params.OptimalityTol = network.infeas
#     m.Params.IntFeasTol = network.infeas
#     m.Params.FeasibilityTol = network.infeas
#     # m.Params.NodefileStart = 0.5
#     m.Params.Presolve = 0
#     if not log: m.Params.OutputFlag = 0
#     if speed: m.Params.NodefileStart = 0.5
#     if threads: m.Threads = 2
#     m.optimize()

#     s = m.Runtime
#     if m.status == GRB.OPTIMAL:
#         chem = m.getObjective().getValue()
#         ys = [m.getVarByName('y[%d]'%j).x for j in network.M]
#         vs = [m.getVarByName('v[%d]'%j).x for j in network.M]
#         soltype = 'Optimal'
#         del_strat = [network.Rxn[i] for i in network.M if ys[i] <.5]

#     elif m.status == GRB.TIME_LIMIT:
#         ys = [m.getVarByName('my[%d]'%j).x for j in network.M]
#         vs = [m.getVarByName('mv[%d]'%j).x for j in network.M]
#         del_strat = [network.Rxn[i] for i in network.M if ys[i] <.5]
#         soltype = 'Time_limit'
    
#     if m.status in (GRB.INFEASIBLE,GRB.INF_OR_UNBD,GRB.UNBOUNDED):
#         # print('Model status: *** INFEASIBLE or UNBOUNDED ***')
#         ys = ['$' for i in network.M]
#         vs = ['~' for i in network.M]
#         # print('Chemical:',vs[network.chemical],sep=' ^ ')
#         # print('Biomass:',vs[network.biomass],sep=' ^ ')
#         del_strat = 'all'

#     print('*'*4,' FINISHED!!! ','*'*4)

#     return  Result(network.Name,del_strat,ys, vs, s,soltype)