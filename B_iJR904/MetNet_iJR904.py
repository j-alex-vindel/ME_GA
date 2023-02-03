import sys
import os

sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Objects')))

from Met_Net import Metabolic_Network
from pymatreader import read_mat

rootname = sys.argv[0]
bacteria = rootname[7:len(rootname)-3]

#========================= Access Data ========================================================================= 
data = read_mat(f"../B_{bacteria}/Data/{bacteria}.mat")

met = data[bacteria]['mets']
rxn = data[bacteria]['rxns']
LB  = data[bacteria]['lb'].tolist()
UB  = data[bacteria]['ub'].tolist()
S   = data[bacteria]['S']


#======================== Biological Assumptions ===============================================================
#======================== Change accordingly ===================================================================
# Comes from iJR904_under_new_ys_reack_myalgo.ipynb
# prespecified amount of glucose uptake 10 mmol/grDW*hr 'EX_glc__D_e' = -10 reaction 

LB[rxn.index('EX_glc__D_e')] = -10
UB[rxn.index('EX_glc__D_e')] = -10

# Unconstrained uptake routes for inorganic phosphate, sulfate and ammonia 
# 'EX_o2_e';'EX_pi_e';'EX_so4_e'; 'EX_nh4_e' index = 184 ; 199 ; 259 ; 169

LB[rxn.index('EX_o2_e')] = 0
LB[rxn.index('EX_pi_e')] = -1000
LB[rxn.index('EX_so4_e')] = -1000 
LB[rxn.index('EX_nh4_e')] = -1000

#Enable secretion routes for acetate, carbon dioxide, ethanol, formate, lactate and succinate
# {'EX_ac_e';'EX_co2_e';'EX_etoh_e';'EX_for_e';'EX_lac__D_e';'EX_succ_e'} change in the upper bound 
# index = 87; 2; 340; 422; 91; 261

sec_routes = ['EX_ac_e','EX_co2_e','EX_etoh_e','EX_for_e','EX_lac__D_e','EX_succ_e']

for i in sec_routes:
    UB[(rxn.index(i))] = 1000

biomass = rxn.index('BIOMASS_Ecoli')
chemical = rxn.index('EX_succ_e')

non_essentials = ['HEX1', 'PGI', 'PFK', 'FBA', 'TPI', 'GAPD','PGK', 'PGM', 'ENO', 'PYK',
'LDH_D', 'PFL','PTAr', 'ACKr','G6PDH2r', 'PGL', 'GND', 'RPI', 'RPE', 'TKT1', 'TALA', 'TKT2', 'FUM',
'FRD2', 'SUCOAS', 'AKGDH','ICDHyr', 'CS', 'MDH','MDH2', 'MDH3', 'ACALD']

ko = [rxn.index(i) for i in non_essentials]



#======================== Metabolic Network Object ==============================================================


MN_ijr904 = Metabolic_Network(S=S,LB=LB,UB=UB,Rxn=rxn,Met=met,KO=ko,Name=bacteria,biomass=biomass,chemical=chemical) 