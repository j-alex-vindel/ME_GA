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
LB[rxn.index('EX_glc__D_e')] = -10
UB[rxn.index('EX_glc__D_e')] = -10

exchange = ['EX_o2_e','EX_pi_e','EX_so4_e','EX_nh4_e']
lbounds = [0,-1000,-1000,-1000]
for index, name in enumerate(exchange):
    LB[rxn.index(name)] = lbounds[index]

secretion = ['EX_ac_e','EX_co2_e','EX_etoh_e','EX_for_e','EX_lac__D_e','EX_succ_e']

ubounds = [1000,1000,1000,1000,1000,1000]
for index, name in enumerate(secretion):
    UB[rxn.index(name)] = ubounds[index]

LB[rxn.index('GLCabcpp')] = -1000
LB[rxn.index('GLCptspp')] = -1000
UB[rxn.index('GLCabcpp')] = 1000
UB[rxn.index('GLCptspp')] = 1000

LB[rxn.index('GLCt2pp')] = 0
UB[rxn.index('GLCt2pp')] = 0

biomass = rxn.index('BIOMASS_Ec_iJO1366_core_53p95M')
chemical = rxn.index('EX_succ_e')

non_essentials = ['GLCabcpp', 'GLCptspp', 'HEX1', 'PGI', 'PFK', 'FBA', 'TPI', 'GAPD','PGK', 'PGM', 'ENO', 'PYK',
'LDH_D', 'PFL', 'ALCD2x', 'PTAr', 'ACKr','G6PDH2r', 'PGL', 'GND', 'RPI', 'RPE', 'TKT1', 'TALA', 'TKT2', 'FUM',
'FRD2', 'SUCOAS', 'AKGDH', 'ACONTa', 'ACONTb', 'ICDHyr', 'CS', 'MDH','MDH2', 'MDH3', 'ACALD']

ko = [rxn.index(i) for i in non_essentials]


#======================== Metabolic Network Object ==============================================================


MN_ijo1366 = Metabolic_Network(S=S,LB=LB,UB=UB,Rxn=rxn,Met=met,KO=ko,Name=bacteria,biomass=biomass,chemical=chemical) 