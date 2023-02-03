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


chemical = str
biomass = str
ko = []

#======================== Metabolic Network Object ==============================================================


MN = Metabolic_Network(S=S,LB=LB,UB=UB,Rxn=rxn,Met=met,KO=ko,Name=bacteria,biomass=biomass,chemical=chemical) 