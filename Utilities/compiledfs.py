import pandas as pd
from os.path import exists



def dfextract(metnet:str=None,generations:int=None,ind:int=None,tgt:int=None):
    
    file = f"../Results/GA{metnet}_M_G{generations}_I{ind}_{tgt}.0.csv"
    
    file_exists = exists(file)
    if file_exists:
        print(f"Reading: \n {file}")
        d = pd.read_csv(f"../Results/GA{metnet}_M_G{generations}_I{ind}_{tgt}.0.csv")
        return d

strains = ['iAF']
tgts = [_ for _ in range(10,100,10)]

gis = [(5,15),(10,20),(20,10),(30,15),(35,35)]

ds = []
for strain in strains:
    for tgt in tgts:
        for numgen,numind in gis:
            df = dfextract(metnet=strain,generations=numgen,ind=numind,tgt=tgt)
            if isinstance(df,pd.DataFrame):
                ds.append(df)

masterdf = pd.concat(ds,ignore_index=True,sort=False)

masterdf.to_csv(f"../Results/MGA_M_{strains[0]}.csv")