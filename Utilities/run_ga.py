import time
import sys
import os

def run_all(*args,sleep:int=None,file:str=None):
    a1,a2,a3,a4,a5 = args
    print(f"Running: {file}.py {a1}")
    os.system(f"python ../Utilities/{file}.py {a1} {a2} {a3} {a4} {a5}")
    time.sleep(sleep)

if __name__ == "__main__":
    file = 'experiment'
    strains = ['iaf']
    tgts = [_ for _ in range(10,100,10)]
    towrite = 'y'
    
    gis = [(5,15),(10,20),(20,10),(30,15),(35,35)]

    for index,strain in enumerate(strains):
        for target in tgts:
            for numgen,numind in gis:
                run_all(strain,numgen,numind,towrite,target,sleep=2,file=file)


    # for index,strain in enumerate(strains):
    #     if index >=1:
    #         break
    #     else:
    #         run_all(strain,numgen,numind,towrite,sleep=50,file=file)
