import time
import sys
import os

def run_all(*args,sleep:int=None,file:str=None):
    a1,a2,a3,a4 = args
    os.system(f"python ../Utilities/{file}.py {a1} {a2} {a3} {a4}")
    time.sleep(sleep)

if __name__ == "__main__":
    file = 'experiment'
    strain = 'ijr'
    numgen = 6
    numind = 5
    towrite = 'y'
    run_all(strain,numgen,numind,towrite,sleep=5,file=file,)
