import sys
import os
import time
import shutil


# Modify so it can take the file's new name

def copy_file_chname(cfile:str=None,bacteria:list=None,nfiles:int=None,t:int=None):

    path = os.getcwd()
    parent_directory = os.path.abspath(os.path.join(path,os.pardir))
    print(f"Parent Directory <{parent_directory}>")
    sfile = f"Utilities/{cfile}"
    source = os.path.abspath(os.path.join(parent_directory,sfile))

    folder = f"B_{bacteria}"
    dest_folder = os.path.join(parent_directory,folder)
    cursizedir = len(os.listdir(dest_folder))
    print(dest_folder)

    ks = [i+1 for i in range(nfiles)]
    
    print(f'From {source} ->')

    for k in ks:
        file = f"run_experiment_p_k{k}.py"
        dest_file = os.path.join(dest_folder,file)
        print(dest_file)
        shutil.copyfile(source,dest_file)

    finsizedir = len(os.listdir(dest_folder))

    print(f"# {finsizedir-cursizedir} Files Added!!")