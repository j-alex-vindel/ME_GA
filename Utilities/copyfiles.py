import sys
import os
import time
import shutil


# Modify so it can take the file's new name

def copy_file_chname(cfile:str=None,bacteria:list=None,nfiles:int=None,t:int=2):

    path = os.getcwd()
    parent_directory = os.path.abspath(os.path.join(path,os.pardir))
    print(f"Parent Dir: <{parent_directory}>")
    sfile = f"Utilities/{cfile}"
    source = os.path.abspath(os.path.join(parent_directory,sfile))

    folder = f"B_{bacteria}"
    dest_folder = os.path.join(parent_directory,folder)
    cursizedir = len(os.listdir(dest_folder))
    print(f"Destiny Dir: {dest_folder}")

    ks = [i for i in range(nfiles)]
    
    print(f'From {source} ->')

    for k in ks:
        if k == 0:
            file = f"{cfile[:7]}{bacteria}.py"
        else:
            file = f"{cfile[:7]}{bacteria}_{k}.py"
        
        dest_file = os.path.join(dest_folder,file)
        print(dest_file)
        shutil.copyfile(source,dest_file)

    finsizedir = len(os.listdir(dest_folder))

    print(f"# {finsizedir-cursizedir} File(s) Added!!")

    time.sleep(t)

if __name__ == "__main__":

    file_to_copy = sys.argv[1]
    bacteria = sys.argv[2]
    n_copies = int(sys.argv[3])
    copy_file_chname(cfile=file_to_copy,bacteria=bacteria,nfiles=n_copies)