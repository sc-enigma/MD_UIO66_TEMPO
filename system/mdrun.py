import os
import time as tm

def run_md(path):
    os.chdir(path)
    os.system('gmx_mpi mdrun -s *.tpr -v')
    while True: 
        tm.sleep(1)
        if 'confout.gro' in os.listdir('.'):
            break
    
path_list = ['']
for path in path_list:
    print(f'gmx_mpi mdrun path')
    run_md(path)
