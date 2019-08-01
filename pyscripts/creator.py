#!/usr/bin/env python3
import fileinput
import subprocess
from os import fdopen, remove
import os

#prots = ["2ci2I","2vik","1wit","1vcc","1tul","1lis","1iibA","1ctf","1cg5B", "1dhn"]
#prots = ["1gvp","1eyvA","1bkrA","1npsA","2chf","1tig","1tit"]
#prots = ["1c9oA", "1kpeA", "256bA", "1elwA"]
#prots = ["1c8cA", "256bA", "1hz6A"]
#prots = ["1acf", "1bgf", "1fna", "1opd", "1rnbA", "1ten", "1who"]
#prots = ["1c9oA", "1tit", "1tig", "2chf", "1npsA", "1bkrA", "1gvp", "1eyvA"]

# Repeticion ultimas mario cde
prots = ["1kpeA", "1wit", "1c9oA", "256bA", "1elwA", "1hz6A"]
# gecco 2015 prots
#prots = ["1tig", "1dtdB", "1sap", "2ci2I", "1c8cA", "1hz6A"]

#prots_faltan_rmsd_factor
#prots = ["1kpeA", "1elwA" , "1tit", "1cg5B", "1lis", "1npsA"]
#prots_faltan_rmsd_factor

#prots = ["1c9oA", "256bA" , "1fna", "1wit", "1ctf", "1eyvA"]

# Solo loops
#prots = ["1eyvA", "1gvp"]

options_name = "options_"

def process_prot(prot_name):
    copy_file_run = "cp run_base.sh run_"+prot_name+".sh"
    copy_file_ini = "cp "+options_name+"base.ini "+options_name+""+prot_name+".ini"
    rpl_prot = "sed -i.bak \\'/^\[Protocol]/,/^\[/{s/^prot*=.*/prot="+prot_name+"/}\\' "+options_name+""+prot_name+".ini"
    #print(rpl_prot)
    process = subprocess.Popen(copy_file_run.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    process = subprocess.Popen(copy_file_ini.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    
    file_in = open(str(options_name+"base.ini"), "r")
    file_out = open(str(options_name+""+prot_name+".ini"), "wt")
    for line in file_in:
        if "protname" in line:
            file_out.write(line.replace("prot=protname", "prot="+prot_name+""))
        else:
            file_out.write(line)
 
    file_in = open(str("run_base.sh"), "r")
    file_out = open(str("run_"+prot_name+".sh"), "wt")
    for line in file_in:
        if "protname" in line:
            if "srun" in line:
                file_out.write(line.replace("srun ./bin/mpi_app "+options_name+"protname.ini", "srun ./bin/mpi_app "+options_name+""+prot_name+".ini"))
            else:
                file_out.write(line.replace("protname", prot_name))
        else:
            file_out.write(line)
            
print("cd .. &")
for prot_name in prots:
    process_prot(prot_name)
    print("sbatch run_"+prot_name+".sh &")
    #print("echo \\"sbatch run_"+prot_name+".sh"\\")    
    #print()


copy_all = "cp run_*.sh ../ ; rm ../run_base.sh"
os.system(copy_all)
copy_all = "cp options_*ini ../ ; rm ../options_base.ini"
os.system(copy_all)

