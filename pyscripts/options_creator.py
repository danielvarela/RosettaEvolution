#!/usr/bin/env python3

import sys
import pandas as pd
import os
import glob
#import matplotlib.pyplot as plt
#import cluster_graph as cg
import numpy as np
import pandas as pd
from sklearn.preprocessing import normalize
from operator import itemgetter
import sys
from sklearn import preprocessing
from sklearn.preprocessing import minmax_scale
from sklearn.preprocessing import MinMaxScaler
import click
from stat import S_ISREG, ST_CTIME, ST_MODE
import time


options_path = "/mnt/netapp2/Store_uni/home/ulc/co/dvm/Code/minifolding/options_dir/"
bash_path = "/mnt/netapp2/Store_uni/home/ulc/co/dvm/Code/minifolding/"

def print_file(f, p, r, d):
    f.writelines("[Protocol]\n")
    f.writelines("prot=%s\n" % p )
    f.writelines("name=MPIMoverDE\n")
    f.writelines("stages=stage2,stage3,stage4\n")
    f.writelines("distance_strategy=euclidean_partial_mario\n")
    f.writelines("init_strategy=init_popul_with_stage\n")
    f.writelines("crossover_strategy=CR\n")
    f.writelines("mutation_strategy=default\n")
    f.writelines("select_parents_strategy=rand\n")
    f.writelines("clean_gen_limit=100000\n")
    f.writelines("fragments_popul_module=10\n")
    f.writelines("map_res=%s\n" % r)
    f.writelines("disturb_degrees=%s\n" % d)
    f.writelines("[Fragments]\n")
    f.writelines("frags_at_popul=true\n")
    f.writelines("strategy_at_population=stage_rosetta_mover\n")
    f.writelines("strategy_at_trials=stage_rosetta_mover\n")
    f.writelines("[PrintPopul]\n")
    f.writelines("output=temp\n")
    f.writelines("gens=-1\n")
    f.writelines("[DE]\n")
    f.writelines("CR=0.99\n")
    f.writelines("F=0.03\n")
    f.writelines("NP=96\n")
    f.writelines("[Extra]\n")
    f.writelines("fitrad=96\n")

def create_files_for(prots, res, dist_deg):
    F_tag = "F003"
    t = len(res) + len(dist_deg)
    created_data = []
    for p in prots:
        for r in res:
            for d in dist_deg:
                file_name = "NewFrags_HybridDE_"+F_tag+"_my_options_res"+str(r)+"_deg"+str(d)+"_"+p+".ini"
                f = open(options_path+file_name,"w")
                f.write("# file "+file_name+"\n")
                f.write("# generated with options_creator.py\n")
                print_file(f, p, r, d)
                d = {"filename" : file_name, "prot": p, "res": r, "deg" : d}
                created_data.append(d)
                f.close()

    return created_data


def print_data_bash_file(f, p, r, d, options_file):
    F_tag = "F003"
    #destiny_exec = {"1cg5B" : "normal" , "1eyvA" :"normal" , "1lis" : "normal"
    #, "1kpeA" : "normal" , "1rnbA" : "shared" , "1tig" : "shared" , "1dhn"
    #                :"shared" , "1bgf" : "shared"}
    # destiny_exec = {"1c8cA" : "normal",  "1elwA" : "normal",  "1c9oA" :
    # "normal",  "256bA" : "normal",  "1ten" : "normal", "1opd" : "normal",
    # "1who" : "normal",  "1fna" : "normal",  "1bgf" : "normal", "1kpeA" :
    # "normal", "1hz6A" : "normal", "1tit" : "normal", "2chf" : "normal", "1npsA"
    # : "normal", "1acf" : "normal", "1bkrA" : "normal", "2ci2I" : "normal",
    # "2vik" : "normal", "1wit" : "normal", "1vcc" : "normal", "1tul" : "normal",
    #                 "1iibA" : "normal", "1ctf" : "normal", "1gvp" : "normal",
    # "1bk2" : "normal"}
    
    destiny_exec = {"n1c8cA" : "shared", "n1elwA":"shared", "n1wit":"shared",
                    "n1kpeA":"shared", "n1opd" :"shared", "n1hz6A" : "shared"}


    # destiny_exec = {"1c8cA" : "normal",  "1elwA" : "normal",  "1c9oA" :
    # "normal",  "256bA" : "normal",  "1ten" : "normal", "1opd" : "normal",
    # "1who" : "normal",  "1fna" : "normal",  "1bgf" : "normal", "1kpeA" :
    # "normal", "1hz6A" : "normal", "1tit" : "normal", "2chf" : "normal", "1npsA"
    # : "normal", "1acf" : "normal", "1bkrA" : "normal", "2ci2I" : "normal",
    # "2vik" : "normal", "1wit" : "normal", "1vcc" : "normal", "1tul" : "normal",
    #                 "1iibA" : "normal", "1ctf" : "normal", "1gvp" : "normal"}


    #destiny_exec = {"1c9oA" : "shared", "1wit" : "shared", "256bA" : "shared",
                    # "dimaio" : "shared", "1elwA" : "shared", "1kpeA" :
                    # "shared", "1hz6A" : "shared"}
 
    #job_name = "DoubleOutDofCryoDE_res"+str(r)+"_dist"+str(d)+"_"+p
    job_name = "NewFrags_HybridDE_"+F_tag+"_"+p
    f.writelines("#!/bin/bash\n")
    f.writelines("#SBATCH -n 12\n")
    f.writelines("#SBATCH --job-name=%s \n" % job_name)
    f.writelines("#SBATCH -a 1-10\n")
    f.writelines("#SBATCH -t 3:40:00\n")
    if (destiny_exec[p] == "shared"):
        f.writelines("#SBATCH --qos=shared\n")
        f.writelines("#SBATCH -p shared\n")
    else:
        f.writelines("#SBATCH -p cola-corta,thinnodes\n")
    f.writelines("#SBATCH -o \"nmpi-%x-%A_%a.out\"\n")
    f.writelines("#SBATCH -e \"error-%x.out\"\n")
    f.writelines("#SBATCH --mem-per-cpu=1GB\n")
    f.writelines("module load intel\n")
    f.writelines("module load impi\n")
    f.writelines("module load boost/1.68.0-python-2.7.15\n")
    f.writelines("export PRINT_ALL=true\n")
    if (destiny_exec[p] == "shared"):
        f.writelines("unset I_MPI_DEVICE\n")
        f.writelines("unset I_MPI_FALLBACK_DEVICE\n")
        f.writelines("export I_MPI_FALLBACK=0\n")
        f.writelines("export I_MPI_FABRICS=shm:tcp\n")
        
    f.writelines("export MPI_ENABLED=true\n")
    exec_options_file = options_path+options_file
    f.writelines("srun ./bin/mpi_app %s\n" % exec_options_file)
    
def create_bash_files(created_data):
    files_list = []
    for e in created_data:
        file_name = "my_bash_run_res"+str(e["res"])+"_deg"+str(e["deg"])+"_"+e["prot"]+".sh"
        f = open(bash_path+file_name,"w")
        print_data_bash_file(f, e["prot"], e["res"], e["deg"], e["filename"])
        f.close()
        files_list.append(file_name)
    return files_list

def make_folder(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def delete_prot_from_name(input_name):
    prots = ["1c9oA", "256bA", "1wit", "dimaio"]
    job_name = input_name
    find_prot = prots[0]
    for p in prots:
        if p in input_name:
            job_name = input_name.replace("_"+p, "")
            find_prot = p
            #print(job_name)
            break

    return (job_name, find_prot)

def found_pdbs(job_name):
    values = job_name.split("_")[1:]
    #print(values)
    res_value = [i for i in values if "res" in i][0]
    deg_value = [i for i in values if "dist" in i or "deg" in i][0]
    deg_value = deg_value.replace("dist","deg")
    file_names = "cart_pose_*"+deg_value+"*"+res_value+"*.pdb"
    #print(file_names)
    files = glob.glob(file_names)
    #print(len(files))
    return len(files)

def find_all():
    file_names = "./nmpi*.out"
    files = glob.glob(file_names)
    dict_execs = {}
    complete_execs = {}
    file_per_job = {}
    prot_per_job = {}
    use_by = time.time() - 30 * 60
    for f in files:
        t_name = f.split("-")[1]
        job_name, find_prot = delete_prot_from_name(t_name)
        if (job_name in dict_execs.keys()):
            dict_execs[job_name] += 1
            if (os.path.getatime(f) > use_by):
                complete_execs[job_name] += 1
            file_per_job[job_name].append(f)
            if not(find_prot in prot_per_job[job_name]):
                prot_per_job[job_name].append(find_prot)
        else:
            dict_execs[job_name] = 1
            file_per_job[job_name] = [f]
            prot_per_job[job_name] = [find_prot]
            if (os.path.getatime(f) > use_by):
                complete_execs[job_name] = 1
            else:
                complete_execs[job_name] = 0
           
    dict_execs = {k:v for k, v in dict_execs.items() if v > 2}
    pdbs_count = {}
    for k, v in dict_execs.items():
        #pdbs_number = found_pdbs(k)
        pdbs_count[k]= 0
        
    for k, v in dict_execs.items():
        print("%s - data: %d/%d - pdbs: %d : " % (k, v, complete_execs[k], pdbs_count[k]),end="")
        for p in prot_per_job[k]:
            print("%s " % p, end="")
        print()
        #make_folder(k+"_data")
    return dict_execs, file_per_job
    # for k, v in dict_execs.items():
    #     print("%s - %d" % (k, v) )

    

def distribute_files(job_name, dict_execs, file_per_job ):
    print(file_per_job[job_name])
    make_folder(job_name+"_data")
    for f in file_per_job[job_name]:
        source = f
        destiny = job_name+"_data/" + f[2:]
        os.rename(source, destiny)
    

def distribute_pdbs(job_name):
    values = job_name.split("_")
    res_value = [i for i in values if "res" in i][0]
    deg_value = [i for i in values if "dist" in i or "deg" in i][0]
    deg_value = deg_value.replace("dist","deg")
    print(res_value)
    print(deg_value)

    prots = ["1wit","1c9oA","256bA", "dimaio"]
    for p in prots:
        #print("./*"+job_name+"*"+p+"*.out")
        file_names = "cart_pose_*"+p+"*"+deg_value+"*"+res_value+"*.pdb"
        distribute_pdb_filenames(file_names, job_name)
        # for f in glob.glob("./*"+job_name+"*"+p+"*.out"):
        #     lines = [line.rstrip('\n') for line in open(f)]
        #     for i, l in enumerate(lines[::-1]):
        #         if ("ind " in l):
        #             pdb_file = l.split(":")[1].strip()
        #             destiny_folder = job_name+"_"+p+"_pdbs" 
        #             make_folder(destiny_folder)
        #             source = pdb_file
        #             destiny = destiny_folder + "/" + pdb_file
        #             if os.path.exists(pdb_file):
        #                 print("mv " +source + " " + destiny)                       
        #                 os.rename(source, destiny)
                        
   
def distribute_pdb_filenames(file_names, job_name):
    print(file_names)
    files = glob.glob(file_names)
    prots = ["1wit","1c9oA","256bA", "dimaio"]
    found_prots = []
    for fname in files:
        find_prot = [prot for prot in prots if prot in fname][0]
        alread_found = True if find_prot in found_prots else False
        destiny_folder = job_name+"_"+find_prot+"_pdbs" 
        if not(alread_found):
            make_folder(destiny_folder)
        source = fname
        destiny = destiny_folder + "/" + fname
        #print("mv " +source + " " + destiny)
        os.rename(source, destiny)
    
@click.command()
@click.option('--v', default=True, help="verbose")
@click.option('--find', default=False, help="find files")
@click.option('--pdb_distribute', default="", help="find files")
def main(v, find, pdb_distribute):
    if (find):
        dict_execs, file_per_job = find_all()
        input_job_name = input(" > job name to distribute ")
        option_to_exec = input(" > exec [data, pdbs] ")
        if (option_to_exec == "data"):
            distribute_pdbs(input_job_name)
            distribute_files(input_job_name, dict_execs, file_per_job)
        else:
            if (option_to_exec == "pdbs"):
                distribute_pdbs(input_job_name)
            else:
                print("wrong exec option")
    else:
        if (len(pdb_distribute) > 0):
            distribute_pdbs(pdb_distribute)           
        else:
            #prots = ["1wit","1c9oA","256bA"]
            #prots = ["1elwA", "1c9oA", "256bA", "1kpeA", "1hz6A", "1wit"]


            # 4 sets enviados
            prots_1 = ["1c8cA", "1elwA", "1c9oA", "256bA", "1ten", "1opd"]
            prots_2 = ["1who", "1fna", "1bgf", "1kpeA","1hz6A", "1tit"]
            prots_3 = ["2chf","1npsA","1acf","1bkrA","2ci2I","2vik","1wit"]
            prots_4 = ["1vcc","1tul","1iibA","1ctf", "1gvp"]
            #prots = prots_1 + prots_2 + prots_3 + prots_4

            # pendientes de enviar a cola-corta (aunque ya están en shared
            # desde ayer)
            #prots = [ "1tit", "1npsA", "1ctf", "1bkrA", "2ci2I" , "1hz6A"
            #          ,"2chf", "1acf" , "2vik" , "1wit", "1vcc" , "1tul" ,
            #          "1iibA" , "1ctf", "1gvp"]


            # prots = [ "1tit", "1npsA", "1ctf", "1bkrA", "2ci2I" , "1hz6A"
            #           ,"2chf", "1acf" , "2vik" , "1wit", "1vcc" , "1tul" ,
            #           "1iibA" , "1ctf", "1gvp"]

            
            # prots = ["2vik"]
            # prots = ["2chf","1npsA","1acf","1bkrA","2vik"]
            # prots = ["1vcc","1tul","1iibA","1gvp"]
            
            #prots = ["1bgf","2chf","1npsA","1bkrA","2vik","1gvp"]
            #prots = ["2chf","1npsA","1bkrA","2vik"]
            # primer set usado (ya estan listas)
            #prots = ["1cg5B" , "1eyvA" , "1lis" , "1kpeA" , "1rnbA" , "1tig" , "1dhn"]

            # new fragments
            prots = ["n1c8cA", "n1elwA", "n1wit", "n1kpeA", "n1opd", "n1hz6A"]

            #prots = ["dimaio"]
            #prots = ["1wit", "256bA"]
            res = [1]
            dist_deg = [10]
            created_data = create_files_for(prots, res, dist_deg)
            #print(created_data)
            bash_list = create_bash_files(created_data)
            for job_file in bash_list:
                if (v is True ):
                    print("* verbose mode activated, use --v False \n")
                    print("prepared for sbatch %s \n" %job_file)
                else:   
                    print("sbatch %s \n" %job_file)
                    os.system("sbatch %s" %job_file)
            
if __name__== "__main__":
    main()
