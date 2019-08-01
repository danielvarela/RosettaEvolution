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
import subprocess
from sklearn import preprocessing
from sklearn.preprocessing import minmax_scale
from sklearn.preprocessing import MinMaxScaler
import click
from stat import S_ISREG, ST_CTIME, ST_MODE

import datetime
import time

path_dir = "/mnt/netapp2/Store_uni/home/ulc/co/dvm/Code/rosetta/demos/public/"
bash_path = "/mnt/netapp2/Store_uni/home/ulc/co/dvm/Code/minifolding/"
copy_path = "/mnt/netapp2/Store_uni/home/ulc/co/dvm/Code/rosetta/demos/public/example_score"
#copy_path  = "/mnt/netapp2/Store_uni/home/ulc/co/dvm/Code/rosetta/demos/public/UP_cryo_score_res2_EXP1_1wit"

def read_processed_files():
    lines = [line.rstrip('\n') for line in open("database_scored_pdbs.csv")]
    job_list = []
    for l in lines:
        j = l.split(",")[0]
        job_list.append(j)

    return job_list

def process_files(pdbs):
    prots = ["1c9oA", "1wit", "256bA"]
    files = glob.glob(pdbs+"/*.pdb")
    job_name = pdbs
    values = job_name.split("_")
    res_value = [i for i in values if "res" in i][0]
    deg_value = [i for i in values if "dist" in i or "deg" in i][0]
    deg_value = deg_value.replace("dist","deg") 
    found_prot = [p for p in prots if p in job_name][0]
    # print(job_name)
    # print("file count : " + str(len(files)))
    # print(res_value)
    # print(deg_value)
    scored = read_score_script_for(pdbs)
    data = {"job_name" : job_name, "count" : len(files), "scored" : scored, "prot": found_prot , "res" : res_value, "deg" : deg_value}
    return data

def make_folder(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def copy_xml_file(directory_fa):
    copy_file = copy_path+"/A_asymm_refine.xml"
    destiny_file = directory_fa+"/A_asymm_refine.xml"
    copy_file_run = "cp "+copy_file+" "+destiny_file
    process = subprocess.Popen(copy_file_run.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

def copy_run_file(directory_fa, prot):
    copy_file = copy_path+"/Score_1wit.sh"
    destiny_file = directory_fa+"/Score_"+prot+".sh"
    copy_file_run = "cp "+copy_file+" "+destiny_file
    process = subprocess.Popen(copy_file_run.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

    
def get_fa_directory(directory):
    directory_fa = directory[2:].replace("pdbs", "FA_score")   
    return path_dir + directory_fa

def read_score_script_for(directory):
    full_directory_fa = get_fa_directory(directory)
    score_file_list = glob.glob(full_directory_fa + "/*.fsc")
    if (len(score_file_list) > 0):
        score_file = score_file_list[0]
        print("score file found for " + directory )
        return True
    else:
        return False
        
def create_score_script_for(directory, data):
    # directory == d["job_name]
    # necesito crear un directorio para crear los ficheros segun la protein
    full_directory_fa = get_fa_directory(directory)
    print("cd " + full_directory_fa) 
    make_folder(full_directory_fa)
    make_folder(full_directory_fa + "/run_scripts")
    copy_xml_file(full_directory_fa)
    copy_run_file(full_directory_fa, data["prot"])
    files = glob.glob(data["job_name"]+"/*.pdb")
    for i, fpdb in enumerate(files):
        script_f = open(full_directory_fa + "/run_scripts/script_file_for_pdb_"+str(i)+".sh", "w+")
        print_score_script_per_file(script_f, bash_path + fpdb[2:], data)

        #os.chmod(bash_path + fpdb[2:], 0o777)
        
    return full_directory_fa

def print_score_script_per_file(f, pdb_path, data):
    prot = data["prot"]
    res_angs_list = ["2A", "5A", "10A"]
    res_angs = "2A"
    for r in res_angs_list:
        if r in data["job_name"]:
            res_angs = r
            break
    res = res_angs.replace("A", "")
    f.writelines("#!/bin/bash \n \n")
    f.writelines("/mnt/netapp2/Store_uni/home/ulc/co/dvm/Code/rosetta/main/source/bin/rosetta_scripts.linuxgccrelease ")
    f.writelines("  -database /mnt/netapp2/Store_uni/home/ulc/co/dvm/Code/rosetta/main/database/ ")
    f.writelines("  -in::file::s %s " % pdb_path)
    native = bash_path + "all_info/info_"+data["prot"]+"/vf_"+data["prot"]+".pdb"
    f.writelines("  -in:file:native %s " % native)
    f.writelines("  -parser::protocol A_asymm_refine.xml ")
    f.writelines("  -ignore_unrecognized_res ")
    prot_map = bash_path +  "all_info/map_"+data["prot"]+"_2A.mrc"
    f.writelines("  -parser::script_vars denswt=30 rms=1.5 reso=%s map=%s testmap=%s " % (res, prot_map, prot_map))
    f.writelines("  -ignore_unrecognized_res ")
    f.writelines("  -edensity::mapreso %s " % res)
    f.writelines("  -default_max_cycles 200 ")
    f.writelines("  -edensity::cryoem_scatterers ")
    f.writelines("  -edensity:sliding_window_wt 2.0 ")
    f.writelines("  -edensity:sliding_window 3 ")
    f.writelines("  -relax:minimize_bond_angles ")
    f.writelines("  -relax:minimize_bond_lengths ")
    f.writelines("  -relax:jump_move true ")
    f.writelines("  -use_bicubic_interpolation ") 
    #f.writelines("  -beta ")
    score_file = "AUTO_score_"+data["deg"]+"_"+data["res"]+"_RELAX_"+data["prot"]+".fsc"
    f.writelines("  -out:file:scorefile %s " % score_file)
    out_pdb = False
    if (out_pdb):
        f.writelines("  -output_secondary_structure true ")
        f.writelines("  -out:pdb true ")
        couple = (data["deg"], data["res"], data["prot"])
        f.writelines("  -out::suffix _$1_%s_%s_%s_RELAX_ " % couple)
    else:
        f.writelines("  -output_secondary_structure false ")
        f.writelines("  -out:pdb false ")
    f.writelines("  -out:level 100 ")
    f.writelines("  -crystal_refine ")
    f.writelines("  -overwrite\n")

    
@click.command()
@click.option('--pdbs', default="", help="pdbs")
def main(pdbs):
    if (len(pdbs) > 0):
        process_files(pdbs)
    else:
        directories = glob.glob("./ExperimentoImita*pdbs/")
        data_list = []
        for d in directories:
            data_list.append(process_files(d))

            
        processed_job_list = read_processed_files()
        for i, d in enumerate(data_list):
            processed = ""
            if (d["job_name"] in processed_job_list):
                processed = "*"
            scored_char = "!YES!" if d["scored"] else "N"
            couple = (i,processed,scored_char, d["job_name"], d["count"], d["prot"], d["res"], d["deg"])
            print("[%d]%s%s \t %s \t %s \t %s \t %s \t %s" % couple)

        select_exp = input(" > select experiment for evaluation ")
        select_exp = int(select_exp)
        fa_directory = create_score_script_for(data_list[select_exp]["job_name"], data_list[select_exp])
        if (len(fa_directory) > 0):
            with open("database_scored_pdbs.csv", "a") as f:
                now_date = datetime.datetime.now().strftime("%Y-%m-%d-%H:%M")
                f.writelines("%s,%s,%s \n" % (data_list[select_exp]["job_name"], now_date, fa_directory ) )
 
            
if __name__== "__main__":
    main()
