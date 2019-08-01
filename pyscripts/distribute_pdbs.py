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

prot = sys.argv[1]

main_path = "/mnt/netapp2/Store_uni/home/ulc/co/dvm/Code/minifolding"
dir_5A = main_path+"/pdbs_cart_1c9oA_EXP5_"+prot+"/"
dir_2A = main_path+"/pdbs_cart_EXP5_res5_"+prot+"/"
dir_10A = main_path+"/pdbs_stage3_EXP3_"+prot+"/"
dir_EXP4 = main_path+"/pdbs_stage3_EXP4_"+prot+"/"
dirpath = "./*"+prot+"*.pdb"
if not (os.path.isdir(dir_2A)):
    os.mkdir(dir_2A)
if not (os.path.isdir(dir_10A)):
    os.mkdir(dir_10A)
if not (os.path.isdir(dir_5A)):
    os.mkdir(dir_5A)
if not (os.path.isdir(dir_EXP4)):
    os.mkdir(dir_EXP4)


#entries = (os.path.join(dirpath, fn) for fn in os.listdir(dirpath))
entries = glob.glob(dirpath)
entries = ((os.stat(path), path) for path in entries)

# leave only regular files, insert creation date
entries = list((stat[ST_CTIME], path)
           for stat, path in entries if S_ISREG(stat[ST_MODE]))

def count_strategy(entries):
    count = 0
    collection_ids = {}
    len_entries = 0
    for d, p in sorted(entries):
        tokens = p.split("_")
        print(tokens[4])
        current_id = tokens[4] 
        if current_id in collection_ids.keys():
            collection_ids[current_id] += 1
        else:
            collection_ids[current_id] = 1
        len_entries += 1
    # df = pd.DataFrame(collection_ids)
    print(collection_ids)
    print("count strategy, found " + str(len_entries))
    print("if 2 : " + str(len_entries / 2) )
    print("if 3 : " + str(len_entries / 3) )
    print("if 4 : " + str(len_entries / 4) )
    number_of_experiments = input("Enter the number of experiments: ")
    limit = len_entries / int(number_of_experiments)
    print(" limit " + str(limit))
    count_of_experiments = 0
    count_processed = 0
    for cdate, path in sorted(entries):
        print("%s %s " % (time.ctime(cdate), os.path.basename(path)))
        if (count_processed < limit):
            if (count_of_experiments == 0):
                print("goes to EXP1")
                os.rename(os.path.basename(path), dir_2A+os.path.basename(path))
            if (count_of_experiments == 1):
                print("goes to EXP2")
                os.rename(os.path.basename(path), dir_5A+os.path.basename(path))
            if (count_of_experiments == 2):
                print("goes to EXP3")
                #os.rename(os.path.basename(path), dir_10A+os.path.basename(path))
            count_processed += 1
        else:
            count_processed = 1
            count_of_experiments += 1
            if (count_of_experiments == 0):
                print("goes to EXP1")
                os.rename(os.path.basename(path), dir_2A+os.path.basename(path))
            if (count_of_experiments == 1):
                print("goes to EXP2")
                os.rename(os.path.basename(path), dir_5A+os.path.basename(path))
            if (count_of_experiments == 2):
                print("goes to EXP3")
                #os.rename(os.path.basename(path), dir_10A+os.path.basename(path))
             
        count += 1
    
def id_strategy(entries):
    count = 0
    previous_id = 0
    count_current_id = 0
    count_current_exp = 1
    for cdate, path in sorted(entries):
        print("%s %s " % (time.ctime(cdate), os.path.basename(path)))
        tokens = path.split("_")
        print(tokens[2])
        current_id = int(tokens[2])
        print(" %s - %s " % (current_id, previous_id))
        if current_id == previous_id :
            if (count_current_exp == 1):
                print("goes to EXP" + str(count_current_exp))
                os.rename(os.path.basename(path), dir_2A+os.path.basename(path))
            if (count_current_exp == 2):
                print("goes to EXP" + str(count_current_exp))
                os.rename(os.path.basename(path), dir_5A+os.path.basename(path))
            if (count_current_exp == 3):
                print("goes to EXP" + str(count_current_exp))
                os.rename(os.path.basename(path), dir_10A+os.path.basename(path))
            if (count_current_exp == 4):
                print("goes to EXP" + str(count_current_exp))
                # os.rename(os.path.basename(path), dir_EXP4+os.path.basename(path))
            
            count_current_id += 1
        else:
            previous_id = current_id
            if (count_current_id > 480):
                count_current_id = 1
                count_current_exp += 1
                print ("init experiment " + str(count_current_exp) )
                #break
            if (count_current_exp == 1):
                print("goeo EXP" + str(count_current_exp))
                os.rename(os.path.basename(path), dir_2A+os.path.basename(path))
            if (count_current_exp == 2):
                print("goes to EXP" + str(count_current_exp))
                os.rename(os.path.basename(path), dir_5A+os.path.basename(path))
            if (count_current_exp == 3):
                print("goes to EXP" + str(count_current_exp))
                os.rename(os.path.basename(path), dir_10A+os.path.basename(path))
            if (count_current_exp == 4):
                print("goes to EXP" + str(count_current_exp))
                # os.rename(os.path.basename(path), dir_EXP4+os.path.basename(path))
            
        count += 1
 
def name_strategy(entries):
    count = 0
    for cdate, path in sorted(entries):
        print("%s %s " % (time.ctime(cdate), os.path.basename(path)))
        if ("res5" in path):
            print("goes to 5A")
            os.rename(os.path.basename(path), dir_5A+os.path.basename(path))
        else:
            if ("res10" in path):
                print("goes to 10A")
                os.rename(os.path.basename(path), dir_10A+os.path.basename(path))
            else:
                print("goes to 2A")
                os.rename(os.path.basename(path), dir_2A+os.path.basename(path))
        count += 1
    
count_strategy(entries)
