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
import shutil
from distutils.dir_util import copy_tree


all_prots = ["1c8cA", "1elwA", "1c9oA", "256bA", "1ten", "1opd", "1who", "1fna", "1rnbA", "1bgf", "1kpeA","1hz6A", "1tit", "1tig", "2chf","1npsA","1acf","1bkrA","2ci2I","2vik","1wit","1vcc","1tul","1lis","1iibA","1ctf","1cg5B", "1dhn", "1eyvA", "1gvp"]

base_path = "/mnt/netapp2/Store_uni/home/ulc/co/dvm/Code/minifolding/"

def find_CR(name):
    tokens = name.split("_")
    for t in tokens:
        if "CR0" in t:
            return t
    return "CR095"

def find_F(name):
    tokens = name.split("_")
    for t in tokens:
        if "F0" in t:
            return t
    return "FX"

def find_prot(name):
    for p in all_prots:
        if p in name:
            return p
    return "PROT"

def make_folder(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

for f in glob.glob(base_path+"nmpi-NoFragsOnly*out"):
    CR = find_CR(f)
    F = find_F(f)
    folder_name = "NoFragsOnlyDE_"+F+"_"+CR+"_RUNS"
    make_folder(base_path+folder_name)
    print(f)
    print(folder_name+"/")
    shutil.copy(f, base_path+folder_name+"/")
 

    
