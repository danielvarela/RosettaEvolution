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
folder = "DE_WithFrags_Largas"

files = glob.glob(bash_path+folder+"/*.out")
error_files = []
for f in files:
    found = False
    lines = [line.rstrip('\n') for line in open(f)]
    for l in lines:
        if "START STAGE stage4" in l:
            found = True
            break
    if found is False:
        print(f)
        print(len(lines))
        error_files.append(f)

for e in error_files:
    print("rm " + e + " && ")
    
