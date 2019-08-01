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

#all_prots = ["1c8cA", "1elwA", "1c9oA", "256bA", "1ten", "1opd", "1who","1fna","1bgf", "1hz6A", "1tit","2chf","1npsA","1acf","1bkrA","2ci2I","2vik","1wit","1vcc","1tul","1iibA","1ctf", "1gvp"]

#all_prots = ["1bgf","2chf","1npsA","1bkrA","2vik","1gvp"]
#prots = ["1cg5B" , "1eyvA" , "1lis" , "1kpeA" , "1rnbA" , "1tig" , "1dhn"]

all_prots = [ "1hz6A", "1elwA", "1c9oA", "256bA", "1ten", "1opd", "1who", "1fna", "1rnbA", "1kpeA", "1tit", "1tig", "2chf", "1npsA", "1acf", "1bkrA", "2vik","1wit","1iibA","1cg5B", "1dhn", "1c8cA", "2ci2I", "1vcc", "1tul", "1lis", "1gvp", "1eyvA", "1ctf", "1bgf"] 
       
found_prots = {}
for p in all_prots:
    files = glob.glob("./nmpi*"+p+"*.out")
    found_prots[p] = 0
    for f in files:
        lines = [s.strip("\n") for s in open(f).readlines()]
        found_stage4 = False
        for l in lines:
            if "START STAGE stage4" in l:
                found_stage4 = True
                break
        if (found_stage4):
            found_prots[p] = found_prots[p] + 1
        else:
            print("mv %s ../remove_wrong_files/ &&" % f)
    #print("prot: %s, files : %s" % (p, found_prots[p] ))

#print(found_prots)
for p in found_prots:
    if found_prots[p] < 5:
        print('"' + p + '"', end=",")
