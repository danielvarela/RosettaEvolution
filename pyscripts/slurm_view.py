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
import subprocess
import datetime

view_all_slurm = "--format=jobid,jobname,elapsed,ncpus,State"
#proc = subprocess.Popen(["sacct", view_all_slurm], stdout=subprocess.PIPE, shell=True)
#(out, err) = proc.communicate()

mylist = []
today_date = datetime.datetime.now()
yesterday_date = datetime.datetime.now() - datetime.timedelta(days=1)
mylist.append(today_date)
today_str = today_date.strftime("%Y-%m-%d-%H:%M")
yesterday_str = yesterday_date.strftime("%Y-%m-%d-%H:%M")
# print(today_str)
# print(yesterday_str)

proc = subprocess.Popen(['sacct',view_all_slurm, "-T", "-S"+yesterday_str, "-E"+today_str+"", "--state=F"],stdout=subprocess.PIPE)
unique_ids = []
while True:
    l = proc.stdout.readline().decode('ascii')
    line = l.strip()
    if ("State" in line):
        print(l, end="")
    if ("---" in line):
        print(l, end="")
    if len(line) > 2: 
        #the real code does filtering here
        job_id = line.split()[0].split("_")[0]
        if (job_id[3:6]).isdigit():
            if not(job_id in unique_ids):
                unique_ids.append(job_id)
                print(l, end="")
    else:
        break
    


