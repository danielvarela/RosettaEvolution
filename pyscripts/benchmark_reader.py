#!/usr/bin/env python

import os
import glob
import re
import mmap
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.preprocessing import normalize
from sklearn.preprocessing import minmax_scale
import string
from operator import itemgetter
import sys
import json


class BenchmarkReader:
    def __init__(self, file_list):
        self.file_list = file_list
        self.open_files = []
        self.total_times = []
        self.selection_times = []
        self.sample_times = []
        self.evaluation_times = []
        self.gen_times = []
        self.process_files(file_list)

    def process_files(self, filenames_list):
        for filename in filenames_list:
            lines = [line.rstrip('\n') for line in open(filename)]
            for i, e in list(enumerate(lines)):
                tokens = e.split();
                if (len(tokens) > 0):
                    if ("TOTAL BENCHMARK time" in e):
                        self.total_times.append(float(tokens[-1]))
                    if ("selection stage time" in e):
                        self.selection_times.append(float(tokens[-1]))
                    if ("sample stage time" in e):
                        self.sample_times.append(float(tokens[-1]))
                    if ("evaluation stage time" in e):
                        self.evaluation_times.append(float(tokens[-1]))
                    if ("gen stage time" in e):
                        self.gen_times.append(float(tokens[-1]))

    def print_raw(self):
        s = pd.Series(self.total_times)
        print("* Data total_times")
        print(s.describe())
        s = pd.Series(self.sample_times)
        print("* Data sample stage")
        print(s.describe())
        s = pd.Series(self.evaluation_times)
        print("* Data evaluation stage")
        print(s.describe())
        s = pd.Series(self.selection_times)
        print("* Data selection stage")
        print(s.describe())
        s = pd.Series(self.gen_times)
        print("* Data gen time")
        print(s.describe())
        
reader = BenchmarkReader(sys.argv[1:])
reader.print_raw()

