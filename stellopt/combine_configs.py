#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 19:16:17 2019

@author: michael
"""

import h5py as h5
import numpy as np
import os

dirc = os.getcwd()
dirc_jobs = os.path.join(dirc, 'done_jobs')

crnt_name = 'input.HSX_test'
neo_name = 'neo_out.HSX_test.00000'

jobs = os.listdir(dirc_jobs)
num = len(jobs)

new_crnt = []
new_eps = []

fails = []
for idx, j in enumerate(jobs):
    job = os.path.join(dirc_jobs, j)
    
    neo_file = os.path.join(job, neo_name)
    crnt_file = os.path.join(job, crnt_name)
    
    try:
        f = open(neo_file, 'r')
        lines = f.readlines()
        f.close()
        
        neo = np.empty(127)
        for ind, line in enumerate(lines):
            line = line.strip()
            line = line.split()
            neo[ind] = float(line[1])
        
        f = open(crnt_file, 'r')
        lines = f.readlines()
        f.close()
        
        line = lines[22].strip()
        line = line.split()
        crnt = np.empty(6)
        for i, l in enumerate(line[3::]):
            crnt[i] = float(l)
            
        new_crnt.append(crnt)
        new_eps.append(neo)
        
    except:
        fails.append(j)
        
new_crnt = np.array(new_crnt)
new_eps = np.array(new_eps)

hf = h5.File('current_data.h5', 'r')

crnt = hf['current profile'][:]
eps = hf['epsilon profile'][:]

hf.close()


all_crnt = np.append(crnt, new_crnt, axis=0)
all_eps = np.append(eps, new_eps, axis=0)

hf_new = h5.File('current_data_new.h5', 'w')

hf_new.create_dataset('current profile', data=all_crnt)
hf_new.create_dataset('epsilon profile', data=all_eps)

hf_new.close()