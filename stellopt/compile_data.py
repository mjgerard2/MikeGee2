#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 12:58:30 2019

@author: michael
"""

import numpy as np
import h5py as hf

import os

dirc = os.getcwd()
dirc_lst = os.listdir(dirc)
new_name = os.path.join(dirc, 'current_data.h5')

jobs = 729
crnt_profile = np.empty((jobs, 6))
eps_profile = np.empty((jobs, 127))

cnt = 0
for job in dirc_lst:
    if job[0:3] == 'job':
        sub = os.path.join(dirc, job)
        name_neo = os.path.join(sub, 'neo_out.HSX_test.00000')
        name_cur = os.path.join(sub, 'input.HSX_test')
        
        f = open(name_cur, 'r')
        lines = f.readlines()
        f.close()
        
        crnt = lines[22]
        crnt = crnt.strip()
        crnt = crnt.split()
        
        f = open(name_neo, 'r')
        lines = f.readlines()
        f.close()
        
        neo = np.empty(127)
        for ind, line in enumerate(lines):
            line = line.strip()
            line = line.split()
            neo[ind] = line[1]
        
        crnt_profile[cnt] = crnt[3::]
        eps_profile[cnt] = neo
        
        cnt+=1

h5 = hf.File(new_name, 'w')
        
h5.create_dataset('current profile', data=crnt_profile)
h5.create_dataset('epsilon profile', data=eps_profile)

h5.close()