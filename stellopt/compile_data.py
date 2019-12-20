#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 12:58:30 2019

@author: michael
"""

import numpy as np
import h5py as hf

import os

cwd = os.getcwd()
name = os.path.join(cwd, 'current_data.h5')

job_dirc = os.path.join(cwd, 'jobs')
dirc_lst = os.listdir(job_dirc)

jobs = len(dirc_lst)

crnt_profile = np.full((jobs, 6), np.nan)
eps_profile = np.full((jobs, 127), np.nan)
for j, job in enumerate(dirc_lst):
    sub = os.path.join(job_dirc, job)
    name_neo = os.path.join(sub, 'neo_out.HSX_test.00000')
    name_cur = os.path.join(sub, 'input.HSX_test')
    
    f = open(name_cur, 'r')
    lines = f.readlines()
    f.close()

    crnt = lines[22]
    crnt = crnt.strip()
    crnt = crnt.split()
        
    try:
        f = open(name_neo, 'r')
        lines = f.readlines()
        f.close()
        
        neo = np.empty(127)
        for ind, line in enumerate(lines):
            line = line.strip()
            line = line.split()
            neo[ind] = line[1]
        
        crnt_profile[j] = crnt[3::]
        eps_profile[j] = neo
        
    except:
        print('{0}, {1}, {2}, {3}, {4}, {5}'.format(*crnt[3::]))
        
crnt_profile = crnt_profile[~np.isnan(crnt_profile)]
eps_profile = eps_profile[~np.isnan(eps_profile)]

crnt_profile = crnt_profile.reshape(int(len(crnt_profile) / 6), 6)
eps_profile = eps_profile.reshape(int(len(eps_profile) / 127), 127)

h5 = hf.File(name, 'r')
        
crnt_old = h5['current profile'][:]
eps_old = h5['epsilon profile'][:]

h5.close()

crnt_new = np.append(crnt_old, crnt_profile).reshape(len(crnt_old)+len(crnt_profile), 6)
eps_new = np.append(eps_old, eps_profile).reshape(len(eps_old)+len(eps_profile), 127)
        
h5 = hf.File(name, 'w')
        
h5.create_dataset('current profile', data=crnt_new)
h5.create_dataset('epsilon profile', data=eps_new)

h5.close()