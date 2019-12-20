#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 13:43:00 2019

@author: michael
"""

import os
import numpy as np

dirc = os.getcwd()
mod_dirc = os.path.join(dirc, 'jobs')
file_name = os.path.join(dirc, 'input.HSX_test')

f = open(file_name, 'r')
lines = f.readlines()
f.close()

log_pts = [-16000, -8000, -4000, 0.0, 4000, 8000, 16000]

ns = 128
sts = len(log_pts)
coils = 6
runs = sts**6

dirc_names = []
cnt=0
for ind1, x1 in enumerate(log_pts):
    for ind2, x2 in enumerate(log_pts):
        for ind3, x3 in enumerate(log_pts):
            for ind4, x4 in enumerate(log_pts):
                for ind5, x5 in enumerate(log_pts):
                    for ind6, x6 in enumerate(log_pts):
                        x = np.array([x1, x2, x3, x4, x5, x6])
                        if any(abs(x) == 4000):
                            ind = sts**5 * ind1 + (sts**4 * ind2 + (sts**3 * ind3 + (sts**2 * ind4 + (sts * ind5 + ind6))))
                                                    
                            dirc_name = 'job_{}'.format(cnt)
                            dirc_names.append(dirc_name)
                            new_dirc = os.path.join(mod_dirc, dirc_name)
                            os.makedirs(new_dirc)
                            
                            new_lines = lines
                            new_lines[19] = '  MGRID_FILE = \'mgrid_hsx_complete.nc\'\n'
                            new_lines[22] = '  EXTCUR =  -10722.0 {0} {1} {2} {3} {4} {5}\n'.format(x1, x2, x3, x4, x5, x6)
                            
                            file_name_new = os.path.join(new_dirc, 'input.HSX_test')
                            
                            f = open(file_name_new, 'w')
                            for l in lines:
                                f.write(l)
                            f.close()
                            
                            cnt+=1
                        else:
                            continue
                          
jobs = len(dirc_names)
runs = np.ceil(jobs/9000)
sets = np.linspace(0, jobs, runs, dtype=int)

for s_ind, s in enumerate(sets[1::]):
    name = 'job_list_{}.txt'.format(s_ind+1)
    list_name = os.path.join(mod_dirc, name)
    f = open(list_name, 'w')  
    for i in np.arange(sets[s_ind], sets[s_ind+1]):
        loc = 'jobs/'+dirc_names[i]
        f.write(loc+'\n')
    f.close()

'''
list_name = os.path.join(dirc, 'job_list.txt')
f = open(list_name, 'w')                
for i in range(len(dirc_names)):
    f.write(dirc_names[i]+'\n')
f.close()


list_name = os.path.join(dirc, 'job_list_1.txt')
f = open(list_name, 'w')                
for i in np.arange(8497, 14896):
    name = 'job_{}'.format(i)
    f.write(name+'\n')
f.close()
'''