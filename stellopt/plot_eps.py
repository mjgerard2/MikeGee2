#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 13:06:27 2019

@author: michael
"""
import os

import numpy as np
import matplotlib.pyplot as plt

dirc = os.getcwd()
dirc_out = os.path.join(dirc, 'output')

cnt = 0

flx_num = 127
inv_flx = 1. / flx_num

dirc_lst = os.listdir(dirc_out)

eps_eff = np.empty((len(dirc_lst), flx_num))
crnt_config = np.empty((len(dirc_lst), 7))

for ind, name in enumerate(dirc_lst):
    dirc_job = os.path.join(dirc_out, name)
    eps_name = os.path.join(dirc_job, 'stellopt.HSX_test')
    
    f = open(eps_name, 'r')
    lines = f.readlines() 
    f.close()
    
    eps = np.empty(127)
    for i in range(4, 131):
        line = lines[i].strip()
        line = line.split()
        eps[i - 4] = float(line[2])
        
    eps_eff[ind] = eps
    
    crnt_name = os.path.join(dirc_job, 'input.HSX_test')
    
    f = open(crnt_name, 'r')
    lines = f.readlines() 
    f.close()
    
    line = lines[22]
    line = line.strip()
    line = line.split()
    
    crnt_config[ind] = np.array(line[2::])

flx = np.linspace(0, 1, flx_num)
plt.plot(flx, eps_eff[0] * 10**4)
    
plt.xlabel('r/a', fontsize=16)
plt.ylabel(r'$\epsilon_{eff} \ (10^{-4})$', fontsize=16)

plt.savefig('eps_eff.png')
plt.close()
