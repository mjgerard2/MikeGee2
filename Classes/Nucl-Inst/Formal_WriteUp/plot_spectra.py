#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 13:41:59 2019

@author: michael
"""

import numpy as np
import matplotlib.pyplot as plt

def plot_spectra(loc, pres):
    f = open(loc, 'r')
    lines = f.readlines()
    
    npts = len(lines) - 12
    data = np.empty((2, npts))
    for i in range(npts):
        line = lines[12+i]
        line = line.strip()
        line = line.split()
        
        data[0][i] = float(line[0])
        data[1][i] = float(line[1])
    
    f.close()
    
    plt.plot(data[0], data[1], label=r'{} $\pm$ {}'.format(pres[0], pres[1]))

loc = ['Wk1_6_1.ASC', 'Wk1_6_5.ASC', 'Wk1_6_8.ASC', 'Wk1_6_10.ASC', 'Wk1_6_13.ASC', 'Wk1_6_15.ASC', 'Wk1_6_20.ASC']
pres = [[1.64, .01],
        [51.2, .1],
        [205, 1],
        [306, 1],
        [461, 1],
        [536, 1],
        [557, 1]]

for ind, l in enumerate(loc):
    plot_spectra(l, pres[ind])

fnt = 14

plt.xlabel('Channel', fontsize=fnt)
plt.ylabel('Counts', fontsize=fnt)
plt.title('Spectra', fontsize=fnt+2)

plt.legend(title='Pressure (Torr)')
#plt.show()
plt.savefig('spectra_samples.png')
plt.close()