#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 16 16:10:16 2019

@author: michael
"""

import numpy as np 
import matplotlib.pyplot as plt

class spectra():
    
    def __init__(self, loc, shield, name):
        f = open(loc, 'r')
        lines = f.readlines()
        
        npts = len(lines) - 12
        data = np.empty((2, npts))
        for i in range(npts):
            line = lines[12+i]
            line = line.strip()
            line = line.split()
            
            data[0][i] = float(line[0])
            data[1][i] = float(line[1]) / 60
        
        f.close()
        
        self.data = data
        self.shield = shield
        self.name = name
        
    def plot_spec(self, fnt=14):
        plt.plot(self.data[0][0:355], self.data[1][0:355], c='k')
        plt.plot(self.data[0][441::], self.data[1][441::], c='k')
        
        plt.plot(self.data[0][355:441], self.data[1][355:441], c='r')
        plt.fill_between(self.data[0][355:441], self.data[1][355:441], [0]*(441-355), color='lightcoral', label='ROI')
        
        plt.xlabel('Channel', fontsize=fnt)
        plt.ylabel('Count Rate (#/s)', fontsize=fnt)
        plt.title(self.shield, fontsize=fnt+2)
        
        plt.legend()
        #plt.show()
        plt.savefig(self.name)
        plt.close()

loc = ['Week_2_Unshielded.ASC', 'Week_2_Al_7.ASC', 'Week_2_Pb_6.ASC', 'Week_2_Background.ASC']
shl = ['No Shield', r'Al Shield: $4.2 \ g/cm^{2}$', r'Pb Shield: $4.568 \ g/cm^{2}$', 'Background']
nm = ['no_shield.png', 'Al_shield.png', 'Pb_shield.png', 'bckgrnd.png']
for ind, l in enumerate(loc):
    absorb = spectra(l, shl[ind], nm[ind])
    absorb.plot_spec()