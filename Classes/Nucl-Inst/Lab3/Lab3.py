#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 13 23:26:53 2019

@author: michael
"""

import matplotlib.pyplot as plt
import scipy.interpolate as spi
import scipy.optimize as spo
import numpy as np

import xlrd 

def line(x, m, b):
    return m*x + b

loc = ("Question_3.xlsx") 

wb = xlrd.open_workbook(loc) 
sheet = wb.sheet_by_index(0) 

npts=10
fnt=16

volt = np.empty(npts)
keys = ['1', '2', '3', '4']
clrs = {'1' : 'k',
        '2' : 'b',
        '3' : 'r',
        '4' : 'g'}
cnts = {}
cnts_err = {}

for key in keys:
    cnts[key] = np.empty(npts)

for i in range(npts):
    volt[i] = sheet.cell_value(i+1, 0)
    for key in keys:
        cnts[key][i] = sheet.cell_value(i+1, int(key))

for key in keys:
    cnts_err[key] = .05 * np.sqrt(cnts[key])
    cnts[key] = .05 * cnts[key]

for key in keys:
    plt.scatter(volt, cnts[key], c=clrs[key])
    plt.errorbar(volt, cnts[key], yerr=cnts_err[key], linestyle='None', c=clrs[key])
    
    popt, pcov = spo.curve_fit(line, volt, cnts[key])
    plt.plot(volt, line(volt, *popt), ls='--', c=clrs[key], label='0.{} V'.format(key))

plt.xlabel('Voltage (kV)', fontsize=fnt)
plt.ylabel(r'Counts $(\# \cdot s^{-1})$', fontsize=fnt)
plt.title(r'$\beta$ Decay', fontsize=fnt+2)

plt.legend(title='Lower Discriminator')
plt.show()
#plt.savefig('Beta_decay.png')
plt.close()