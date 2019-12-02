#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 19 10:19:56 2019

@author: michael
"""
import matplotlib.pyplot as plt
import numpy as np

import xlrd 

loc = ("B_3_data.xlsx") 

wb = xlrd.open_workbook(loc) 
sheet = wb.sheet_by_index(0) 

npts = 25

gain = np.empty(npts)
cnts = np.empty(npts)
err = np.empty(npts)
for i in range(npts):
    gain[i] = sheet.cell_value(i+1, 2)
    cnts[i] = .05*sheet.cell_value(i+1, 3)
    err[i] = .05*sheet.cell_value(i+1, 4)
        
fnt=16

plt.scatter(gain, cnts, s=20)
plt.errorbar(gain, cnts, yerr=err)

plt.xlabel('Gain', fontsize=fnt)
plt.ylabel(r'Counts ($\# / s$)', fontsize=fnt)
plt.title('Lower Discriminator: {}'.format(0.3), fontsize=fnt+2)

plt.show()
#plt.savefig('fig_B3.png')
plt.close()