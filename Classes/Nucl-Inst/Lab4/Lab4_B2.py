#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 19 10:53:16 2019

@author: michael
"""

import matplotlib.pyplot as plt
import numpy as np

import xlrd 

loc = ("B2_data.xlsx")

wb = xlrd.open_workbook(loc)
sheet = wb.sheet_by_index(0)

npts = 17

hv = np.empty(npts)
cnts = np.empty(npts)
err = np.empty(npts)
for i in range(npts):
    hv[i] = sheet.cell_value(i+6, 0)
    cnts[i] = .05*sheet.cell_value(i+6, 1)
    err[i] = .05*sheet.cell_value(i+6, 2)
    
fnt=16

plt.scatter(hv, cnts, s=20)
plt.errorbar(hv, cnts, yerr=err)

plt.xlabel('Applied Voltage (kV)', fontsize=fnt)
plt.ylabel(r'Counts ($\# / s$)', fontsize=fnt)
plt.title('Lower Discriminator: {}'.format(0.3), fontsize=fnt+2)

#plt.show()
plt.savefig('fig_B2.png')
plt.close()
