#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 15:22:35 2019

@author: michael
"""

import matplotlib.pyplot as plt
import numpy as np

import xlrd 

loc = ("Wk2_2.xlsx")

wb = xlrd.open_workbook(loc)
sheet = wb.sheet_by_index(0)

npts = 8

hv = np.empty(npts)

cnts_408_Abs = np.empty(npts)
cnts_408_NoAbs = np.empty(npts)

cnts_400_Abs = np.empty(npts)
cnts_400_NoAbs = np.empty(npts)

for i in range(npts):
    hv[i] = sheet.cell_value(i+7, 0)
    cnts_408_Abs[i] = sheet.cell_value(i+7, 1)
    cnts_408_NoAbs[i] = sheet.cell_value(i+7, 2)
    
    cnts_400_Abs[i] = sheet.cell_value(i+17, 1)
    cnts_400_NoAbs[i] = sheet.cell_value(i+17, 2)
    
cnts_err_408_Abs = np.sqrt(cnts_408_Abs) * .05
cnts_408_Abs = cnts_408_Abs * .05

cnts_err_408_NoAbs = np.sqrt(cnts_408_NoAbs) * .05
cnts_408_NoAbs = cnts_408_NoAbs * 0.05

cnts_408 = cnts_408_NoAbs - cnts_408_Abs
cnts_err_408 = np.sqrt(cnts_err_408_Abs**2 + cnts_err_408_NoAbs**2)

cnts_err_400_Abs = np.sqrt(cnts_400_Abs) * .05
cnts_400_Abs = cnts_400_Abs * .05

cnts_err_400_NoAbs = np.sqrt(cnts_400_NoAbs) * .05
cnts_400_NoAbs = cnts_400_NoAbs * .05

cnts_400 = cnts_400_NoAbs - cnts_400_Abs
cnts_err_400 = np.sqrt(cnts_err_400_Abs**2 + cnts_err_400_NoAbs**2)
'''
plt.scatter(hv, cnts_408_NoAbs, label='Without Absorber')
plt.errorbar(hv, cnts_408_NoAbs, yerr=cnts_err_408_NoAbs)

plt.scatter(hv, cnts_408_Abs, label='With Absorber')
plt.errorbar(hv, cnts_408_Abs, yerr=cnts_err_408_Abs)

plt.scatter(hv, cnts_408, label='Difference')
plt.errorbar(hv, cnts_408, yerr=cnts_err_408, ls='--')
'''
plt.scatter(hv, cnts_400_NoAbs, label='Without Absorber')
plt.errorbar(hv, cnts_400_NoAbs, yerr=cnts_err_400_NoAbs)

plt.scatter(hv, cnts_400_Abs, label='With Absorber')
plt.errorbar(hv, cnts_400_Abs, yerr=cnts_err_400_Abs)

plt.scatter(hv, cnts_400, label='Difference')
plt.errorbar(hv, cnts_400, yerr=cnts_err_400, ls='--')

plt.xlabel('Channel', fontsize=14)
plt.ylabel('Counts (#/s)', fontsize=14)
plt.title(r'$^{198}Au$ with 400 Detector', fontsize=16)


plt.legend(loc='upper left')
#plt.show()
plt.savefig('Wk2_2_1.png')
plt.close()