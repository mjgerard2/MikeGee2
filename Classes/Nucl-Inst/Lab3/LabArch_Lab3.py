#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 21:27:36 2019

@author: michael
"""
import matplotlib.pyplot as plt
import scipy.interpolate as spi
import scipy.optimize as spo
import numpy as np

import xlrd 

def func(x, a, b, c, d, e, f, g):
    return (a*x*x*x + b*x*x + c*x + d) / (e*x*x + f*x +g)

def func_hg(x, a, b, c):
    return a * np.exp(b * x) + c

def line(x, a, b):
    return a * x + b
  
# Give the location of the file 
loc = ("SCA_count_1.xlsx") 
loc1 = ('parafin_thinkness.xlsx')
  
# To open Workbook 
wb = xlrd.open_workbook(loc) 
sheet = wb.sheet_by_index(0) 

wb1 = xlrd.open_workbook(loc1)
sheet1 = wb1.sheet_by_index(0)

fnt = 16

h_scl = 1000./1024  
npts = 15
hv = np.empty(npts)
cnts = np.empty(npts)
hght = np.empty(npts)
for i in range(npts):
    hv[i] = sheet.cell_value(i+1, 0)
    cnts[i] = sheet.cell_value(i+1, 1) 
    hght[i] = sheet.cell_value(i+1, 2) * h_scl
    
thck = np.empty(6)
p_cnt = np.empty(6)
for i in range(6):
    thck[i] = sheet1.cell_value(i+1, 0)
    p_cnt[i] = sheet1.cell_value(i+1, 1)
    
p_cnt_err = .05 * np.sqrt(p_cnt)
p_cnt = .05 * p_cnt

popt, pcov = spo.curve_fit(line, thck, p_cnt)

plt.scatter(thck, p_cnt, s=10)
plt.errorbar(thck, p_cnt_err, linestyle="None", elinewidth=1)

plt.plot(thck, line(thck, *popt), ls='--', c='b')

plt.ylim(250, 375)
plt.ylabel(r'Neutron Count $(\# \cdot s^{-1})$', fontsize=fnt)
plt.xlabel('Parafin Thickness (~.5 in)', fontsize=fnt)
plt.title('Neutron Thermalization', fontsize=fnt+2)
#plt.show()
plt.savefig('parafin_thck.png')
plt.close()
    
'''
cnts_err = np.sqrt(cnts)
plt.scatter(hv, .01*cnts, s=5)
plt.errorbar(hv, .01*cnts, yerr=.01*cnts_err, linestyle="None", elinewidth=1)

cnt_plt = .01*np.delete(cnts, [1, 3])
hv_plt = np.delete(hv, [1, 3])

#cnt_intp = spi.interp1d(hv_plt, cnt_plt)
popt, pcov = spo.curve_fit(func, hv_plt, cnt_plt)
hv_int = np.linspace(np.min(hv_plt), np.max(hv_plt), 100)

plt.plot(hv_plt[10:13], func(hv_plt[10:13], *popt), ls='--', c='orange', label='Rising')
plt.plot(hv_plt[7:11], func(hv_plt[7:11], *popt), ls='--', c='g', label='Transition')
plt.plot(hv_plt[0:8], func(hv_plt[0:8], *popt), ls='--', c='r', label='Plateau')


plt.ylabel(r'Counts ($\# \cdot s^{-1}$)', fontsize=fnt)
plt.xlabel('Voltage (kV)', fontsize=fnt)

plt.title('SCA Counter', fontsize=fnt+2)

plt.legend(loc='lower left')
plt.savefig('Counting_Curve.png')
#plt.show()
plt.close()

popt, pcov = spo.curve_fit(func_hg, hv, hght)

plt.scatter(hv, hght)
plt.plot(hv_int, func_hg(hv_int, *popt), ls='--', c='b')

plt.ylabel('Peak Height', fontsize=fnt)
plt.xlabel('Voltage (kV)', fontsize=fnt)

plt.title('MCA Peaks', fontsize=fnt+2)

#plt.show()
plt.savefig('Peak_curve.png')
plt.close()
'''