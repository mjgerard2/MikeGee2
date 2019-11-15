#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 13:48:58 2019

@author: michael
"""

import xlrd

import numpy as np
import matplotlib.pyplot as plt

from matplotlib import gridspec
from scipy.stats import linregress

def line(x, m, b):
    return m * x + b

loc = ('Wk1_MetaData.xlsx')

wb = xlrd.open_workbook(loc)
sheet = wb.sheet_by_index(0)
sheet.cell_value(0, 0)

mV_Torr = 760. / 7440
npts = sheet.nrows -1

### data[0:5] = (pressure, peak: position, width, counts, error in counts)
data = np.empty((5, npts))
for i in range(npts):
    I = i+1
    data[0][i] = sheet.cell_value(I, 0)
    data[1][i] = sheet.cell_value(I, 3)
    data[2][i] = sheet.cell_value(I, 4)
    data[3][i] = sheet.cell_value(I, 5)
    data[4][i] = sheet.cell_value(I, 6)
    
pres_err = np.array([.1, .1, .1, .1, 1, 1, 1, 1, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10])
pres_err = pres_err * mV_Torr
data[0] = mV_Torr * data[0]

arg_sort = np.argsort(data[0])
sort_data = np.empty(data.shape)
for ind, i in enumerate(arg_sort):
    sort_data[0:5, ind] = data[0:5, i]

beg1 = 0
end1 = 15
m1, b1, r_value, p_value, std_err = linregress(sort_data[0][beg1:end1], sort_data[3][beg1:end1])

beg2 = end1
end2 = npts
m2, b2, r_value, p_value, std_err = linregress(sort_data[0][beg2:end2], sort_data[3][beg2:end2])

fnt = 12

fig = plt.figure()
gs = gridspec.GridSpec(4, 2)

ax1 = fig.add_subplot(gs[0:2, 0:2])
ax2 = fig.add_subplot(gs[2, 0:2], sharex=ax1)
ax3 = fig.add_subplot(gs[3, 0:2], sharex=ax2)

ax1.set_title(r'$\alpha$ Emitter Peak Data', fontsize=fnt+2)
ax1.errorbar(sort_data[0], sort_data[3], xerr=pres_err, yerr=sort_data[4], c='b')

ax1.plot(sort_data[0][beg1:end1+1], line(sort_data[0][beg1:end1+1], m1, b1), c='b', ls='--')
ax1.plot(sort_data[0][beg2:end2], line(sort_data[0][beg2:end2], m2, b2), c='b', ls='--')
ax1.set_ylabel('Counts', fontsize=fnt)

ax2.plot(sort_data[0], sort_data[1], c='g')
ax2.scatter(sort_data[0], sort_data[1], s=15, c='g')
ax2.set_ylabel('Channel', fontsize=fnt)

ax3.plot(sort_data[0], sort_data[2], c='r')
ax3.scatter(sort_data[0], sort_data[2], s=15, c='r')
ax3.set_xlabel('Pressure (Torr)', fontsize=fnt)
ax3.set_ylabel('FWHM', fontsize=fnt)

plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)

#plt.show()
plt.savefig('alf_peak.png')
plt.close()