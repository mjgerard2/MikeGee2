# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import xlrd

import numpy as np
import matplotlib.pyplot as plt

file_name = 'Count_Data.xlsx'
file = xlrd.open_workbook(file_name)
file_sheet = file.sheet_by_index(0)

Ypts = 5
Xpts = 7

data = np.empty((Ypts, Xpts))

for y in range(Ypts):
    for x in range(Xpts):
        data[y, x] = file_sheet.cell_value(x+1, y)

fnt=16

fig = plt.figure()
ax1 = plt.subplot(211)
ax2 = plt.subplot(212)
ax3 = ax2.twinx()

ax1.errorbar(data[0], .1*data[1], yerr=.1*np.sqrt(data[1]), xerr=0.01, label='MCA')
ax1.errorbar(data[0], .1*data[2], yerr=.1*np.sqrt(data[2]), xerr=0.01, label='Counter')

lns1 = ax2.errorbar(data[0], data[3], yerr=5, xerr=0.01, color='k')
lns2 = ax3.errorbar(data[0], data[4], yerr=0.01, xerr=0.01, color='r')

ax1.set_ylabel(r'Counts (#/s)', fontsize=fnt)

ax2.set_xlabel('Applied Voltage (kV)', fontsize=fnt)
ax2.set_ylabel('Channel', fontsize=fnt, color='k')
        
ax3.set_ylabel('FWHM', fontsize=fnt, color='r')
ax3.tick_params(axis='y', colors='r')
        
ax1.legend()

plt.tight_layout()
#plt.show()
plt.savefig('count_plot.png')
plt.close()