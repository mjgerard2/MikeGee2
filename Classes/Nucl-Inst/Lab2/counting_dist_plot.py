# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 15:44:45 2019

@author: mjgerard2
"""

import xlrd

import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import norm

def gaus(x, avg, std):
    return (1. / (std * np.sqrt(2*np.pi))) * np.exp(-.5*((x - avg) / std)**2)

file_name = 'counting_distribution.xlsx'
file = xlrd.open_workbook(file_name)
file_sheet = file.sheet_by_index(0)

Ypts = 3
Xpts = 50

data = np.empty((Ypts, Xpts))

for y in range(Ypts):
    for x in range(Xpts):
        data[y, x] = file_sheet.cell_value(x+1, y)

fnt=16

#avg, std = norm.fit(data[1])

avg = np.mean(data[1])
std = np.sqrt( (1. / (Xpts-1)) * np.sum((data[1] - avg)**2) )
chi_sqr = (1. / avg) * np.sum((data[1] - avg)**2)

print('Mean = %.3f' % avg)
print('Variance = %.3f' % (std))
print('Chi Squared = %.3f' % chi_sqr)

cnt_min = np.min(data[1])
cnt_max = np.max(data[1])

plt.hist(data[1], bins=25, range=(cnt_min, cnt_max), density=True, alpha=0.6, color='g')

x_min, x_max = plt.xlim()
x_dom = np.linspace(x_min, x_max, 100)
p_dom = norm.pdf(x_dom, avg, std)

plt.plot(x_dom, p_dom, c='k', linewidth=2)

plt.xlabel('Counts', fontsize=fnt)
plt.ylabel('Prob. Den.', fontsize=fnt)

plt.title('Mean = %.0f , std = %.1f' % (avg, std), fontsize=fnt)

#plt.show()
plt.savefig('hist.png')
plt.close()
