#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 14:35:54 2019

@author: michael
"""

import xlrd
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import linregress

def line(x, m, b):
    return m * x + b

class data():
    
    def __init__(self, sheet, shld):
        self.sheet = sheet
        self.shld = shld
        
        data = np.empty((3, 12))
        
        data[0][0] = 0.
        data[1][0] = sheet.cell_value(1, 6)
        data[2][0] = sheet.cell_value(1, 7)
        
        for i in range(11):
            data[0][i+1] = sheet.cell_value(5+i, 0)
            data[1][i+1] = sheet.cell_value(5+i, 1)
            data[2][i+1] = sheet.cell_value(5+i, 2)
            
        self.data = data
        
    def plot_linReg(self):
        x, y, y_err = self.data[0]*.001, np.log10(self.data[1]), (self.data[2]/self.data[1])
        self.m, self.b, r_value, p_value, std_err = linregress(x, y)
        
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        
        ax.errorbar(x, y, yerr=y_err, ls='None', c='k')
        ax.scatter(x, y, s=25, c='k')
        ax.plot(x, line(x, self.m, self.b), ls='--', c='k', label=r'Slope = %.3f $\pm$ %.3f' % (self.m, std_err))
        
        fnt = 14
        
        ax.set_xlabel('Pressure (kTorr)', fontsize=fnt)
        ax.set_ylabel(r'$log_{10}(\#)$', fontsize=fnt)
        
        ax.set_title(self.shld+' Shield', fontsize=fnt+2)
        
        plt.legend()
        #plt.show()
        plt.savefig(self.shld+'_shield.png')
        plt.close()
        
    

loc = ('Range_Energy_Week_2.xlsx')

wb = xlrd.open_workbook(loc)

Al_sheet = data(wb.sheet_by_index(0), 'Aluminum')
Al_sheet.plot_linReg()

Pb_sheet = data(wb.sheet_by_index(1), 'Lead')
Pb_sheet.plot_linReg()