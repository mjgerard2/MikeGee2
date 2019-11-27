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
        data[1][0] = sheet.cell_value(1, 6) / 60.
        data[2][0] = sheet.cell_value(1, 7) / 60.
        
        for i in range(11):
            data[0][i+1] = .001 * sheet.cell_value(5+i, 0)
            data[1][i+1] = sheet.cell_value(5+i, 1) / 60.
            data[2][i+1] = sheet.cell_value(5+i, 2) / 60.
        
        if shld == 'Aluminum':
            bk_pts = 8
        else:
            bk_pts = 9
        bck_data = np.empty((3, bk_pts))
        
        for i in range(bk_pts):
            bck_data[0][i] = .001 * sheet.cell_value(5+i, 7)
            bck_data[1][i] = sheet.cell_value(5+i, 8) / 60.
            bck_data[2][i] = sheet.cell_value(5+i, 9) / 60.
            
        m, b, r_val, p_val, std = linregress(bck_data[0], bck_data[1])
        self.back = line(data[0], m, b)
        self.back_err = np.mean(bck_data[2])
        self.data = data
        
    def plot_cntRate(self, fnt=14):
        x, y = self.data[0], self.data[1]
        
        sub_err = np.sqrt(self.data[2]**2 + self.back_err**2)
        
        #plt.scatter(x, y, c='k')
        plt.errorbar(x, y, yerr=self.data[2], c='k', label='With Background')
        
        #plt.scatter(x, y-y_back, c='r')
        plt.errorbar(x, y-self.back, yerr=sub_err, c='r', label='Background Subtracted')
        
        plt.xlabel('Thickness $(g/cm^{2})$', fontsize=fnt)
        plt.ylabel('Counts $(\#/s^{-1})$', fontsize=fnt)
        plt.title(self.shld, fontsize=fnt+2)
        
        plt.legend()
        plt.grid()
        #plt.show()
        plt.savefig(self.shld+'_bckg.png')
        plt.close()
        
    def plot_linReg(self):
        x, y, y_err = self.data[0], np.log(self.data[1] - self.back), (np.sqrt(self.data[2]**2 + self.back_err**2)/self.data[1])
        self.m, self.b, r_value, p_value, std_err = linregress(x, y)
        
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        
        ax.errorbar(x, y, yerr=y_err, ls='None', c='k')
        ax.scatter(x, y, s=25, c='k')
        if self.shld == 'Aluminum':
            ax.plot(x, line(x, self.m, self.b), ls='--', c='k', label=r'$\mu_{Al}$'+r' = (%.3f $\pm$ %.3f)' % (-self.m, std_err) +r' $cm^{2}/g$')
        else:
            ax.plot(x, line(x, self.m, self.b), ls='--', c='k', label=r'$\mu_{Pb}$'+r' = (%.3f $\pm$ %.3f)' % (-self.m, std_err) +r' $cm^{2}/g$')

        
        fnt = 14
        
        ax.set_xlabel('Thickness $(g/cm^{2})$', fontsize=fnt)
        ax.set_ylabel(r'$\ln(\#/s)$', fontsize=fnt)
        
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