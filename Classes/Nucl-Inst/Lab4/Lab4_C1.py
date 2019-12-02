#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 19 11:38:38 2019

@author: michael
"""

import matplotlib.pyplot as plt
import numpy as np
import xlrd 

import scipy.optimize as spo

def lineFit(x, m, b):
    return m*x + b

def parabFit(x, a, b, c):
    return a*x*x + b*x + c

def gausFit(x, a, x0, sig):
    return a * np.exp( - (x-x0)**2 / (2*sig**2) )


class spectra():
    
    def __init__(self, loc, npts, spc, xrd=False):
        self.spc = spc
        
        if xrd:
            wb = xlrd.open_workbook(loc)
            sheet = wb.sheet_by_index(0)
            
            chan = np.empty(npts)
            sign = np.empty(npts)
            
            for i in range(npts):
                chan[i] = sheet.cell_value(i, 0)
                sign[i] = sheet.cell_value(i, 1)
                
            self.chan = chan
            self.sign_err = np.sqrt(sign) *.01
            self.sign = sign * .01
            
        else:
            f = open(loc, 'r')
            lines = f.readlines()
            
            live = lines[4]
            live = live.strip()
            live = live.split()
            elt = float(live[-1])
    
            chan = np.empty(npts)
            sign = np.empty(npts)
            
            for l, line in enumerate(lines[12::]):
                line = line.strip()
                line = line.split() 
                chan[l] = float(line[0])
                sign[l] = float(line[1]) 
            
            self.chan = chan
            self.sign_err = np.sqrt(sign) / elt
            self.sign = sign / elt
            
            f.close()
            
    def fullSpec(self, save=False, name='None'):
        plt.scatter(self.chan, self.sign, s=2, c='k')
        
        plt.xlabel('Channel', fontsize=14)
        plt.ylabel('Counts (#/s)', fontsize=14)
        plt.title(self.spc+' Spectra', fontsize=16)
        
        if save:
            plt.savefig(name)
        else:
            plt.show()
        plt.close()
        
    def left_min(self, lft, rgt, wind=100, plot=True):
        chan = self.chan[lft:rgt]
        sign = self.sign[lft:rgt]
        
        popt, pcov = spo.curve_fit(parabFit, chan, sign)
        
        self.chan_min_lft = -.5 * (popt[1] / popt[0])
        self.sign_min_lft = parabFit(self.chan_min_lft, *popt)
        
        if plot:
            plt.scatter(self.chan[lft-wind:lft], self.sign[lft-wind:lft], c='k', s=5)
            plt.scatter(self.chan[rgt:rgt+wind], self.sign[rgt:rgt+wind], c='k', s=5)
        
            plt.scatter(chan, sign, c='r', s=5)
            plt.plot(chan, parabFit(chan, *popt), c='r', ls='--')
            
            plt.scatter([self.chan_min_lft], [self.sign_min_lft], c='k', marker='X', s=40, label='Min: (%.0f, %.2f)' % (self.chan_min_lft, self.sign_min_lft))
            
            plt.xlabel('Channel', fontsize=14)
            plt.ylabel('Counts (#/s)', fontsize=14)
            plt.title('Left Min', fontsize=16)
            
            plt.legend()            
            plt.show()
            plt.close()
            
    def right_min(self, lft, rgt, wind=100, plot=True):
        chan = self.chan[lft:rgt]
        sign = self.sign[lft:rgt]
        
        popt, pcov = spo.curve_fit(parabFit, chan, sign)
        
        self.chan_min_rgt = -.5 * (popt[1] / popt[0])
        self.sign_min_rgt = parabFit(self.chan_min_rgt, *popt)
        
        if plot:
            plt.scatter(self.chan[lft-wind:lft], self.sign[lft-wind:lft], c='k', s=5)
            plt.scatter(self.chan[rgt:rgt+wind], self.sign[rgt:rgt+wind], c='k', s=5)
        
            plt.scatter(chan, sign, c='r', s=5)
            plt.plot(chan, parabFit(chan, *popt), c='r', ls='--')
            
            plt.scatter([self.chan_min_rgt], [self.sign_min_rgt], c='k', marker='X', s=40, label='Min: (%.0f, %.2f)' % (self.chan_min_rgt, self.sign_min_rgt))
            
            plt.xlabel('Channel', fontsize=14)
            plt.ylabel('Counts (#/s)', fontsize=14)
            plt.title('Right Min', fontsize=16)
            
            plt.legend()            
            plt.show()
            plt.close()
            
    def fullSpec_pk(self, plot=True, save=False, name='None'):
        chan_lft = self.chan[self.chan <= self.chan_min_lft]
        chan_rgt = self.chan[self.chan >= self.chan_min_rgt]
        chan_pk = self.chan[(self.chan >= self.chan_min_lft) & (self.chan <= self.chan_min_rgt)]
        
        sign_lft = self.sign[self.chan <= self.chan_min_lft]
        sign_rgt = self.sign[self.chan >= self.chan_min_rgt]
        sign_pk = self.sign[(self.chan >= self.chan_min_lft) & (self.chan <= self.chan_min_rgt)]
        
        plin, pcv = spo.curve_fit(lineFit, [self.chan_min_lft, self.chan_min_rgt], [self.sign_min_lft, self.sign_min_rgt])
        bckg = lineFit(chan_pk, *plin)
        
        sign_bckg = sign_pk - bckg
        a0 = np.max(sign_bckg)
        x0 = .5 * (self.chan_min_lft + self.chan_min_rgt)
        sig0 = .5 * x0
        popt, pcov = spo.curve_fit(gausFit, chan_pk, sign_bckg, p0=[a0, x0, sig0])
        
        a, xo, sig = int(popt[0]), int(popt[1]), int(popt[2])
        FWHM_lft = xo - 1.177410023*sig
        FWHM_rgt = xo + 1.177410023*sig
        FWHM_plt = np.array([FWHM_lft, FWHM_rgt])
        
        if plot:
            plt.scatter(chan_lft, sign_lft, s=2, c='k')            
            plt.scatter(chan_rgt, sign_rgt, s=2, c='k')
            
            plt.scatter(chan_pk, sign_pk, s=2, c='r', label='Photopeak')
            plt.plot(chan_pk, bckg, c='r')
            plt.plot([xo, xo], [0, a]+bckg[chan_pk==xo], c='r', label='Pulse Height')
            
            plt.plot(FWHM_plt, gausFit(FWHM_plt, *popt)+lineFit(FWHM_plt, *plin), c='b', label='FWHM')
            
            plt.xlabel('Channel', fontsize=14)
            plt.ylabel('Counts (#/s)', fontsize=14)
            plt.title(self.spc+' Spectra', fontsize=16)
            
            plt.legend(loc='upper right')   
            if save:
                plt.savefig(name)
            else:
                plt.show()
            plt.close()
            
    def peak(self, plot=True, save=False, name='None'):
        chan = self.chan[(self.chan >= self.chan_min_lft) & (self.chan <= self.chan_min_rgt)]
        sign = self.sign[(self.chan >= self.chan_min_lft) & (self.chan <= self.chan_min_rgt)]
        sign_err = self.sign_err[(self.chan >= self.chan_min_lft) & (self.chan <= self.chan_min_rgt)]
        
        plin, pcv = spo.curve_fit(lineFit, [self.chan_min_lft, self.chan_min_rgt], [self.sign_min_lft, self.sign_min_rgt])
        bckg = lineFit(chan, *plin)
        
        sign_bckg = sign - bckg
        a0 = np.max(sign_bckg)
        x0 = .5 * (self.chan_min_lft + self.chan_min_rgt)
        sig0 = .25 * (self.chan_min_rgt - self.chan_min_lft)
        popt, pcov = spo.curve_fit(gausFit, chan, sign_bckg, p0=[a0, x0, sig0])
        gaus = gausFit(chan, *popt)
                
        a, xo, sig = int(popt[0]), int(popt[1]), int(popt[2])
        FWHM_lft = xo - 1.177410023*sig
        FWHM_rgt = xo + 1.177410023*sig
        FWHM_plt = np.array([FWHM_lft, FWHM_rgt])
        FWHM = 2.354820045*sig

        if plot:
            print('Resolution: %.2f' % (FWHM/xo))
            
            plt.scatter(chan, sign, c='r', s=5)
            plt.errorbar(chan, sign, yerr=sign_err, linestyle='None', c='r', elinewidth=.5)
            plt.plot(chan, gaus+bckg, c='r', ls='--', label='Gaussian Fit')
            
            plt.plot(chan, bckg, c='r')
            plt.plot([xo, xo], [0, a]+bckg[chan==xo], c='r' , label='Pulse Height')
            
            plt.plot(FWHM_plt, gausFit(FWHM_plt, *popt)+lineFit(FWHM_plt, *plin), c='b', label='FWHM')
            
            plt.xlabel('Channel', fontsize=14)
            plt.ylabel('Counts (#/s)', fontsize=14)
            plt.title(self.spc+' Peak Resolution: %.2f' % (FWHM/xo), fontsize=16)
            
            plt.legend(loc='upper left')   
            if save:
                plt.savefig(name)
            else:
                plt.show()
            plt.close()
        

loc = 'C_3_Cs.ASC'
npts=2048
spec = spectra(loc, npts, spc=r'$^{137}_{55}Cs$')

spec.fullSpec(save=True, name='Cs_Full.png')

#spec.left_min(1210, 1250, plot=True)
#spec.right_min(1350, 1400, plot=True)

#spec.fullSpec_pk(plot=True, save=False, name='fig_C1.png')
#spec.peak(plot=True, save=False, name='fig_C2.png')