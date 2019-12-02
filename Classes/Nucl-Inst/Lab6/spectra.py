#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 19 11:38:38 2019

@author: michael
"""

import matplotlib.pyplot as plt
import matplotlib.lines as mlines

import numpy as np
import xlrd 

import scipy.optimize as spo

def lineFit(x, m, b):
    return m*x + b

def parabFit(x, a, b, c):
    return a*x*x + b*x + c

def gausFit(x, a, x0, sig):
    return a * np.exp( - (x-x0)**2 / (2*sig**2) )

def extNum(line, idx):
    line = line.strip()
    line = line.split()
    return float(line[idx])


class spectra():
    
    def __init__(self, loc, spc, xrd=False):
        self.spc = spc
        
        f = open(loc, 'r')
        lines = f.readlines()
        
        self.elt = extNum(lines[4], 3)
        npts = 6000#int(extNum(lines[5], 2))
        
        self.fit_a = extNum(lines[6], 5)
        self.fit_b = extNum(lines[7], 5)
        self.fit_c = extNum(lines[8], 5)
        
        chan = np.empty(npts)
        sign = np.empty(npts)
        
        for l, line in enumerate(lines[12:12+npts]):
            line = line.strip()
            line = line.split() 
            chan[l] = float(line[0])
            sign[l] = float(line[1]) 
                        
        self.chan = chan
        self.sign_err = np.sqrt(sign) / self.elt
        self.sign = sign / self.elt
        
        f.close()
            
    def fullSpec(self, save=False, name='None'):
        #plt.scatter(self.chan, self.sign, s=2, c='k')
        #plt.plot(lineFit(self.chan, self.fit_b, self.fit_a), self.sign, c='k')
        plt.plot(self.chan, self.sign, c='k')
        
        plt.xlabel('Energy (keV)', fontsize=14)
        plt.ylabel('Counts (#/s)', fontsize=14)
        plt.title(self.spc, fontsize=16)
        
        if save:
            plt.savefig(name)
        else:
            plt.show()
            
        plt.close()
        
    def peak_id(self, peaks, save=False, name='None'):
        fnt=14
        eng = parabFit(self.chan, self.fit_c, self.fit_b, self.fit_a)
        
        plt.plot(eng, self.sign, c='k')
        
        for p in peaks:
            plt.plot([p[0], p[0]], [0, p[1]], ls='--', label=str(p[0])+' keV')
        
        plt.xlabel('Energy (keV)', fontsize=fnt)
        plt.ylabel('Counts (#/s)', fontsize=fnt)
        plt.title(self.spc, fontsize=fnt+2)
        
        plt.grid()
        plt.legend()
        if save:
            plt.savefig(name)
        else:
            plt.show()
        plt.close()
        
    def calib(self, peaks, save=False, name='None'):
        fig, axs = plt.subplots(2, 1, sharex=True)
        
        axs[0].plot(self.chan, self.sign, c='k')
        
        axs[0].set_ylabel('Counts (#/s)', fontsize=14)
        axs[0].set_title(self.spc+' Calibration', fontsize=16)
        axs[1].set_xlabel('Channel', fontsize=14)
        axs[1].set_ylabel('Energy (keV)', fontsize=14)
        
        for key in peaks:
            data = peaks[key]
            if data[1] == 'pp':
                clr = 'c'
            elif data[1] == 'Na':
                clr = 'r'
            elif data[1] == 'err':
                clr = 'b'

            axs[0].plot([data[0], data[0]], [0, self.sign[data[0]]], ls='--', color=clr)
            
        blue_line = mlines.Line2D([], [], color='c', ls='--', label=r'$^{60}Co$ Peaks')
        red_line = mlines.Line2D([], [], color='r', ls='--', label=r'$^{133}Ba$ Peaks')
        
        axs[0].legend(handles=[blue_line, red_line])
        
        axs[1].plot(self.chan, parabFit(self.chan, self.fit_c, self.fit_b, self.fit_a), color='k', label='Quadratic Fit')
        axs[1].plot(self.chan, lineFit(self.chan, self.fit_b, self.fit_a), ls='--', color='r', label='Linear Fit')
        axs[1].legend()
        
        plt.grid()
        plt.tight_layout()
        if save:
            plt.savefig(name)
        else:
            plt.show()
            
        
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
        

loc = 'montana_soil.ASC'
spec = spectra(loc, spc='Montana Soil')

spec.fullSpec(save=False, name='samp_2.png')

# Sample 1 (Gallium)
'''
peaks = [[411, .34],
         [600, .15],
         [629, .56],
         [833, 1.26],
         [894, .15],
         [1050, .08],
         [2209, .15]]
'''
# Sample 2 (Silver)
'''
peaks = [[446, 1.15],
         [620, .6],
         [657, 12],
         [677, 1.4],
         [686, .9],
         [706, 2],
         [743, .7],
         [763, 2.6],
         [817, .8],
         [884, 6.5],
         [937, 3],
         [1385, 1.6],
         [1477, .4],
         [1506, .8],
         [1544, .4]]
'''
# Sample 3 (Iron)
'''
peaks = [[411, .133],
         [1098, .1],
         [1294, .067],
         [1463, 0.05]]
'''
# Sample 4 (Aluminum)
'''
peaks = [[411, .129],
         [1335, .045],
         [1461, .056]]
'''
# Sample 5 (Arsenic)
'''
peaks = [[558, 73],
         [656, 8],
         [1216, 3],
         [1228, 2]]
'''
# Sample 6 (Gold)
'''
peaks = [[70, 70],
         [411, 980]]
'''
# Sample 7 (Nickel)
'''
peaks = [[74, .14],
         [411, .11],
         [510, .24],
         [810, .655]]
'''
#spec.peak_id(peaks, save=True, name='samp_1.png')

#spec.left_min(1210, 1250, plot=True)
#spec.right_min(1350, 1400, plot=True)

#spec.fullSpec_pk(plot=True, save=False, name='fig_C1.png')
#spec.peak(plot=True, save=False, name='fig_C2.png')