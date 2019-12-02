#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 13:25:51 2019

@author: michael
"""

import matplotlib.pyplot as plt
import numpy as np

import scipy.optimize as spo

class detector():
    
    def __init__(self, spc, loc_400_Abs, loc_400_NoAbs, loc_408_Abs, loc_408_NoAbs, npts):
        self.spc = spc
        
        self.chan_400_Abs, self.sign_400_Abs, self.sign_err_400_Abs = self.readFile(loc_400_Abs, npts)
        self.chan_400_NoAbs, self.sign_400_NoAbs, self.sign_err_400_NoAbs = self.readFile(loc_400_NoAbs, npts)
        
        self.chan_408_Abs, self.sign_408_Abs, self.sign_err_408_Abs = self.readFile(loc_408_Abs, npts)
        self.chan_408_NoAbs, self.sign_408_NoAbs, self.sign_err_408_NoAbs = self.readFile(loc_408_NoAbs, npts)
        
    def readFile(self, loc, npts):
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
        
        sign_err = np.sqrt(sign) / elt
        sign = sign / elt
        
        f.close()
        
        return chan, sign, sign_err
        
    def subplot_400_together(self, lft_bnd, rgt_bnd, save=False, name='None'):
        chan_Abs = self.chan_400_Abs[lft_bnd:rgt_bnd]
        sign_Abs = self.sign_400_Abs[lft_bnd:rgt_bnd]
        sign_err_Abs = self.sign_err_400_Abs[lft_bnd:rgt_bnd]
        
        chan_NoAbs = self.chan_400_NoAbs[lft_bnd:rgt_bnd]
        sign_NoAbs = self.sign_400_NoAbs[lft_bnd:rgt_bnd]
        sign_err_NoAbs = self.sign_err_400_Abs[lft_bnd:rgt_bnd]
                
        plt.scatter(chan_Abs, sign_Abs, c='k', s=5, label='With Absorber')
        plt.errorbar(chan_Abs, sign_Abs, yerr=sign_err_Abs, c='k', elinewidth=1)#, linestyle='None')
        
        plt.scatter(chan_NoAbs, sign_NoAbs, c='r', s=5, label='Without Absorber')
        plt.errorbar(chan_NoAbs, sign_NoAbs, yerr=sign_err_NoAbs, c='r', elinewidth=1)#, linestyle='None')
        
        plt.xlabel('Channel', fontsize=14)
        plt.ylabel('Counts (#/s)', fontsize=14)
        plt.title('400 Detector '+self.spc, fontsize=16)
        
        plt.tight_layout()
        plt.legend()
        if save:
            plt.savefig(name)
        else:
            plt.show()
        plt.close()
        
    def subplot_408_together(self, lft_bnd, rgt_bnd, save=False, name='None'):
        chan_Abs = self.chan_408_Abs[lft_bnd:rgt_bnd]
        sign_Abs = self.sign_408_Abs[lft_bnd:rgt_bnd]
        sign_err_Abs = self.sign_err_408_Abs[lft_bnd:rgt_bnd]
        
        chan_NoAbs = self.chan_408_NoAbs[lft_bnd:rgt_bnd]
        sign_NoAbs = self.sign_408_NoAbs[lft_bnd:rgt_bnd]
        sign_err_NoAbs = self.sign_err_408_Abs[lft_bnd:rgt_bnd]
                
        plt.scatter(chan_Abs, sign_Abs, c='k', s=5, label='With Absorber')
        plt.errorbar(chan_Abs, sign_Abs, yerr=sign_err_Abs, c='k', elinewidth=1)#, linestyle='None')
        
        plt.scatter(chan_NoAbs, sign_NoAbs, c='r', s=5, label='Without Absorber')
        plt.errorbar(chan_NoAbs, sign_NoAbs, yerr=sign_err_NoAbs, c='r', elinewidth=1)#, linestyle='None')
        
        plt.xlabel('Channel', fontsize=14)
        plt.ylabel('Counts (#/s)', fontsize=14)
        plt.title('408 Detector '+self.spc, fontsize=16)
        
        plt.tight_layout()        
        plt.legend()
        if save:
            plt.savefig(name)
        else:
            plt.show()
        plt.close()
        
    def subplot_diff(self, lft_bnd, rgt_bnd, save=False, name='None'):
        chan = self.chan_400_Abs[lft_bnd:rgt_bnd]
        sign_400_Abs = self.sign_400_Abs[lft_bnd:rgt_bnd]
        sign_err_400_Abs = self.sign_err_400_Abs[lft_bnd:rgt_bnd]
        
        sign_400_NoAbs = self.sign_400_NoAbs[lft_bnd:rgt_bnd]
        sign_err_400_NoAbs = self.sign_err_400_Abs[lft_bnd:rgt_bnd]
        
        sign_400 = sign_400_NoAbs - sign_400_Abs
        sign_err_400 = np.sqrt(sign_err_400_NoAbs**2 + sign_err_400_Abs**2)
        
        chan = self.chan_408_Abs[lft_bnd:rgt_bnd]
        sign_408_Abs = self.sign_408_Abs[lft_bnd:rgt_bnd]
        sign_err_408_Abs = self.sign_err_408_Abs[lft_bnd:rgt_bnd]
        
        sign_408_NoAbs = self.sign_408_NoAbs[lft_bnd:rgt_bnd]
        sign_err_408_NoAbs = self.sign_err_408_Abs[lft_bnd:rgt_bnd]
        
        sign_408 = sign_408_NoAbs - sign_408_Abs
        sign_err_408 = np.sqrt(sign_err_408_NoAbs**2 + sign_err_408_Abs**2)
        
        plt.scatter(chan, sign_400, c='k', s=5, label='Detector 400')
        plt.errorbar(chan, sign_400, yerr=sign_err_400, c='k', elinewidth=.5)#, linestyle='None')
        
        plt.scatter(chan, sign_408, c='r', s=5, label='Detector 408')
        plt.errorbar(chan, sign_408, yerr=sign_err_408, c='r', elinewidth=.5)#, linestyle='None')
        
        plt.xlabel('Channel', fontsize=14)
        plt.ylabel('Counts (#/s)', fontsize=14)
        plt.title(self.spc+' Spectra', fontsize=16)
        
        plt.legend()
        if save:
            plt.savefig(name)
        else:
            plt.show()
        plt.close()
        
npts=2048

loc_400_Abs = 'Wk2_400_Co_WithShield.ASC'
loc_400_NoAbs = 'Wk2_400_Co_NoShield.ASC'

loc_408_Abs = 'Wk2_408_Co_WithShield.ASC'
loc_408_NoAbs = 'Wk2_408_Co_NoShield.ASC'

spc = r'$^{60}_{27}Co$'

det = detector(spc, loc_400_Abs, loc_400_NoAbs, loc_408_Abs, loc_408_NoAbs, npts)
#det.subplot_diff(0, 200, save=True, name='Wk2_3_Co_408.png')
det.subplot_400_together(0, 200, save=True, name='Wk2_400_Co')
det.subplot_408_together(0, 200, save=True, name='Wk2_408_Co')