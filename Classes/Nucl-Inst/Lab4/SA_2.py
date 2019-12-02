#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 26 22:16:40 2019

@author: michael
"""

import numpy as np
import matplotlib.pyplot as plt

import scipy.optimize as spo

def lineFit(x, m, b):
    return m*x + b

res = np.array([.1, .07, .05, .11])
chan = np.array([679, 1151, 1304, 418])
enrg = np.array([.662, 1.173, 1.331, .412])
'''
popt, pcov = spo.curve_fit(lineFit, enrg, chan)

plt.scatter(enrg[0], chan[0], c='r', label=r'$^{137}_{55}Cs$')
plt.scatter(enrg[3], chan[3], c='g', label=r'$^{198}Au$')
plt.scatter(enrg[1:3], chan[1:3], c='c', label=r'$^{60}_{27}Co$')

plt.plot(enrg, lineFit(enrg, *popt), ls='--', c='b', label='Slope = %.0f' % popt[0])

plt.xlabel('Gamma Energy (MeV)', fontsize=14)
plt.ylabel('Channel', fontsize=14)
#plt.title(': %.2f' % (FWHM/xo), fontsize=16)

plt.legend()
#plt.show()
plt.savefig('SA_2.png')
plt.close()
'''
popt, pcov = spo.curve_fit(lineFit, enrg, res)

plt.scatter(enrg[0], res[0], c='r', label=r'$^{137}_{55}Cs$')
plt.scatter(enrg[3], res[3], c='g', label=r'$^{198}Au$')
plt.scatter(enrg[1:3], res[1:3], c='c', label=r'$^{60}_{27}Co$')

plt.plot(enrg, lineFit(enrg, *popt), ls='--', c='b', label='Slope = %.3f' % popt[0])

plt.xlabel('Gamma Energy (MeV)', fontsize=14)
plt.ylabel('Resolution', fontsize=14)
#plt.title(': %.2f' % (FWHM/xo), fontsize=16)

plt.legend()
#plt.show()
plt.savefig('SA_3.png')
plt.close()