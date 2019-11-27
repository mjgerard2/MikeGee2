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
from scipy.interpolate import interp1d
from scipy.stats import linregress

def grad_err(f, x):
    err = np.empty(f.shape)
    
    err[0] = .5 * (f[1] - f[0])
    err[-1] = .5  * (f[-1] - f[-2])
    
    for ind, val in enumerate(f[1:-1]):
        x1 = [x[ind-1], x[ind], x[ind+1]]
        y1 = [f[ind-1], f[ind], f[ind+1]]
        m, b, r, s, std = linregress(x1, y1)
        err[ind] = .5 * m * (x[ind+1] - x[ind-1])
    
    return err        
    
def line(x, m, b):
    return m * x + b

def find_nearest(array, val):
    return np.argmin((array - val)**2)

loc = ('Wk1_MetaData.xlsx')

wb = xlrd.open_workbook(loc)
sheet = wb.sheet_by_index(0)
sheet.cell_value(0, 0)

mV_Torr = 760. / 7440
eff_dist = 100 * 0.0000726698
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
    
#pres_err = np.array([.1, .1, .1, .1, 1, 1, 1, 1, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10])
pres_err = np.array([2]*npts)
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

Deff = eff_dist*sort_data[0]
Crat = 0.00833333 * sort_data[3]

Cerr = 0.00833333*sort_data[4]
Derr = eff_dist*np.hypot(pres_err, .6)

D_mark = 3.66
 
indx = find_nearest(Deff, D_mark)
knee_Deff = Deff[indx]
knee_Crat = Crat[indx]

lin_C = Crat[Deff < D_mark]
lin_C_avg = np.mean(lin_C)
lin_C_std = np.std(lin_C)
lin_C_err = np.sum(Cerr[Deff < D_mark]**2) / len(lin_C)

stop_C = Crat[Deff > D_mark]
stop_D = Deff[Deff > D_mark]
ms, bs, r_value, p_value, std_err = linregress(stop_D, stop_C)

stop_C_err = np.sum(Cerr[Deff > D_mark]**2) / len(stop_C)

stop_Ceff = .5*lin_C_avg
stop_Ceff_err = np.hypot(lin_C_std, lin_C_err)

#Cerr_sum = np.sum(Cerr**2) / len(Cerr)
#Cerr_tot = np.hypot(Cerr_sum, stop_Ceff_err)

stop_Deff = (stop_Ceff - bs) / ms
stop_Deff_err =  np.hypot((2*stop_Ceff_err) / ms, stop_C_err) #abs(Cerr_tot / ms)

en_scl = 5.486 / sort_data[1][0]
avg_eng = sort_data[1] * en_scl

mEn, bEn, rEn, sEn, std = linregress(sort_data[0][0:10], sort_data[1][0:10])
#en_err = abs(mEn * en_scl * sort_data[0][0])
en_err = round(abs(en_scl * (bEn - sort_data[1][0])), 2)

en_err = np.hypot(en_scl, en_err)

avg_eng_intp = interp1d(Deff, avg_eng, kind='cubic')
Deff_new = np.linspace(Deff[0], Deff[-1], 1000)
Eng_new = avg_eng_intp(Deff_new)
stop_pwr = -np.gradient(Eng_new, Deff_new[1]-Deff_new[0])

MR_ind = np.argmax(stop_pwr)
stp_MR = stop_pwr[MR_ind]
eng_MR = Eng_new[MR_ind]
Deff_MR = Deff_new[MR_ind]

stp_pwr_intp = interp1d(Eng_new, stop_pwr)

stp_pwr_data = stp_pwr_intp(avg_eng)
stp_err = grad_err(stp_pwr_data, avg_eng)

fnt = 14

### Plot Peak Chan and FWHM ###

fig = plt.figure()
gs = gridspec.GridSpec(4, 2)

ax1 = fig.add_subplot(gs[0:2, 0:2])
ax2 = fig.add_subplot(gs[2:4, 0:2], sharex=ax1)

ax1.errorbar(sort_data[0], sort_data[1], xerr=pres_err, yerr=1, c='k')
ax1.scatter(sort_data[0], sort_data[1], s=15, c='k')
ax1.set_ylabel('Channel', fontsize=fnt)
ax1.set_title('Features of Alpha Decay Peak', fontsize=fnt+2)

ax2.errorbar(sort_data[0], sort_data[2], xerr=pres_err, yerr=1, c='r')
ax2.scatter(sort_data[0], sort_data[2], s=15, c='r')
ax2.set_xlabel('Pressure (Torr)', fontsize=fnt)
ax2.set_ylabel('FWHM (Ch#)', fontsize=fnt)

plt.setp(ax1.get_xticklabels(), visible=False)

ax1.grid()
ax2.grid()
#plt.show()
plt.savefig('alf_peak_feat.png')
plt.close()

### Plot Average Energy ###

plt.errorbar(Deff, avg_eng, yerr=en_err, color='k')
plt.scatter(Deff, avg_eng, s=25, color='k')
plt.scatter([Deff_MR], [eng_MR], c='r', s=150, marker='*', label='( %.2f cm, %.2f MeV)' % (Deff_MR, eng_MR))

plt.xlabel(r'$d_{eff}$ (cm)', fontsize=fnt)
plt.ylabel('Mean Energy (MeV)', fontsize=fnt)
plt.title('Alpha Particle Energy', fontsize=fnt+2)

plt.legend()
plt.grid(b=True)
#plt.show()
plt.savefig('alph_energy.png')
plt.close()

### Plot Stopping Power ###

plt.plot(Eng_new, stop_pwr, color='k')
plt.errorbar(avg_eng, stp_pwr_intp(avg_eng), yerr=stp_err, color='k', marker='o', linestyle='None')
plt.scatter([eng_MR], [stp_MR], c='r', s=150, marker='*', label='( %.2f MeV, %.2f MeV/cm)' % (eng_MR, stp_MR))

plt.xlabel(r'Mean Energy $(MeV)$', fontsize=fnt)
plt.ylabel(r'$\langle - \frac{dE}{dx} \rangle \ (MeV/cm)$', fontsize=fnt)
plt.title('Stopping Power', fontsize=fnt+2)

plt.legend()
plt.grid()
#plt.show()
plt.savefig('stop_pwr.png')
plt.close()

### Plot Peak Counts ###

plt.plot([0, 4], [stop_Ceff, stop_Ceff], c='b', ls='--', label=r'$50\%$ Count Rate = '+r'(%.1f $\pm$ %.1f) $s^{-1}$' % (stop_Ceff, stop_Ceff_err))
plt.plot([stop_Deff, stop_Deff], [0, 10], c='r', ls='--', label=r'Mean Range = '+r'(%.2f $\pm$ %.2f) $cm$' % (stop_Deff, stop_Deff_err))
plt.scatter([Deff_MR], [9.1], c='r', s=150, marker='*', label=r'$d_{eff}^{*}$')

plt.fill_between([0, 4], [stop_Ceff+stop_Ceff_err, stop_Ceff+stop_Ceff_err], [stop_Ceff-stop_Ceff_err, stop_Ceff-stop_Ceff_err], color='lightsteelblue')
plt.fill_between([stop_Deff-stop_Deff_err, stop_Deff+stop_Deff_err], [10, 10], [0, 0], color='salmon')
plt.errorbar(Deff, Crat, xerr=Derr, yerr=Cerr, c='k')

plt.ylabel('Count Rate (#/s)', fontsize=fnt)
plt.xlabel(r'$d_{eff}$ (cm)', fontsize=fnt)
plt.title('Alpha Decay Counts', fontsize=fnt+2)

plt.grid(b=True)
plt.legend()
#plt.show()
plt.savefig('alf_peak_cnts.png')
plt.close()
