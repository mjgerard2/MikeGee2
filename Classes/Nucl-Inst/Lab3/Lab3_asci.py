#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 21:41:10 2019

@author: michael
"""
import matplotlib.pyplot as plt
import numpy as np

loc = 'lab3_spectra1.ASC'

f = open(loc, 'r')

lines = f.readlines()

npts=1024
chan = np.empty(npts)
sign = np.empty(npts)

for l, line in enumerate(lines[12::]):
    line = line.strip()
    line = line.split() 
    chan[l] = float(line[0])
    sign[l] = float(line[1])

f.close()

fnt=16
scl = 1. / 180
sign = scl * sign

chan_we = chan[17:125]
sign_we = sign[17:125]

chan_Li = chan[124:230]
sign_Li = sign[124:230]

chan_Li_ = chan[230:270]
sign_Li_ = sign[230:270]

plt.plot(chan[0:17], sign[0:17], c='k')
plt.plot(chan_we, sign_we, label='Wall Effect')
#plt.plot(chan[60:125], sign[60:125], label='Wall 2')
plt.plot(chan_Li, sign_Li, label=r'$^{7}_{3}Li^{*}$')
plt.plot(chan_Li_, sign_Li_, label=r'$^{7}_{3}Li$')
plt.plot(chan[270:400], sign[270:400], c='k')

plt.xlabel('Channel', fontsize=fnt)
plt.ylabel(r'Counts ($\# /s}$)', fontsize=fnt)
plt.title('Spectrum', fontsize=fnt+2)

plt.legend()
#plt.show()
plt.savefig('Spectra_1.png')
plt.close()

print(np.max(sign_Li_) / np.max(sign_Li))
print(np.trapz(sign_Li_) / np.trapz(sign_Li))

'''
loc = 'lab3wk2.ASC'

f = open(loc, 'r')

lines = f.readlines()

npts=1024
chan = np.empty(npts)
sign = np.empty(npts)

for l, line in enumerate(lines[12::]):
    line = line.strip()
    line = line.split() 
    chan[l] = float(line[0])
    sign[l] = float(line[1])

f.close()

fnt=16
scl = 1. / 20
sign = scl * sign

plt.plot(chan, sign)

plt.xlabel('Channel', fontsize=fnt)
plt.ylabel(r'Counts ($\# /s}$)', fontsize=fnt)
plt.title('Spectrum', fontsize=fnt+2)

#plt.legend()
#plt.show()
plt.savefig('Spectra_2.png')
plt.close()
'''