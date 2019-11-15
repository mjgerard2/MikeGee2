# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 09:33:39 2019

@author: mjgerard2
"""
import numpy as np

def DebyeLen(n, T, q=1.6022e-19, k=1.3807e-23, ep=8.8542e-12):
    return np.sqrt( (k * T * ep) / (n * q * q) )

def PlasPram(n, d):
    return n * d**3

def PlasFreq(n, m, q=1.6022e-19, ep=8.8542e-12):
    return np.sqrt( (n * q * q) / (ep * m) )

def VelTherm(m, T, k=1.3807e-23):
    return np.sqrt( (T * k) / m )

n_e = 1e32
n_i = n_e

T_e = 10**7.2
T_i = T_e


m_e = 9.1094e-31
m_p = 1.6726e-27

Z = 1
m_i = Z * m_p

d_e = DebyeLen(n_e, T_e)
d_i = DebyeLen(n_i, T_i)
 
N_e = PlasPram(n_e, d_e)
N_i = PlasPram(n_i, d_i)

Om_e = PlasFreq(n_e, m_e)
Om_i = PlasFreq(n_i, m_i)

V_e = VelTherm(m_e, T_e)
V_i = VelTherm(m_i, T_i)

print('Elec Debye Length: {0}\n'
      'Ion Debye Length: {1}\n\n'
      'Elec Parameter: {2}\n'
      'Ion Parameter: {3}\n\n'
      'Elec Frequency: {4}\n'
      'Ion Frequency: {5}\n\n'
      'Elec Velocity: {6}\n'
      'Ion Velocity: {7}'.format(d_e, d_i, N_e, N_i, Om_e, Om_i, V_e, V_i))