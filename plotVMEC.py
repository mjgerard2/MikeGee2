#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 13:44:55 2019

@author: michael
"""
import os

import h5py as hf
import numpy as np

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt

def find_nearest(array, value):
    return (np.abs(array - value)).argmin()

class s_func_cyln():
    
    def __init__(self, name):
        h5 = hf.File(name, 'r')
        
        self.s_surf = h5['s surf'][:]
        
        self.r_dom = h5['r dom'][:]
        self.p_dom = h5['p dom'][:]
        self.z_dom = h5['z dom'][:]
        
        self.r_max = np.max(self.r_dom)
        self.z_min = np.min(self.z_dom)
        self.z_max = np.max(self.z_dom)
                
        h5.close()
        
    def plot_torCS(self, loc, p):
        """ Plot S surfaces along toroidal cross section and save as png
        
        Parameters
        -----------------------------------------------------------------------
        loc : str
            location where figure will be saved
            
        p : float
            toroidal angle in radians
        """
        fnt=14
        p_ind = find_nearest(self.p_dom, p)
        
        fig, axs = plt.subplots(1)
        s_map = axs.pcolormesh(self.r_dom, self.z_dom, self.s_surf[0::, p_ind, 0::].T, vmin=0, vmax=1, cmap=plt.get_cmap('tab20c'))
        
        cbar = fig.colorbar(s_map, ax=axs)
        cbar.ax.set_ylabel(r'$\sqrt{s}$', fontsize=fnt, rotation=0)
        
        axs.set_xlabel('r (m)', fontsize=fnt)
        axs.set_ylabel('z (m)', fontsize=fnt)
        
        axs.set_title(r'$\theta$ = %.1f$^{\circ}$' % (self.p_dom[p_ind]*57.29577951), fontsize=fnt+2)
        
        plt.savefig( os.path.join( loc, 'S_Map_interp_{}.png'.format( int( self.p_dom[p_ind]*57.29577951 ) ) ) )
        plt.close()
        
    def plot_fluxSurf_3D(self, loc):
        """ Plot S sourface contours along entire toroidal domain and save as png
        
        Parameters
        -----------------------------------------------------------------------
        loc : str
            location where figure will be saved
            
        """
        fig = plt.figure(figsize=(16.,4.))
        ax = fig.add_subplot(111, projection = '3d')
        
        
        plt.ioff()
        fig1 = plt.figure()
        _ax = fig1.add_subplot(111)
        plt.ion()
        
        #plot every n-th element in the array
        for ii in range(4):
            for i in range(self.p_dom.shape[0]):
                alph = 0.5
                
                #get the contours of different s levels in (r,z) 
                cont = _ax.contour(self.r_dom, self.z_dom, self.s_surf[:,i,:].T, levels = [.99], linewidths = .1)#[0.1, 0.3, 0.5, 0.7, 0.95], linewidths = 0.1)  
                ver = cont.allsegs[0][0]
                
                #transform from (r,z,p) to (x,y,z)
                x = []
                y = []
                z = []
                for j in range(len(ver[:,0])):
                    _x = ver[j,0]*np.cos(self.p_dom[i]+ii*0.5*np.pi)
                    _y = ver[j,0]*np.sin(self.p_dom[i]+ii*0.5*np.pi)
                    _z = ver[j,1]
                    x.append(_x)
                    y.append(_y)
                    z.append(_z)
                ax.plot3D(x,y,z,color = 'k', lw = 1., alpha = alph)  
                
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        ax.set_zlabel('z [m]')
        ax.set_xlim([-self.r_max, self.r_max])
        ax.set_ylim([-self.r_max, self.r_max])
        ax.set_zlim([self.z_min, self.z_max])
        
        ax.view_init(azim=30, elev=60)
        
        ### Lance - plot Lorentz fource data here ###
        
        fig.savefig(os.path.join(loc, 'LCFS_3D.png'))
        plt.close()        
        
dirc = os.getcwd()
data_dirc = os.path.join(dirc, 'Data')
plot_dirc = os.path.join(dirc, 'Figures')

file_name = os.path.join(data_dirc, 's_surf_polar.h5')

Sval = s_func_cyln(file_name)
Sval.plot_fluxSurf_3D(plot_dirc)