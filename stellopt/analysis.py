#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 13:53:55 2019

@author: michael
"""

from matplotlib import animation as an  
import matplotlib.pyplot as plt

import numpy as np
import h5py as hf

import os

class data():
    
    def __init__(self, loc):
        
        f = hf.File(loc, 'r')
        
        self.crt_profile = f['current profile'][:]
        self.eps_profile = f['epsilon profile'][:]
        
        self.jobs = self.crt_profile.shape[0]
        self.ns = self.eps_profile.shape[1]
        self.s_dom = np.linspace(0, 1, self.ns)
        
        f.close()
        
    def eps_compare(self, init=0):
        """ Takes the indexed epsilon effective profile and orders all other 
        epsilon profiles by their similarity to the initial one.
        
        Parameters
        --------------------------------------------------------
        init : int
            index for the initial epsilon effective profile
        """
        eps_init = self.eps_profile[init]
        crt_init = self.crt_profile[init]
        
        self.eps_ordered = np.empty((self.jobs, self.ns))
        self.crt_ordered = np.empty((self.jobs, 6))
        self.idx_ordered = np.empty(self.jobs)
        
        indx = np.arange(self.jobs)
        indx = np.delete(indx, init)
        
        self.eps_ordered[0] = eps_init
        self.crt_ordered[0] = crt_init
        self.idx_ordered[0] = init
        
        for c in np.arange(1, self.jobs):
            compare = np.empty(len(indx))
            for j, i in enumerate(indx):
                eps_chk = np.sum( (eps_init - self.eps_profile[i])**2 )
                compare[j] = eps_chk
            
            j_min = np.argmin(compare)
            i_min = indx[j_min]
                        
            self.eps_ordered[c] = self.eps_profile[i_min]
            self.crt_ordered[c] = self.crt_profile[i_min]
            self.idx_ordered[c] = i_min
                        
            indx = np.delete(indx, j_min)
            
    def color_map(self):
        """ Generate color map of epsilon effective profiles, ordered 
        by their similarity to the chosen initial epsilon profile. 
        """
        try:
            self.eps_ordered
        except:
            self.eps_compare()
            
        fnt = 14
        
        job_dom = np.arange(self.jobs)
        ns_dom = np.linspace(0, 1, self.ns)
        
        fig, axs = plt.subplots(1)
        s_map = axs.pcolormesh(ns_dom, job_dom, np.log10(self.eps_ordered))
        
        cbar = fig.colorbar(s_map, ax=axs)
        cbar.ax.set_ylabel(r'$\log\left(\epsilon_{eff}\right)$', fontsize=fnt, rotation=0)
        
        axs.set_xlabel(r'$\frac{r}{a}$', fontsize=fnt)
        axs.set_ylabel('Job', fontsize=fnt)
        
        plt.tight_layout()
        #plt.show()
        plt.savefig('eps_eff_scan.png')
        plt.close()
        
    def plot_profile(self, ax, idx):
        """ Plot individual epsilon effective profile
        
        Parameters
        --------------------------------------------------------
        ax : obj
            pyplot axis object
            
        idx : int
            index of epsilon effective profile, taken from self.eps_profile
        """
        ax.plot(self.s_dom, self.eps_ordered[idx])
        
    def find_config(self, idx):
        """ Finds current configuration for a particular job in the ordered
        job list.
        
        Parameters
        --------------------------------------------------------
        idx : int
            Job index to identify
        """
        try:
            self.idx_ordered
        except:
            self.eps_compare()
            
        indx = int(self.idx_ordered[idx])
        self.eps_compare(init=indx)
        
    def animate_search(self):
        job_dom = np.arange(self.jobs)

        fnt = 14
        fig, axs = plt.subplots(2, 1, sharex=True)
        
        plt.subplots_adjust(hspace=0.05, right=0.8)
        #plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
        
        axs[0].set_ylabel(r'$\epsilon_{eff}$', fontsize=fnt)
        
        s_map = axs[1].pcolormesh(self.s_dom, job_dom, np.log10(self.eps_ordered))
        cbar_ax = fig.add_axes([.83, 0.15, 0.03, 0.7])
        cbar = fig.colorbar(s_map, cax=cbar_ax)
        cbar.ax.set_ylabel(r'$\log\left(\epsilon_{eff}\right)$', fontsize=fnt)
        
        axs[1].set_xlabel(r'$\frac{r}{a}$', fontsize=fnt)
        axs[1].set_ylabel('Job', fontsize=fnt)
        
        artist = []
        for e, eps in enumerate(self.eps_ordered[0:100]):
            ax_1, = axs[1].plot([0, 1], [e, e], color='r')
            ax_2, = axs[0].plot(self.s_dom, np.log10(eps), color='r')
            
            text = 'Job {0}: {1}, {2}, {3}, {4}, {5}, {6}'.format(e, *self.crt_ordered[e])
            texts1 = axs[0].text(0.5, 1.1, text, bbox={'facecolor':'white', 'alpha':1, 'pad':5}, 
                 transform=axs[0].transAxes, ha='center', fontsize=fnt-2, animated=True)
                    
            artist.append([ax_1, ax_2, texts1])
                
        ani = an.ArtistAnimation(fig, artist)
        writer = an.FFMpegWriter(fps=1, codec='h264')
        ani.save('animated_search_slow.mp4', writer=writer)
        
        plt.close()
        
    def check_zero(self):
        """ Check if any epsilon effective values are equal to zero. 
        Then print the index of profile with a zero value and how many
        zeros there are.
        """
        for e, eps in enumerate(self.eps_profile):
            if any(eps == 0):
                eps_thresh = eps[eps == 0]
                print('{0}: {1}\n'.format(e, len(eps_thresh)))
            
dirc = os.getcwd()
name = os.path.join(dirc, 'current_data.h5')

datum = data(name)

datum.eps_compare()
datum.find_config(276)
datum.find_config(1)
#datum.animate_search()
datum.color_map()

#look = datum.eps_ordered