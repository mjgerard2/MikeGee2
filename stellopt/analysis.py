#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 13:53:55 2019

@author: michael
"""

from matplotlib import animation as an  
from matplotlib.colors import LogNorm
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
    
    def find_crnt_profile(self, crnt_pro):
        """ Find index of a particular current profile.
        
       Parameters
        --------------------------------------------------------
        crnt_pro : array
            Current profile whose index we want to identify
        """
        self.base_idx = np.argmin( np.sum((self.crt_profile - crnt_pro)**2, axis=1) )
        
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
            
        eps_min = np.min(self.eps_ordered[0::, 0:26])
        eps_max = np.max(self.eps_ordered[0::, 0:26])
            
        fnt = 14
        
        job_dom = np.arange(self.jobs)
        ns_dom = np.linspace(0, 1, self.ns)
        
        fig, axs = plt.subplots(1)
        s_map = axs.pcolormesh(ns_dom, job_dom, self.eps_ordered, norm=LogNorm(vmin=eps_min, vmax=eps_max))
        
        cbar = fig.colorbar(s_map, ax=axs)
        cbar.ax.set_ylabel(r'$\epsilon_{eff}^{3/2}$', fontsize=fnt)
        
        axs.set_xlabel(r'$\frac{r}{a}$', fontsize=fnt)
        axs.set_ylabel('Job', fontsize=fnt)
        
        plt.tight_layout()
        #plt.show()
        plt.savefig('eps_eff_scan.png')
        plt.close()
        
    def eps_extrema(self):
        """ Find the index of the maximum and minimum epsilon effective values 
        for r/a < 0.2.
        """
        ind_min = np.argmin(self.eps_profile[0::, 0:26])
        ind_max = np.argmax(self.eps_profile[0::, 0:26])
        
        self.eps_min = int((ind_min - (ind_min % 25)) / 25)
        self.eps_max = int((ind_max - (ind_max % 25)) / 25)
                
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
        
    def animate_search(self, name, slow=False):
        """ Creates and animated stream through the ordered epsilon effective 
        profiles, which helps with a manuel search of favorable profiles.
        
        Parameters
        --------------------------------------------------------
        name : str
            name of animation to generate
        """        
        eps_min = np.min(self.eps_ordered[0::, 0:26])
        eps_max = np.max(self.eps_ordered[0::, 0:26])
        
        eps_argmin = np.argmin(self.eps_ordered[0::, 0:26]) % 26
        eps_argmax = np.argmax(self.eps_ordered[0::, 0:26]) % 26
        
        fnt = 14
        fig, axs = plt.subplots(2, 1, sharex=True)
        
        plt.subplots_adjust(hspace=0.05, right=0.8)
        
        axs[0].plot(self.s_dom, self.eps_ordered[0], ls='--', c='k', label='QHS')
        axs[0].plot(self.s_dom[eps_argmax], eps_max, '^', markersize=7, markerfacecolor='None', markeredgecolor='b', label=r'$Max\left(\epsilon_{eff}^{3/2}\right)$')
        axs[0].plot(self.s_dom, self.eps_profile[self.eps_min], ls='--', c='b')
        axs[0].plot(self.s_dom[eps_argmin], eps_min, 's', markersize=7, markerfacecolor='None', markeredgecolor='g', label=r'$Min\left(\epsilon_{eff}^{3/2}\right)$')
        axs[0].plot(self.s_dom, self.eps_profile[self.eps_max], ls='--', c='g')

        axs[0].set_yscale('log')
        axs[0].set_ylabel(r'$\epsilon_{eff}^{3/2}$', fontsize=fnt)
        axs[0].legend(loc='upper right')
        
        s_map = axs[1].pcolormesh(self.s_dom, np.arange(self.jobs), self.eps_ordered, norm=LogNorm(vmin=eps_min, vmax=eps_max))
        cbar_ax = fig.add_axes([.83, 0.15, 0.03, 0.7])
        cbar = fig.colorbar(s_map, cax=cbar_ax)
        cbar.ax.set_ylabel(r'$\epsilon_{eff}^{3/2}$', fontsize=fnt)
        
        axs[1].set_xlabel(r'$\frac{r}{a}$', fontsize=fnt)
        axs[1].set_ylabel('Job', fontsize=fnt)
        
        artist = []
        for e, eps in enumerate(self.eps_ordered):
            ax_1, = axs[1].plot([0, 1], [e, e], color='r')
            ax_2, = axs[0].plot(self.s_dom, eps, color='r')
            
            crnts = 1e-3 * self.crt_ordered[e]
            text = 'Job {0}: ({1}, {2}, {3}, {4}, {5}, {6})'.format(e, *crnts)+r'$\cdot 10^{3}$'
            texts1 = axs[0].text(0.5, 1.1, text, bbox={'facecolor':'white', 'alpha':1, 'pad':5}, 
                 transform=axs[0].transAxes, ha='center', fontsize=fnt-2, animated=True)
                    
            artist.append([ax_1, ax_2, texts1])
        
        if slow:
            ani = an.ArtistAnimation(fig, artist[0:100])
            writer = an.FFMpegWriter(fps=1, codec='h264')
            ani.save(name+'_slow.mp4', writer=writer)
            
            plt.close()
            
        else:
            ani = an.ArtistAnimation(fig, artist)
            writer = an.FFMpegWriter(fps=25, codec='h264')
            ani.save(name+'.mp4', writer=writer)
            
            plt.close()
        
    def crnt_profile_prams(self):
        """ Calculate the mirror value, inversion value and entropy of all
        current configurations, as well as average and std of the epsilon
        effective profiles.
        
        Very Much so In Devolopment
        """
        eps_noEnt = np.zeros(self.ns)
        eps_pram = {}     
        for idx, crt in enumerate(self.crt_ordered):
            crt_evn = np.array([crt[0], crt[2], crt[4]])
            crt_odd = np.array([crt[5], crt[3], crt[1]])
            
            inv = int( np.sum( (crt_evn - crt_odd)**2 ) / 256e6 )
            mrr = int( np.sum( (crt[0:3] - crt[3::])**2 ) / 256e6 )
            ent = len(np.unique(crt))-1
            
            if ent == 0:
                eps_noEnt = eps_noEnt + self.eps_ordered[idx]
            else:
                key = '{0}-{1}-{2}'.format(ent, mrr, inv)
                if key in eps_pram:
                    eps_pram[key] = np.append(eps_pram[key], np.log10(self.eps_ordered[idx][0:26]))
                else:
                    eps_pram[key] = np.log10(self.eps_ordered[idx][0:26])
                
        eps_noEnt = eps_noEnt / 3
        eps_avg = np.ones((2, 13, 13))
        eps_std = np.ones((2, 13, 13))
        for et in range(2):
            for mr in range(13):
                for iv in range(13):
                    key = '{0}-{1}-{2}'.format(et+1, mr, iv)
                    if key in eps_pram:
                        eps_avg[et,mr,iv] = np.mean(eps_pram[key])
                        eps_std[et,mr,iv] = np.std(eps_pram[key])
        
        
        ### Color Map ###
        fnt = 14
        fig, axs = plt.subplots()
        
        axs.set_xlabel('Mirror Value', fontsize=fnt)
        axs.set_ylabel('Inversion Value', fontsize=fnt)

        s_map = axs.pcolormesh(np.arange(13), np.arange(13), eps_avg[1], vmin=-3, vmax=-2)
        cbar = fig.colorbar(s_map, ax=axs)
        #cbar.ax.set_ylabel(r'$\log\left(\epsilon_{eff}\right)$', fontsize=fnt)
        
        plt.show()
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
name = os.path.join(dirc, 'current_data_new.h5')

datum = data(name)
datum.eps_extrema()

compare = 'Min'

if compare == 'QHS':
    datum.find_crnt_profile(np.zeros(6))
    datum.eps_compare(init=datum.base_idx)
    
    datum.animate_search('QHS_compare')
    datum.animate_search('QHS_compare', slow=True)
    
if compare == 'Min':
    datum.eps_compare(init=datum.eps_min)
    
    datum.animate_search('Min_compare')
    datum.animate_search('Min_compare', slow=True)
    