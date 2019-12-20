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
        eps_init_inv = 1. /  eps_init
        crt_init = self.crt_profile[init]
        
        self.eps_ordered = np.empty((self.jobs, self.ns))
        self.crt_ordered = np.empty((self.jobs, 6))
        self.idx_ordered = np.empty(self.jobs)
        
        indx = np.arange(self.jobs)
        indx = np.delete(indx, init)
        
        self.eps_ordered[0] = eps_init
        self.crt_ordered[0] = crt_init
        self.idx_ordered[0] = init
        
        compare = np.empty((len(indx), 2))
        for j, i in enumerate(indx):
            eps_chk = np.mean( np.abs(self.eps_profile[i] - eps_init) * eps_init_inv )
            compare[j] = np.array([eps_chk, i])
        
        compare = np.array( sorted(compare, key=lambda x: x[0]) )
        
        for c, com in enumerate(compare):
            ind = int(com[1])
            
            self.eps_ordered[c] = self.eps_profile[ind]
            self.crt_ordered[c] = self.crt_profile[ind]
            self.idx_ordered[c] = ind
                                    
    def save_order(self, name):
        """ Save a set of ordered profiles.
        
        Parameters
        --------------------------------------------------------
        name : str
            name of hdf5 file generated
        """
        try:
            self.crt_ordered
        except:
            self.eps_compare()
        
        f = hf.File(name+'.h5', 'w')
        
        f.create_dataset('current profile', data=self.crt_ordered)
        f.create_dataset('epsilon profile', data=self.eps_ordered)
        
        f.close()
            
    def color_map(self):
        """ Generate color map of epsilon effective profiles, ordered 
        by their similarity to the chosen initial epsilon profile. 
        """
        try:
            self.eps_ordered
        except:
            self.eps_compare()
            
        eps_min = np.min(self.eps_ordered[np.nonzero(self.eps_ordered[0::, 0:26])])
        eps_max = np.max(self.eps_ordered[0::, 0:26])
            
        fnt = 14
        
        job_dom = np.arange(self.jobs)
        ns_dom = np.linspace(0, 1, self.ns)
        
        fig, axs = plt.subplots(1)
        s_map = axs.pcolormesh(ns_dom, job_dom, self.eps_ordered, norm=LogNorm(vmin=eps_min, vmax=eps_max), cmap='jet')
        
        cbar = fig.colorbar(s_map, ax=axs)
        cbar.ax.set_ylabel(r'$\epsilon_{eff}^{3/2}$', fontsize=fnt)
        
        axs.set_xlabel(r'$\frac{r}{a}$', fontsize=fnt)
        axs.set_ylabel('Job', fontsize=fnt)
        
        plt.tight_layout()
        #plt.show()
        plt.savefig('eps_eff_scan.png')
        plt.close()
        
    def plot_10(self, name):
        """ Plot the first 10 profiles in the list.
        
        Parameters
        --------------------------------------------------------
        name : str
            Name of figure generated
        """
        self.find_crnt_profile(np.zeros(6))
        qhs_profile = self.eps_profile[self.base_idx]
        
        fnt=16
        
        plt.yscale('log')
        plt.xlabel('r/a', fontsize=fnt)
        plt.ylabel(r'$\epsilon_{eff}^{3/2}$', fontsize=fnt)
        
        plt.plot(self.s_dom, qhs_profile, ls='--', c='k', label='QHS')
        
        for i in range(8):
            idx = self.jobs - (1 + i)
            crnts = -6.25e-6 * self.crt_profile[idx]
            plt.plot(self.s_dom, self.eps_profile[idx], label='(%.2f, %.2f, %.2f, %.2f, %.2f, %.2f)' % (crnts[0], crnts[1], crnts[2], crnts[3], crnts[4], crnts[5]))
        
        plt.legend(fontsize=fnt-5)
        plt.tight_layout()
        plt.savefig(name+'.png')
        plt.close()
        
    def eps_extrema(self):
        """ Find the index of the maximum and minimum epsilon effective values 
        for r/a < 0.2.
        """
        ind_min = np.argmin(self.eps_profile[0::, 0:26])
        ind_max = np.argmax(self.eps_profile[0::, 0:26])
                
        self.eps_min = int((ind_min - (ind_min % 26)) / 26)
        self.eps_max = int((ind_max - (ind_max % 26)) / 26)
                
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
        
    def qhs_difference(self, name):
        """ Plot the percent difference for all configurations agains QHS.
        
        Parameters
        --------------------------------------------------------
        name : str
            name of plot generated
        """
        self.find_crnt_profile(np.zeros(6))
        qhs_profile = self.eps_profile[self.base_idx][0:26]
        
        compare = np.empty(self.jobs)
        for ind, eps in enumerate(self.eps_profile[0::, 0:26]):
            compare[ind] = np.mean((eps - qhs_profile) / qhs_profile)
        
        qhs_diff = np.sort(compare)
        qhs_diff = qhs_diff[qhs_diff < 1]
        x_dom = np.arange(len(qhs_diff))
                
        fnt = 16
        
        plt.xlabel('Configuration', fontsize=fnt)
        plt.ylabel('% Difference', fontsize=fnt)
        plt.plot([0, np.max(x_dom)], [0, 0], ls='--', label='QHS')
        plt.plot(x_dom, qhs_diff)
        
        plt.legend(loc='upper left', fontsize=fnt)
        plt.tight_layout()
        plt.savefig(name+'.png')
        plt.close()
        
    def still_search(self, name):
        """ Creates still image of animation.
        
        Parameters
        --------------------------------------------------------
        name : str
            Name of figure generated
        """
        eps_min = np.min(self.eps_profile[np.nonzero(self.eps_profile[0::, 0:26])])
        eps_max = np.max(self.eps_profile[0::, 0:26])
        
        fnt = 14
        fig, axs = plt.subplots(2, 1, sharex=True)
        
        plt.subplots_adjust(hspace=0.05, right=0.8)
        
        self.find_crnt_profile(np.zeros(6))
        qhs_profile = self.eps_profile[self.base_idx]
        axs[1].plot([0,1], [self.base_idx, self.base_idx], c='k', ls='--')
        
        find = [[16000, 16000, 16000, -16000, -16000, -16000],
                [-16000, -16000, -16000, 16000, 16000, 16000],
                [16000, -16000, -16000, -16000, 16000, 16000],
                [-16000, 16000, 16000, 16000, -16000, -16000]]
        
        labs = ['Mirror Inv: (-.1, -.1, -.1, .1, .1, .1)',
                'Mirror: (.1, .1, .1, -.1, -.1, -.1)',
                'Flip14: (-.1, .1, .1, .1, -.1, -.1)',
                'Flip14 Inv: (.1, -.1, -.1, -.1, .1, .1)'                
                ]
        
        cols = ['b',
                'orange',
                'g',
                'c']
        
        for idx, f in enumerate(find):
            self.find_crnt_profile(f)
            axs[0].plot(self.s_dom, self.eps_profile[self.base_idx], label=labs[idx], c=cols[idx])
            axs[1].plot([0, 1], [self.base_idx, self.base_idx], c=cols[idx])
        
        axs[0].plot(self.s_dom, qhs_profile, ls='--', c='k', label='QHS')
        axs[0].plot(self.s_dom, self.eps_profile[self.eps_max], ls='--', c='r', label=r'$Max\left(\epsilon_{eff}^{3/2}\right)$')

        axs[0].set_yscale('log')
        axs[0].set_ylabel(r'$\epsilon_{eff}^{3/2}$', fontsize=fnt)
        leg = axs[0].legend()
        leg.set_bbox_to_anchor([1, -.1])
                
        s_map = axs[1].pcolormesh(self.s_dom, np.arange(self.jobs), self.eps_profile, norm=LogNorm(vmin=eps_min, vmax=eps_max), cmap='jet')
        cbar_ax = fig.add_axes([.83, 0.15, 0.03, 0.7])
        cbar = fig.colorbar(s_map, cax=cbar_ax)
        cbar.ax.set_ylabel(r'$\epsilon_{eff}^{3/2}$', fontsize=fnt)
        
        axs[1].set_xlabel(r'$\frac{r}{a}$', fontsize=fnt)
        axs[1].set_ylabel('Configuration', fontsize=fnt)
        
        axs[1].set_zorder(-1)
    
        plt.savefig(name+'.png')
        
    def animate_search(self, name, slow=False):
        """ Creates and animated stream through the ordered epsilon effective 
        profiles, which helps with a manuel search of favorable profiles.
        
        Parameters
        --------------------------------------------------------
        name : str
            name of animation to generate
        """        
        eps_min = np.min(self.eps_profile[0::, 0:26])
        eps_max = np.max(self.eps_profile[0::, 0:26])
        
        eps_argmin = np.argmin(self.eps_profile[0::, 0:26]) % 26
        eps_argmax = np.argmax(self.eps_profile[0::, 0:26]) % 26
        
        fnt = 14
        fig, axs = plt.subplots(2, 1, sharex=True)
        
        plt.subplots_adjust(hspace=0.05, right=0.8)
        
        self.find_crnt_profile(np.zeros(6))
        qhs_profile = self.eps_profile[self.base_idx]
        
        axs[0].plot(self.s_dom, qhs_profile, ls='--', c='k', label='QHS')
        axs[0].plot(self.s_dom[eps_argmax], eps_max, '^', markersize=7, markerfacecolor='None', markeredgecolor='b')
        axs[0].plot(self.s_dom, self.eps_profile[self.eps_max], ls='--', c='b', label=r'$Max\left(\epsilon_{eff}^{3/2}\right)$')
        axs[0].plot(self.s_dom[eps_argmin], eps_min, 's', markersize=7, markerfacecolor='None', markeredgecolor='g')
        axs[0].plot(self.s_dom, self.eps_profile[self.eps_min], ls='--', c='g', label=r'$Min\left(\epsilon_{eff}^{3/2}\right)$')

        axs[0].set_yscale('log')
        axs[0].set_ylabel(r'$\epsilon_{eff}^{3/2}$', fontsize=fnt)
        axs[0].legend(loc='upper right')
        
        s_map = axs[1].pcolormesh(self.s_dom, np.arange(self.jobs), self.eps_profile, norm=LogNorm(vmin=eps_min, vmax=eps_max))
        cbar_ax = fig.add_axes([.83, 0.15, 0.03, 0.7])
        cbar = fig.colorbar(s_map, cax=cbar_ax)
        cbar.ax.set_ylabel(r'$\epsilon_{eff}^{3/2}$', fontsize=fnt)
        
        axs[1].set_xlabel(r'$\frac{r}{a}$', fontsize=fnt)
        axs[1].set_ylabel('Configuration', fontsize=fnt)
        
        artist = []
        
        if slow:
            profs = self.eps_profile[0:100]
        else:
            profs = self.eps_profile
            
        for e, eps in enumerate(profs):
            ax_1, = axs[1].plot([0, 1], [e, e], color='r')
            ax_2, = axs[0].plot(self.s_dom, eps, color='r')
            
            crnts = -6.25e-6 * self.crt_profile[e]
            text = 'Aux. Coils: (%.2f, %.2f, %.2f, %.2f, %.2f, %.2f)' % (crnts[0], crnts[1], crnts[2], crnts[3], crnts[4], crnts[5])
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
        eps_pram = {}     
        for idx, crt in enumerate(self.crt_profile):
            crt_evn = np.array([crt[0], crt[2], crt[4]])
            crt_odd = np.array([crt[5], crt[3], crt[1]])
            
            inv = int(np.sum((crt_evn - crt_odd)**2 / 16e6))
            mrr = int(np.sum((crt[0:3] - crt[3::])**2 / 16e6))
            
            key = '{0}-{1}'.format(mrr, inv)
            if key in eps_pram:
                eps_pram[key] = np.append(eps_pram[key], self.eps_profile[idx])#[0:26])
            else:
                eps_pram[key] = self.eps_profile[idx]#[0:26]
        '''      
        clr = ['tab:blue', 'tab:green', 'tab:red', 'tab:orange', 'tab:cyan', 'tab:purple', 'tab:olive']
        for et in range(6):
            plt.yscale('log')
            eps_avg = []
            eps_err = []
            dom = []
            for mriv in range(192):
                key = '{0}-{1}-{1}'.format(et+1, mriv)
                if key in eps_pram:
                    eps_avg = np.append(eps_avg, np.mean(eps_pram[key]))
                    eps_err = np.append(eps_err, np.std(eps_pram[key]))
                    dom = np.append(dom, mriv)
            
            plt.scatter(dom, eps_avg, c=clr[et], s=3)
            
        plt.savefig('Center_line.png')
        plt.close()
        '''
        eps_pro = self.eps_profile[0]
        eps_pro_inv = 1. / self.eps_profile[0]
        
        eps_avg = np.ones((192, 192))
        eps_std = np.ones((192, 192))
        for mr in range(192):
            for iv in range(192):
                key = '{0}-{1}'.format(mr, iv)
                if key in eps_pram:
                    js = int( len(eps_pram[key]) / 127 )
                    eps_pros = eps_pram[key].reshape(js, 127)
                    eps_diff = np.empty(js)
                    for i in range(js):
                        eps_diff[i] = np.mean( np.abs(eps_pro - eps_pros[i]) * eps_pro_inv * 100) #np.mean(eps_pram[key])
                    
                    eps_avg[mr,iv] = np.mean(eps_diff)
                    eps_std[mr,iv] = np.std(eps_diff)
                        
        ### Color Map ###
        fnt = 14
        fig, axs = plt.subplots()
        '''
        axs.set_xlabel('Mirror Value', fontsize=fnt)
        axs.set_ylabel('Inversion Value', fontsize=fnt)
        
        s_map = axs.pcolormesh(np.arange(192), np.arange(192), eps_avg, norm=LogNorm(vmin=np.min(eps_avg), vmax=np.max(eps_avg)), cmap='jet')
        cbar = fig.colorbar(s_map, ax=axs)
        cbar.ax.set_ylabel(r'$\%$ Diff w/ QHS', fontsize=fnt)
        
        #plt.show()
        plt.savefig('QHS_diff.png')
        plt.close()
        
        
        
        fig, axs = plt.subplots()
        
        axs.set_xlabel('Mirror Value', fontsize=fnt)
        axs.set_ylabel('Inversion Value', fontsize=fnt)
        
        s_map = axs.pcolormesh(np.arange(192), np.arange(192), eps_std, norm=LogNorm(vmin=np.min(eps_std), vmax=np.max(eps_std)), cmap='jet')
        cbar = fig.colorbar(s_map, ax=axs)
        cbar.ax.set_ylabel(r'STD of $\%$ Diff w/ QHS', fontsize=fnt)
        
        #plt.show()
        plt.savefig('QHS_diff_std.png')
        plt.close()
        '''
        eps_diff = np.ones((192, 192))
        for i in range(192):
            for j in range(192):
                if eps_avg[i,j] != 0:
                    eps_diff[i,j] = (eps_avg[i,j] - eps_std[i,j]) / eps_avg[i,j]**2
        
        fig, axs = plt.subplots()
        
        axs.set_xlabel('Mirror Value', fontsize=fnt)
        axs.set_ylabel('Inversion Value', fontsize=fnt)
        
        s_map = axs.pcolormesh(np.arange(192), np.arange(192), eps_diff, vmin=np.min(eps_diff), vmax=np.max(eps_diff), cmap='jet')
        cbar = fig.colorbar(s_map, ax=axs)
        cbar.ax.set_ylabel(r'Diff of STD of $\%$ Diff w/ Min', fontsize=fnt)
        
        #plt.show()
        plt.savefig('Min_diff_diff.png')
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
name_og = os.path.join(dirc, 'current_data.h5')
name_min = os.path.join(dirc, 'EpsOrder_MIN.h5')
name_qhs = os.path.join(dirc, 'EpsOrder_QHS.h5')

datum = data(name_og)

#datum.plot_10('plot_10')
'''
print('Extrema Ordering')
datum.eps_extrema()
datum.eps_compare(init=datum.eps_min)
datum.save_order('EpsOrder_MIN')

print('QHS Ordering')
datum.find_crnt_profile(np.zeros(6))
datum.eps_compare(init=datum.base_idx)
datum.save_order('EpsOrder_QHS')
'''
compare = 'Min'

if compare == 'QHS':
    anal_qhs = data(name_qhs)
    #anal_qhs.eps_extrema()
    #anal_qhs.animate_search('QHS_compare')
    #anal_qhs.animate_search('QHS_compare', slow=True)
    #anal_qhs.color_map()
    anal_qhs.crnt_profile_prams()
    
    
elif compare == 'Min':
    anal_min = data(name_min)
    #anal_min.eps_extrema()
    #anal_min.animate_search('Min_compare')
    #anal_min.animate_search('Min_compare', slow=True)
    #anal_min.qhs_difference('qhs_diff')
    #anal_min.plot_10('plot_10_max')
    #anal_min.still_search('Sig_Configs')
    anal_min.crnt_profile_prams()
