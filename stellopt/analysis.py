#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 13:53:55 2019

@author: michael
"""

from matplotlib import animation as an  
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D

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
            idx = i#self.jobs - (1 + i)
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
        
    def qhs_extrema(self):
        """ Find the index of the epsilon effective profile whose percent comaprison
        with the qhs configuration is the smallest.
        """
        self.find_crnt_profile(np.zeros(6))
        qhs_profile = self.eps_profile[self.base_idx]
        qhs_profile_inv = 1. / qhs_profile
        
        compare = np.empty((self.jobs, 2))
        for ind, eps in enumerate(self.eps_profile):
            compare[ind] = [np.mean((eps - qhs_profile) * qhs_profile_inv), ind]
            
        compare = np.array( sorted(compare, key=lambda x: x[0]) )
        
        self.qhs_min = int(compare[0,1])
                
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
        qhs_profile = self.eps_profile[self.base_idx]#[0:26]
        qhs_profile_inv = 1. / qhs_profile
        
        compare = np.empty(self.jobs)
        for ind, eps in enumerate(self.eps_profile):#[0::, 0:26]):
            compare[ind] = np.mean((eps - qhs_profile) * qhs_profile_inv)
        
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
        
    def crnt_profile_prams_3D(self, pts=1000):
        tst = np.empty(pts)
        inv = np.empty(pts)
        mrr = np.empty(pts)
        
        for idx, crt in enumerate(self.crt_profile[0:pts]):
            inv_evn = np.array([crt[0], crt[2], crt[4]])
            inv_odd = np.array([crt[5], crt[1], crt[3]])
            
            tst_evn = np.array([crt[1], crt[4], crt[5]])
            tst_odd = np.array([crt[0], crt[2], crt[4]])
            
            tst[idx] = int(np.sum((tst_evn - tst_odd)**2 / 16e6))
            inv[idx] = int(np.sum((inv_evn - inv_odd)**2 / 16e6))
            mrr[idx] = int(np.sum((crt[0:3] - crt[3::])**2 / 16e6))
            
        fnt=14
        cm = plt.cm.get_cmap('jet')
        fig = plt.figure()
        axs = fig.add_subplot(111, projection='3d')
        
        axs.set_xlim(0, 200)
        axs.set_ylim(0, 200)
        axs.set_zlim(0, 200)
        
        axs.set_xlabel('Mirror Value', fontsize=fnt)
        axs.set_ylabel('Inversion Value', fontsize=fnt)
        
        #axs.set_yscale('log')
        #axs.set_xscale('log')
        
        mrr_inv = {}
        for i in range(pts):
            key = '{0} {1} {2}'.format(mrr[i], inv[i], tst[i])
            if key in mrr_inv:
                mrr_inv[key]+=1
            else:
                mrr_inv[key]=1
        
        mrr_seg = []
        inv_seg = []
        tst_seg = []
        cnts = []
        for key in mrr_inv:
            key_spt = key.split()
            mrr_seg.append( float(key_spt[0]) )
            inv_seg.append( float(key_spt[1]) )
            tst_seg.append( float(key_spt[2]) )
            
            cnts.append(float(mrr_inv[key]))
        
        ax = axs.scatter(mrr_seg, inv_seg, tst_seg, c=cnts, s=1, marker='X', norm=LogNorm(vmin=1, vmax=68), cmap=cm)
        fig.colorbar(ax)
        
        axs.view_init(0, 270)
        
        if pts==self.jobs:
            plt.savefig('crnt_pram.png')
        else:
            plt.savefig('crnt_pram3D_{}.png'.format(pts))
        plt.close()
    
    def crnt_profile_prams(self, pts=1000):
        inv = np.empty(pts)
        mrr = np.empty(pts)
        
        for idx, crt in enumerate(self.crt_profile[0:pts]):
            inv_evn = np.array([crt[0], crt[2], crt[4]])
            inv_odd = np.array([crt[5], crt[1], crt[3]])
            
            inv[idx] = int(np.sum((inv_evn - inv_odd)**2 / 16e6))
            mrr[idx] = int(np.sum((crt[0:3] - crt[3::])**2 / 16e6))
            
        fnt=14
        cm = plt.cm.get_cmap('jet')
        fig, axs = plt.subplots(1, 1)
        
        axs.set_xlim(0, 200)
        axs.set_ylim(0, 200)
        
        axs.set_xlabel('Mirror Value', fontsize=fnt)
        axs.set_ylabel('Inversion Value', fontsize=fnt)
        
        #axs.set_yscale('log')
        #axs.set_xscale('log')
        
        mrr_inv = {}
        for i in range(pts):
            key = '{0} {1}'.format(mrr[i], inv[i])
            if key in mrr_inv:
                mrr_inv[key]+=1
            else:
                mrr_inv[key]=1
        
        mrr_seg = np.array([])
        inv_seg = np.array([])
        cnts = np.array([])
        for key in mrr_inv:
            key_spt = key.split()            
            mrr_seg = np.append(mrr_seg, float(key_spt[0]))
            inv_seg = np.append(inv_seg, float(key_spt[1]))

            cnts = np.append(cnts, float(mrr_inv[key]))
        
        net_cnts_inv = 1. / np.sum(cnts)
        mrr_cen = net_cnts_inv * np.sum(mrr_seg * cnts)
        inv_cen = net_cnts_inv * np.sum(inv_seg * cnts)
                
        ax = axs.scatter(mrr_seg, inv_seg, c=cnts, s=1, norm=LogNorm(vmin=1, vmax=1144), cmap=cm)
        axs.scatter(36.05569212584262, 36.04722952659123, c='k', s=50, marker='X')
        axs.scatter(mrr_cen, inv_cen, c='tab:red', s=50, marker='X')
        fig.colorbar(ax)
                
        if pts==self.jobs:
            plt.savefig('crnt_pram.png')
        else:
            plt.savefig('crnt_pram_{}.png'.format(pts))
        plt.close()
    
    def animate_crnt_profile_prams(self, pts=1000):
        inv = np.empty(pts)
        mrr = np.empty(pts)
        
        for idx, crt in enumerate(self.crt_profile[0:pts]):
            crt_evn = np.array([crt[0], crt[2], crt[4]])
            crt_odd = np.array([crt[5], crt[1], crt[3]])
            
            inv[idx] = int(np.sum((crt_evn - crt_odd)**2 / 16e6))
            mrr[idx] = int(np.sum((crt[0:3] - crt[3::])**2 / 16e6))
            
        fnt=14
        cm = plt.cm.get_cmap('jet')
        fig, axs = plt.subplots(1, 1)
        
        axs.set_xlim(0, 200)
        axs.set_ylim(0, 200)
        
        axs.set_xlabel('Mirror Value', fontsize=fnt)
        axs.set_ylabel('Inversion Value', fontsize=fnt)
        
        #axs.set_yscale('log')
        #axs.set_xscale('log')
        
        artist = []
        segs = np.linspace(1, pts, 750, dtype=int)
        for s in segs:
            mrr_inv = {}
            for i in range(s):
                key = '{0} {1}'.format(mrr[i], inv[i])
                if key in mrr_inv:
                    mrr_inv[key]+=1
                else:
                    mrr_inv[key]=1
            
            mrr_seg = np.array([])
            inv_seg = np.array([])
            cnts = np.array([])
            for key in mrr_inv:
                key_spt = key.split()
                mrr_seg = np.append(mrr_seg, float(key_spt[0]))
                inv_seg = np.append(inv_seg, float(key_spt[1]))
    
                cnts = np.append(cnts, float(mrr_inv[key]))
                    
            net_cnts_inv = 1. / np.sum(cnts)
            mrr_cen = net_cnts_inv * np.sum(mrr_seg * cnts)
            inv_cen = net_cnts_inv * np.sum(inv_seg * cnts)
                
            ax = axs.scatter(mrr_seg, inv_seg, c=cnts, s=1, marker='X', norm=LogNorm(vmin=1, vmax=1144), cmap=cm)
            pt = axs.scatter(mrr_cen, inv_cen, c='tab:red', s=50, marker='X')
            artist.append([ax, pt,])
            
        fig.colorbar(ax)
        axs.scatter(36.05569212584262, 36.04722952659123, c='k', s=50, marker='X')
        
        ani = an.ArtistAnimation(fig, artist)
        writer = an.FFMpegWriter(fps=25, codec='h264')
        if pts == self.jobs:
            ani.save('crnt_pram.mp4', writer=writer)
        else:
            ani.save('crnt_pram_{}.mp4'.format(pts), writer=writer)
        
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
name_dif = os.path.join(dirc, 'EpsOrder_DIFF.h5')
'''
datum = data(name_min)

datum.plot_10('plot_10_dif')

print('QHS Lowest Difference')
datum.qhs_extrema()
datum.eps_compare(init=datum.qhs_min)
datum.save_order('EpsOrder_DIFF')


print('Extrema Ordering')
datum.eps_extrema()
datum.eps_compare(init=datum.eps_min)
datum.save_order('EpsOrder_MIN')


print('QHS Ordering')
datum.find_crnt_profile(np.zeros(6))
datum.eps_compare(init=datum.base_idx)
datum.save_order('EpsOrder_QHS')
'''
compare = 'Dif'

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

elif compare == 'Dif':
    anal_dif = data(name_dif)
    
    anal_dif.animate_crnt_profile_prams(pts=3706)#anal_dif.jobs)
    #for i in range(37):
    #    anal_dif.crnt_profile_prams(pts=100*(i+1))
    
    #anal_dif.animate_crnt_profile_prams(pts=3706)