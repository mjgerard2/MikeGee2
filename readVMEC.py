# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 14:38:30 2019

@author: mjgerard2
"""
import os

import h5py as hf
import numpy as np

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt

from netCDF4 import Dataset 
from matplotlib.lines import Line2D   

from scipy.interpolate import griddata
from scipy.interpolate import interp1d

def find_nearest(arr, val):
    """ Returns index of the element in an array closest to a given value.
    
    Parameters
    ------------------------------------------------------------------------
    arr : array
        Array to be searched
        
    val : float
        Value whose closest element's index in array you want returned
    """
    return (np.abs(arr - val)).argmin()

def save_data(name, data):
    """ Save data set as hdf5 file.
    
    Parameters
    ------------------------------------------------------------------------
    name : str
        Location and name where file will be saved
        
    data : dict
        Data to be saved
    """
    h5 = hf.File(name, 'w')
    
    for key in data:
        h5.create_dataset(key, data=data[key])
    
    h5.close()

class readVMEC:
    
    def __init__(self, file_name):
        rootgrp = Dataset(file_name, 'r')
        
        self.name = file_name
        
        self.nfp = rootgrp['/nfp'][0]
        self.mpol = rootgrp['/mpol'][0]
        self.ntor = rootgrp['/ntor'][0]
        
        self.xm = rootgrp['/xm'][:]
        self.xn = rootgrp['/xn'][:]
        self.md = len(self.xm)
        
        self.xm_nyq = rootgrp['/xm_nyq'][:]
        self.xn_nyq = rootgrp['/xn_nyq'][:]
        self.md_nyq = len(self.xm_nyq)
        
        self.rmnc = rootgrp['/rmnc'][:,:]
        self.zmns = rootgrp['/zmns'][:,:]
        self.lmns = rootgrp['/lmns'][:,:]
        
        self.bumnc = rootgrp['/bsubumnc'][:,:]
        self.bvmnc = rootgrp['/bsubvmnc'][:,:]
        self.bsmns = rootgrp['/bsubsmns'][:,:]
        self.bmnc = rootgrp['/bmnc'][:,:]
                    
        self.ns = rootgrp['/ns'][0]
        
        rootgrp.close()
        
        mult = 40
        self.u_num = self.mpol * mult
        self.v_num = self.ntor * mult
        
        self.s_dom = np.linspace(0, 1, self.ns)
        self.u_dom = np.linspace(0, 2*np.pi, self.u_num)
        self.v_dom = np.linspace(0, .5*np.pi, self.v_num)
        
    def polarCoord(self):
        """ Produces the cylindrical coordinates of magnetic flux surfaces, 
        indexed by VMEC flux coordinates [s,v,u].
        """
        pol, tor = np.meshgrid(self.u_dom, self.v_dom)
        
        pol_xm = np.dot(self.xm.reshape(self.md, 1), pol.reshape(1, self.v_num * self.u_num))
        tor_xn = np.dot(self.xn.reshape(self.md, 1), tor.reshape(1, self.v_num * self.u_num))
        
        cos_pol = np.cos(pol_xm)
        cos_tor = np.cos(tor_xn)

        sin_pol = np.sin(pol_xm)
        sin_tor = np.sin(tor_xn)
        
        self.cos_mu_nv = cos_pol*cos_tor + sin_pol*sin_tor
        self.sin_mu_nv = sin_pol*cos_tor - sin_tor*cos_pol
        
        self.R_coord = np.dot(self.rmnc, self.cos_mu_nv).reshape(self.ns, self.v_num, self.u_num)
        self.Z_coord = np.dot(self.zmns, self.sin_mu_nv).reshape(self.ns, self.v_num, self.u_num)
        
        self.dRu_coord = np.dot(-self.rmnc*self.xm, self.sin_mu_nv).reshape(self.ns, self.v_num, self.u_num)
        self.dRv_coord = np.dot(self.rmnc*self.xn, self.sin_mu_nv).reshape(self.ns, self.v_num, self.u_num) 
        
        self.dZu_coord = np.dot(self.zmns*self.xm, self.cos_mu_nv).reshape(self.ns, self.v_num, self.u_num)
        self.dZv_coord = np.dot(-self.zmns*self.xn, self.cos_mu_nv).reshape(self.ns, self.v_num, self.u_num) 
        
        ds = self.s_dom[1] - self.s_dom[0]
        dr = np.gradient(self.rmnc, ds, axis=0)
        dz = np.gradient(self.zmns, ds, axis=0)
        
        self.dRs_coord = np.dot(dr, self.cos_mu_nv).reshape(self.ns, self.v_num, self.u_num)
        self.dZs_coord = np.dot(dz, self.sin_mu_nv).reshape(self.ns, self.v_num, self.u_num)
        
    def polarCoord_Bfield(self):
        """ Produces the magnetic field in real space in cylindrical coordinates, 
        indexed by VMEC flux coordinates [s,v,u].        
        """
        try:
            self.dRu_coord
        except:
            self.polarCoord()

        Bs_coord = np.dot(self.bsmns, self.sin_mu_nv).reshape(self.ns, self.v_num, self.u_num)
        Bu_coord = np.dot(self.bumnc, self.cos_mu_nv).reshape(self.ns, self.v_num, self.u_num)
        Bv_coord = np.dot(self.bvmnc, self.cos_mu_nv).reshape(self.ns, self.v_num, self.u_num)
        
        B_norm = 1. / (self.dRs_coord * self.dZu_coord - self.dRu_coord * self.dZs_coord)
        Br_coord = (self.dZu_coord * Bs_coord - self.dZs_coord * Bu_coord) * B_norm 
        Bp_coord = ( ( (Bs_coord * (self.dRu_coord*self.dZv_coord - self.dRv_coord*self.dZu_coord) + Bu_coord * (self.dRv_coord*self.dZs_coord - self.dRs_coord*self.dZv_coord)) * B_norm ) + Bv_coord) / self.R_coord
        Bz_coord = (self.dRs_coord * Bu_coord - self.dRu_coord * Bs_coord) * B_norm
        
        self.B_field = np.stack((Br_coord, Bp_coord, Bz_coord), axis=3)
        
    def polarCoord_Bmod(self):
        """ Calculate mod B in real space, indexed by VMEC flux 
        coordintes [s,v,u].
        """
        try:
            self.cos_mu_nv
        except:
            self.polarCoord()
            
        self.B_mod = np.dot(self.bmnc, self.cos_mu_nv).reshape(self.ns, self.v_num, self.u_num)
        
    def fluxSurf(self, save=False, loc=os.getcwd(), r_res=0.001, z_res=0.001):
        """ Produce a 3D array of flux surface s values, indexed in cylindrical coordinates [r,p,z]
        
        Parameters
        ------------------------------------------------------------------------
        r_res : float
            r coordinate resolution in meters
            
        z_res : array
            z coordinate resolution in meters
        """
        try:
            self.R_coord
        except:
            self.polarCoord()
            
        r_min = np.min(self.R_coord)
        r_max = np.max(self.R_coord)
        r_num = int((r_max - r_min) / r_res)
        self.r_dom = np.linspace(r_min, r_max, r_num)
        
        z_min = np.min(self.Z_coord)
        z_max = np.max(self.Z_coord)
        z_num = int((z_max - z_min) / z_res)
        self.z_dom = np.linspace(z_min, z_max, z_num)
        
        r_grid, z_grid = np.meshgrid(self.z_dom, self.r_dom)
        
        self.s_surf = np.full((r_num, self.v_num, z_num), np.nan)
        
        for v_ind, v in enumerate(self.v_dom):
            r_tor = self.R_coord[0::, v_ind, 0::]
            z_tor = self.Z_coord[0::, v_ind, 0::]
            
            s_temp = np.full((r_num, z_num), np.nan)
            for s_ind, s in enumerate(self.s_dom):
                for u_ind, u in enumerate(self.u_dom):
                    r = r_tor[s_ind, u_ind]
                    z = z_tor[s_ind, u_ind]
                    
                    r_ind = find_nearest(self.r_dom, r)
                    z_ind = find_nearest(self.z_dom, z)
                    
                    s_temp[r_ind, z_ind] = s
                    
            s_temp = np.ma.masked_invalid(s_temp)
            s_temp_inv = ~s_temp.mask
            
            r_grid_mask = r_grid[s_temp_inv]
            z_grid_mask = z_grid[s_temp_inv]
            s_temp_mask = s_temp[s_temp_inv]
            
            s_temp_intp = griddata((z_grid_mask, r_grid_mask), s_temp_mask, (z_grid, r_grid))
            self.s_surf[0::, v_ind, 0::] = s_temp_intp
            
        if save:
            name = os.path.join(loc, 'flux_surface_HSX.h5')
            
            data = {'s surf' : self.s_surf,
                    'r dom' : self.r_dom,
                    'p dom' : self.v_dom,
                    'z dom' : self.z_dom}
            
            save_data(name, data)
        
    def cyldCoord_Bfield(self, save=False, loc=os.getcwd(), r_res=0.001, z_res=0.001):
        """ Interpolates the magnetic field in real space in cylindrical coordinates, 
        indexed by the cylindrical domain [r,p,z].        
        
        Parameters
        ------------------------------------------------------------------------
        r_res : float
            r coordinate resolution in meters
            
        z_res : array
            z coordinate resolution in meters
        """
        try:
            self.B_field
        except:
            self.polarCoord_Bfield()
        
        r_min = np.min(self.R_coord)
        r_max = np.max(self.R_coord)
        r_num = int((r_max - r_min) / r_res)
        self.r_dom = np.linspace(r_min, r_max, r_num)
        
        z_min = np.min(self.Z_coord)
        z_max = np.max(self.Z_coord)
        z_num = int((z_max - z_min) / z_res)
        self.z_dom = np.linspace(z_min, z_max, z_num)

        self.B_cyld = np.full((r_num, self.v_num, z_num, 3), np.nan)
        for s in range(1, self.ns):            
            for v in range(self.v_num):
                for u in range(self.u_num):
                    r = self.R_coord[s,v,u]
                    z = self.Z_coord[s,v,u]
                    
                    r_ind = find_nearest(self.r_dom, r)
                    z_ind = find_nearest(self.z_dom, z)
                    
                    self.B_cyld[r_ind, v, z_ind] = self.B_field[s,v,u]
        
        z_grd, r_grd = np.meshgrid(self.z_dom, self.r_dom)
        for v in range(self.v_num):                        
            Br = self.B_cyld[0::, v, 0::, 0]
            Bv = self.B_cyld[0::, v, 0::, 1]
            Bz = self.B_cyld[0::, v, 0::, 2]
            
            Br = np.ma.masked_invalid(Br)
            Br_inv = ~Br.mask
            
            r_grd_mask = r_grd[Br_inv]
            z_grd_mask = z_grd[Br_inv]
            
            Br_mask = Br[Br_inv]
            Bv_mask = Bv[Br_inv]
            Bz_mask = Bz[Br_inv]
            
            Br = griddata((r_grd_mask, z_grd_mask), Br_mask, (r_grd, z_grd))    
            Bv = griddata((r_grd_mask, z_grd_mask), Bv_mask, (r_grd, z_grd))    
            Bz = griddata((r_grd_mask, z_grd_mask), Bz_mask, (r_grd, z_grd))  
            
            self.B_cyld[0::, v, 0::] = np.stack((Br, Bv, Bz), axis=2)
        
        if save:
            name = os.path.join(loc, 'Bfield_HSX.h5')
            
            data = {'B field': self.B_cyld,
                    'r dom' : self.r_dom,
                    'p dom' : self.v_dom,
                    'z dom' : self.z_dom}
            
            save_data(name, data)
        
    def plot_fluxSurf(self, loc, v):
        """ Plot interpolated s surfaces from fluxSurf() at a given toroidal cross section
        
        Parameters
        ------------------------------------------------------------------------
        loc : str
            Directory in which png will be stored
            
        v : float
            Toroidal angle of cross section
        """
        try:
            self.s_surf
        except:
            self.fluxSurf()
        
        fnt=14
        v_ind = find_nearest(self.v_dom, v)
        
        fig, axs = plt.subplots(1)
        s_map = axs.pcolormesh(self.r_dom, self.z_dom, self.s_surf[0::, v_ind, 0::].T, vmin=0, vmax=1, cmap=plt.get_cmap('tab20c'))
        
        cbar = fig.colorbar(s_map, ax=axs)
        cbar.ax.set_ylabel(r'$\sqrt{s}$', fontsize=fnt, rotation=0)
        
        axs.set_xlabel('r (m)', fontsize=fnt)
        axs.set_ylabel('z (m)', fontsize=fnt)
        
        axs.set_title(r'$\theta$ = %.1f$^{\circ}$' % (self.v_dom[v_ind]*57.29577951), fontsize=fnt+2)
        
        plt.savefig( os.path.join( loc, 'S_Map_V_{}.png'.format( int( self.v_dom[v_ind]*57.29577951 ) ) ) )
        plt.close()
        
    def plot_fluxCont(self, dirc, tor=0):
        """ Plot a selection of flux surfaces in real space across a toroidal plane
        
        Parameters
        ------------------------------------------------------------------------
        dirc : str
            location of png image produced
            
        tor : float
            angle in radians of toroidal cross section to plot
        """
        try:
            self.R_coord
            self.Z_coord
        except:
            self.polarCoord()
            
        tor_ind = find_nearest(self.v_dom, tor)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        fnt = 16
        
        clrs = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
        
        f_dom = np.linspace(1, self.ns-1, 7, dtype=int) 
        for ind, f in enumerate(f_dom):
            line = Line2D(self.R_coord[f, tor_ind, 0::], self.Z_coord[f, tor_ind, 0::], color=clrs[ind], label='%.2f' % self.s_dom[f])
            ax.add_line(line)
        
        ax.set_xlim(np.min(self.R_coord), np.max(self.R_coord))
        ax.set_ylim(np.min(self.Z_coord), np.max(self.Z_coord))
        
        ax.set_xlabel('R', fontsize=fnt)
        ax.set_ylabel('Z', fontsize=fnt)
        ax.set_title(r'$\phi$ = %.1f$^{\circ}$' % (tor*(180/np.pi)), fontsize=fnt+2)
        
        plt.legend()
        plt.savefig( os.path.join( dirc, 'flux_tor_{}.png'.format( int( tor * (180/np.pi) ) ) ) )
        plt.close()
        
    def paraPlot_B(self, loc, flx=126, mod=True, R=False, P=False, Z=False):
        """ Make Parametric plot of the magnetic field along a flux surface.
        
        Parameters
        ------------------------------------------------------------------------
        loc : str
            Directory in which png will be saved
            
        flx : int
            The index of the flux surface to plot
        """
        if mod:
            try:
                self.B_mod
            except:
                self.polarCoord_Bmod()
    
            B_min = np.min(self.B_mod[np.nonzero(self.B_mod)])
            B_max = np.max(self.B_mod)
            B_show = self.B_mod[flx, 0::, 0::].T
            
            title='Bmod_para_{}.png'.format(flx+1)
            label=r'$|\mathbf{B}|$'
        else:
            try:
                self.B_field()
            except:
                self.polarCoord_Bfield()
            if R:
                B_min = self.B_field[flx,0::,0::,0].min()
                B_max = self.B_field[flx,0::,0::,0].max()
                B_show = self.B_field[flx, 0::, 0::, 0].T
                
                title='Br_para_{}.png'.format(flx+1)
                label=r'$B_{r}$'
            
            elif P:
                B_min = self.B_field[flx,0::,0::,1].min()
                B_max = self.B_field[flx,0::,0::,1].max()
                B_show = self.B_field[flx,0::,0::,1].T
                
                title='Bp_para_{}.png'.format(flx+1)
                label=r'$B_{\phi}$'
                
            elif Z:
                B_min = self.B_field[flx,0::,0::,2].min()
                B_max = self.B_field[flx,0::,0::,2].max()
                B_show = self.B_field[flx,0::,0::,2].T
                
                title='Bz_para_{}.png'.format(flx+1)
                label=r'$B_{z}$'

        fnt = 14
        
        fig, axs = plt.subplots(1)
        s_map = axs.pcolormesh(self.v_dom, self.u_dom, B_show, vmin=B_min, vmax=B_max, cmap=plt.get_cmap('jet'))
        
        cbar = fig.colorbar(s_map, ax=axs)
        cbar.ax.set_ylabel(label, fontsize=fnt, rotation=0)
        
        axs.set_xlabel(r'Tor', fontsize=fnt)
        axs.set_ylabel(r'Pol', fontsize=fnt)
        
        axs.set_title(r'Flux Surface: {0} of 127'.format(flx+1), fontsize=fnt+2)
        
        plt.savefig( os.path.join( loc, title ) )
        plt.close()
        
    def plot_3D_B(self, flx=126, mod=True, R=False, P=False, Z=False):
        """ Plot a 3D flux surface showing mod B.
        
        Parameters
        ------------------------------------------------------------------------
        flx : int
            Indexed flux surface to plot
        """
        if mod:
            try:
                self.B_mod
            except:
                self.polarCoord_Bmod()    
            B_min = np.min(self.B_mod[np.nonzero(self.B_mod)])
            B_max = np.max(self.B_mod)
            B_scale = self.B_mod[flx].T
            B_rescale = (B_scale - B_min) / (B_max - B_min)
            title=r'$\|B\|$'
            
        else:
            try:
                self.B_field
            except:
                self.polarCoord_Bfield()
            if R:
                B_min = self.B_field[flx,0::,0::,0].min()
                B_max = self.B_field[flx,0::,0::,0].max()
                B_scale = self.B_field[flx,0::,0::,0].T
                B_rescale = (B_scale - B_min) / (B_max - B_min)
                title=r'$B_{r}$'
            elif P:
                B_min = self.B_field[flx,0::,0::,1].min()
                B_max = self.B_field[flx,0::,0::,1].max()
                B_scale = self.B_field[flx,0::,0::,1].T
                B_rescale = (B_scale - B_min) / (B_max - B_min)
                title=r'$B_{\phi}$'
            elif Z:
                B_min = self.B_field[flx,0::,0::,2].min()
                B_max = self.B_field[flx,0::,0::,2].max()
                B_scale = self.B_field[flx,0::,0::,2].T
                B_rescale = (B_scale - B_min) / (B_max - B_min)
                title=r'$B_{z}$'
        
        fnt=14
        fig = plt.figure()

        v_2dom, u_2dom = np.meshgrid(self.v_dom, self.u_dom)
        
        X = self.R_coord[flx].T * np.cos(v_2dom)
        Y = self.R_coord[flx].T * np.sin(v_2dom)
        
        fig.patch.set_facecolor('white')
        ax = fig.gca(projection='3d')
        ax.plot_surface(X, Y, self.Z_coord[flx].T, facecolors = cm.jet(B_rescale), rstride=1, cstride=1, antialiased=False)
        ax.auto_scale_xyz([X.min(), X.max()], [X.min(), X.max()], [X.min(), X.max()])
        
        m = cm.ScalarMappable(cmap=cm.jet)
        m.set_array(B_scale)
        cbar = plt.colorbar(m)
        cbar.ax.set_ylabel(title, fontsize=fnt+4)
        
        ax.set_xlabel('X (m)', fontsize=fnt)
        ax.set_ylabel('Y (m)', fontsize=fnt)
        ax.set_zlabel('Z (m)', fontsize=fnt)
        
        plt.show()
        
    def plot_B(self, loc, v, mod=True, R=False, P=False, Z=False):
        """ Plots poloidal cross sections of the magnetic field.
        
        Parameters
        ------------------------------------------------------------------------
        loc : str
            Directory in which png will be saved
            
        v : int
            Toroidal angle of cross section
        """
        fnt=14
        v_ind = find_nearest(self.v_dom, v)

        if mod:
            try:
                self.B_mod
            except:
                self.polarCoord_Bmod()    
            B_min = np.min(self.B_mod[np.nonzero(self.B_mod)])
            B_max = np.max(self.B_mod)
            B_show = self.B_mod[0::, v_ind, 0::]
        else:
            try:
                self.B_field
            except:
                self.polarCoord_Bfield()
            if R:
                B_min = np.min(self.B_field[0::,0::,0::,0])
                B_max = np.max(self.B_field[0::,0::,0::,0])
                B_show = self.B_field[0::, v_ind, 0::, 0]
            elif P:
                B_min = np.min(self.B_field[np.nonzero(self.B_field[0::,0::,0::,1])])
                B_max = np.max(self.B_field[0::,0::,0::,1])
                B_show = self.B_field[0::, v_ind, 0::, 1]
            elif Z:
                B_min = np.min(self.B_field[0::,0::,0::,2])
                B_max = np.max(self.B_field[0::,0::,0::,2])
                B_show = self.B_field[0::, v_ind, 0::, 2]
        
        fig, axs = plt.subplots(1)
        s_map = axs.pcolormesh(self.R_coord[0::, v_ind, 0::], self.Z_coord[0::, v_ind, 0::], B_show, vmin=B_min, vmax=B_max, cmap=plt.get_cmap('jet'))
        
        cbar = fig.colorbar(s_map, ax=axs)
        cbar.ax.set_ylabel(r'$|\mathbf{B}|$', fontsize=fnt, rotation=0)
        
        axs.set_xlabel('r (m)', fontsize=fnt)
        axs.set_ylabel('z (m)', fontsize=fnt)
        
        axs.set_title(r'$\theta$ = %.1f$^{\circ}$' % (self.v_dom[v_ind]*57.29577951), fontsize=fnt+2)
        
        plt.savefig( os.path.join( loc, 'B_field_V_{}.png'.format( int( self.v_dom[v_ind]*57.29577951 ) ) ) )
        plt.close()
        
    def plot_poloidal(self, loc, v):
        fnt=14
        v_ind = find_nearest(self.v_dom, v)
        
        R_ma = np.array([self.R_coord[0,v_ind,0], self.Z_coord[0,v_ind,0]])
        r_ma = np.linalg.norm(R_ma)
        
        R_lcfs = np.array([self.R_coord[-1,v_ind,0], self.Z_coord[-1,v_ind,0]])
        r_lcfs = np.linalg.norm(R_lcfs)
        
        R_pol = R_lcfs - R_ma
        r_pol = np.linalg.norm(R_pol)
        
        beta = (np.pi - np.arccos( ((r_ma**2 + r_pol**2) - r_lcfs**2) / (2 * r_ma * r_pol) ) ) * (180 / np.pi)
        
        plt.plot(self.R_coord[-1, v_ind, 0::], self.Z_coord[-1, v_ind, 0::])
        plt.plot([0,R_ma[0]], [0,R_ma[1]])
        plt.plot([R_ma[0],R_lcfs[0]], [R_ma[1],R_lcfs[1]], ls='--')
        plt.plot(self.R_coord[0::,v_ind,0], self.Z_coord[0::,v_ind,0])
        
        plt.xlim([self.R_coord.min(), self.R_coord.max()])   
        plt.ylim([self.Z_coord.min(), self.Z_coord.max()])   
        plt.title(r'$\beta = %s $' % int(beta))
        
        plt.savefig( os.path.join(loc, 'pol_{}.png'.format( int( self.v_dom[v_ind]*57.29577951 ) ) ) )
        plt.close()
        
    def fit_func(self, loc, v):
        fnt=14
        v_ind = find_nearest(self.v_dom, v)
        
        # Fit S (beg)
        r = self.R_coord[0::,v_ind,0] - self.R_coord[0,v_ind,0]
        z = self.Z_coord[0::,v_ind,0] - self.Z_coord[0,v_ind,0]
        R = np.stack((r,z), axis=1)
        R_norm = np.linalg.norm(R, axis=1)
        s = (R_norm / R_norm[-1])**2
        # Fit S (end)
        
        R_vals_1 = self.R_coord[30,v_ind,0::]#[self.R_coord[30,v_ind,0], self.R_coord[30,v_ind,10], self.R_coord[30,v_ind,20], self.R_coord[30,v_ind,30], self.R_coord[30,v_ind,40], self.R_coord[30,v_ind,50], self.R_coord[30,v_ind,60], self.R_coord[30,v_ind,70], self.R_coord[30,v_ind,80], self.R_coord[30,v_ind,90], self.R_coord[30,v_ind,100]]
        Z_vals_1 = self.Z_coord[30,v_ind,0::]#[self.Z_coord[30,v_ind,0], self.Z_coord[30,v_ind,10], self.Z_coord[30,v_ind,20], self.Z_coord[30,v_ind,30], self.Z_coord[30,v_ind,40], self.Z_coord[30,v_ind,50], self.Z_coord[30,v_ind,60], self.Z_coord[30,v_ind,70], self.Z_coord[30,v_ind,80], self.Z_coord[30,v_ind,90], self.Z_coord[30,v_ind,100]]

        R_vals_2 = self.R_coord[90,v_ind,0::]#[self.R_coord[90,v_ind,0], self.R_coord[90,v_ind,10], self.R_coord[90,v_ind,20], self.R_coord[90,v_ind,30], self.R_coord[90,v_ind,40], self.R_coord[90,v_ind,50], self.R_coord[90,v_ind,60], self.R_coord[90,v_ind,70], self.R_coord[90,v_ind,80], self.R_coord[90,v_ind,90], self.R_coord[90,v_ind,100]]
        Z_vals_2 = self.Z_coord[90,v_ind,0::]#[self.Z_coord[90,v_ind,0], self.Z_coord[90,v_ind,10], self.Z_coord[90,v_ind,20], self.Z_coord[90,v_ind,30], self.Z_coord[90,v_ind,40], self.Z_coord[90,v_ind,50], self.Z_coord[90,v_ind,60], self.Z_coord[90,v_ind,70], self.Z_coord[90,v_ind,80], self.Z_coord[90,v_ind,90], self.Z_coord[90,v_ind,100]]
        
        R_samp = {'Samp 1' : np.stack((R_vals_1, Z_vals_1), axis=1),
                  'Samp 2' : np.stack((R_vals_2, Z_vals_2), axis=1)}
        flux = {'Samp 1' : 30,
                'Samp 2' : 90}
        R_ma = np.array([self.R_coord[0,v_ind,0], self.Z_coord[0,v_ind,0]])
        r_ma = np.linalg.norm(R_ma)
        
        R_lcfs = np.array([self.R_coord[-1,v_ind,0], self.Z_coord[-1,v_ind,0]])
        r_lcfs = np.linalg.norm(R_lcfs)
        
        R_pol = R_lcfs - R_ma
        r_pol = np.linalg.norm(R_pol)
        
        beta = (np.pi - np.arccos( ((r_ma**2 + r_pol**2) - r_lcfs**2) / (2 * r_ma * r_pol) ) ) * (180 / np.pi)
        
        if abs(beta) > 1e-5:
            R_orth = R_ma - R_pol * (np.dot(R_ma, R_pol) / r_pol**2)
            r_orth = R_orth / np.linalg.norm(R_orth)
            r_strt = R_pol / r_pol
            #print('GS: %.5f, %.5f' % (np.dot(r_orth, r_strt), beta))
        else:
            r_orth = np.array([0,1])
            r_strt = np.array([1,0])
            #print('Arb: %.5f, %.5f' % (np.dot(r_orth, r_strt), beta))
        
        R_apprx = {}
        for key in R_samp:
            R_sample = R_samp[key]
            R_guess = np.empty(R_sample.shape)
            
            flx = flux[key]
            for ind, R in enumerate(R_sample):
                '''
                R_vec = R - R_ma
                R_vec_norm = np.linalg.norm(R_vec)
                r_vec = R_vec / R_vec_norm
                u_guess = np.arccos(np.dot(r_vec, r_strt))
                
                if key == 'Samp 2':
                    plt.plot([R_ma[0], R_ma[0]+.1*r_vec[0]], [R_ma[1], R_ma[1]+.1*r_vec[1]], c='tab:green')
                '''
                '''
                if np.dot(r_vec, r_orth) < 0:
                    u_guess = 2*np.pi - u_guess
                '''
                a = r_pol
                b = np.linalg.norm(R - R_ma)
                c = np.linalg.norm(R - R_lcfs)
                
                u_guess = np.arccos( ((a**2 + b**2) - c**2) / (2*a*b) )
                
                u_apprx = find_nearest(self.u_dom, u_guess)
                '''
                R_max = np.hypot(self.R_coord[-1,v_ind,u_apprx], self.Z_coord[-1,v_ind,u_apprx])
                s_guess = (R_vec_norm / R_max)**2
                s_apprx = find_nearest(self.s_dom, s_guess)
                '''
                R_guess[ind] = np.array([self.R_coord[flx, v_ind, u_apprx], self.Z_coord[flx, v_ind, u_apprx]])
                    
            R_apprx[key] = R_guess
        
        #plt.plot(self.s_dom, self.s_dom, label='exact')
        #plt.plot(self.s_dom, s, ls='--', label='approx.')
        
        plt.plot(self.R_coord[-1,v_ind,0::], self.Z_coord[-1,v_ind,0::])
        
        #plt.plot(R_vals_1, Z_vals_1, color='tab:red')
        #plt.plot(R_vals_2, Z_vals_2, color='tab:green')
        
        plt.scatter(R_vals_1, Z_vals_1, color='tab:red', s=10, facecolor='white')
        plt.scatter(R_vals_2, Z_vals_2, color='tab:green', s=10, facecolor='white')
        
        #plt.scatter(R_apprx['Samp 1'][0::,0], R_apprx['Samp 1'][0::,1], color='tab:red', s=9)
        #plt.scatter(R_apprx['Samp 2'][0::,0], R_apprx['Samp 2'][0::,1], color='tab:green', s=9)
        
        plt.plot([R_ma[0], R_ma[0]+.1*r_strt[0]], [R_ma[1], R_ma[1]+.1*r_strt[1]])
        #plt.plot([R_ma[0], R_ma[0]+.1*r_orth[0]], [R_ma[1], R_ma[1]+.1*r_orth[1]])
        plt.plot(self.R_coord[0::,v_ind,0], self.Z_coord[0::,v_ind,0], ls='--')
        
        plt.xlim(self.R_coord.min(), self.R_coord.max())
        z_lim = (self.R_coord.max() - self.R_coord.min()) * .5
        plt.ylim(-z_lim, z_lim)
        #plt.legend()
        plt.savefig( os.path.join(loc, 'pol_{}.png'.format( int( self.v_dom[v_ind]*57.29577951 ) ) ) )
        plt.close()
        
    def arc_leng(self, loc, v, flx):
        v_ind = find_nearest(self.v_dom, v)
        
        R_ma = np.array([self.R_coord[0,v_ind,0], self.Z_coord[0,v_ind,0]])
        
        R_vals = self.R_coord[flx,v_ind,0::]
        Z_vals = self.Z_coord[flx,v_ind,0::]
        R_vecs = np.stack((R_vals, Z_vals), axis=1)
        
        dRu = self.dRu_coord[flx,v_ind,0::]
        
        r_norm = np.empty(self.u_num-1)
        d_theta = np.empty(self.u_num-1)
        for ind, R in enumerate(R_vecs[1:self.u_num]):
            r_1 = R_vecs[ind] - R_ma
            r_2 = R - R_ma
            
            r1_norm = np.linalg.norm(r_1)
            r2_norm = np.linalg.norm(r_2)
            
            d_theta[ind] = np.arccos( np.dot(r_1, r_2) / (r1_norm*r2_norm) )
            r_norm[ind] = np.mean([dRu[ind+1], dRu[ind]])
            
        arc = r_norm / d_theta
        plt.plot(range(len(arc)), arc)
        
        plt.savefig( os.path.join(loc, 'pol_{}.png'.format( int( self.v_dom[v_ind]*57.29577951 ) ) ) )
        plt.close()
        
dirc_name = os.getcwd()
file_name = os.path.join(dirc_name, 'wout_HSX_test_opt0.nc')
fig_dirc = os.path.join(dirc_name, 'Figures')

vmec_data = readVMEC(file_name)
vmec_data.polarCoord()

v_dom = np.linspace(0, .5*np.pi, 15)
for v in v_dom:
    vmec_data.arc_leng(fig_dirc, v, 126)
