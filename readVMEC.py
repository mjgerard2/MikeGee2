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

def find_nearest(array, value):
    return (np.abs(array - value)).argmin()

class readVMEC:
    
    def __init__(self, file_name):
        rootgrp = Dataset(file_name, 'r')
            
        self.nfp = rootgrp['/nfp'][0]
        self.mpol = rootgrp['/mpol'][0]
        self.ntor = rootgrp['/ntor'][0]
        
        self.xm = rootgrp['/xm'][:]
        self.xn = rootgrp['xn'][:]
        self.md = len(self.xm)
        
        self.xm_nyq = rootgrp['/xm_nyq'][:]
        self.xn_nyq = rootgrp['/xn_nyq'][:]
        self.md_nyq = len(self.xm_nyq)
        
        self.rmnc = rootgrp['/rmnc'][:,:]
        self.zmns = rootgrp['/zmns'][:,:]
        self.lmns = rootgrp['/lmns'][:,:]
        
        self.bumnc = rootgrp['/bsupumnc'][:,:]
        self.bvmnc = rootgrp['/bsupvmnc'][:,:]
                    
        self.ns = rootgrp['/ns'][0]
        
        rootgrp.close()
        
        self.u_num = self.mpol * 4
        self.v_num = self.ntor 
        
        self.s_dom = np.linspace(0, 1, self.ns)
        self.u_dom = np.linspace(0, 2*np.pi, self.u_num)
        self.v_dom = np.linspace(0, .5*np.pi, self.v_num)
        
    def polarCoord(self):
        """
        Produces the cylindrical coordinates of magnetic flux surfaces, 
        indexed by VMEC flux coordinates [s,v,u].
        """
        pol, tor = np.meshgrid(self.u_dom, self.v_dom)
        
        pol_xm = np.dot(self.xm.reshape(self.md, 1), pol.reshape(1, self.v_num * self.u_num))
        tor_xn = np.dot(self.xn.reshape(self.md, 1), tor.reshape(1, self.v_num * self.u_num))
        
        cos_pol = np.cos(pol_xm)
        cos_tor = np.cos(tor_xn)

        sin_pol = np.sin(pol_xm)
        sin_tor = np.sin(tor_xn)
        
        cos_mu_nv = cos_pol*cos_tor + sin_pol*sin_tor
        sin_mu_nv = sin_pol*cos_tor - sin_tor*cos_pol
        
        self.R_coord = np.dot(self.rmnc, cos_mu_nv).reshape(self.ns, self.v_num, self.u_num)
        self.Z_coord = np.dot(self.zmns, sin_mu_nv).reshape(self.ns, self.v_num, self.u_num)
        self.L_coord = np.dot(self.lmns, sin_mu_nv).reshape(self.ns, self.v_num, self.u_num)
        
        self.dRu_coord = np.dot(-self.rmnc*self.xm, sin_mu_nv).reshape(self.ns, self.v_num, self.u_num)
        self.dRv_coord = np.dot(self.rmnc*self.xn, sin_mu_nv).reshape(self.ns, self.v_num, self.u_num)
        
        self.dZu_coord = np.dot(self.zmns*self.xm, cos_mu_nv).reshape(self.ns, self.v_num, self.u_num)
        self.dZv_coord = np.dot(-self.zmns*self.xn, cos_mu_nv).reshape(self.ns, self.v_num, self.u_num)
        
    def fluxSurf(self, r_res=0.001, z_res=0.001):
        """
        Produce a 3D array of flux surface s values, indexed in cylindrical coordinates [r, theta, z]
        
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
                    
                    s_temp[r_ind, z_ind] = np.sqrt(s)
                    
            s_temp = np.ma.masked_invalid(s_temp)
            s_temp_inv = ~s_temp.mask
            
            r_grid_mask = r_grid[s_temp_inv]
            z_grid_mask = z_grid[s_temp_inv]
            s_temp_mask = s_temp[s_temp_inv]
            
            s_temp_intp = griddata((z_grid_mask, r_grid_mask), s_temp_mask, (z_grid, r_grid))
            self.s_surf[0::, v_ind, 0::] = s_temp_intp
            
    def fluxSurf_Polar(self, r_res=0.001, p_res=1., z_res=0.001):
        """
        Produce a 3D array of flux surface s values, indexed in polar coordinates [r,p,z]
        
        Parameters
        ------------------------------------------------------------------------
        r_res : float
            coordinate resolution in meters
            
        p_res : float
            coordinate resolution in degrees
            
        z_res : float
            coordinate resolution in meters
        """
        try:
            self.R_coord
        except:
            self.polarCoord()
            
        r_min = np.min(self.R_coord)
        r_max = np.max(self.R_coord)
        r_num = int((r_max - r_min) / r_res)
        self.r_dom = np.linspace(r_min, r_max, r_num)
        
        p_min = 0
        p_max = .5*np.pi
        p_num = int(((p_max - p_min)*180) / (p_res*np.pi))
        self.p_dom = np.linspace(p_min, p_max, p_num)
        
        z_min = np.min(self.Z_coord)
        z_max = np.max(self.Z_coord)
        z_num = int((z_max - z_min) / z_res)
        self.z_dom = np.linspace(z_min, z_max, z_num)
        
        p_grid, r_grid, z_grid = np.meshgrid(self.p_dom, self.r_dom, self.z_dom)
        
        s_surf = np.full((r_num, p_num, z_num), np.nan)
        for s_ind, s in enumerate(self.s_dom):
            for v_ind, v in enumerate(self.v_dom):
                p_ind = find_nearest(self.p_dom, v)
                for u_ind, u in enumerate(self.u_dom):
                    r = self.R_coord[s_ind, v_ind, u_ind]
                    z = self.Z_coord[s_ind, v_ind, u_ind]
                    
                    r_ind = find_nearest(self.r_dom, r)
                    z_ind = find_nearest(self.z_dom, z)
                    
                    s_surf[r_ind, p_ind, z_ind] = np.sqrt(s)
        
        s_surf = np.ma.masked_invalid(s_surf)
        s_surf_inv = ~s_surf.mask
        
        r_grid_mask = r_grid[s_surf_inv]
        p_grid_mask = p_grid[s_surf_inv]
        z_grid_mask = z_grid[s_surf_inv]
        s_surf_mask = s_surf[s_surf_inv]
                
        s_surf_polar = griddata((r_grid_mask, p_grid_mask, z_grid_mask), s_surf_mask, (r_grid, p_grid, z_grid))        
        
        for p_ind, p in enumerate(self.p_dom):
            self.R_coord[:,p,:]
            
    def fluxSurf_Cartesian(self, res=0.005):
        """
        Produce a 3D array of flux surface s values, indexed in cartesian coordinates [x,y,z]
        
        Parameters
        ------------------------------------------------------------------------
        res : float
            coordinate resolution in meters
        """
        try:
            self.R_coord
        except:
            self.polarCoord()
            
        x_min = 0
        x_max = np.max(self.R_coord)
        x_num = int((x_max - x_min) / res)
        self.x_dom = np.linspace(x_min, x_max, x_num)
        self.y_dom = self.x_dom
        
        z_min = np.min(self.Z_coord)
        z_max = np.max(self.Z_coord)
        z_num = int((z_max - z_min) / res)
        self.z_dom = np.linspace(z_min, z_max, z_num)
        
        y_grid, x_grid, z_grid = np.meshgrid(self.y_dom, self.x_dom, self.z_dom)
        
        s_surf = np.full((x_num, x_num, z_num), np.nan)
        for s_ind, s in enumerate(self.s_dom):
            for v_ind, v in enumerate(self.v_dom):
                for u_ind, u in enumerate(self.u_dom):
                    r = self.R_coord[s_ind, v_ind, u_ind]
                    z = self.Z_coord[s_ind, v_ind, u_ind]
                    
                    x = r * np.cos(v)
                    y = r * np.sin(v)
                    
                    x_ind = find_nearest(self.x_dom, x)
                    y_ind = find_nearest(self.y_dom, y)
                    z_ind = find_nearest(self.z_dom, z)
                    
                    s_surf[x_ind, y_ind, z_ind] = np.sqrt(s)
                    
        s_surf = np.ma.masked_invalid(s_surf)
        s_surf_inv = ~s_surf.mask
        
        x_grid_mask = x_grid[s_surf_inv]
        y_grid_mask = y_grid[s_surf_inv]
        z_grid_mask = z_grid[s_surf_inv]
        s_surf_mask = s_surf[s_surf_inv]

        self.s_surf_cart = griddata((x_grid_mask, y_grid_mask, z_grid_mask), s_surf_mask, (x_grid, y_grid, z_grid))
        
    def fluxCoord_Bfield(self):
        """
        Produces the magnetic field in real space in cylindrical coordinates, 
        indexed by VMEC flux coordinates [s,v,u].        
        """
        try:
            self.dRu_coord
            self.dZu_coord
        except:
            self.polarCoord()
            
        pol, tor = np.meshgrid(self.u_dom, self.v_dom)
        
        pol_xm = np.dot(self.xm_nyq.reshape(self.md_nyq, 1), pol.reshape(1, self.v_num * self.u_num))
        tor_xn = np.dot(self.xn_nyq.reshape(self.md_nyq, 1), tor.reshape(1, self.v_num * self.u_num))
        
        cos_mu_nv = np.cos(pol_xm)*np.cos(tor_xn) + np.sin(pol_xm)*np.sin(tor_xn)
        
        Bu_coord = np.dot(self.bumnc, cos_mu_nv).reshape(self.ns, self.v_num, self.u_num)
        Bv_coord = np.dot(self.bvmnc, cos_mu_nv).reshape(self.ns, self.v_num, self.u_num)
        
        Br_coord = Bu_coord * self.dRu_coord
        Bz_coord = Bu_coord * self.dZu_coord
        self.B_flux = np.stack((Br_coord, Bv_coord, Bz_coord), axis=3)
        
    def cylndCoord_Bfield(self, r_res=0.005, z_res=0.005):
        """
        Interpolates the magnetic field in real space in cylindrical coordinates, 
        indexed by the cylindrical domain [r,theta,z].        
        
        Parameters
        ------------------------------------------------------------------------
        r_res : float
            r coordinate resolution in meters
            
        z_res : array
            z coordinate resolution in meters
        """
        try:
            self.B_flux
        except:
            self.fluxCoord_Bfield()
        
        r_min = np.min(self.R_coord)
        r_max = np.max(self.R_coord)
        r_num = int((r_max - r_min) / r_res)
        self.r_dom = np.linspace(r_min, r_max, r_num)
        
        z_min = np.min(self.Z_coord)
        z_max = np.max(self.Z_coord)
        z_num = int((z_max - z_min) / z_res)
        self.z_dom = np.linspace(z_min, z_max, z_num)
        
        self.B_cyld = np.empty((r_num, self.v_num, z_num, 3))
        
        for v in range(self.v_num):
            r = self.R_coord[0,v,0]
            z = self.Z_coord[0,v,0]
            
            r_ind = find_nearest(self.r_dom, r)
            z_ind = find_nearest(self.z_dom, z)
            
            self.B_cyld[r_ind, v, z_ind] = self.B_flux[0,v,0]
        
        for s in range(1, self.ns):            
            for v in range(self.v_num):
                for u in range(self.u_num):
                    r = self.R_coord[s,v,u]
                    z = self.Z_coord[s,v,u]
                    
                    r_ind = find_nearest(self.r_dom, r)
                    z_ind = find_nearest(self.z_dom, z)
                    
                    self.B_cyld[r_ind, v, z_ind] = self.B_flux[s,v,u]
        
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
    
    def plot_fluxSurf(self, dirc, tor=0):
        """
        Plots a selection of flux surfaces in real space across a toroidal plane
        
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
        
    def plot_S_Surf(self, loc, v):
        """
        Produces png of the interpolated s surfaces from fluxSurf() at a given toroidal cross section
        
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
        
    def save_Sval_Polar(self, loc):
        """
        Save an hdf5 file of the 3D array of flux surfaces indexed in polar coordinates [r,p,z]
        
        Parameters
        ------------------------------------------------------------------------
        loc : str
            Directory in which hdf5 file is saved
        """
        try:
            self.s_surf_polar
        except:
            self.fluxSurf_Polar()
        
        name = os.path.join(loc, 's_surf_polar.h5')
        h5 = hf.File(name, 'w')
        
        h5.create_dataset('s surf', data=self.s_surf_polar)
        
        h5.create_dataset('r dom', data=self.r_dom)
        h5.create_dataset('p dom', data=self.p_dom)
        h5.create_dataset('z dom', data=self.z_dom)
        
        h5.close()
        
    def save_Sval_Cartesian(self, loc):
        """
        Save an hdf5 file of the 3D array of flux surfaces indexed in cartesian coordinates [x,y,z]
        
        Parameters
        ------------------------------------------------------------------------
        loc : str
            Directory in which hdf5 file is saved
        """
        try:
            self.s_surf_cart
        except:
            self.fluxSurf_Cartesian()
        
        name = os.path.join(loc, 's_surf_cart.h5')
        h5 = hf.File(name, 'w')
        
        h5.create_dataset('s surf', data=self.s_surf_cart)
        
        h5.create_dataset('x dom', data=self.x_dom)
        h5.create_dataset('y dom', data=self.y_dom)
        h5.create_dataset('z dom', data=self.z_dom)
        
        h5.close()
        
    def save_Bfield_Cyldrical(self, loc):
        """
        Save an hdf5 file of the B field, indexed in cylindrical coordinates [r,p,z]
        
        Parameters
        ------------------------------------------------------------------------
        loc : str
            Directory in which hdf5 file is saved
        """
        try:
            self.B_cyld
        except:
            self.cylndCoord_Bfield()
            
        name = os.path.join(loc, 'B_field.h5')
        h5 = hf.File(name, 'w')
        
        h5.create_dataset('B field', data=self.B_cyld)
        
        h5.create_dataset('r dom', data=self.r_dom)
        h5.create_dataset('z dom', data=self.z_dom)
        h5.create_dataset('p dom', data=self.v_dom)
        
        h5.close()

dirc_name = os.getcwd()
file_name = os.path.join(dirc_name, 'wout_QHS_Rstart_1_513_32polmodes_18x24_axis_v2.nc')

vmec_data = readVMEC(file_name)
vmec_data.save_Bfield_Cyldrical(dirc_name)