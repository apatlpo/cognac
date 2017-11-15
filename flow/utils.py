# encoding: utf-8
'''
Useful functions to analyze roms simulations
'''


import sys
from glob import glob

from netCDF4 import Dataset

import numpy as np
import warnings




def minmax(x,xname):
    print(xname+' min= %e, max=%e'%(x.min(),x.max()))
    
    

#--------------------------------------------------------------------------------------------

class grid(object):
    
    def __init__(self,datadir='/home/datawork-lops-osi/jgula/NESED/'):
        
        self._datadir = datadir
        
        self._load_hgrid()
        self._load_vgrid()
        
        
    def _load_hgrid(self):
        ''' load horizontal grid variables
        '''
        # search for files
        grd_file = glob(self._datadir+'*grd.nc')
        if len(grd_file)==0:
            print('No grid file found in'+self._datadir)
            sys.exit()
        else:
            print('%i grid file found, uses: %s'%(len(grd_file),grd_file[0]))
            grd_file = grd_file[0]
        # open and load variables
        nc = Dataset(grd_file,'r')
        self._nch = nc
        self.lon_rho = nc['lon_rho'][:]
        self.lat_rho = nc['lat_rho'][:]
        self.h = nc['h'][:]
        #
        self.Lp = self.lon_rho.shape[1]
        self.Mp = self.lon_rho.shape[0]
        #
        print('dim lon_rho: %i  %i' %self.lon_rho.shape)
        minmax(self.lon_rho,'lon_rho')
        minmax(self.lat_rho,'lat_rho')
        minmax(self.h,'h')
        #
        self.hextent = [self.lon_rho.min(), self.lon_rho.max(), \
                        self.lat_rho.min(), self.lat_rho.max()]


    def _load_vgrid(self):
        ''' load vertical grid variables
        '''
        # search for files with 
        files = glob(self._datadir+'*.nc')
        if len(files)==0:
            print('No nc file found in'+self._datadir)
            sys.exit()
        for file in files:
            nc = Dataset(file,'r')
            if 'sc_w' in nc.ncattrs():
                print('vertical grid parameters found in %s'%(file))
                break
            else:
                nc.close()
        self.hc = nc.getncattr('hc')
        self.sc_w = nc.getncattr('sc_w')
        self.Cs_w = nc.getncattr('Cs_w')
        self.sc_r = nc.getncattr('sc_r')
        self.Cs_r = nc.getncattr('Cs_r')
        #self.z_r = z_r()

        
    def get_z(self,zeta,h,sc,cs):
        ''' compute vertical coordinates
            zeta should have the size of the final output
            vertical coordinate should be first
        '''
        z0 = (self.hc * sc + h * cs) / (self.hc + h)
        z = zeta + (zeta + h) * z0
        return z

        
#--------------------------------------------------------------------------------------------
def get_soundc(t,s,z,lon,lat):
    ''' compute sound velocity
    '''
    import gsw
    p = gsw.p_from_z(z,lat.mean())
    SA = gsw.SA_from_SP(s,p, lon, lat)
    CT = gsw.CT_from_pt(SA,t)
    c = gsw.sound_speed(s,t,p)
    # inputs are: SA (absolute salinity) and CT (conservative temperature)
    return c
    
    
    