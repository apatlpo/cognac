import numpy as np
import pandas as pd
import xarray as xr
from netCDF4 import Dataset

from scipy.interpolate import interp1d

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import hvplot.pandas  # noqa
import hvplot

import gsw


#----------------------------------- core --------------------------------------

#
class waterp():
    ''' Data holder for a water column based on either:
     - in situ temperature, salinity, pressure profiles
     - the World Ocean Atlas (WOA) climatology, see:
        https://www.nodc.noaa.gov/OC5/woa18/

    Should base this on pandas !!

    Parameters
    ----------
    pressure: np.ndarray, optional
        pressure in dbar
    temperature: np.ndarray, optional
        in situ temperature in degC
    salinity: np.ndarray, optional
        salinity in PSU
    lon: float, optional
        longitude of the selected location
    lat: float, optional
        latitude of the selected location

    '''

    def __init__(self,
                 pressure=None,
                 temperature=None,
                 salinity=None,
                 lon=None,
                 lat=None,
                 name=None,
                 interp_kwargs=None,
                 **kwargs
                 ):

        self._pts, self._woa = False, False
        self._interp_kwargs = {'kind': 'quadratic',
                               'fill_value': None,
                               'bounds_error': False,
                               }
        if interp_kwargs:
            self._interp_kwargs.update(**interp_kwargs)

        args = [pressure, temperature, salinity, lon, lat]
        if all([a is not None for a in args]):
            self._load_from_pts(pressure, temperature, salinity,
                                lon, lat, name)
        elif all([lon ,lat]):
            self._load_from_woa(lon,lat,name)
        else:
            print('Inputs missing')

        # init isopycnal displacement and velocity
        self.eta = 0.
        self.detadt = 0.


    def _load_from_pts(self, pressure, temperature, salinity, lon, lat,
                       name):
        self._pts=True
        #
        if all([isinstance(i, float)
                for i in [pressure, temperature, salinity, lon, lat]]):
            #pressure = np.array([-1, 1e4])
            #temperature = np.array([temperature, temperature])
            #salinity = np.array([salinity, salinity])
            #lon = np.array([lon,lon])
            #lat = np.array([lat,lat])
            pressure = [-1, 1e4]
            temperature = [temperature, temperature]
            salinity = [salinity, salinity]
            lon = [lon,lon]
            lat = [lat,lat]
        #
        self.lon, self.lat = lon, lat
        #
        self.p = pressure
        self.z = gsw.z_from_p(self.p, self.lat)
        #
        self.temp, self.s = temperature, salinity
        #
        self._update_eos()
        #
        if name is None:
            self.name = 'Provided water profile at lon=%.0f, lat=%.0f'%(self.lon,self.lat)

    def _load_from_woa(self, lon, lat, name):
        self._woa=True
        #
        self._tfile = 'woa18_A5B7_t00_01.nc'
        nc = Dataset(self._tfile,'r')
        #
        glon = nc.variables['lon'][:]
        glat = nc.variables['lat'][:]
        ilon = np.argmin(np.abs(lon-glon))
        ilat = np.argmin(np.abs(lat-glat))
        self.lon = glon[ilon]
        self.lat = glat[ilat]
        #
        self.z = -nc.variables['depth'][:].data
        self.p = gsw.p_from_z(self.z,self.lat)
        #
        self.temp = nc.variables['t_an'][0,:,ilat,ilon]
        nc.close()
        #
        self._sfile = 'woa18_A5B7_s00_01.nc'
        nc = Dataset(self._sfile,'r')
        self.s = nc.variables['s_an'][0,:,ilat,ilon]
        nc.close()
        #
        self._update_eos()
        #
        if name is None:
            self.name = 'WOA water profile at lon=%.0f, lat=%.0f'%(self.lon,self.lat)

    def _update_eos(self):
        # derive absolute salinity and conservative temperature
        self.SA = gsw.SA_from_SP(self.s, self.p, self.lon, self.lat)
        self.CT = gsw.CT_from_t(self.SA, self.temp, self.p)
        # derive N2
        self.N2, self.p_mid = gsw.Nsquared(self.SA, self.CT, self.p, lat=self.lat)
        self.z_mid = gsw.z_from_p(self.p_mid, self.lat)

    def _get_df(self, z=None):
        if z is None:
            z = self.z
        return pd.DataFrame(data={'temperature': self.get_temp(z),
                                  'salinity': self.get_s(z),
                                  'rho': self.get_rho(z),
                                  },
                            index=z,
                            )

    def show_on_map(self, zdiff=None):
        if self._woa:
            ds = xr.open_dataset(self._tfile, decode_times=False).squeeze()
            #
            if zdiff is not None:
                toplt = ds['t_an'].sel(depth=zdiff[0]) \
                            - ds['t_an'].sel(depth=zdiff[1])
                title = 't(%dm) - t(%dm) [degC]'%(zdiff[0],zdiff[1])
                print('t(0m) - t(500m) global maximum = %.1f degC'%toplt.max())
            else:
                toplt = ds['t_an'].sel(depth=0)
                title = 'sea surface temperature [degC]'
            #
            crs=ccrs.PlateCarree()
            plt.figure(figsize=(10, 5))
            ax = plt.axes(projection=crs)
            hdl = ax.pcolormesh(ds.lon, ds.lat, toplt, transform=crs,
                                cmap=plt.get_cmap('CMRmap_r'))
            ax.plot(self.lon,self.lat, '*', markersize=10,
                    markerfacecolor='CadetBlue', markeredgecolor='w',
                    transform=crs)
            ax.coastlines(resolution='110m')
            ax.gridlines()
            plt.colorbar(hdl, ax=ax)
            ax.set_title(title)
            plt.show()
        else:
            print('No map to show')

    def _plot_matplotlib(self):
        plt.figure(figsize=(9,5))
        ax = plt.subplot(131)
        ax.plot(self.get_temp(self.z),self.z,'k')
        ax.set_ylabel('z [m]')
        ax.set_title('in situ temperature [degC]')
        plt.grid()
        #
        ax = plt.subplot(132)
        ax.plot(self.get_s(self.z),self.z,'k')
        ax.set_yticklabels([])
        ax.set_title('practical salinity [psu]')
        plt.grid()
        #
        ax = plt.subplot(133)
        ax.plot(self.get_rho(self.z),self.z,'k')
        ax.set_yticklabels([])
        ax.set_title('density [kg/m3]')
        plt.grid()

    def _plot_hv(self):
        df = self._get_df()
        return df.hvplot(grid=True,
                         subplots=True,
                         shared_axes=False,
                         width=200,
                         invert=True,
                         )

    def plot(self):
        print(self)
        return self._plot_hv()

    def __repr__(self):
        return self.name

    def get_temp(self,z):
        ''' get in situ temperature
        '''
        #return self.interp(self.temp, z)
        SA = self.interp(self.SA, z-self.eta)
        CT = self.interp(self.CT, z-self.eta)
        p = self.get_p(z)
        return gsw.conversions.t_from_CT(SA,CT,p)

    def get_s(self, z):
        ''' get practical salinity
        '''
        #return self.interp(self.s, z-self.eta)
        SA = self.interp(self.SA, z-self.eta)
        CT = self.interp(self.CT, z-self.eta)
        p = self.get_p(z)
        return gsw.conversions.SP_from_SA(SA, p, self.lon, self.lat)

    def get_p(self, z):
        ''' get pressure
        '''
        return self.interp(self.p, z)

    def get_theta(self, z):
        ''' get potential temperature
        '''
        SA = self.interp(self.SA, z-self.eta)
        CT = self.interp(self.CT, z-self.eta)
        return gsw.conversions.pt_from_CT(SA,CT)

    def get_rho(self, z, ignore_temp=False):
        p = self.get_p(z)
        SA = self.interp(self.SA, z-self.eta)
        CT = self.interp(self.CT, z-self.eta)
        if ignore_temp:
            CT[:]=self.CT[0]
            print('Uses a uniform conservative temperature in water density computation, CT= %.1f degC' %self.CT[0])
        return gsw.density.rho(SA, CT, p)

    def get_N2(self, z):
        return self.interp(self.N2, z-self.eta, z_in=self.z_mid)

    def get_compressibility(self, z):
        ' returns compressibility in 1/dbar'
        p = self.get_p(z)
        SA = self.interp(self.SA, z-self.eta)
        CT = self.interp(self.CT, z-self.eta)
        return gsw.density.kappa(SA, CT, p)*1e4

    def get_alpha(self, z):
        ' returns thermal expansion in 1/degC'
        p = self.get_p(z)
        SA = self.interp(self.SA, z-self.eta)
        CT = self.interp(self.CT, z-self.eta)
        return gsw.density.alpha(SA, CT, p)

    def get_beta(self, z):
        ' returns thermal expansion in kg/g'
        p = self.get_p(z)
        SA = self.interp(self.SA, z-self.eta)
        CT = self.interp(self.CT, z-self.eta)
        return gsw.density.beta(SA, CT, p)

    def get_adiabatic_lapse_rate(self, z):
        ' returns adiabatic lapse rate in degC/dbar'
        p = self.get_p(z)
        SA = self.interp(self.SA, z-self.eta)
        CT = self.interp(self.CT, z-self.eta)
        return gsw.adiabatic_lapse_rate_from_CT(SA, CT, p)*1e4

    def update_eta(self, eta, t):
        ''' Update isopycnal diplacement and water velocity given a function
        for isopycnal displacement and time

        Parameters
        ----------
        eta: func
            Isopycnal as a function time
        t: float
            Time in seconds
        '''
        self.eta = eta(t)
        self.detadt = (eta(t+.1)-eta(t-.1))/.2

    # utils
    def interp(self, v_in, z_out, z_in=None):
        if z_in is None:
            z_in = self.z
        _z = z_in[~(np.isnan(z_in)|np.isnan(v_in))]
        _v = v_in[~(np.isnan(z_in)|np.isnan(v_in))]
        return interp1d(_z,
                        _v,
                        **self._interp_kwargs
                        )(z_out)
