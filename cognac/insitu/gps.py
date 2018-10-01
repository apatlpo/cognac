#
# ------------------------- GPS data -----------------------------------
#

import pandas as pd
import matplotlib.pyplot as plt
from  matplotlib.dates import date2num, datetime, num2date
import pynmea2

from .utils import ll_lim_default

class gps_data(object):
    ''' Data container for gps data
    '''
    def __init__(self, lon=[], lat=[], time=[]):
        self.d = pd.DataFrame()
        #self.d = xr.Dataset({'lon': (['time'], lon), 'lat': (['time'], lat)},
        #                    coords = {'time': time})

    def __getitem__(self, item):
        return self.d[item]

    def __add__(self, other):
        self.d = self.d.append(other.d)
        return self

    def add(self, lon, lat, time, sort=False):
        if not isinstance(lon,list): lon=[lon]
        if not isinstance(lat,list): lat=[lat]
        if not isinstance(time,list): time=[time]
        d = xr.Dataset({'lon': (['time'], lon), 'lat': (['time'], lat)},
                       coords = {'time': time})
        self.d.update(d, inplace=True)
        if sort:
            self.d = self.d.sortby('time')

    def trim(self, t0, t1):
        ''' select data between t0 and t1 '''
        self.d = self.d[t0:t1]

    def clean(self, d):
        '''Use deployment to clean gps data'''
        self.trim(d.start.time, d.end.time)

    def plot(self, fac, label='', linestyle='-', lw=2., t0=None, t1=None, ll_lim=None, \
              **kwargs):
        fig, ax, crs = fac
        lon, lat = self.d['lon'], self.d['lat']
        if t0 is not None:
            lon = lon[t0:]
            lat = lat[t0:]
        if t1 is not None:
            lon = lon[:t1]
            lat = lat[:t1]
        if ll_lim is None:
            ll_lim = ll_lim_default
        ax.plot(lon, lat, lw=lw, linestyle=linestyle, transform=crs, **kwargs)
        if 'marker' in kwargs:
            del kwargs['marker']
        ax.scatter(lon[0], lat[0], 10, marker='o', transform=crs, **kwargs)
        ax.scatter(lon[-1], lat[-1], 10, marker='*', transform=crs, **kwargs)
        xoffset = 0.01 * (ll_lim[1] - ll_lim[0])
        yoffset = 0.01 * (ll_lim[3] - ll_lim[2])
        plt.text(lon[-1]+xoffset, lat[-1]-yoffset, label, fontsize=9, transform=crs, **kwargs)

    def plot_chunks(self, map, linestyle='-', **kwargs):
        # needs update
        from itertools import groupby
        long = [list(v) for k, v in groupby(self.lon, np.isfinite) if k]
        latg = [list(v) for k, v in groupby(self.lat, np.isfinite) if k]
        for i, (lon, lat) in enumerate(zip(long, latg)):
            x, y = map(np.array(lon), np.array(lat))
            map.plot(x, y, linestyle=linestyle, **kwargs)
            map.scatter(x[0], y[0], 10, marker='o', **kwargs)
            map.scatter(x[-1], y[-1], 10, marker='*', **kwargs)
            xoffset = 0.01 * (map.xmax - map.xmin)
            yoffset = 0.01 * (map.ymax - map.ymin)
            plt.text(x[-1] + xoffset, y[-1] - yoffset, 'd'+str(i+1), fontsize=9, **kwargs)

    def plot_scatter(self, map, label='', markersize=10, **kwargs):
        # needs update
        x, y = map(np.array(self.lon), np.array(self.lat))
        map.scatter(x, y, markersize, **kwargs)
        if 'marker' in kwargs:
            del kwargs['marker']
        map.scatter(x[0], y[0], 10, marker='o', **kwargs)
        map.scatter(x[-1], y[-1], 10, marker='*', **kwargs)
        if 'edgecolor' in kwargs:
            del kwargs['edgecolor']
        xoffset = 0.01 * (map.xmax - map.xmin)
        yoffset = 0.01 * (map.ymax - map.ymin)
        plt.text(x[-1]+xoffset, y[-1]-yoffset, label, fontsize=9, **kwargs)

    def fill_with_d(self, d):
        # needs update
        # convert to array
        lon = np.array(self.lon)
        lat = np.array(self.lat)
        time = np.array(self.time)
        # add first point
        i = np.where(time >= d.start.time)[0][0]
        if np.isnan(lon[i]):
            lon[i] = d.start.lon
            lat[i] = d.start.lat
        # add last point
        i = np.where(time <= d.end.time)[0][-1]
        if np.isnan(lon[i]):
            lon[i] = d.end.lon
            lat[i] = d.end.lat
        # fill in gaps
        inonan = np.where( (time >= d.start.time) & (time <= d.end.time) & (~np.isnan(lon)))
        inan = np.where( (time >= d.start.time) & (time <= d.end.time) & (np.isnan(lon)) )
        lon[inan] = interp1d(time[inonan], lon[inonan], kind='slinear',
                                  bounds_error=False)(time[inan])
        lat[inan] = interp1d(time[inonan], lat[inonan], kind='slinear',
                                  bounds_error=False)(time[inan])
        # back to list
        self.lon = lon.tolist()
        self.lat = lat.tolist()
        # store array of interpolated values
        self.i_filled = inan

    def fill_with_deps(self, dep):
        '''wrapper for lists of deployments'''
        # needs update
        for d in dep:
            self.fill_with_d(d)


def read_gps_tois(file, verbose=False):

    # init gps container
    gps = gps_data()

    if isinstance(file, list):

        for f in file:
            gps = gps + read_gps_tois(f)

    else:

        print('Reads ' + file)

        gpsfile = pynmea2.NMEAFile(file)

        data = pd.DataFrame()
        for d in gpsfile:
            if verbose:
                print(d)
            time = datetime.datetime.combine(d.datestamp, d.timestamp)
            data = data.append({'lon': d.longitude, 'lat': d.latitude, 'time': time}, ignore_index=True)
        #
        gps.d = data.set_index('time')

    return gps

def interp_gps(time, gps):
    '''Interpolate lists of gps onto a given timeline
    '''
    nan = [np.NaN for t in time]
    gps_out = gps_data(lon=nan, lat=nan, time=time)
    lon_out = np.array(gps_out.lon)
    lat_out = np.array(gps_out.lat)
    #
    for i, lgps in enumerate(gps):
        # interp
        lon_i = interp1d(lgps.time, lgps.lon, kind='slinear', bounds_error=False)(time)
        lat_i = interp1d(lgps.time, lgps.lat, kind='slinear', bounds_error=False)(time)
        # fill target array
        ig = np.where( (time>=lgps.time[0]) & (time<=lgps.time[-1]) )
        lon_out[ig] = lon_i[ig]
        lat_out[ig] = lat_i[ig]
    #
    gps_out.lon = lon_out.tolist()
    gps_out.lat = lat_out.tolist()
    return gps_out
