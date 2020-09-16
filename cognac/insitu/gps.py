#
# ------------------------- GPS data -----------------------------------
#

import pandas as pd
import pynmea2
import pickle
import copy

import matplotlib.pyplot as plt
from  matplotlib.dates import date2num, datetime, num2date

from bokeh.io import output_notebook, show
from bokeh.layouts import gridplot
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.plotting import figure

from .utils import get_distance

gps_attrs = ['d']

class gps(object):
    ''' Data container for gps data
    '''
    def __init__(self, lon=[], lat=[], time=[]):
        self.d = pd.DataFrame()

    def __getitem__(self, item):
        if item=='time':
            # time or depth?
            return self.d.index
        else:
            return self.d[item]

    def __add__(self, other):
        self.d = self.d.append(other.d)
        return self

    def __bool__(self):
        return not self.empty()

    def empty(self):
        return self.d.empty

    def add(self, lon, lat, time, sort=False):
        if not isinstance(lon,list): lon=[lon]
        if not isinstance(lat,list): lat=[lat]
        if not isinstance(time,list): time=[time]
        #d = xr.Dataset({'lon': (['time'], lon), 'lat': (['time'], lat)},
        #               coords = {'time': time})
        d = pd.DataFrame({'lon': lon, 'lat': lat}, index=time)
        d.index.rename('time', inplace=True)
        self.d = self.d.append(d)
        if sort:
            self.d = self.d.sortby('time')

    def trim(self, t0, t1, inplace=True):
        ''' select data between t0 and t1 '''
        if inplace:
            # beware of exact time indexing:
            # https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#exact-indexing
            #self.d = self.d[t0:t1]
            time = self.d.index
            self.d = self.d.loc[(time > t0) & (time <= t1)]
        else:
            gp = copy.deepcopy(self)
            gp.trim(t0,t1)
            return gp

    def clean(self, d, inplace=True):
        '''Use deployment to clean gps data'''
        if inplace:
            self.trim(d.start.time, d.end.time)
        else:
            return self.trim(d.start.time, d.end.time, inplace=False)

    def sort(self):
        self.d.sort_index(inplace=True)

    def resample(self, rule, inplace=False, **kwargs):
        if inplace:
            self.d = self.d.resample(rule, **kwargs).mean()
        else:
            gp = copy.deepcopy(self)
            gp.resample(rule, inplace=True, **kwargs)
            return gp

    #
    def compute_velocity(self):
        if not self:
            return
        dl = get_distance(self.d['lon'] , self.d['lat'],
                          self.d['lon'].shift(periods=1),
                          self.d['lat'].shift(periods=1))
        dt = pd.Series(self.d.index).diff().dt.total_seconds()
        dt.index = dl.index
        v = (dl/dt).rename('velocity')
        if 'velocity' in self.d:
            self.d.update(v)
        else:
            self.d = pd.concat([self.d,v], axis=1)

    #
    def to_pickle(self, file):
        dictout = {key: getattr(self,key) for key in gps_attrs}
        pickle.dump( dictout, open( file, 'wb' ) )
        print('Data store to '+file)

    def _read_pickle(self, file):
        p = pickle.load( open( file, 'rb' ) )
        for key in gps_attrs:
            setattr(self, key, p[key])

    #
    def plot(self, fac=None, ax=None, label='', linestyle='-', lw=2.,
             t0=None, t1=None, ll_lim=None,
             **kwargs):
        lon, lat = self.d['lon'], self.d['lat']
        if t0 is not None:
            lon = lon[t0:]
            lat = lat[t0:]
        if t1 is not None:
            lon = lon[:t1]
            lat = lat[:t1]
        if fac:
            ax = fac[1]
        else:
            ax = plt.subplot(111)
        if ll_lim is None:
            ll_lim = ax.get_extent()
        ax.plot(lon, lat, lw=lw, linestyle=linestyle, **kwargs)
        if 'marker' in kwargs:
            del kwargs['marker']
        ax.scatter(lon[0], lat[0], 10, marker='o', **kwargs)
        ax.scatter(lon[-1], lat[-1], 10, marker='*', **kwargs)
        xoffset = 0.01 * (ll_lim[1] - ll_lim[0])
        yoffset = 0.01 * (ll_lim[3] - ll_lim[2])
        ax.text(lon[-1]+xoffset, lat[-1]-yoffset, label,
                fontsize=9, **kwargs)

    def plot_bk(self):

        if 'velocity' not in self.d:
            self.compute_velocity()

        output_notebook()
        TOOLS = 'pan,wheel_zoom,box_zoom,reset,help'

        _d = self.d

        # create a new plot and add a renderer
        s1 = figure(tools=TOOLS, plot_width=300, plot_height=300, title=None,
                      x_axis_type='datetime')
        s1.line('time', 'lon', source=_d)
        s1.add_tools(HoverTool(
            tooltips=[('time','@time{%T}'),('longitude','@{lon}{%0.3f}'),],
            formatters={'time': 'datetime','longitude' : 'printf',},
            mode='vline'
            ))
        #
        s2 = figure(tools=TOOLS, plot_width=300, plot_height=300, title=None,
                      x_axis_type='datetime', x_range=s1.x_range)
        s2.line('time', 'velocity', source=_d)
        s2.add_tools(HoverTool(
            tooltips=[('time','@time{%T}'),('velocity','@{velocity}{%0.3f}'),],
            formatters={'time': 'datetime','velocity' : 'printf',},
            mode='vline'
            ))

        p = gridplot([[s1, s2]])
        show(p)

    def plot_chunks(self, map, linestyle='-', **kwargs):
        # !! needs update
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
        # !! needs update
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
        # !! needs update
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
        # !! needs update
        for d in dep:
            self.fill_with_d(d)


def read_gps_alees(file, verbose=False):

    # init gps container
    gp = gps()
    if isinstance(file, list):
        for f in file:
            gp = gp + read_gps_alees(f, verbose=verbose)
    else:
        print('Reads ' + file)
        gpsfile = pynmea2.NMEAFile(file)
        data = pd.DataFrame()
        for d in gpsfile:
            if verbose:
                print(d)
            if all([d.datestamp, d.timestamp]):
                time = datetime.datetime.combine(d.datestamp, d.timestamp)
                data = data.append({'lon': d.longitude,
                                    'lat': d.latitude,
                                    'time': time
                                    },
                                    ignore_index=True)
        #
        gp.d = data.set_index('time')

    return gp

def interp_gps(time, gp):
    '''Interpolate lists of gps onto a given timeline
    '''
    nan = [np.NaN for t in time]
    gp_out = gps(lon=nan, lat=nan, time=time)
    lon_out = np.array(gp_out.lon)
    lat_out = np.array(gp_out.lat)
    #
    for i, lgp in enumerate(gp):
        # interp
        lon_i = interp1d(lgps.time, lgp.lon, kind='slinear', bounds_error=False)(time)
        lat_i = interp1d(lgps.time, lgp.lat, kind='slinear', bounds_error=False)(time)
        # fill target array
        ig = np.where( (time>=lgp.time[0]) & (time<=lgp.time[-1]) )
        lon_out[ig] = lon_i[ig]
        lat_out[ig] = lat_i[ig]
    #
    gp_out.lon = lon_out.tolist()
    gp_out.lat = lat_out.tolist()
    return gp_out
