#
# ------------------------- GPS data -----------------------------------
#

import numpy as np
import pandas as pd
import xarray as xr

import pynmea2
import pickle
import copy

import matplotlib.pyplot as plt
from  matplotlib.dates import date2num, datetime, num2date
from matplotlib.colors import cnames

import folium

from bokeh.io import output_notebook, show
from bokeh.layouts import gridplot
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.plotting import figure

from .utils import get_distance, print_degmin
from .events import event, deployment

gps_attrs = ['label', 'd']

class gps(object):
    ''' Data container for gps data
    '''
    def __init__(self,
                 lon=None,
                 lat=None,
                 time=None,
                 df=None,
                 file=None,
                 label='gps',
                 ):
        self.label = label
        #
        if file is not None:
            if file.split('.')[-1]=='p':
                # pickle file
                self._read_pickle(file)
            elif file.split('.')[-1]=='nc':
                # netcdf file
                self._read_nc(file)
            else:
                assert False, 'You need to pass either a pickle or netcdf file'
            return
        #
        if (df is not None
            and isinstance(df, pd.DataFrame)
            and "lon" in df
            and "lat" in df
            ):
            self.d = df
            return
        #
        self.d = pd.DataFrame()
        if lon is not None and lat is not None and time is not None:
            for _lon, _lat, _time in zip(lon, lat, time):
                self.add(_lon, _lat, _time)

    def __getitem__(self, item):
        if item=='time':
            # time or depth?
            return self.d.index
        elif isinstance(item, str) and item in self.d.columns:
            return self.d[item]
        else:
            # date
            return self.d.loc[item]

    def print_position(self, time):
        lon, lat = self[time][['lon', 'lat']]
        print(print_degmin(lat)+' / '+print_degmin(lon))

    def __add__(self, other):
        self.d = self.d.append(other.d)
        return self

    def __bool__(self):
        return not self.empty()

    def __repr__(self):
        return 'cognac.insitu.gps.gps({})'.format(str(self))

    def __str__(self):
        return self.label+" - {} points".format(len(self.d.index))

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

    def trim(self, t0, t1, inplace=False):
        ''' select data between t0 and t1 '''
        if inplace:
            # beware of exact time indexing:
            # https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#exact-indexing
            #self.d = self.d[t0:t1]
            time = self.d.index
            self.d = self.d.loc[(time >= t0) & (time <= t1)]
        else:
            gp = copy.deepcopy(self)
            gp.trim(t0,t1, inplace=True)
            return gp

    def clean(self, d, inplace=False):
        '''Use deployment to clean gps data'''
        if inplace:
            self.trim(d.start.time, d.end.time)
        else:
            gp = copy.deepcopy(self)
            gp.trim(d.start.time, d.end.time, inplace=True)
            return gp

    def sort(self):
        """ sort by time, inplace method
        """
        self.d.sort_index(inplace=True)

    def resample(self,
                 rule,
                 inplace=False,
                 interpolate=False,
                 **kwargs):
        ''' temporal resampling
        https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.resample.html

        Parameters
        ----------
        rule: DateOffset, Timedelta or str
            Passed to pandas.DataFrame.resample, examples:
                - '10T': 10 minutes
                - '10S': 10 seconds
        inplace: boolean, optional
            turn inplace resampling on, default is False
        interpolate: boolean, optional
            turn on interpolation for upsampling
        kwargs:
            passed to resample
        '''
        if inplace:
            _d = (self
                  .d
                  .resample(rule,
                            **kwargs,
                            )
                  ).mean()
            if interpolate:
                _d = _d.interpolate(method='linear')
            self.d = _d
        else:
            gp = copy.deepcopy(self)
            gp.resample(rule,
                        inplace=True,
                        interpolate=interpolate,
                        **kwargs
                        )
            return gp

    def reindex(self, time, inplace=False):
        """ Reindex time line based on another time line

        Parameters
        ----------
        time: pd.Series, tuple
            Series of timestamps or tuple containing start and end time
        inplace: boolean, optional
            Turn on in place reindexing, default to False
        """
        if inplace:
            if isinstance(time, deployment):
                time = (time.start.time, time.end.time)
            if isinstance(time, tuple):
                dt = self.d.reset_index()['time'].diff().median()
                time = pd.date_range(time[0], time[1], freq=dt).rename('time')
            self.d = self.d.reindex(time)
        else:
            gp = copy.deepcopy(self)
            gp.reindex(time, inplace=True)
            return gp

    def reindex_like(self, other, inplace=False):
        """ Reindex time line based on another gps object

        Parameters
        ----------
        other: gps obj
            gps object with time line to reindex on
        inplace: boolean, optional
            Turn on in place reindexing, default to False
        """
        if inplace:
            self.d = self.d.reindex_like(other.d)
        else:
            gp = copy.deepcopy(self)
            gp.reindex_like(other, inplace=True)
            return gp

    def fillna(self, other=None, interp=None, inplace=False):
        """ fill NaNs with another gps object and/or interpolate existing data
        """
        if not inplace:
            gp = copy.deepcopy(self)
            gp.fillna(other=other, interp=interp, inplace=True)
            return gp
        if other is not None:
            self.d = self.d.fillna(other.d)
        if interp is not None:
            if isinstance(interp, dict):
                kwargs=interp
            else:
                kwargs=dict(method='linear')
            self.d = self.d.interpolate(**kwargs)

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

    def join(self, other, rule):
        ''' Resample and join 2 gps objects
        Returns a dataframe
        '''
        assert self.label is not other.label, 'labels must differ'
        # update velocity fields
        self.compute_velocity(inplace=True)
        other.compute_velocity(inplace=True)
        # resample
        d_resampled = self.resample(rule).d
        other_resampled = other.resample(rule).d
        d = d_resampled.join(other_resampled,
                           how='outer',
                           lsuffix='_'+self.label,
                           rsuffix='_'+other.label,
                           )
        # add separation and its time rate of change
        s = (get_distance(d['lon_'+self.label],
                          d['lat_'+self.label],
                          d['lon_'+other.label],
                          d['lat_'+other.label],
                          )
             .rename('separation')
             )
        ds = s.diff()
        dt = pd.Series(d.index).diff().dt.total_seconds()
        dt.index = ds.index
        v = (ds/dt).rename('separation_velocity')
        d = pd.concat([d, s, v], axis=1)
        # reorder columns
        d = d.reindex(sorted(d.columns), axis=1)
        return d

    #
    def compute_velocity(self, inplace=False):
        """ update velocity
        """
        if not self:
            return
        d = self.d
        # remove duplicates
        d = d[~d.index.duplicated(keep='first')]
        #
        lon, lat = d['lon'], d['lat']
        lon_p = d['lon'].shift(periods=1)
        lat_p = d['lat'].shift(periods=1)
        dl = get_distance(lon, lat, lon_p, lat_p)
        dl_x = get_distance(lon, lat, lon_p, lat)
        dl_x = dl_x * np.sign(lon_p - lon)
        dl_y = get_distance(lon, lat, lon, lat_p)
        dl_y = dl_y * np.sign(lat_p - lat)
        #
        dt = pd.Series(d.index).diff().dt.total_seconds()
        dt.index = dl.index
        V = [(dl/dt).rename('velocity'),
             (dl_x/dt).rename('u'),
             (dl_y/dt).rename('v'),
             ]
        if inplace:
            for v in V:
                self.d[v.name] = v
        else:
            _d = d.drop(columns=[v.name for v in V], errors="ignore")
            return pd.concat([_d,] + V, axis=1)

    ## I/O

    def to_pickle(self, file):
        dictout = {key: getattr(self,key) for key in gps_attrs}
        pickle.dump( dictout, open( file, 'wb' ) )
        print('Data store to '+file)

    def _read_pickle(self, file):
        p = pickle.load( open( file, 'rb' ) )
        for key in gps_attrs:
            setattr(self, key, p[key])

    #
    def to_nc(self, file, **kwargs):
        ds = self.d.to_xarray()
        ds.attrs['label'] = self.label
        ds.to_netcdf(file, **kwargs)
        print('Data store to '+file)

    def _read_nc(self, file):
        ds = xr.open_dataset(file)
        self.d = ds.to_dataframe()
        self.label = ds.label


    ## plot

    def plot(self,
             fac=None,
             label='',
             linestyle='-',
             lw=2.,
             t0=None,
             t1=None,
             bounds=None,
             offset=0.02,
             start_marker='o',
             end_marker='*',
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
            if bounds is None:
                bounds = ax.get_extent()
        else:
            ax = plt.subplot(111)
            if bounds is None:
                bounds = [min(lon)-.1, max(lon)+.1, min(lat)-.1, max(lat)+.1]
        ax.plot(lon, lat, lw=lw, linestyle=linestyle, **kwargs)
        if 'marker' in kwargs:
            del kwargs['marker']
        if start_marker:
            ax.scatter(lon[0], lat[0], 20, marker=start_marker, **kwargs)
        if end_marker:
            ax.scatter(lon[-1], lat[-1], 20, marker=end_marker, **kwargs)
        xoffset = offset * (bounds[1] - bounds[0])
        yoffset = offset * (bounds[3] - bounds[2])
        ax.text(lon[-1]+xoffset, lat[-1]-yoffset, label,
                fontsize=9, **kwargs)

    def plot_folium(self,
                    m,
                    label=None,
                    rule='10T',
                    color='black',
                    ):
        """ add trajectory to folium map
        """
        if label is None:
            label = self.label
        # resample data
        d = self.resample(rule).d
        #
        for t, row in d.iterrows():
            folium.Circle((row['lat'], row['lon']),
                          tooltip=str(t),
                          popup=folium.Popup(label+'<br>'+str(t),
                                             max_width=150,
                                             sticky=True),
                          radius=1e2,
                          color=cnames[color],
                          fill=True,
                         ).add_to(m)
        return m

    def plot_bk(self, unit=None, rule=None):

        if 'velocity' not in self.d:
            self.compute_velocity(inplace=True)

        output_notebook()
        TOOLS = 'pan,wheel_zoom,box_zoom,reset,help'

        # line specs
        lw = 3
        c = 'black'

        if rule is not None:
            d = self.resample(rule).d
        else:
            d = self.d

        def _add_start_end(s, y):
            #_y = y.iloc[y.index.get_loc(_d.start.time), method='nearest')]
            if unit is not None:
                for _d in unit:
                    s.line(x=[_d.start.time, _d.start.time],
                           y=[y.min(), y.max()],
                           color='cadetblue', line_width=2)
                    s.line(x=[_d.end.time, _d.end.time],
                           y=[y.min(), y.max()],
                           color='salmon', line_width=2)

        # create a new plot and add a renderer
        s1 = figure(tools=TOOLS,
                    plot_width=300, plot_height=300,
                    title='longitude',
                    x_axis_type='datetime')
        s1.line('time', 'lon', source=d, line_width=lw, color=c)
        s1.add_tools(HoverTool(
            tooltips=[('Time','@time{%F %T}'),('longitude','@{lon}{0.000f}'),],
            formatters={'@time': 'datetime','@lon' : 'printf',},
            mode='vline'
            ))
        _add_start_end(s1, d['lon'])
        #
        s2 = figure(tools=TOOLS,
                    plot_width=300, plot_height=300,
                    title='latitude',
                    x_axis_type='datetime',
                    x_range=s1.x_range
                    )
        s2.line('time', 'lat', source=d, line_width=lw, color=c)
        s2.add_tools(HoverTool(
            tooltips=[('Time','@time{%F %T}'),('latitude','@{lat}{0.000f}'),],
            formatters={'@time': 'datetime','@lat' : 'printf',},
            mode='vline'
            ))
        _add_start_end(s2, d['lat'])
        #
        s3 = figure(tools=TOOLS,
                    plot_width=300, plot_height=300,
                    title='speed',
                    x_axis_type='datetime',
                    x_range=s1.x_range
                    )
        s3.line('time', 'velocity', source=d, line_width=lw, color=c)
        s3.add_tools(HoverTool(
            tooltips=[('Time','@time{%F %T}'),('Velocity','@{velocity}{%0.000f}'),],
            formatters={'@time': 'datetime','@velocity' : 'printf',},
            mode='vline'
            ))
        _add_start_end(s3, d['velocity'])

        p = gridplot([[s1, s2, s3]])
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


def read_gps_alees(file, label='gps', verbose=False):

    # init gps container
    gp = gps(label=label)
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

def read_gps_lops(file, label='gps', verbose=False):

    # init gps container
    gp = gps(label=label)
    if isinstance(file, list):
        for f in file:
            gp = gp + read_gps_lops(f, verbose=verbose)
    else:
        print('Reads ' + file)
        #gpsfile = pynmea2.NMEAFile(file)
        gpsfile = open(file)
        data = pd.DataFrame()
        t = None
        #for s in gpsfile:
        for s in gpsfile.readlines():
            if verbose:
                print(s)
            #  if len(s.fields)>0 and hasattr(s, 'is_valid') and s.is_valid:
            d = parse_nmea_sentence(s, t=t)
            if d:
                data = data.append(d, ignore_index=True)
                t = d['time']
        #
        gp.d = data.set_index('time')

    return gp

def parse_nmea_sentence(line, t=None):

    try:
        s = pynmea2.parse(line)
    except:
        return None

    time = None

    if 'GPRMC' in s.identifier():
        #print(line)
        #print(s.datestamp, s.timestamp, type(s.datestamp), type(s.timestamp))
        if (len(s.fields)>0
            and isinstance(s.datestamp, datetime.date)
            and isinstance(s.timestamp, datetime.time)
            ):
            try:
                time = datetime.datetime.combine(s.datestamp, s.timestamp)
            except:
                print(line)
    elif 'GPGGA' in s.identifier():
        if t:
            time = datetime.datetime(t.year,
                                     t.month,
                                     t.day,
                                     int(s.data[0][:2]),
                                     int(s.data[0][2:4]),
                                     int(s.data[0][4:6]),
                                     )
            if time<t:
                time+=datetime.timedelta(days=1)
        else:
            time = t

    if (
        time is not None
        and hasattr(s, 'longitude')
        and hasattr(s, 'latitude')
    ):
        return dict(lon=s.longitude,
                    lat=s.latitude,
                    time=time
                    )
    else:
        return None

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
