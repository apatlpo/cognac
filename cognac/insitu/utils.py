
import os, sys, pickle, glob
import csv, yaml
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
#from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
from cartopy.io import shapereader
from  matplotlib.dates import date2num, datetime, num2date
import gsw
from scipy.interpolate import interp1d
from netCDF4 import Dataset

# gps data
#import pynmea2


#
# ------------------------- Event/Deployment objects -----------------------------------
#

#_ll_lim_default = [6.4, 6.6, 42.92, 43.2]
_ll_lim_default = [6., 6.6, 42.7, 43.2]

def dec2sec(dec):
    # return coordinates with seconds
    idec = np.trunc(dec)
    sec = np.abs(dec - idec) * 60.
    return [idec, sec]

class event(object):

    def __init__(self, label=None, logline=None, coord_min=True):

        # split string
        l=logline.split()

        # label
        self.label = label

        # time information
        # assumes: 02/09/2016 05:35:00 7 17.124 43 19.866
        self.time = datetime.datetime(int(l[0][6:10]), int(l[0][3:5]), int(l[0][0:2]),
                                      int(l[1][0:2]), int(l[1][3:5]), int(l[1][6:8]))
        #self.time = date2num(self.dtime)

        if len(l)>2:
            # lon, lat data
            if coord_min:
                self.lon = float(l[2]) + float(l[3])/60.
                self.lat = float(l[4]) + float(l[5])/60.
            else:
                self.lon = float(l[2])
                self.lat = float(l[3])

    def __str__(self):
        print('-')
        print('Event label: '+str(self.label)) #+'\n')
        print('Time: ')
        print(self.time)
        lon=self.lon
        lat=self.lat
        print('Lon:'+str(lon)+' = '+str(dec2sec(lon)[0])+'deg '+str(dec2sec(lon)[1]))
        print('Lat:' + str(lat) + ' = ' + str(dec2sec(lat)[0]) + 'deg ' + str(dec2sec(lat)[1]))
        return ''


class deployment(object):

    def __init__(self, label=None, start=None, end=None, loglines=None):
        self.label=label
        if loglines is None:
            self.start=start
            self.end=end
        else:
            self.start=event(label='start', logline=loglines[0])
            self.end = event(label='end', logline=loglines[1])

    def __str__(self):
        print('--')
        print('Deployment label: '+self.label)
        print('Start:')
        print(self.start)
        print('End:')
        print(self.end)
        return ''

    def plot_time(self, axis=None, y0=0., dy=0.5, **kwargs):
        t0 = self.start.time
        t1 = self.end.time
        rect = Rectangle((t0, y0-dy/2.), t1-t0, dy, **kwargs)
        axis.add_patch(rect)
        return

    def plot_map(self, map, line=True, label=True, **kwargs):
        #
        x0, y0 = map(self.start.lon, self.start.lat)
        map.scatter(x0, y0, 10, marker='o', **kwargs)
        #
        x1, y1 = map(self.end.lon, self.end.lat)
        map.scatter(x1, y1, 10, marker='*', **kwargs)
        #
        if line:
            map.plot([x0,x1],[y0,y1],'-', **kwargs)
        if label:
            xoffset = 0.02 * (map.xmax - map.xmin)
            yoffset = 0.02 * (map.ymax - map.ymin)
            plt.text(x0+xoffset, y0-yoffset, self.label, fontsize=9, **kwargs)
        return

class objdict(object):
    ''' Dict like object that treats path and color as attributes
    '''
    def __init__(self, *args,**kwargs):
        self._dict= dict(*args,**kwargs)

    def __getitem__(self, key):
        #if key not in ['path', 'color']:
        #    val = self._dict[key]
        #return val
        return self._dict[key]

    def __setitem__(self, key, value):
        #if key not in ['path', 'color']:
        #    self._dict[key] = value
        #else:
        #    setattr(self, key, value)
        self._dict[key] = value

    def __iter__(self):
        for key, value in self._dict.items():
            if key not in ['path', 'color']:
                yield value

    def __str__(self):
        return self._dict.__str__()

class campaign(object):
    ''' Campaign object, gathers deployments information
    '''

    def __init__(self, file):

        with open(file, 'r') as stream:
            cp = yaml.load(stream)

        default_attr = {'name': 'unknown',
                        'lon_lim': None, 'lat_lim': None,
                        'path': './'}
        for key, value in default_attr.items():
            if key in cp:
                setattr(self, key, cp[key])
            else:
                setattr(self, key, value)

        if 'pathp' in cp:
            self.pathp = cp['pathp']
        else:
            self.pathp = self.path+'datap/'

        self._units = {}
        for i, idep in cp['units'].items():
            self._units[i] = objdict(path=self.path)
            for d, value in idep.items():
                if d[0] is not '_':
                    self._units[i][d] = deployment(label=d,loglines=value)
                elif d == '_path':
                    self._units[i]['path'] = os.path.join(self.path,value)
                else:
                    self._units[i][d[1:]] = value

    def __getitem__(self,item):
        if item in self._units:
            return self._units[item]
        else:
            return None

    def load_data(self):
        ''' load processed data
        '''
        pass


#
# ------------------------- source data -----------------------------------
#

# need to move to an independent file

class emission_data(object):
    ''' Data container for emission data
    '''
    def __init__(self, time=None, sound=None, lon=None, lat=None):
        if time is None:
            self.time = []
            self.sound = []
            self.lon = []
            self.lat = []
        else:
            self.time = time
            self.sound = sound
            self.lon = lon
            self.lat = lat

    def add(self, time, sound, lon, lat, sort=False):
        self.time+=time
        self.sound+=sound
        self.lon+=lon
        self.lat+=lat
        #self.time.append(time)
        #self.sound.append(sound)
        #self.lon.append(lon)
        #self.lat.append(lat)
        if sort:
            i=np.argsort(time)
            self.time=self.time[i]
            self.sound=self.sound[i]
            self.lon=self.lon[i]
            self.lat=self.lat[i]

    def trim(self, t0, t1):
        ''' select data between t0 and t1 '''
        self.sound = [gsound for i, gsound in enumerate(self.sound) if self.time[i] > t0 and self.time[i] < t1]
        self.lon = [glon for i, glon in enumerate(self.lon) if self.time[i] > t0 and self.time[i] < t1]
        self.lat = [glat for i, glat in enumerate(self.lat) if self.time[i] > t0 and self.time[i] < t1]
        ltime = [gtime for i, gtime in enumerate(self.time) if self.time[i] > t0 and self.time[i] < t1]
        self.time=ltime

    def clean(self, d):
        '''Use deployment to clean gps data'''
        self.trim(d.start.time, d.end.time)




#
# ------------------------- Maps and Metrics -----------------------------------
#


def get_distance(lon1 , lat1 , lon2 , lat2):
    ''' wrapper around distance calculator
    '''
    if isinstance(lon1, list):
        lon1 = np.array(lon1)
        lon2 = np.array(lon2)
        lat1 = np.array(lat1)
        lat2 = np.array(lat2)
    #return _distance_on_spherical_earth(lon1 , lat1 , lon2 , lat2)
    #return
    earth_radius = 6373.e3
    d2r = np.pi/180.
    #
    return earth_radius * np.arccos(np.sin(d2r*lat1)*np.sin(d2r*lat2) \
                                  +np.cos(d2r*lat1)*np.cos(d2r*lat2)*np.cos(d2r*(lon2-lon1)))


#def _distance_on_spherical_earth(lon1 , lat1 , lon2 , lat2) :
#    earth_radius = 6373.e3
#    d2r=np.pi/180.
#    return earth_radius * np.arccos(np.sin(d2r*lat1)*np.sin(d2r*lat2) \
#                                  +np.cos(d2r*lat1)*np.cos(d2r*lat2)*np.cos(d2r*(lon2-lon1)))

def lstr(l):
    ''' Print lon/lat, deg + minutes decimales
    '''
    return '%d deg %.5f' %(int(l), (l-int(l))*60.)


def plot_map(fig=None, coast='med', figsize=(10, 10), ll_lim=None, cp=None):
    crs = ccrs.PlateCarree()
    #
    if fig is None:
        fig = plt.figure(figsize=figsize)
    else:
        fig.clf()

    if cp is not None:
        ll_lim = cp.lon_lim+cp.lat_lim
    elif ll_lim is None:
        ll_lim = _ll_lim_default

    ax = fig.add_subplot(111, projection=crs)
    ax.set_extent(ll_lim, crs=crs)
    gl = ax.gridlines(crs=crs, draw_labels=True, linewidth=2, color='k',
                      alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    #
    if coast in ['10m', '50m', '110m']:
        ax.coastlines(resolution=coast, color='k')
    elif coast in ['auto', 'coarse', 'low', 'intermediate', 'high', 'full']:
        shpfile = shapereader.gshhs('h')
        shp = shapereader.Reader(shpfile)
        ax.add_geometries(
            shp.geometries(), crs, edgecolor='black', facecolor='none')
    elif coast is 'med':
        # conda install -c conda-forge gdal
        # ogr2ogr -f "ESRI Shapefile" med_coast.shp /Users/aponte/.local/share/cartopy/shapefiles/gshhs/h/GSHHS_h_L1.shp -clipsrc 6 7 42.5 43.5
        shp = shapereader.Reader(os.getenv('HOME')+'/data/OSM/med/med_coast')
        for record, geometry in zip(shp.records(), shp.geometries()):
            ax.add_geometries([geometry], crs, facecolor='None', edgecolor='black')
    elif coast is 'med_high':
        # conda install -c conda-forge gdal
        # ogr2ogr -f "ESRI Shapefile" med_high_coast.shp ../coastlines/lines.shp -clipsrc 6 7 42.5 43.5
        shp = shapereader.Reader(os.getenv('HOME')+'/data/OSM/med/med_high_coast')
        for record, geometry in zip(shp.records(), shp.geometries()):
            ax.add_geometries([geometry], crs, facecolor='None', edgecolor='black')

    return [fig, ax, crs]

def plot_bathy(fac):
    fig, ax, crs = fac
    ### GEBCO bathymetry
    bathy = os.getenv('HOME') + \
            '/data/bathy/RN-1994_1473924981206/GEBCO_2014_2D_5.625_42.0419_8.8046_44.2142.nc'
    ds = xr.open_dataset(bathy)
    cs = ax.contour(ds.lon, ds.lat, ds.elevation, [-2000., -1000., -500., -200., -100.],
                    linestyles='-', colors='black', linewidths=0.5, )
    plt.clabel(cs, cs.levels, inline=True, fmt='%.0f', fontsize=9)


#
# ------------------------- EOS wrappers -----------------------------------
#


def dens0(S,T,P):
    # ctd file contain in situ temperature and practical salinity
    # CT= conservative temperature
    # SA= absolute salinity
    # p is sea pressure (dbar)
    SA = gsw.SA_from_SP(S, P, 6., 42.)
    CT = gsw.CT_from_t(SA, T, P)
    # sound speed
    #C = gsw.sound_speed(SA, CT, P)
    return gsw.sigma0(SA,CT)+1000.


def alphabeta(S,T,P):
    SA = gsw.SA_from_SP(S, P, 6., 42.)
    CT = gsw.CT_from_t(SA, T, P)
    return gsw.alpha(SA, CT,P), gsw.beta(SA, CT,P)


#
# ------------------------- Time plotting -----------------------------------
#


def get_time_ticks():

    # for time plotting purposes
    t0 = datetime.datetime(2016, 9, 2, 0)
    t1 = datetime.datetime(2016, 9, 4, 12)

    ### plot time line
    # assign date locator / formatter to the x-axis to get proper labels
    # dloc = mdates.DayLocator()
    # dform = mdates.DateFormatter('%Y-%m-%d')
    # hloc = mdates.HourLocator(interval=6)
    # hform = mdates.DateFormatter('%H:%M')
    # locator = AutoDateLocator(minticks=3)
    # formatter = AutoDateFormatter(locator)
    # formatter = DateFormatter('%Y-%m-%d %H:%M:%S')
    # plt.gcf().axes[0].xaxis.set_major_formatter(formatter)
    # ax.xaxis.set_minor_locator(hloc)
    # ax.xaxis.set_minor_formatter(hform)
    # ax.xaxis.set_major_locator(hloc)
    # ax.xaxis.set_major_formatter(dform)
    # fig.autofmt_xdate()
    tck = [];
    tck_label = [];
    tdays = []
    t = t0
    while t < t1 + datetime.timedelta(hours=6):
        tck.append(t)
        if t.hour is 12:
            tck_label.append(t.strftime('%Y-%m-%d'))
        else:
            tck_label.append('')
        t += datetime.timedelta(hours=6)
        if t.hour is 0:
            tdays.append(date2num(t))

    return t0, t1, tck, tck_label, tdays



def get_time_ticks_zoomed(iday):

    # compute ticks
    t0 = datetime.datetime(2016, 9, 2 + iday, 0)
    t1 = datetime.datetime(2016, 9, 3 + iday, 0)

    # major tick
    tck = [];
    tck_label = [];
    t = t0
    while t < t1 + datetime.timedelta(hours=1):
        tck.append(t)
        if np.mod(t.hour, 3) == 0:
            tck_label.append(t.strftime('%H'))
        else:
            tck_label.append('')
        t += datetime.timedelta(hours=1)

    # minor tick
    tck_minor = [];
    t = t0
    while t < t1 + datetime.timedelta(hours=1):
        tck_minor.append(t)
        t += datetime.timedelta(minutes=15)

    # string of day
    day_str = t0.strftime('%Y-%m-%d')

    return t0, t1, tck, tck_label, day_str, tck_minor
