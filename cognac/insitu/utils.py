
import os, sys, pickle, glob
import csv, yaml
import numpy as np
import xarray as xr
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from  matplotlib.dates import date2num, datetime, num2date
from matplotlib.colors import cnames
#from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
from cartopy.io import shapereader

import folium
from folium.plugins import MeasureControl, MousePosition
import geojson

import gsw
from netCDF4 import Dataset

#
# ------------------------- Event/Deployment objects -----------------------------------
#

#_ll_lim_default = [6.4, 6.6, 42.92, 43.2]
_ll_lim_default = [6., 6.6, 42.7, 43.2]

class event(object):

    def __init__(self, label=None, logline=None, coord_min=True):

        # split string
        l=logline.split()

        # label
        self.label = label

        # time information
        # assumes: 02/09/2016 05:35:00 7 17.124 43 19.866
        self.time = datetime.datetime(int(l[0].split('/')[2]),
                                      int(l[0].split('/')[1]),
                                      int(l[0].split('/')[0]),
                                      int(l[1].split(':')[0]),
                                      int(l[1].split(':')[1]),
                                      int(l[1].split(':')[2]))

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

    def plot_time(self, ax, y0=0., dy=0.5, **kwargs):
        #t0 = self.start.time
        #t1 = self.end.time
        #rect = Rectangle((t0, y0-dy/2.), t1-t0, dy, **kwargs)
        #ax.add_patch(rect)
        #ax.plot([t0,t1],[y0,y0],lw=20)
        #ax.plot([self.start.time,self.end.time],[y0,y0],lw=20)
        pass

    def plot_on_map(self, ax, line=True, label=False, yshift=1, s=30,
                    **kwargs):
        #
        x0, y0 = self.start.lon, self.start.lat
        x1, y1 = self.end.lon, self.end.lat
        #
        ax.scatter(x0, y0, s, marker='o', **kwargs)
        ax.scatter(x1, y1, s, marker='*', **kwargs)
        #
        if line:
            ax.plot([x0, x1], [y0, y1], '-', **kwargs)
        if label:
            xb, yb = ax.get_xbound(), ax.get_ybound()
            xoffset = 0.02 * (xb[1]-xb[0])
            yoffset = 0.02 * (yb[1]-yb[0]) * yshift
            if type(label) is not str:
                label = self.label
            ax.text(x0+xoffset, y0+yoffset, label, fontsize=10)
        return

class objdict(object):
    ''' Dict like object that treats some parameters (e.g. path and color)
    as attributes
    '''
    def __init__(self, *args,**kwargs):
        self._dict = dict(*args,**kwargs)
        self._skip = ['path']

    def __getitem__(self, key):
        return self._dict[key]

    def __setitem__(self, key, value):
        self._dict[key] = value

    def __iter__(self):
        for key, value in self._dict.items():
            if key not in self._skip:
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

        if self.lon_lim and self.lat_lim:
            self.lon_mid = (self.lon_lim[0]+self.lon_lim[1])*.5
            self.lat_mid = (self.lat_lim[0]+self.lat_lim[1])*.5

        if 'pathp' in cp:
            self.pathp = cp['pathp']
        else:
            self.pathp = self.path+'datap/'

        self._units = {}
        for i, idep in cp['units'].items():
            self._units[i] = objdict(path=self.path)
            for d, value in idep.items():
                if d[0]!='_':
                    self._units[i][d] = deployment(label=d,loglines=value)
                elif d=='_path':
                    self._units[i]['path'] = os.path.join(self.path,value)
                else:
                    self._units[i][d[1:]] = value
                    self._units[i]._skip.append(d[1:])

    def __repr__(self):
        return self.name

    def __getitem__(self, item):
        if item in self._units:
            return self._units[item]
        else:
            return None

    def __iter__(self):
        for key, value in self._units.items():
            yield value

    def items(self):
        for key, value in self._units.items():
            yield key, value

    def add_legend(self, ax, labels=None, **kwargs):
        """ Add legend for units on an axis
        """
        from matplotlib.lines import Line2D
        if labels is None:
            labels = list(self._units)
        custom_lines = [Line2D([0], [0], color=self[label]['color'], lw=4)
                        for label in labels
                        ]
        ax.legend(custom_lines, labels, **kwargs)

    def map(self,
            width='100%',
            height='100%',
            tiles='Cartodb Positron',
            ignore_labels=[],
            ):
        ''' Plot overview map with folium

        Parameters:
        ----------

        tiles: str
            tiles used, see `folium.Map?``
                - "OpenStreetMap"
                - "Mapbox Bright" (Limited levels of zoom for free tiles)
                - "Mapbox Control Room" (Limited levels of zoom for free tiles)
                - "Stamen" (Terrain, Toner, and Watercolor)
                - "Cloudmade" (Must pass API key)
                - "Mapbox" (Must pass API key)
                - "CartoDB" (positron and dark_matter)

        '''

        m = folium.Map(location=[self.lat_mid, self.lon_mid],
                       width=width,
                       height=height,
                       zoom_start=11,
                       tiles=tiles,
                      )

        # bathymetric contours
        contours_geojson = load_bathy_contours()
        tooltip = folium.GeoJsonTooltip(fields=['title'],
                                        aliases=['depth'],
                                        )
        popup = folium.GeoJsonPopup(fields=['title'],
                                    aliases=['depth'],
                                    )
        #colorscale = branca.colormap.linear.Greys_03.scale(levels[-1],levels[0])
        def style_func(feature):
            return {'color':   feature['properties']['stroke'], #colorscale(feature['properties']['level-value']),
                    'weight':  3, #x['properties']['stroke-width'],
                    #'fillColor': x['properties']['fill'],
                    'opacity': 1.,
                    #'popup': feature['properties']['title'],
                   }
        folium.GeoJson(contours_geojson,
                       name='geojson',
                       style_function=style_func,
                       tooltip=tooltip,
                       popup=popup,
                       ).add_to(m)

        # campaign details
        for uname, u in self.items():
            if uname not in ignore_labels:
                for d in u:
                    folium.Polygon([(d.start.lat, d.start.lon),
                                    (d.end.lat, d.end.lon)
                                    ],
                                   tooltip=uname+' '+d.label+'<br>'
                                            +str(d.start.time)+'<br>'
                                            +str(d.end.time),
                                   color=cnames[u['color']],
                                   dash_array='10 20',
                                   opacity=.5
                                  ).add_to(m)
                    folium.Circle((d.start.lat, d.start.lon),
                                  tooltip=uname+' '+d.label+'<br>'
                                            +str(d.start.time),
                                  radius=2*1e2,
                                  color=cnames[u['color']],
                                 ).add_to(m)
                    folium.Circle((d.end.lat, d.end.lon),
                                  tooltip=uname+' '+d.label+'<br>'
                                            +str(d.end.time),
                                  radius=1e2,
                                  color=cnames[u['color']],
                                 ).add_to(m)

        # useful plugins
        MeasureControl().add_to(m)
        MousePosition().add_to(m)

        return m

    def load_data(self):
        ''' load processed data
        '''
        pass

#
# ------------------------- Maps and Metrics -----------------------------------
#

def dec2sec(dec):
    # return coordinates with seconds
    idec = np.trunc(dec)
    sec = np.abs(dec - idec) * 60.
    return [idec, sec]

def ll_degmin(l):
    ''' Print lon/lat, deg + minutes decimales
    '''
    return '%d deg %.5f' %(int(l), (l-int(l))*60.)

def get_distance(lon1 , lat1 , lon2 , lat2):
    ''' wrapper around distance calculator in meters
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

def plot_map(fig=None, coast='med', figsize=(10, 10), ll_lim=None, cp=None):
    crs = ccrs.PlateCarree()
    #
    if fig is None:
        fig = plt.figure(figsize=figsize)
    else:
        fig.clf()

    if ll_lim is None:
        if cp is not None:
            ll_lim = cp.lon_lim+cp.lat_lim
        else:
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
    elif coast=='med':
        # conda install -c conda-forge gdal
        # ogr2ogr -f "ESRI Shapefile" med_coast.shp /Users/aponte/.local/share/cartopy/shapefiles/gshhs/h/GSHHS_h_L1.shp -clipsrc 5 42. 7 44.
        shp = shapereader.Reader(os.getenv('HOME')+'/data/OSM/med/med_coast.shp')
        for record, geometry in zip(shp.records(), shp.geometries()):
            ax.add_geometries([geometry], crs, facecolor='None', edgecolor='black')
    elif coast=='med_high':
        # conda install -c conda-forge gdal
        # ogr2ogr -f "ESRI Shapefile" med_high_coast.shp ../coastlines/lines.shp -clipsrc 6 7 42.5 43.5
        shp = shapereader.Reader(os.getenv('HOME')+'/data/OSM/med/med_high_coast.shp')
        for record, geometry in zip(shp.records(), shp.geometries()):
            ax.add_geometries([geometry], crs, facecolor='None', edgecolor='black')

    return [fig, ax, crs]


# GEBCO bathymetry
_bathy_file = os.getenv('HOME') + '/data/bathy/' \
        'gebco1/gebco_2020_n44.001617431640625_s41.867523193359375_w4.61151123046875_e8.206787109375.nc'
_bathy_dir = '/'.join(_bathy_file.split('/')[:-1])

def plot_bathy(fac):
    fig, ax, crs = fac
    #bfile = 'gebco0/GEBCO_2014_2D_5.625_42.0419_8.8046_44.2142.nc'
    ds = xr.open_dataset(_bathy_file)
    cs = ax.contour(ds.lon, ds.lat, ds.elevation, [-2000., -1000., -500., -200., -100.],
                    linestyles='-', colors='black', linewidths=0.5, )
    plt.clabel(cs, cs.levels, inline=True, fmt='%.0f', fontsize=9)

def store_bathy_contours(contour_file='contours.geojson',
                         levels=[0, 100, 500, 1000, 2000, 3000],
                         ):
    """ Store bathymetric contours as a geojson
    The geojson may be used for folium plots
    """
    # Create contour data lon_range, lat_range, Z
    depth = -xr.open_dataset(_bathy_file)['elevation']
    contours = depth.plot.contour(levels=levels, cmap='gray_r')

    # Convert matplotlib contour to geojson
    from geojsoncontour import contour_to_geojson
    contours_geojson = contour_to_geojson(
                            contour=contours,
                            geojson_filepath=os.path.join(bathy_dir,
                                                          contour_file),
                            ndigits=3,
                            unit='m',
                        )
def load_bathy_contours(contour_file='contours.geojson'):
    ''' load bathymetric contours as geojson
    '''
    with open(os.path.join(_bathy_dir,contour_file), 'r') as f:
        contours = geojson.load(f)
    return contours

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
        if t.hour==12:
            tck_label.append(t.strftime('%Y-%m-%d'))
        else:
            tck_label.append('')
        t += datetime.timedelta(hours=6)
        if t.hour==0:
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
