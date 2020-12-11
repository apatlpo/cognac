
import os, sys, pickle
import csv
import numpy as np
import xarray as xr

import matplotlib.pyplot as plt
from  matplotlib.dates import date2num, datetime
from matplotlib.colors import cnames

#from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
from cartopy.io import shapereader

import folium
import geojson

import gsw
from netCDF4 import Dataset

#
# ------------------------- Maps and Metrics -----------------------------------
#

def dec2sec(dec):
    # return coordinates with seconds
    idec = np.trunc(dec)
    sec = np.abs(dec - idec) * 60.
    return [idec, sec]

def ll_dec(deg, min):
    """ converts lon or lat in deg, min to decimal
    """
    return deg + min/60.

def ll_degmin(l):
    """ converts lon or lat in decimal to deg, min

    Parameters
    ----------
    l: float

    Return
    ------
    deg, min

    """
    return int(l), (l-int(l))*60.

def print_degmin(l):
    ''' Print lon/lat, deg + minutes decimales
    '''
    dm = ll_degmin(l)
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

def plot_bathy(fac, levels=[-2000., -1000., -500., -200., -100.]):
    fig, ax, crs = fac
    #bfile = 'gebco0/GEBCO_2014_2D_5.625_42.0419_8.8046_44.2142.nc'
    ds = xr.open_dataset(_bathy_file)
    cs = ax.contour(ds.lon, ds.lat, ds.elevation, levels,
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

def plot_track_folium(df,
                m,
                label='gps',
                rule=None,
                color='black',
                radius=1e2,
                line=False,
                ):
    """ add trajectory to folium map
    """
    # resample data
    if rule:
        df = df.resample(rule).mean()
    # drop NaNs
    df = df.dropna()
    #
    for t, row in df.iterrows():
        folium.Circle((row['lat'], row['lon']),
                      tooltip=str(t),
                      popup=folium.Popup(label+'<br>'+str(t),
                                         max_width=150,
                                         sticky=True),
                      radius=radius,
                      color=cnames[color],
                      fill=True,
                     ).add_to(m)
    if line:
        folium.Polygon([(row['lat'], row['lon'])
                        for t, row in df.iterrows()
                        ],
                       color=cnames[color],
                       opacity=.5
                      ).add_to(m)
    return m

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
