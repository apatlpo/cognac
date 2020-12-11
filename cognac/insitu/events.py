#
# ------------------------- Event/Deployment objects -----------------------------------
#
import os
from glob import glob
import yaml

import pandas as pd
import xarray as xr

import matplotlib.pyplot as plt
from  matplotlib.dates import date2num, datetime
import matplotlib.dates as mdates
from matplotlib.patches import Rectangle
from matplotlib.colors import cnames

import folium
from folium.plugins import MeasureControl, MousePosition

from .utils import load_bathy_contours, dec2sec
from .gps import gps as gps_obj

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

    def __contains__(self, item):
        return item in self._dict

    def __getitem__(self, key):
        return self._dict[key]

    def __setitem__(self, key, value):
        self._dict[key] = value

    def __iter__(self):
        for key, value in self._dict.items():
            if key not in self._skip:
                yield value

    def __str__(self):
        #return {key: value for key, value in self._dict.items() if key not in self._skip}.__str__()
        return self._dict.__str__()

class campaign(object):
    ''' Campaign object, gathers deployments information
    '''

    def __init__(self, file):

        with open(file, 'r') as stream:
            cp = yaml.load(stream)

        default_attr = {'name': 'unknown',
                        'lon_lim': None, 'lat_lim': None,
                        'start': None, 'end': None,
                        'path': './'
                        }
        for key, value in default_attr.items():
            if key in cp:
                setattr(self, key, cp[key])
            else:
                setattr(self, key, value)

        if self.lon_lim and self.lat_lim:
            self.lon_mid = (self.lon_lim[0]+self.lon_lim[1])*.5
            self.lat_mid = (self.lat_lim[0]+self.lat_lim[1])*.5

        if self.start:
            self.start = pd.Timestamp(self.start)
        if self.end:
            self.end = pd.Timestamp(self.end)

        if 'pathp' in cp:
            self.pathp = cp['pathp']
        else:
            self.pathp = self.path+'datap/'

        self._units = {}
        for i, unit in cp['units'].items():
            self._units[i] = objdict(path=self.path)
            for d, value in unit['deployments'].items():
                self._units[i][d] = deployment(label=d,loglines=value)
            for d, value in unit.items():
                if d=='path':
                    self._units[i]['path'] = os.path.join(self.path, value)
                elif d!='deployments':
                    self._units[i][d] = value
                    self._units[i]._skip.append(d)

    def __repr__(self):
        return self.name

    def __getitem__(self, item):
        if item in self._units:
            return self._units[item]
        else:
            return None

    def __iter__(self):
        #for key, value in self._units.items():
        #    yield value
        for key in self._units:
            yield key

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
            ignore=[],
            **kwargs,
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

        if ignore=='all':
            ignore=self._units

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
            if uname not in ignore:
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

        fmtr_lon = "function(dec) {var min= (dec-Math.round(dec))*60; " \
                    +"direction = (dec < 0) ? 'W' : 'E'; " \
                    +"return L.Util.formatNum(dec, 0) + direction + L.Util.formatNum(min, 2);};"
        fmtr_lat = "function(dec) {var min= (dec-Math.round(dec))*60; " \
                    +"direction = (dec < 0) ? 'S' : 'N'; " \
                    +"return L.Util.formatNum(dec, 0) + direction + L.Util.formatNum(min, 2);};"
        MousePosition(lat_formatter=fmtr_lon, lng_formatter=fmtr_lat).add_to(m)

        return m

    def timeline(self, height=.3, legend_loc=3):
        """ Plot the campaign deployment timeline
        """

        fig = plt.figure(figsize=(15,5))
        ax=fig.add_subplot(111)

        y=0
        yticks, yticks_labels = [], []
        starts, ends = [], []
        for uname, u in self.items():
            for d in u:
                start = mdates.date2num(d.start.time)
                end = mdates.date2num(d.end.time)
                rect = Rectangle((start, y-height/2.), end-start, height, color=u['color'])
                ax.add_patch(rect)
                starts.append(start)
                ends.append(end)
            yticks.append(y)
            yticks_labels.append(uname)
            y+=-1

        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks_labels)
        self.add_legend(ax, loc=legend_loc)

        # assign date locator / formatter to the x-axis to get proper labels
        locator = mdates.AutoDateLocator(minticks=3)
        formatter = mdates.AutoDateFormatter(locator)
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(formatter)

        # set the limits
        delta_time = max(ends) - min(starts)
        plt.xlim([min(starts)-delta_time*.05, max(ends)+delta_time*.05])
        plt.ylim([y+1-2*height,2*height])

    def timeline_old(self):
        """ older version of timeline
        """
        fig = plt.figure(figsize=(15,5))
        ax=fig.add_subplot(111)

        y=0
        yticks, yticks_labels = [], []
        for uname, u in self.items():
            for d in u:
                ax.plot([d.start.time,d.end.time],[y,y], lw=4, color=u['color'])
            yticks.append(y)
            yticks_labels.append(uname)
            y+=1

        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks_labels)
        self.add_legend(ax, loc=2)

    def load(self, item, unit=None):
        """ load processed data files

        Returns
        -------
        iridium: dict
            {'unit0': {'deployment0': data, ...}}
        """

        if item=='ship':
            return xr.open_dataset(self.pathp+'ship.nc').to_dataframe()

        data_files = self.get_pfiles(item=item)
        _files = [f.split('/')[-1] for f in data_files]

        units = set(f.split('_')[0] for f in _files)
        if unit:
            if isinstance(unit,str):
                units = [unit]
            else:
                units = [u for u in units if u in unit]

        data = {}
        for u in units:
            _files = [f.split('/')[-1]
                      for f in self.get_pfiles(unit=u,
                                               item=item,
                                               )
                      ]
            deployments = [f.split('_')[2].split('.')[0] for f in _files]
            data[u] = {d: gps_obj(file=self.get_pfiles(unit=u,
                                                       item=item,
                                                       deployment=d,
                                                       )
                              )
                       for d in deployments
                       }
        return data

    def get_pfiles(self,
                   unit='*',
                   item='*',
                   deployment='*',
                   extension='nc',
                   ):
        """ Return processed data files
        Parameters
        ----------
        unit, item, d, extention: str, optional
            Defaults: '*', '*', '*', 'nc'
            Typical file path: self.pathp+unit+'_'+item+'_'+deployment+'.'+extension

        """
        if any([_=='*' for _ in [unit, item, deployment]]):
            return glob(self.pathp+unit+'_'+item+'_'+deployment+'.'+extension)
        else:
            return self.pathp+unit+'_'+item+'_'+deployment+'.'+extension
