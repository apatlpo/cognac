#
# ------------------------- RBR data -----------------------------------
#

from glob import glob
import copy

import pandas as pd
import xarray as xr

from bokeh.io import output_notebook, show
from bokeh.layouts import gridplot
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.plotting import figure

import pickle

rbr_attrs = ['d', 'id', 'file', 'dt']

class rbr(object):
    """ Contains DST inclino data
    """
    def __init__(self, file, id=None, **kwargs):
        if file is not None:
            if '*' in file:
                files_match = glob(file)
                assert len(files_match)==1, \
                    'Invalid file name or multiple file match'
                file = files_match[0]
            self.file = file
            extension = file.split('.')[-1]
            if extension=='txt':
                self._file_striped = self.file.rstrip('.txt')
                self.d = self._read_txt(file, **kwargs)
                self.dt = (self.d.index[1]-self.d.index[0]).total_seconds()
            elif extension=='p':
                self._file_striped = self.file.rstrip('.p')
                self._read_pickle(file)
                return
            elif extension=='nc':
                self._file_striped = self.file.rstrip('.nc')
                self._read_nc(file, **kwargs)
                return
        else:
            print('You need to provide a file name')
        self.id = id

    def __str__(self):
        return 'rbr '+str(self.id)+self.d.__str__()

    def __getitem__(self, item):
        if item=='time':
            # time or depth?
            return self.d.index
        else:
            return getattr(self.d, item)

    def _read_txt(self, file, **kwargs):
        #version = 'standard'
        version = kwargs['version'] if 'version' in kwargs else 'standard'
        if version=='lopstech':
            d = pd.read_table(file,
                              names=['temperature','pressure','sea_pressure','depth'],
                              skiprows=[0], sep=',', decimal='.', index_col=0, parse_dates=[0])
            # dtype={'temp': np.float64} is ignored if passed to read_table
            d.temperature = d.temperature.astype(float)
            #d = d.set_index('time')
            d.index.name='time'
        elif version=='standard':
            d = pd.read_csv(file,
                            comment='//',
                            header=1,
                            sep='\s+|\t+',
                            parse_dates=[[3,4]],
                            infer_datetime_format=True,
                            index_col=False,
                            names=['Cruise', 'Station', 'Type',
                                   'yyyy-mm-ddT', 'hh:mm:ss.sss',
                                   'Longitude [degrees_east]', 'Latitude [degrees_north]',
                                   'Bot. Depth [m]','depth', 'temperature',
                                   'pressure', 'Sea pressure [dbar]',
                                  ]
                            )
            d = (d
                 .rename(columns={'yyyy-mm-ddT_hh:mm:ss.sss':'time'})
                 .set_index('time')
                 )
            d = d[['temperature','depth']]
        return d

    def trim(self, t0=None, t1=None, d=None, inplace=False):
        ''' select data between t0 and t1 '''
        if inplace:
            if any([t0, t1]):
                self.d = self.d[t0:t1]
            elif 'deployment' in str(type(d)):
                self.d = self.d[d.start.time:d.end.time]
        else:
            out = copy.deepcopy(self)
            if any([t0, t1]):
                out.d = out.d[t0:t1]
            elif 'deployment' in str(type(d)):
                out.d = out.d[d.start.time:d.end.time]
            return out

    #
    def to_pickle(self, file=None):
        dictout = {key: getattr(self,key) for key in rbr_attrs}
        if file is None:
            file = self._file_striped+'.p'
        if '.p' not in file:
            file = file+'.p'
        pickle.dump( dictout, open( file, 'wb' ) )
        print('Data store to '+file)

    def _read_pickle(self, file):
        p = pickle.load( open( file, 'rb' ) )
        for key in rbr_attrs:
            setattr(self, key, p[key])

    #
    def to_nc(self, file, **kwargs):
        ds = self.d.to_xarray()
        for key in rbr_attrs:
            if key != 'd':
                ds.attrs[key] = getattr(self,key)
        ds.to_netcdf(file, **kwargs)
        print('Data store to '+file)

    def _read_nc(self, file):
        ds = xr.open_dataset(file)
        self.d = ds.to_dataframe()
        for key in rbr_attrs:
            if key != 'd':
                setattr(self, key, getattr(ds,key))

    #
    def plot(self, figsize=(15,10), **kwargs):
        return self.d.plot(subplots=True,
                           figsize=figsize,
                           **kwargs
                           )

    def plot_bk(self, rule=None):

        output_notebook()
        TOOLS = 'pan,wheel_zoom,box_zoom,reset,help'

        # line specs
        lw = 3
        c = 'black'

        # subsample and compute speed
        if rule is not None:
            # 1T = 1 minute, 30S = 30 seconds
            d = self.d.resample(rule).mean()
        else:
            d = self.d

        # create a new plot and add a renderer
        s1 = figure(tools=TOOLS,
                    plot_width=600, plot_height=300,
                    title='depth',
                    x_axis_type='datetime')
        s1.line('time', 'depth', source=d, line_width=lw, color=c)
        s1.add_tools(HoverTool(
            tooltips=[('Time','@time{%F %T}'),('depth','@{depth}{0.00f}'),],
            formatters={'@time': 'datetime','@depth' : 'printf',},
            mode='vline'
            ))
        #_add_start_end(s1, d['lon'])
        #
        s2 = figure(tools=TOOLS,
                    plot_width=600, plot_height=300,
                    title='temperature',
                    x_axis_type='datetime',
                    x_range=s1.x_range
                    )
        s2.line('time', 'temperature', source=d, line_width=lw, color=c)
        s2.add_tools(HoverTool(
            tooltips=[('Time','@time{%F %T}'),('temperature','@{temperature}{0.00f}'),],
            formatters={'@time': 'datetime','@temperature' : 'printf',},
            mode='vline'
            ))
        #_add_start_end(s2, d['lat'])

        p = gridplot([[s1], [s2]])
        show(p)
