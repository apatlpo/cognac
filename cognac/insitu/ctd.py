# !!! This should be a more general routine for ctd type data
#
# ------------------------- ctd data -----------------------------------
#

import pandas as pd
import numpy as np

import gsw

import matplotlib.pyplot as plt

from bokeh.io import output_notebook, show
from bokeh.layouts import gridplot
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.plotting import figure

import pickle

ctd_attrs = ['d_time', 'd_depth', 'file', 'start_date', 'dt']

class ctd(object):
    """ Contains CTD cast data
    """
    def __init__(self, file, **kwargs):
        if file is not None:
            self.file = file
            if '.cnv' in file:
                self._file_striped = self.file.rstrip('.cnv')
                self.d_time = self._read_cnv(file, **kwargs)
                self._read_cnv_header(file)
            elif '.p' in file:
                self._file_striped = self.file.rstrip('.p')
                self._read_pickle(file)
        else:
            print('You need to provide a file name')

    def __str__(self):
        if hasattr(self, 'd_depth'):
            return self.d_depth.__str__()
        else:
            return self.d_time.__str__()

    def __getitem__(self, item):
        if item=='time':
            # time or depth?
            return self.d_time.index
        if hasattr(self, 'd_depth'):
            _d = self.d_depth
        else:
            _d = self.d_time
        return _d[item]

    # IO
    def _read_cnv(self, file, **kwargs):
        if 'columns' in kwargs:
            columns = kwargs['columns']
        else:
            columns = ['sample','pressure','temperature',
                       'salinity','conductivity','flag'
                       ]
        #
        d = pd.read_table(file,
                          names=columns,
                          comment='#', delim_whitespace=True,
                          index_col=0,
                          )
        d = d[~d.index.str.contains("\*")]
        for c in columns:
            if c!='sample':
                d[c] = d[c].astype(np.float64)
        #dtypes = {c: np.float64 if }
        #dtypes = {'pressure': np.float64, 'temperature': np.float64,
        #          'salinity': np.float64, 'conductivity': np.float64,
        #          'flag': np.float64}
        #for key, value in dtypes.items():
        #        setattr(d, key, getattr(d, key).astype(value))
        d.index = d.index.astype(int)
        #encoding='iso-8859-1' # decimal='.',
        # dtype={'temp': np.float64} is ignored if passed to read_table
        #d.pressure = d.pressure.astype(float)
        #d.temperature = d.temperature.astype(float)
        #d.index.name='time'
        return d

    def _read_cnv_header(self, file):
        with open(file) as f:
            for line in f:
                if line.strip()[0] in ['#','*']:
                    # extract useful information
                    if 'start_time' in line:
                        date = line.split('=')[1].split('[')[0].strip()
                        self.start_date = pd.datetime.strptime(date, '%b %d %Y %H:%M:%S')
                    elif 'interval' in line:
                        self.dt = float(line.split(':')[1].strip())
        # compute time line
        if hasattr(self,'start_date') and hasattr(self,'dt'):
            self.d_time['time'] = self.start_date \
                + pd.to_timedelta(self.d_time.index * self.dt, unit='s')
            self.d_time.set_index('time', inplace=True)

    #
    def to_pickle(self, file):
        dictout = {key: getattr(self,key) for key in ctd_attrs}
        pickle.dump( dictout, open( file, 'wb' ) )
        print('Data store to '+file)

    def _read_pickle(self, file):
        p = pickle.load( open( file, 'rb' ) )
        for key in ctd_attrs:
            setattr(self, key, p[key])

    # cleaning and resampling
    def resample(self, rule, inplace=True, **kwargs):
        # should add option to compute in place or not
        _d = self.d_time.resample(rule=rule, **kwargs).mean()
        if inplace:
            self.d_time = _d
            self._update_time_info()
        else:
            return _d

    def _update_time_info(self):
        self.dt = (self.d_time.index[1]-self.d_time.index[0]).total_seconds()
        self.start_date = self.d_time.index[0]

    def clean_and_depth_bin(self,
                            bin_size=1,
                            ):
        ''' select descent and bin by depth
        '''
        # compute speed of descent
        _d = self.d_time
        dpdt = _d.pressure.diff()/self.dt
        threshold = .2 # dbar/s
        self.d_time = _d[dpdt>threshold]
        self._update_time_info()
        #
        p = np.round(_d.pressure/bin_size)
        self.d_depth = _d.groupby(by=p).mean()
        del self.d_depth['pressure']
        self.d_depth = self.d_depth[self.d_depth.index>=0]
        #
        self._update_eos()

    # EOS
    def _update_eos(self):
        assert hasattr(self, 'd_depth'), \
                'You need to run clean_and_depth_bin first'
        d = self.d_depth
        d['SA'] = gsw.SA_from_SP(d.salinity, d.index, d.longitude, d.latitude)
        d['CT'] = gsw.CT_from_t(d.SA, d.temperature, d.index)
        d['sound_speed'] = gsw.sound_speed(d.SA, d.CT, d.index)

    # plotting
    def plot_tseries(self):
        d = self.d_time[['pressure','temperature','salinity']]
        d.plot(lw=1,
               subplots=True,
               grid=True,
               layout=(2,3),
               figsize=(10,10)
               )

    def plot_depth(self,
                   variables=['temperature','salinity', 'sound_speed'],
                   figsize=(15,7),
                   grid=True,
                   **kwargs
                   ):

        assert hasattr(self, 'd_depth'), \
                'You need to run clean_and_depth_bin first'

        d = self.d_depth[variables]

        Nx = 4
        Ny = int(np.ceil(((d).shape[1]-1)/Nx))
        i=1
        plt.figure(figsize=(20,7))
        for name, series in d.iteritems():
            if not name in ['pressure', 'z']:
                ax = plt.subplot(Ny, Nx, i)
                ax.plot(series, -d.index, **kwargs)
                ax.set_xlabel(name)
                if grid:
                    ax.grid()
                i+=1

    def plot_bk(self):

        output_notebook()
        TOOLS = 'pan,wheel_zoom,box_zoom,reset,help'

        _d = self.d_depth
        _d['z'] = -_d.index # approx

        # create a new plot and add a renderer
        s1 = figure(tools=TOOLS, plot_width=300, plot_height=300, title=None)
        s1.line('temperature', 'z', source=_d)
        s1.add_tools(HoverTool(
            tooltips=[('z','@z{0.0f}'),('temperature','@{temperature}{0.0000f}'),],
            formatters={'@z': 'printf','@temperature' : 'printf',},
            mode='hline'
            ))

        s2 = figure(tools=TOOLS, plot_width=300, plot_height=300, title=None,
                    y_range=s1.y_range)
        s2.line('salinity', 'z', source=_d)
        s2.add_tools(HoverTool(
            tooltips=[('z','@z{0.0f}'),('salinity','@{salinity}{0.0000f}'),],
            formatters={'@z': 'printf','@salinity' : 'printf',},
            mode='hline'
            ))

        p = gridplot([[s1, s2]])
        show(p)
