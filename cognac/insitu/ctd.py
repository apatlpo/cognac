# !!! This should be a more general routine for ctd type data
#
# ------------------------- RBR data -----------------------------------
#

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

from bokeh.io import output_notebook, show
from bokeh.layouts import gridplot
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.plotting import figure

import pickle

ctd_attrs = ['d', 'file', 'start_date', 'dt']

class ctd(object):
    """ Contains CTD cast data
    """
    def __init__(self, file, **kwargs):
        if file is not None:
            self.file = file
            if '.cnv' in file:
                self._file_striped = self.file.rstrip('.cnv')
                self.d = self._read_cnv(file, **kwargs)
                self._read_cnv_header(file)
            elif '.p' in file:
                self._file_striped = self.file.rstrip('.p')
                self._read_p(file)
        else:
            print('You need to provide a file name')

    def __str__(self):
        return self.d.__str__()

    def __getitem__(self, item):
        if item is 'time':
            # time or depth?
            return self.d.index
        else:
            return getattr(self.d, item)

    # IO
    def _read_cnv(self, file, **kwargs):
        d = pd.read_table(file,
                             names=['sample','pressure','temperature',
                                    'salinity','conductivity','flag'],
                             comment='#', delim_whitespace=True,
                             index_col=0, **kwargs)
        d = d[~d.index.str.contains("\*")]
        dtypes = {'pressure': np.float64, 'temperature': np.float64,
                  'salinity': np.float64, 'conductivity': np.float64,
                  'flag': np.float64}
        for key, value in dtypes.items():
            setattr(d, key, getattr(d, key).astype(value))
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
            self.d['time'] = self.start_date \
                + pd.to_timedelta(self.d.index * self.dt, unit='s')
            self.d.set_index('time', inplace=True)

    #
    def to_pickle(self, file=None):
        dictout = {key: getattr(self,key) for key in ctd_attrs}
        if file is None:
            file = self._file_striped+'.p'
        pickle.dump( dictout, open( file, 'wb' ) )
        print('Data store to '+file)

    def _read_p(self, file):
        p = pickle.load( open( file, 'rb' ) )
        for key in ctd_attrs:
            setattr(self, key, p[key])

    #
    def resample(self, *args, **kwargs):
        # should add option to compute in place or not
        self.d = self.d.resample(*args, **kwargs).mean()
        self._update_time_info()

    def _update_time_info(self):
        self.dt = (self.d.index[1]-self.d.index[0]).total_seconds()
        self.start_date = self.d.index[0]

    def clean_and_depthbin(self, dp=1, plot=False):
        ''' select descent and bin by depth
        '''
        # compute speed of descent
        dpdt = self.d.pressure.diff()/self.dt
        threshold = .2 # dbar/s
        self.d = self.d[dpdt>threshold]
        self._update_time_info()
        #
        p = np.round(self.d.pressure/dp)
        self.d = self.d.groupby(by=p).mean()
        #
        if plot:
            dpdt.plot()
            dpdt[dpdt>threshold].plot(color='orange')

    #
    def plot(self, **kwargs):
        #self.d.plot(**kwargs)
        Nx = 4
        Ny = int(np.ceil(((self.d).shape[1]-1)/Nx))
        i=1
        plt.figure(figsize=(15,7))
        for name, series in self.d.iteritems():
            if not name in ['pressure', 'z']:
                ax = plt.subplot(Ny, Nx, i)
                ax.plot(series, -self.d.index)
                ax.set_xlabel(name)
                i+=1
        plt.show()

    def plot_bk(self):

        output_notebook()
        TOOLS = 'pan,wheel_zoom,box_zoom,reset,help'

        _d = self.d
        _d['z'] = -_d.index # approx

        # create a new plot and add a renderer
        s1 = figure(tools=TOOLS, plot_width=300, plot_height=300, title=None)
        s1.line('temperature', 'z', source=_d)
        s1.add_tools(HoverTool(
            tooltips=[('z','@z{%0.1f}'),('temperature','@{temperature}{%0.4f}'),],
            formatters={'z': 'printf','temperature' : 'printf',},
            mode='hline'
            ))

        s2 = figure(tools=TOOLS, plot_width=300, plot_height=300, title=None,
                    y_range=s1.y_range)
        s2.line('salinity', 'z', source=_d)
        s2.add_tools(HoverTool(
            tooltips=[('z','@z{%0.1f}'),('salinity','@{salinity}{%0.4f}'),],
            formatters={'z': 'printf','salinity' : 'printf',},
            mode='hline'
            ))

        p = gridplot([[s1, s2]])
        show(p)
