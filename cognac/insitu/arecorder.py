
# timeline -> files

import os
#import csv
import pandas as pd

arec_attrs = ['map', 'path']

#from ..acoustic.thinkdsp import *

class acoustic_recorder(object):

    def __init__(self, path, recorder):
        self.path = path
        self.recorder = recorder
        if recorder == 'logger_head':
            self._load_logger_head(path)

    def _load_logger_head(self, path):
        # load CSV file
        skip = lambda x: 'filename' in x
        df = pd.read_csv(path+'LOG.CSV',
                         names=['file_name','ID','gain', 'Voltage', 'Version'],
                         comment='f')
        df = df.set_index('file_name')
        # list files actually available:
        dfiles = pd.DataFrame()
        for dirname, dirnames, filenames in os.walk(path):
            if dirname != path:
                dfiles = dfiles.append(pd.DataFrame({'file_path': [os.path.join(dirname, f) for f in filenames if '.wav' in f],
                                                     'file_name': [f for f in filenames if '.wav' in f]}), ignore_index=True)
        dfiles = dfiles.set_index('file_name')
        #
        df = pd.concat([dfiles, df], join='inner', axis=1)
        time = pd.to_datetime(df.reset_index()['file_name'],
                              format='%Y%m%dT%H%M%S.wav').rename('time').to_frame()
        df = pd.concat([df.reset_index(), time], axis=1).set_index('time')

        # build timeline and file mapping
        self.df = df.sort_index()
        #

    def map(self, t, dt=None):
        _df = self.df.reset_index()
        tp = _df.shift(1)
        if type(t) is list:
            if t[0] < _df['time'].iloc[1]:
                out = _df.loc[ (_df['time']<=t[1]) ]
            else:
                tm = _df.shift(-1)
                out = _df.loc[ (tm['time']>=t[0]) & (_df['time']<=t[1]) ]
            return out.set_index('time')
        elif dt is not None:
            return self.map([t-pd.Timedelta(dt)/2,t+pd.Timedelta(dt)/2])
