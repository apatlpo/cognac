
# timeline -> files

import os
#import csv
import numpy as np
import pandas as pd

arec_attrs = ['map', 'path']

from acoustics import Signal
#from ..acoustic.thinkdsp import *

def join(s1, s2):
    ''' join two signals
    '''
    assert s1.fs == s2.fs
    return Signal(np.concatenate([s1,s2]), fs=s1.fs)

class acoustic_recorder(object):

    def __init__(self, path, recorder):
        self.path = path
        self.recorder = recorder
        if recorder == 'logger_head':
            self._load_logger_head(path)

    def __getitem__(self,t):
        df = self._map(t).reset_index()
        s = None
        for f in df['file_path']:
            if s is None:
                s = Signal.from_wav(f)
            else:
                s = join(s, Signal.from_wav(f))
        # further triming
        i0 = int( (t.start - df['time'].iloc[0]).seconds*s.fs )
        i1 = int( (df['time'].iloc[-1] - t.stop).seconds*s.fs )
        s = Signal(s[i0:i1], s.fs)
        t = df['time'].iloc[0] + pd.Timedelta( i0/s.fs ,unit='seconds')
        return s, t

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

    def _map(self, t, dt=None):
        _df = self.df.reset_index()
        tp = _df.shift(1)
        if type(t) is slice:
            if t.start < _df['time'].iloc[1]:
                out = _df.loc[ (_df['time']<=t.stop) ]
            else:
                tm = _df.shift(-1)
                out = _df.loc[ (tm['time']>=t.start) & (_df['time']<=t.stop) ]
            return out.set_index('time')
        elif dt is not None:
            return self.map([t-pd.Timedelta(dt)/2,t+pd.Timedelta(dt)/2])
