#
# ------------------------- RBR data -----------------------------------
#

import pandas as pd
import pickle

# containment and delegation

rbr_attrs = ['d', 'file', 'dt']

class rbr(object):
    """ Contains DST inclino data
    """
    def __init__(self, file=None, **kwargs):
        if file is not None:
            self.file = file
            if '.txt' in file:
                self._file_striped = self.file.rstrip('.txt')
                self.d = self._read_txt(file, **kwargs)
                self.dt = (self.d.index[1]-self.d.index[0]).total_seconds()
            elif '.p' in file:
                self._file_striped = self.file.rstrip('.p')
                self._read_pickle(file)
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

    def _read_txt(self, file):
        d = pd.read_table(file,
                          names=['temperature','pressure','sea_pressure','depth'],
                          skiprows=[0], sep=',', decimal='.', index_col=0, parse_dates=[0])
        # dtype={'temp': np.float64} is ignored if passed to read_table
        d.temperature = d.temperature.astype(float)
        #d = d.set_index('time')
        d.index.name='time'
        return d

    def trim(self, t0=None, t1=None, d=None):
        ''' select data between t0 and t1 '''
        if any([t0, t1]):
            self.d = self.d[t0:t1]
        elif 'deployment' in str(type(d)):
            self.d = self.d[d.start.time:d.end.time]
        return self

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
