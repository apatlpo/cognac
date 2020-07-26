#
# ------------------------- inclino data -----------------------------------
#

import pandas as pd
import pickle
import copy

# containment and delegation
inclino_attrs = ['d', 'id',]

class inclino(object):
    """ Contains DST inclino data
    """
    def __init__(self, file, id, **kwargs):
        if file is not None:
            self.d = self._read(file, **kwargs)
        else:
            print('You need to provide a file name')
        self.id = id

    def __str__(self):
        return self.d.__str__()

    def __getitem__(self, item):
        if item=='time':
            # time or depth?
            return self.d.index
        else:
            return getattr(self.d, item)

    def _read(self, file, **kwargs):
        d = pd.read_table(file,
                             names=['sample','temperature','depth','tilt_x',
                                    'tilt_y','tilt_z','EAL'],
                             comment='#', sep='\t', decimal=',', index_col=1,
                             parse_dates=[1],
                             encoding='iso-8859-1', **kwargs)
        # dtype={'temp': np.float64} is ignored if passed to read_table
        d.temperature = d.temperature.astype(float)
        d.index.name='time'
        return d

    def trim(self, t0=None, t1=None, d=None, inplace=True):
        ''' select data between t0 and t1 '''
        if inplace:
            i = self
            if any([t0, t1]):
                self.d = self.d[t0:t1]
            elif 'deployment' in str(type(d)):
                self.d = self.d[d.start.time:d.end.time]
        else:
            i = copy.deepcopy(self)
            i.trim(t0=t0, t1=t1, d=d)
            return i

    #
    def to_pickle(self, file):
        dictout = {key: getattr(self,key) for key in inclino_attrs}
        pickle.dump( dictout, open( file, 'wb' ) )
        print('Data store to '+file)

    def _read_pickle(self, file):
        p = pickle.load( open( file, 'rb' ) )
        for key in ctd_attrs:
            setattr(self, key, p[key])
