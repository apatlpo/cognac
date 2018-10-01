#
# ------------------------- inclino data -----------------------------------
#

import pandas as pd

# containment and delegation

class inclino(object):
    """ Contains DST inclino data
    """
    def __init__(self, file, **kwargs):
        if file is not None:
            self.d = self._read(file, **kwargs)
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

    def trim(self, t0=None, t1=None, d=None):
        ''' select data between t0 and t1 '''
        if any([t0, t1]):
            self.d = self.d[t0:t1]
        elif 'deployment' in str(type(d)):
            self.d = self.d[d.start.time:d.end.time]
        return self
