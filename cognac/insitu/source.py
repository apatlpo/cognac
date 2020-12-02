
from glob import glob
import re
import numpy as np
import pandas as pd

# cognac data and tools
from .gps import *
from .arecorder import *

source_attrs = ['label', 'gps', 'emission']

def load_emission_sequence(path):
    files = sorted(glob(path+'son*.wav'),
                   key=lambda x: int(re.match('\D*(\d+)', x.split('/')[-1]).group(1)))
    sequence = [Signal.from_wav(f) for f in files]
    #if len(sequence)==
    fs = set([s.fs for s in sequence])
    return sequence, fs

def load_source_files(path,
                      label='source',
                      ignore=[],
                      verbose=-1):
    if path:
        files = sorted([os.path.join(w[0],f) for w in os.walk(path)
                        for f in w[2] if ('.txt' in f) and ('mission_' in f)]
                       )
        files = [f for f in files if f not in ignore]
    sdata = source_rtsys(label=label)
    for f in files:
        sdata += source_rtsys(f, verbose=verbose)
    return sdata

class source_rtsys(object):
    ''' Data container for rtsys acoustical source log
    '''
    def __init__(self,
                 file=None,
                 clean=True,
                 label='source',
                 verbose=-1):
        self.label = label
        #
        if file is not None:
            if file.split('.')[-1]=='p':
                # pickle file
                self._read_pickle(file)
                return
            elif file.split('.')[-1]=='nc':
                # netcdf file
                self._read_nc(file)
                return
            else:
                self.gps, self.emission = read_log_file(file, label, verbose)
        else:
            return
        if clean:
            self.clean()

    def __bool__(self):
        if hasattr(self,'gps') and hasattr(self,'emission'):
            if self.gps and self.emission:
                return True
        return False

    def __add__(self, other):
        if hasattr(self, 'gps') and hasattr(self, 'emission'):
            self.gps.d = self.gps.d.append(other.gps.d)
            self.emission.d = self.emission.d.append(other.emission.d)
        else:
            self.gps = other.gps
            self.emission = other.emission
        return self

    def trim(self, t0, t1, inplace=True):
        ''' select data between t0 and t1 '''
        if inplace:
            self.gps.trim(t0,t1)
            self.emission.trim(t0,t1)
        else:
            selfc = copy.deepcopy(self)
            selfc.trim(t0,t1)
            return selfc

    def clean_deployment(self, d, inplace=True):
        '''Use deployment to clean gps data'''
        if inplace:
            self.trim(d.start.time, d.end.time)
        else:
            return self.trim(d.start.time, d.end.time, inplace=False)

    def clean(self):
        """ wrapper around sort, drop_duplicates, compute_velocity
        """
        self.sort()
        self.drop_duplicates()
        self.gps.compute_velocity()

    def sort(self):
        self.gps.sort()
        self.emission.sort()

    def drop_duplicates(self):
        if not self:
            return
        d = (self.gps.d.reset_index().drop_duplicates(subset='time')
                 .set_index('time'))
        self.gps.d = d
        d = (self.emission.d.reset_index().drop_duplicates(subset='time')
                 .set_index('time'))
        self.emission.d = d

    #
    def to_pickle(self, file):
        dictout = {key: getattr(self,key) for key in source_attrs}
        pickle.dump( dictout, open( file, 'wb' ) )
        print('Data store to '+file)

    def _read_pickle(self, file):
        p = pickle.load( open( file, 'rb' ) )
        for key in source_attrs:
            setattr(self, key, p[key])

    #
    def to_nc(self, file):
        self.gps.to_nc(file.rstrip('.nc')+'_gps.nc')
        self.emission.to_nc(file.rstrip('.nc')+'_emission.nc')

    def _read_nc(self, file):
        self.gps = gps(file=file.rstrip('.nc')+'_gps.nc')
        self.emission = gps(file=file.rstrip('.nc')+'_emission.nc')
        self.label = self.gps.label

class emissions(gps):
    ''' Data container for emission data
    '''

    def add(self, time, sound, lon, lat, sort=False):
        if not isinstance(lon,list): lon=[lon]
        if not isinstance(lat,list): lat=[lat]
        if not isinstance(time,list): time=[time]
        if not isinstance(sound,list): sound=[sound]
        d = pd.DataFrame({'lon': lon, 'lat': lat, 'sound': sound},
                         index = time)
        d.index.rename('time', inplace=True)
        self.d = self.d.append(d)
        if sort:
            self.d = self.d.sortby('time')

def read_log_file(file, label, verbose):

    print('Reads '+file)

    # init gps container
    gp = gps(label=label)

    # open and read data
    fptr = open(file, 'r')
    lines = fptr.readlines()
    fptr.close()

    # init variables
    gps_sync_start = -1
    gps_sync_stop = -1
    pps_delta_sync = -1
    current_idx = -1
    sync_timer = 0
    cycle_time = -1
    max_idx = -1
    idx_son=-1
    #t_nan = datetime.datetime(1, 1, 1, 0, 0, 0)
    t_nan = pd.NaT
    pps_sync_date = t_nan

    # init final arrays
    edata = emissions()
    e_time=[]; e_sound=[]; e_lon=[]; e_lat=[]

    for i, line in enumerate(lines):
        line = line.replace('\n', '').replace('\r', '')
        if "PPS :: Configured next TX" in line:
            # (Occurs once early in the file)
            #if verbose: print Fore.WHITE, i, line
            pass
        if "PPS :: sync" in line:
            # (doesn't mean there is valid position)
            date = line.split(' ')[3].split('T')[1].replace('Z', '').split(':')
            # date is the number of seconds since midnight
            date = int(date[0]) * 3600 + int(date[1]) * 60 + int(date[2])
            if pps_delta_sync == -1:
                # skips first sync
                #print Fore.CYAN, i, line
                pass
            else:
                # cycle_time = time between pps_sync, typiquement 3s
                ## print Fore.CYAN, i, line, " delta=%d"%(date-pps_delta_sync)
                if cycle_time == -1:
                    # use first cycle_time as reference
                    cycle_time = date - pps_delta_sync
                else:
                    if (date - pps_delta_sync != cycle_time and current_idx != max_idx):
                        # indicate when the time between synchro is the not initial one
                        #print Fore.LIGHTRED_EX, i, line, " delta=%d, idx:%d" % (date - pps_delta_sync, current_idx)
                        pass
                    else:
                        #print Fore.CYAN, i, line, " delta=%d" % (date - pps_delta_sync)
                        pass
                # sync_timer is the number of sync achieved since last emission
                sync_timer += 1
                if verbose>1: print('sync_timer='+str(sync_timer))
            pps_delta_sync = date
            # store date
            year = int(line.split(' ')[3].split('T')[0].split('-')[0])
            month = int(line.split(' ')[3].split('T')[0].split('-')[1])
            day = int(line.split(' ')[3].split('T')[0].split('-')[2])
            h = int(line.split(' ')[3].split('T')[1].replace('Z', '').split(':')[0])
            m = int(line.split(' ')[3].split('T')[1].replace('Z', '').split(':')[1])
            s = int(line.split(' ')[3].split('T')[1].replace('Z', '').split(':')[2])
            #pps_sync_date = date2num(datetime.datetime(year, month, day, h, m, s))
            #pps_sync_date = datetime.datetime(year, month, day, h, m, s)
            pps_sync_date = pd.Timestamp(year, month, day, h, m, s)
        if "PPS :: idle" in line:
            # debug: once sync is done, checks that sync is activated
            #print Fore.GREEN, i, line
            pass
        if "GPS :: $GNRMC" in line:
            # GPS data
            if verbose>0:
                print(line)
            if line.split(',')[2] == 'V' and gps_sync_start == -1:
                # no coordinates available
                if gps_sync_stop == -1 and gps_sync_start == -1:
                    #
                    time = float(line.split(',')[1])
                    h = int( time / 10000 )
                    m = int( (time - h * 10000) / 100 )
                    s = int( time - h * 10000 - m * 100 )
                    gps_sync_start = h * 3600 + m * 60 + s
                    if verbose>1: print(i, "(GPS sync start)", line)
                    #
                    if line.split(',')[9] is not '':
                        date = float(line.split(',')[9])
                        day = int( date/1e4 )
                        month = int( (date-day*1e4)/1e2 )
                        year = 2000+int( date-day*1e4-month*1e2 )
                    #
                if gps_sync_stop != -1:
                    if verbose>1: print(i, "(GPS sync lost!)", line)
            if line.split(',')[2] == 'A':
                # coordinates available
                time = float(line.split(',')[1])
                h = int(time / 10000)
                m = int( (time - h * 10000) / 100 )
                s = int( time - h * 10000 - m * 100 )
                #
                date = float(line.split(',')[9])
                day = int(date / 1e4)
                month = int((date - day * 1e4) / 1e2)
                year = 2000+int(date - day * 1e4 - month * 1e2)
                #
                lon = float(line.split(',')[5])
                lond = np.floor(lon / 100)
                lon = lond + (lon - lond * 100) / 60.
                #
                lat = float(line.split(',')[3])
                latd = np.floor(lat / 100)
                lat = latd + (lat - latd * 100) / 60.
                # store data
                #time = date2num(datetime.datetime(year, month, day, h, m, s))
                time = datetime.datetime(year, month, day, h, m, s)

                gp.add(lon, lat, time)
                #
                if gps_sync_stop == -1:
                    gps_sync_stop = h * 3600 + m * 60 + s
                    if verbose>1: print(i, "(GPS sync done)", line, "(Took %.2f seconds)" % ( \
                    gps_sync_stop - gps_sync_start))
                #
        if "DSP :: Transmission done" in line:
            # this is when the sound is produced
            # checks that the sound a synchronisation occurred since last emission
            if sync_timer != 1 and cycle_time != -1:
                if verbose>1: print(i, line,
                "(ERROR detected wrong cycle time: %ds)" % (cycle_time * sync_timer))
            sync_timer = 0
        if "WAV :: Reading" in line:
            # store file number
            idx_son = int( line.split('/')[-1].replace('son','').replace('.wav','') )
        if "WAV :: sent repondeur_idx" in line:
            # This piece of codes checks which sounds is emitted
            idx_data = line.split(' ')[-1].split('/')
            idx = int(idx_data[0])
            max_idx = int(idx_data[1])
            if current_idx == -1:
                # first emissions
                if idx != 0:
                    # first sound should be 0
                    if verbose>1: print(i, line, "(Error idx should be 0)")
            else:
                # following emissions
                good_value = (current_idx + 1) % (max_idx + 1)
                if idx != good_value and verbose>-1:
                    if verbose>1: print(i, line, "(Error idx should be %d)" % good_value)
            current_idx = idx
            # store sound and time
            e_time.append(pps_sync_date)
            e_sound.append(idx)
            if idx != idx_son and verbose>0:
                print(i, line, idx_son, 'Error repondeur and wav reading line do not match')

    # find coordinates corresponding to emission time
    if gp.d.size>0:
        for t in e_time:
            l = gp['lon'].loc[gp['time'] == t ].values
            if l.size>0:
                e_lon.append( l[0] )
            else:
                e_lon.append( np.NaN )
            #
            l = gp['lat'].loc[gp['time'] == t ].values
            if l.size>0:
                e_lat.append( l[0] )
            else:
                e_lat.append( np.NaN )

        # fill in emission data container
        edata.add(e_time, e_sound, e_lon, e_lat)
        edata.d = edata.d.loc[edata['time']!=t_nan]

    return gp, edata
