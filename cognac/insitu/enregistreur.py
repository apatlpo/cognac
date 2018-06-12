
import os, sys, pickle, glob
import csv
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import copy

import datetime

import pandas as pd
import pynmea2

# cognac data and tools
from data import *
from utils import *


# flag controlling the production or not of figures
prtfig=True



#
# ------------------------- File in/out -----------------------------------
#



def read_logger_file(file):

    # init gps container
    gps = gps_data()

    if isinstance(file, list):

        for f in file:
            gps = gps + read_logger_file(f)

    else:

        print('Reads ' + file)

        if True:

            #np.genfromtxt(file, usecols=(), delimiter=',')

            gpsfile = pynmea2.NMEAFile(file)

            data = pd.DataFrame()
            for d in gpsfile:
                time = datetime.datetime.combine(d.datestamp, d.timestamp)
                data = data.append({'lon': d.longitude, 'lat': d.latitude, 'time': time}, ignore_index=True)
            #
            gps.d = data.set_index('time')

        else:

            # open and read data
            fptr = open(file, 'r')
            lines = fptr.readlines()
            fptr.close()

            for i, line in enumerate(lines):
                line = line.replace('\n', '').replace('\r', '')
                if "$GPRMC" in line:
                    if line.split(',')[2] == 'A':
                        # coordinates available
                        #print line.split(',')[1]
                        time = float(line.split(',')[1])
                        h = int(time / 10000)
                        m = int( (time - h * 10000) / 100 )
                        s =  time - h * 10000 - m * 100
                        #
                        #print line.split(',')[9]
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
                        #time = date2num(datetime.datetime(year, month, day, h, m, int(s)))
                        time = datetime.datetime(year, month, day, h, m, int(s))
                        gps.add(lon, lat, time)

    return gps


#
# ------------------------- Main code -----------------------------------
#


if __name__ == "__main__":

    #
    # -------- read gps logger
    #

    # data directory and filename
    #
    data_dir = '../data/enregistreur_gps/'

    gps_files = sorted(glob.glob(data_dir+'*.DAT')) #[:2]
    gps_all = read_logger_file(gps_files)

    #
    # clean gps data with deployment log
    gps = []
    for i,lgps in enumerate(rec):
        gps.append(copy.deepcopy(gps_all))
        gps[i].clean(rec[i])

    # plot map
    if prtfig:
        map = plot_map()
        for i,lgps in enumerate(gps):
            lgps.plot(map, label='d'+str(i+1), color=c_enregistreur)

        #plt.ion();
        plt.show()
        #plt.savefig('figs/enregistreur_gps_logger.pdf')

    sys.exit()


    #
    # plot lon/lat time series
    #
    # for time plotting purposes
    t0 = datetime.datetime(2016, 9, 2, 0)
    t1 = datetime.datetime(2016, 9, 4, 12)
    tck = [];
    tck_label = [];
    tdays = []
    t = t0
    while t < t1 + datetime.timedelta(hours=6):
        tck.append(t)
        if t.hour is 12:
            tck_label.append(t.strftime('%Y-%m-%d'))
        else:
            tck_label.append('')
        t += datetime.timedelta(hours=6)
        if t.hour is 0:
            tdays.append(date2num(t))

    if prtfig:
        fig = plt.figure(figsize=(10, 10))

        # lon
        ax = fig.add_subplot(211)

        for i, lgps in enumerate(gps_logger):
            ax.plot(lgps.time,lgps.lon, color=c_enregistreur)

        for i, lgps in enumerate(gps_spot):
            #ax.plot(lgps.time,lgps.lon, linestyle=':', color='black')
            ax.scatter(lgps.time, lgps.lon, 10, marker='.', color='black')

        for i, lgps in enumerate(gps_argos):
            ax.plot(lgps.time,lgps.lon, linestyle=' ', color='grey', marker='+')

        for i, ld in enumerate(rec):
            ax.plot(ld.start.time,ld.start.lon, color=c_enregistreur, marker='o')
            ax.plot(ld.end.time, ld.end.lon, color=c_enregistreur, marker='*')

        plt.xlim([t0, t1])
        ax.set_xticks(tck)
        ax.set_xticklabels(tck_label)
        [plt.axvline(x=td, color='k', lw=1) for td in tdays]
        # [plt.axvline(x=d.start.time, color=c_enregistreur, lw=1) for d in rec]
        # [plt.axvline(x=d.end.time, color=c_enregistreur, lw=1) for d in rec]
        plt.grid()

        # lat
        ax = fig.add_subplot(212)

        for i, lgps in enumerate(gps_logger):
            ax.plot(lgps.time, lgps.lat, color=c_enregistreur)

        for i, lgps in enumerate(gps_spot):
            #ax.plot(lgps.time, lgps.lat, linestyle=':', color='black')
            ax.scatter(lgps.time, lgps.lat, 10, marker='.', color='black')

        for i, lgps in enumerate(gps_argos):
            ax.plot(lgps.time, lgps.lat, linestyle=' ', color='grey', marker='+')

        for i, ld in enumerate(rec):
            ax.plot(ld.start.time, ld.start.lat, color=c_enregistreur, marker='o')
            ax.plot(ld.end.time, ld.end.lat, color=c_enregistreur, marker='*')

        plt.xlim([t0, t1])
        ax.set_xticks(tck)
        ax.set_xticklabels(tck_label)
        [plt.axvline(x=td, color='k', lw=1) for td in tdays]
        #[plt.axvline(x=d.start.time, color=c_enregistreur, lw=1) for d in rec]
        #[plt.axvline(x=d.end.time, color=c_enregistreur, lw=1) for d in rec]
        plt.grid()

        plt.savefig('figs/enregistreur_gps_lonlat.pdf')


    #
    # interpolate gps data onto a global timeline and interpolate with ship log
    #
    gps_logger_interp = interp_gps(time, gps_logger)
    gps_spot_interp = interp_gps(time, gps_spot)
    gps_argos_interp = interp_gps(time, gps_argos)
    # fills start and end of deployments with interpolation from log
    gps_spot_interp_filled = copy.deepcopy(gps_spot_interp)
    gps_spot_interp_filled.fill_with_deps(rec)

    fig = plt.figure(figsize=(10, 10))

    # lon
    ax = fig.add_subplot(211)

    for i, lgps in enumerate(gps_logger):
        ax.plot(lgps.time, lgps.lon, linewidth=2, color=c_enregistreur)

    for i, ld in enumerate(rec):
        ax.plot(ld.start.time, ld.start.lon, color=c_enregistreur, marker='o')
        ax.plot(ld.end.time, ld.end.lon, color=c_enregistreur, marker='*')

    for i, lgps in enumerate(gps_spot):
        ax.scatter(lgps.time, lgps.lon, 10, marker='.', color='grey')

    for i, lgps in enumerate(gps_argos):
        ax.plot(lgps.time, lgps.lon, linestyle=' ', color='grey', marker='+')

    ax.plot(gps_spot_interp_filled.time, gps_spot_interp_filled.lon, linestyle='-', color='black')


    plt.xlim([t0, t1])
    ax.set_xticks(tck)
    ax.set_xticklabels(tck_label)
    [plt.axvline(x=td, color='k', lw=1) for td in tdays]
    plt.grid()

    # lat
    ax = fig.add_subplot(212)

    for i, lgps in enumerate(gps_logger):
        ax.plot(lgps.time, lgps.lat, linewidth=2, color=c_enregistreur)

    for i, ld in enumerate(rec):
        ax.plot(ld.start.time, ld.start.lat, color=c_enregistreur, marker='o')
        ax.plot(ld.end.time, ld.end.lat, color=c_enregistreur, marker='*')

    for i, lgps in enumerate(gps_spot):
        ax.scatter(lgps.time, lgps.lat, 10, marker='.', color='grey')

    for i, lgps in enumerate(gps_argos):
        ax.plot(lgps.time, lgps.lat, linestyle=' ', color='grey', marker='+')

    ax.plot(gps_spot_interp_filled.time, gps_spot_interp_filled.lat, linestyle='-', color='black')

    plt.xlim([t0, t1])
    ax.set_xticks(tck)
    ax.set_xticklabels(tck_label)
    [plt.axvline(x=td, color='k', lw=1) for td in tdays]
    plt.grid()

    plt.savefig('figs/enregistreur_gps_lonlat_interp.pdf')


    #
    # Distance between interpolated trajectory and actual data
    #
    gps = gps_spot_interp_filled
    distance_logger = get_distance(gps_logger_interp.lon , gps_logger_interp.lat ,
                                                  gps.lon , gps.lat)

    distance_argos = get_distance(gps_argos_interp.lon, gps_argos_interp.lat,
                                                  gps.lon, gps.lat)

    fig = plt.figure(figsize=(10, 10))

    ax = fig.add_subplot(211)
    ax.plot(gps.time, distance_logger, color='k', linestyle='-', linewidth=2)
    ax.set_title('Distance between final track and logger')
    plt.xlim([t0, t1])
    ax.set_ylabel('[m]')
    ax.set_xticks(tck)
    ax.set_xticklabels(tck_label)
    [plt.axvline(x=td, color='k', lw=1) for td in tdays]
    plt.grid()

    ax = fig.add_subplot(212)
    ax.plot(gps.time, distance_argos, color='k', linestyle='-', linewidth=2)
    ax.set_title('Distance between final track and argos positions')
    plt.xlim([t0, t1])
    ax.set_ylabel('[m]')
    ax.set_xticks(tck)
    ax.set_xticklabels(tck_label)
    [plt.axvline(x=td, color='k', lw=1) for td in tdays]
    plt.grid()


    plt.savefig('figs/enregistreur_gps_lonlat_interp_dist.pdf')


    #
    # store final gps track of enregistreur
    #
    pickle.dump(gps, open('pydata/enregistreur.gps', 'wb'))


    sys.exit()