# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 10:13:48 2018

@author: mhamon
"""

import os, sys
from glob import glob
import struct
import datetime as dt
#from datetime import datetime
import time
import imaplib
import email

import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

import warnings
warnings.filterwarnings('ignore')

# cognac data and tools
from .gps import *

def Interrogation_Planeurs(mailbox='INBOX', sbd_dir='../data/iridium'):
    #
    try:
        server = imaplib.IMAP4_SSL('domicile.ifremer.fr', 993)
        print('connected to server')
    except imaplib.IMAP4.error:
        print('pb login serveur')
        return
    server.login('planeurs_mp1', 'viadomicile29,caplane4')
    #print(server.list())
    server.select(mailbox=mailbox, readonly=True)
    typ, data = server.search(None, '(FROM "sbdservice@sbd.iridium.com")')
    #print(data)
    for num in data[0].split():
        typ, data = server.fetch(num, '(RFC822)')
        email_body = data[0][1].decode('utf-8') #récupère le contenu du mail
        #print(email_body)
        mail = email.message_from_string(email_body)
        #Vérifie si il y a une pièce jointe
        if mail.get_content_maintype() != 'multipart':
            continue
        IMEI = mail["Subject"].split(": ")[1] #On récupère le numéro IMEI depuis le sujet
        #Détermine le dossier de destination des pièces jointes
        #  detach_dir = "N:\GESTION_LOCALE_ISO9001\I_Instrumentation\Balise_Iridium" + "\\" + IMEI
        #print(IMEI)
        detach_dir = os.path.join(sbd_dir, IMEI)
        #print(detach_dir)
        if not os.path.exists(detach_dir) :
            os.mkdir(detach_dir)
        #récupère la pièce jointe
        # we use walk to create a generator so we can iterate on the parts and forget about the recursive headache
        for part in mail.walk():
            # multipart are just containers, so we skip them
            if part.get_content_maintype() == 'multipart':
                continue
            # is this part an attachment ?
            if part.get_content_type() == "text/plain" :
                continue
            filename = part.get_filename()
            # if there is no filename, skip the part
            if not filename:
               continue
            att_path = os.path.join(detach_dir, filename)
            #print(att_path)
            #Si le fichier existe déjà, passe au suivant
            if not os.path.isfile(att_path) :
                # finally write the stuff
                print('write '+att_path)
                fp = open(att_path, 'wb')
                fp.write(part.get_payload(decode=True))
                fp.close()

    #Ferme la connexion
    server.close()
    server.logout()
    print('data updated from server')



def analyse_sbd(IMEI, sbd_dir='../data/iridium'):
    #
    # init gps container
    gp = gps()
    #print('search sbd files in '+sbd_dir)
    if not (os.path.exists(sbd_dir)):
        return
    files = sorted(glob(os.path.join(sbd_dir,'%d/*.sbd'%IMEI)))
    #print(files)
    #
    Liste_points=[]
    f_synthese=open(os.path.join(sbd_dir,'synthese_'+str(IMEI)+'.txt'),'w')
    for fichier in files:
        f=open(fichier,"rb")
        Data_SBD=f.read()
        f.close()
        if len(Data_SBD)<3:
            continue
        s=struct.Struct(">BBB")
        Header_Data=s.unpack(Data_SBD[0:3])
        CMD_TYPE=Header_Data[0]
        CMD_SUB_TYPE=Header_Data[1]
        SEQ_NUM=Header_Data[2]
        if CMD_TYPE==0x01 and CMD_SUB_TYPE==0xfd and SEQ_NUM==0x00:
            Nb_T_POSITION_REPORT=int((len(Data_SBD)-3)/16)
            for i in range(Nb_T_POSITION_REPORT):
                s=struct.Struct(">IBIIBBB")
                T_POSITION_REPORT_DATA=s.unpack(Data_SBD[i*16+3:i*16+19])
                GPS_TIME=T_POSITION_REPORT_DATA[0] \
                         +time.mktime(dt.datetime.strptime("06/01/1980 00:00:00",
                                "%d/%m/%Y %H:%M:%S").timetuple())
                RECORD_TYPE=(T_POSITION_REPORT_DATA[1]&0x80)>>7
                RESERVED=(T_POSITION_REPORT_DATA[1]&0x40)>>6
                NUMBER_SATS_VISIBLE=((T_POSITION_REPORT_DATA[1]&0x38)>>3) +3
                HDOP=(T_POSITION_REPORT_DATA[1]&0x07)*0.5
                GPS_LATITUDE=T_POSITION_REPORT_DATA[2]*0.000001-90
                GPS_LONGITUDE=T_POSITION_REPORT_DATA[3]*0.000001-180
                SPEED=T_POSITION_REPORT_DATA[4]
                GPS_FIX_ACCURACY=T_POSITION_REPORT_DATA[4]
                TIME_TO_FIX=T_POSITION_REPORT_DATA[5]
                if RECORD_TYPE==1:
                    ti = dt.datetime.fromtimestamp(GPS_TIME)
                    #Liste_points.append((GPS_LONGITUDE, GPS_LATITUDE, time))
                    gp.add(GPS_LONGITUDE, GPS_LATITUDE, ti)
                    f_synthese.write( ti.strftime('%d/%m/%Y %H:%M:%S')+ "\t"
                                      +str(GPS_LATITUDE)
                                      + "\t"+str(GPS_LONGITUDE)+'\r\n')
    f_synthese.close()
    #
    #if Liste_points:
    #    x = [t[0] for t in Liste_points]
    #    y = [t[1] for t in Liste_points]
    #    t = [t[2] for t in Liste_points]
    #else:
    #    x, y, t = None, None, None

    #return x, y, t
    return gp

def lstr(l):
    return '%d deg %.5f' %(int(l), (l-int(l))*60.)

def plot_map(fig=None, coast_resolution='10m', figsize=(10, 10)):
    #
    if fig is None:
        fig = plt.figure(figsize=figsize)
    else:
        fig.clf()
    ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
    ax.set_extent(ll_lim, crs=ccrs.PlateCarree())
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='k',
                      alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    ax.coastlines(resolution=coast_resolution, color='k')

    return fig, ax


def plot_bathy(ax):
    ### GEBCO bathymetry
    dpath = '/Users/aponte/Current_projects/cognac/campagnes_techno/cognac_pilote/manip_europe/bathy/RN-1994_1473924981206'
    ds = xr.open_dataset(dpath+'/GEBCO_2014_2D_5.625_42.0419_8.8046_44.2142.nc')
    ds = ds.sel(lon=slice(ll_lim[0], ll_lim[1]), lat=slice(ll_lim[2], ll_lim[3]))
    cs = ds['elevation'].plot.contour(ax=ax, levels=[-2000., -1000., -100., -50.], linestyles='-',
                                 colors='black', linewidths=0.5, transform=ccrs.PlateCarree())
    plt.clabel(cs, cs.levels, inline=True, fontsize=10)




#
ll_lim = [6.4, 6.6, 42.92, 43.2]

def main():

    fig=None

    while True:


        for maildir in ['INBOX', 'BALISES_NOVATECH']:
            Interrogation_Planeurs(mailbox=maildir)

        fig, ax = plot_map(fig=fig)

        IMEIS = [300434060873240, 300434062298540, 300434060657120]
        instrums = ['source', 'recepteur', 'vmp']
        h = []
        for instrum, imei in zip(instrums, IMEIS):
            x, y, t = Analyse_SBD(imei)
            if x is not None:
                #print( ' %d : x = %.3f , y = %.3f, z = %s' %(imei, x[-1], y[-1], \
                #print(' %d : lon = %s , lat = %s, t = %s' %(imei, \
                print(' %s : lon = %s , lat = %s, t = %s' % (instrum, \
                                                             lstr(x[-1]), lstr(y[-1]), \
                                                             t[-1].strftime("%d/%m/%Y %H:%M:%S")) )
                #
                #label = '%d %s %s' %(imei, lstr(x[-1]), lstr(y[-1]))
                label = '%s %s %s' %(instrum, lstr(x[-1]), lstr(y[-1]))
                ax.plot(x,y,'-', transform=ccrs.PlateCarree())
                lh = ax.plot(x[-1], y[-1], 'o', transform=ccrs.PlateCarree(), label=label)
                h.append(lh)
        plt.legend()
        plot_bathy(ax)
        plt.pause(10.)
        print('\n')
        time.sleep(10.)

    plt.show(block=False)


if __name__ == "__main__":
    main()
