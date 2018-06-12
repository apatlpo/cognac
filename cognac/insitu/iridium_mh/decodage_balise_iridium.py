# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 10:13:48 2018

@author: mhamon
"""

import struct
import datetime
import time
import os
import imaplib
import email
import sys
from PyQt4.QtGui import *
from PyQt4.QtCore import *
import  PyQt4.Qwt5     as Qwt

import ihm_balises_iridiumUI

class MainDialog(QDialog,ihm_balises_iridiumUI.Ui_Dialog):
    def __init__(self,parent=None):
        super(MainDialog,self).__init__(parent)
        self.setupUi(self)
        pixmap=QPixmap('logo_ifremer.png')
        self.label_logo_ifremer.setPixmap(pixmap)
        pixmap=QPixmap('logo_lops.png')
        self.label_logo_lops.setPixmap(pixmap)
        self.timer=QTimer()
        self.connect(self.timer, SIGNAL('timeout()'), self.on_timer)
        self.create_plot()
#        self.spinBox.setValue(1)
        periode=int(self.spinBox.value()*60000)
        self.timer.start(periode)

    def create_plot(self):
        self.qwtPlot.setCanvasBackground(Qt.white)
        self.qwtPlot.setAxisTitle(self.qwtPlot.xBottom, 'Longitude')
        self.qwtPlot.setAxisTitle(self.qwtPlot.yLeft, 'Latitude')       
        grid = Qwt.QwtPlotGrid()
        grid.attach(self.qwtPlot)
        self.qwtPlot.replot()
        self.curve = [None]*3
        pen = [ QPen(QColor('red')) ,QPen(QColor('black')), QPen(QColor('blue')) ]
        for i in range(3):
            self.curve[i] =  Qwt.QwtPlotCurve('')
            self.curve[i].setRenderHint(Qwt.QwtPlotItem.RenderAntialiased)
            pen[i].setWidth(2)
            self.curve[i].setPen(pen[i])
            self.curve[i].attach(self.qwtPlot)

    def on_pushbutton_Rafraichir_Clicked(self):
        self.on_timer()

    def on_pushbutton_Parcourir_Clicked(self):
        self.lineEdit_Rpertoire_Travail.setText(QFileDialog.getExistingDirectory())

    def on_checkbox_Clicked(self):
        self.Analyse_SBD(300434060873240)
        self.Analyse_SBD(300434062298540)
        self.Analyse_SBD(300434060657120)   #VMP    

    def on_timer(self):
        periode=int(self.spinBox.value()*60000)
        self.timer.start(periode)
        self.Interrogation_Planeurs()
        self.on_checkbox_Clicked()

    def Analyse_SBD(self,IMEI):
        repertoire=str(self.lineEdit_Rpertoire_Travail.text())+'/'+str(IMEI)
        print repertoire
        if not (os.path.exists(repertoire)):
            return
        dirs=os.listdir(repertoire)
        dirs.sort()
        Liste_points=[]
        f_synthese=open(repertoire+'\\'+'synthese_'+str(IMEI)+'.txt','w')        
        for fichier in dirs: 
            if fichier.split('.')[-1]=='sbd':
                f=open(repertoire+'\\'+fichier,"rb")
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
                    Nb_T_POSITION_REPORT=(len(Data_SBD)-3)/16
                    for i in range(Nb_T_POSITION_REPORT):
                        s=struct.Struct(">IBIIBBB")
                        T_POSITION_REPORT_DATA=s.unpack(Data_SBD[i*16+3:i*16+19])
                        GPS_TIME=T_POSITION_REPORT_DATA[0]+time.mktime(datetime.datetime.strptime("06/01/1980 00:00:00", "%d/%m/%Y %H:%M:%S").timetuple())
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
                            Liste_points.append((GPS_LONGITUDE,GPS_LATITUDE))
                            f_synthese.write( datetime.datetime.fromtimestamp(GPS_TIME).strftime('%d/%m/%Y %H:%M:%S')+ "\t"+str(GPS_LATITUDE) + "\t"+str(GPS_LONGITUDE)+'\r\n')
        f_synthese.close()
        x = [t[0] for t in Liste_points]
        y = [t[1] for t in Liste_points]

        if IMEI==300434060873240:
            if  not(self.checkBox_300434060873240.isChecked()):
                x=[]
                y=[]
            self.curve[0].setData(x,y)  
            self.lineEdit_300434060873240.setText(datetime.datetime.fromtimestamp(GPS_TIME).strftime('%d/%m/%Y %H:%M:%S')+ "   Lat : "+str(GPS_LATITUDE) + "   Lon : "+str(GPS_LONGITUDE))

        if IMEI==300434062298540:
            if not(self.checkBox_300434062298540.isChecked()):
                x=[]
                y=[]
            self.curve[1].setData(x,y)
            self.lineEdit_300434062298540.setText(datetime.datetime.fromtimestamp(GPS_TIME).strftime('%d/%m/%Y %H:%M:%S')+ "   Lat : "+str(GPS_LATITUDE) + "   Lon : "+str(GPS_LONGITUDE))

        if IMEI==300434060657120:
            if not(self.checkBox_VMP.isChecked()):
                x=[]
                y=[]
            self.curve[2].setData(x,y)
            self.lineEdit_VMP.setText(datetime.datetime.fromtimestamp(GPS_TIME).strftime('%d/%m/%Y %H:%M:%S')+ "   Lat : "+str(GPS_LATITUDE) + "   Lon : "+str(GPS_LONGITUDE))

        self.qwtPlot.replot()


    def Interrogation_Planeurs(self):
        try:
            server = imaplib.IMAP4_SSL('domicile.ifremer.fr', 993)
            self.label_Date_MAJ_Serveur.setText(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        except imaplib.IMAP4.error:
            self.label_Date_MAJ_Serveur.setText("pb login serveur")
            return
        server.login('planeurs_mp1', 'viadomicile29,caplane4')
        server.select()
        typ, data = server.search(None,  '(FROM "sbdservice@sbd.iridium.com")')
        for num in data[0].split():
            typ, data = server.fetch(num, '(RFC822)')
            email_body = data[0][1] #récupère le contenu du mail
            mail = email.message_from_string(email_body)
            #Vérifie si il y a une pièce jointe
            if mail.get_content_maintype() != 'multipart':	
                continue   			
            IMEI = mail["Subject"].split(": ")[1] #On récupère le numéro IMEI depuis le sujet
            #Détermine le dossier de destination des pièces jointes
#            detach_dir = "N:\GESTION_LOCALE_ISO9001\I_Instrumentation\Balise_Iridium" + "\\" + IMEI 
            detach_dir =str(self.lineEdit_Rpertoire_Travail.text())+ "\\" + IMEI 
            print detach_dir
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
                print att_path
                #Si le fichier existe déjà, passe au suivant
                if os.path.isfile(att_path) :
                    continue
               # finally write the stuff
                fp = open(att_path, 'wb')
                fp.write(part.get_payload(decode=True))
                fp.close()
            if 'OK'==server.copy(num,'BALISES_NOVATECH')[0]:
                        # remove it from inbox
                server.store(num,'+FLAGS','\\Deleted')

        #Ferme la connexion
        server.close()
        server.logout()


def main():
    app=QApplication(sys.argv)
    form=MainDialog()
    form.show()
    app.exec_()


if __name__ == "__main__":
    main()








