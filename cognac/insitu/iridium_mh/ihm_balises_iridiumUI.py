# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ihm_balises_iridium.ui'
#
# Created: Fri Apr 06 13:44:02 2018
#      by: PyQt4 UI code generator 4.9.6
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName(_fromUtf8("Dialog"))
        Dialog.resize(874, 486)
        self.groupBox = QtGui.QGroupBox(Dialog)
        self.groupBox.setGeometry(QtCore.QRect(10, 10, 551, 111))
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.checkBox_300434060873240 = QtGui.QCheckBox(self.groupBox)
        self.checkBox_300434060873240.setGeometry(QtCore.QRect(10, 20, 161, 17))
        self.checkBox_300434060873240.setChecked(True)
        self.checkBox_300434060873240.setObjectName(_fromUtf8("checkBox_300434060873240"))
        self.checkBox_300434062298540 = QtGui.QCheckBox(self.groupBox)
        self.checkBox_300434062298540.setGeometry(QtCore.QRect(10, 50, 171, 17))
        self.checkBox_300434062298540.setChecked(True)
        self.checkBox_300434062298540.setObjectName(_fromUtf8("checkBox_300434062298540"))
        self.checkBox_VMP = QtGui.QCheckBox(self.groupBox)
        self.checkBox_VMP.setGeometry(QtCore.QRect(10, 80, 171, 17))
        self.checkBox_VMP.setChecked(True)
        self.checkBox_VMP.setObjectName(_fromUtf8("checkBox_VMP"))
        self.lineEdit_300434060873240 = QtGui.QLineEdit(self.groupBox)
        self.lineEdit_300434060873240.setGeometry(QtCore.QRect(190, 20, 351, 20))
        self.lineEdit_300434060873240.setObjectName(_fromUtf8("lineEdit_300434060873240"))
        self.lineEdit_300434062298540 = QtGui.QLineEdit(self.groupBox)
        self.lineEdit_300434062298540.setGeometry(QtCore.QRect(190, 50, 351, 20))
        self.lineEdit_300434062298540.setObjectName(_fromUtf8("lineEdit_300434062298540"))
        self.lineEdit_VMP = QtGui.QLineEdit(self.groupBox)
        self.lineEdit_VMP.setGeometry(QtCore.QRect(190, 80, 351, 20))
        self.lineEdit_VMP.setObjectName(_fromUtf8("lineEdit_VMP"))
        self.qwtPlot = Qwt5.QwtPlot(Dialog)
        self.qwtPlot.setGeometry(QtCore.QRect(10, 180, 851, 291))
        self.qwtPlot.setObjectName(_fromUtf8("qwtPlot"))
        self.groupBox_2 = QtGui.QGroupBox(Dialog)
        self.groupBox_2.setGeometry(QtCore.QRect(570, 10, 191, 91))
        self.groupBox_2.setObjectName(_fromUtf8("groupBox_2"))
        self.spinBox = QtGui.QSpinBox(self.groupBox_2)
        self.spinBox.setGeometry(QtCore.QRect(20, 20, 42, 22))
        self.spinBox.setProperty("value", 10)
        self.spinBox.setObjectName(_fromUtf8("spinBox"))
        self.label = QtGui.QLabel(self.groupBox_2)
        self.label.setGeometry(QtCore.QRect(20, 50, 181, 16))
        self.label.setObjectName(_fromUtf8("label"))
        self.label_Date_MAJ_Serveur = QtGui.QLabel(self.groupBox_2)
        self.label_Date_MAJ_Serveur.setGeometry(QtCore.QRect(20, 70, 151, 16))
        self.label_Date_MAJ_Serveur.setAlignment(QtCore.Qt.AlignCenter)
        self.label_Date_MAJ_Serveur.setObjectName(_fromUtf8("label_Date_MAJ_Serveur"))
        self.pushButton_Rafraichir = QtGui.QPushButton(self.groupBox_2)
        self.pushButton_Rafraichir.setGeometry(QtCore.QRect(80, 20, 75, 23))
        self.pushButton_Rafraichir.setObjectName(_fromUtf8("pushButton_Rafraichir"))
        self.label_logo_ifremer = QtGui.QLabel(Dialog)
        self.label_logo_ifremer.setGeometry(QtCore.QRect(640, 110, 221, 61))
        self.label_logo_ifremer.setObjectName(_fromUtf8("label_logo_ifremer"))
        self.label_logo_lops = QtGui.QLabel(Dialog)
        self.label_logo_lops.setGeometry(QtCore.QRect(780, 10, 91, 71))
        self.label_logo_lops.setObjectName(_fromUtf8("label_logo_lops"))
        self.groupBox_3 = QtGui.QGroupBox(Dialog)
        self.groupBox_3.setGeometry(QtCore.QRect(10, 120, 541, 51))
        self.groupBox_3.setObjectName(_fromUtf8("groupBox_3"))
        self.lineEdit_Rpertoire_Travail = QtGui.QLineEdit(self.groupBox_3)
        self.lineEdit_Rpertoire_Travail.setGeometry(QtCore.QRect(10, 20, 421, 20))
        self.lineEdit_Rpertoire_Travail.setObjectName(_fromUtf8("lineEdit_Rpertoire_Travail"))
        self.pushButton_Parcourir = QtGui.QPushButton(self.groupBox_3)
        self.pushButton_Parcourir.setGeometry(QtCore.QRect(450, 20, 75, 23))
        self.pushButton_Parcourir.setObjectName(_fromUtf8("pushButton_Parcourir"))

        self.retranslateUi(Dialog)
        QtCore.QObject.connect(self.pushButton_Rafraichir, QtCore.SIGNAL(_fromUtf8("clicked()")), Dialog.on_pushbutton_Rafraichir_Clicked)
        QtCore.QObject.connect(self.checkBox_300434060873240, QtCore.SIGNAL(_fromUtf8("clicked()")), Dialog.on_checkbox_Clicked)
        QtCore.QObject.connect(self.checkBox_300434062298540, QtCore.SIGNAL(_fromUtf8("clicked()")), Dialog.on_checkbox_Clicked)
        QtCore.QObject.connect(self.checkBox_VMP, QtCore.SIGNAL(_fromUtf8("clicked()")), Dialog.on_checkbox_Clicked)
        QtCore.QObject.connect(self.pushButton_Parcourir, QtCore.SIGNAL(_fromUtf8("clicked()")), Dialog.on_pushbutton_Parcourir_Clicked)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(_translate("Dialog", "Decodage balises IRIDIUM NOVATECH", None))
        self.groupBox.setTitle(_translate("Dialog", "Selection IMEI balises                           Dernier point recu", None))
        self.checkBox_300434060873240.setText(_translate("Dialog", "300434060873240 (Rouge)", None))
        self.checkBox_300434062298540.setText(_translate("Dialog", "300434062298540 (Noir)", None))
        self.checkBox_VMP.setText(_translate("Dialog", "300434060657120 VMP (Bleu)", None))
        self.groupBox_2.setTitle(_translate("Dialog", "Periode interrogation serveur (mn)", None))
        self.label.setText(_translate("Dialog", "Derniere interrogation serveur :", None))
        self.label_Date_MAJ_Serveur.setText(_translate("Dialog", "##/##/#### ##:##:##", None))
        self.pushButton_Rafraichir.setText(_translate("Dialog", "Rafraichir", None))
        self.label_logo_ifremer.setText(_translate("Dialog", "TextLabel", None))
        self.label_logo_lops.setText(_translate("Dialog", "TextLabel", None))
        self.groupBox_3.setTitle(_translate("Dialog", "Repertoire de travail", None))
        self.lineEdit_Rpertoire_Travail.setText(_translate("Dialog", "c:\\", None))
        self.pushButton_Parcourir.setText(_translate("Dialog", "Parcourir", None))

from PyQt4 import Qwt5

if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    Dialog = QtGui.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())

