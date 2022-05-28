import sys
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QDialog, QApplication
from poprawne3 import *
import numpy as np
from math import sqrt, sin, cos, radians,degrees
import math

class MyForm(QDialog):
    def __init__(self):
        super().__init__() 
        self.ui = Ui_Dialog()
        self.ui.setupUi(self) 
        self.setWindowIcon(QtGui.QIcon('ikona.jpg'))
        
        if self.ui.radioButton_2.isChecked:
            self.a = 6378137.0 
            self.b = 6356752.31424518

        elif self.ui.radioButton.isChecked:
            self.a = 6378137.0
            self.b = 6356752.31414036

        self.flat = (self.a - self.b) / self.a
        self.ecc = sqrt(2 * self.flat - self.flat ** 2) 
        self.ecc2 = (2 * self.flat - self.flat ** 2)
    
        if self.ui.radioButton_3.isChecked:
            self.LAM0 = 15
            self.s = 5
        
        if self.ui.radioButton_4.isChecked:
            self.LAM0 = 18
            self.s = 6
            
        if self.ui.radioButton_5.isChecked:
            self.LAM0 = 21
            self.s = 7
        if self.ui.radioButton_6.isChecked:
            self.LAM0 = 24
            self.s = 8

            
        self.ui.pushButton.clicked.connect(self.u2000)
        self.ui.pushButton_2.clicked.connect(self.u1992) 
        self.ui.pushButton_3.clicked.connect(self.N)

    
        self.show()



    def u2000(self):
        m=0.999923
        
        if len(self.ui.lineEdit.text())!=0:
            fi=float(self.ui.lineEdit.text())
        else:
            fi=0
                
        if len(self.ui.lineEdit_2.text())!=0:
            lam=float(self.ui.lineEdit_2.text())
        else:
            lam=0
        
        fi=radians(fi)
        lam=radians(lam)
        N = self.a/math.sqrt(1-self.ecc2*math.sin(fi)**2)
        t = np.tan(fi)
        e_2 = self.ecc2/(1-self.ecc2)
        n2 = e_2 * (np.cos(fi))**2    

        LAM0 = math.radians(self.LAM0)
        l = lam - LAM0
        A0 = 1 - (self.ecc2/4) - ((3*(self.ecc2**2))/64) - ((5*(self.ecc2**3))/256)   
        A2 = (3/8) * (self.ecc2 + ((self.ecc2**2)/4) + ((15 * (self.ecc2**3))/128))
        A4 = (15/256) * (self.ecc2**2 + ((3*(self.ecc2**3))/4))
        A6 = (35 * (self.ecc2**3))/3072 
        sig = self.a * ((A0*fi) - (A2*np.sin(2*fi)) + (A4*np.sin(4*fi)) - (A6*np.sin(6*fi))) 
        x = sig + ((l**2)/2) * N *np.sin(fi) * np.cos(fi) * (1 + ((l**2)/12) * ((math.cos(fi))**2) * (5 - t**2 + 9*n2 + 4*(n2**2)) + ((l**4)/360) * ((math.cos(fi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
        y = l * (N*math.cos(fi)) * (1 + ((((l**2)/6) * (math.cos(fi))**2) * (1-t**2+n2)) +  (((l**4)/(120)) * (math.cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))
        x00 = round(m * x, 3) 
        y00 = round(m * y + (self.s*1000000) + 500000, 3)
        self.ui.label_3.setText('wynik ' + str(x00) + ', ' + str(y00))


    

        
    def u1992(self):
        if len(self.ui.lineEdit.text())!=0:
            phi=float(self.ui.lineEdit.text())
        else:
            phi=0
                
        if len(self.ui.lineEdit_2.text())!=0:
            lam=float(self.ui.lineEdit_2.text())
        else:
            lam=0
            
            
        phi = radians(phi)
        lam = radians(lam)
        L0 = radians(19)
        e2_ = (self.a**2 - self.b**2)/(self.b**2)
        eta2 = e2_ * cos(phi)**2
        t = math.tan(phi)
        l = lam - L0
        N = self.a/(1-self.ecc2*(sin(phi))**2)**(0.5)
        A0 = 1 - (self.ecc2/4) - ((3*(self.ecc2**2))/64) - ((5*(self.ecc2**3))/256)   
        A2 = (3/8) * (self.ecc2 + ((self.ecc2**2)/4) + ((15 * (self.ecc2**3))/128))
        A4 = (15/256) * (self.ecc2**2 + ((3*(self.ecc2**3))/4))
        A6 = (35 * (self.ecc2**3))/3072 
        sig = self.a * ((A0*phi) - (A2*sin(2*phi)) + (A4*sin(4*phi)) - (A6*sin(6*phi)))
        x = sig + (l**2)/2 * N * sin(phi) * cos(phi) * ( 1+(((l**2)/12) * (cos(phi)**2) * (5 - t**2 + 9*eta2 + 4*(eta2**2)) ) + (((l**4)/360) * (cos(phi)**4) * (61 - 58*(t**2) + (t**4) + 270*eta2 - 330*(eta2)*(t**2)) ) )
        y = l * N * cos(phi) * ( 1+(((l**2)/6) * (cos(phi)**2) * (1 - t**2 + eta2) ) + (((l**4)/120) * (cos(phi)**4) * (5 - 18*(t**2) + (t**4) + 14*eta2 - 58*(eta2)*(t**2)) ) )
        m0 = 0.9993
        x92 = round(m0*x - 5300000, 3)
        y92 = round(m0*y + 500000,3)
        self.ui.label_4.setText('wynik: ' + str(x92) + ', ' + str(y92))
            
            
            
    def N(self):
        if len(self.ui.lineEdit.text())!=0:
            phi=float(self.ui.lineEdit.text())
        else:
            phi=0
                
            
        Np=self.a / np.sqrt((1 - self.ecc2 * (np.sin(phi)) ** 2))
        self.ui.label_5.setText('wynik: ' + str(round(Np,3)))
            
    
    

            
    
        
   
        
        
        
        
if __name__=="__main__":
    app = QApplication(sys.argv)
    w = MyForm() 
    w.show() 
    sys.exit(app.exec_())

