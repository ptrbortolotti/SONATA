#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
ZetCode PyQt5 tutorial 

This example shows an icon
in the titlebar of the window.

author: Jan Bodnar
website: zetcode.com 
last edited: January 2015
"""

import sys
from PyQt5.QtWidgets import (QWidget, QToolTip, QAction, qApp, QMainWindow,
    QPushButton, QApplication, QMessageBox)
from PyQt5.QtGui import QFont,QIcon    
from PyQt5.QtCore import QCoreApplication


class Example(QMainWindow):
    
    def __init__(self):
        super(Example,self).__init__()
        self.initUI()
        
        
    def initUI(self):
        
              
        QToolTip.setFont(QFont('SansSerif', 8))
        
        self.setToolTip('This is a <b>QWidget</b> widget')
        
        btn = QPushButton('Button', self)
        btn.setToolTip('This is a <b>QPushButton</b> widget')
        btn.resize(btn.sizeHint())
        btn.move(5, 25)       
    
        qbtn = QPushButton('Quit', self)
        qbtn.resize(qbtn.sizeHint())
        qbtn.move(5, 55)       
        
        exitAction = QAction(QIcon('logo.png'), 'Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.triggered.connect(qApp.quit)
        
        self.toolbar = self.addToolBar('Exit')
        self.toolbar.addAction(exitAction)
        
        self.setGeometry(300, 300, 300, 220)
        self.setWindowTitle('SONATA')
        self.setWindowIcon(QIcon('logo.png'))        
        self.show()
        

if __name__ == '__main__':
    
    app = QApplication(sys.argv)
    ex = Example()
    sys.exit(app.exec_()) 