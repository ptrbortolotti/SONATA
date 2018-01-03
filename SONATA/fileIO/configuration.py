# -*- coding: utf-8 -*-
"""
Created on Wed Jan 03 13:24:23 2018

@author: TPflumm
"""

from SONATA.topo.utils import allunique
from SONATA.fileIO.readinput import read_segment, read_TXTrowSTR, read_INTrowSTR,\
                                    read_BOOLrowSTR, read_FLOATrowSTR, read_layup


#================ SECTION-CONFIGURATION OBJECT ========================================

class Configuration(object):
    
    def __init__(self, filename=None):
        self.SETUP_scale_factor = 1
        self.SETUP_Theta = 0

        
        self.WEB_ID = []
        self.WEB_Pos1 = []
        self.WEB_Pos2 = []
        self.SEG_ID = []    
        self.SEG_CoreMaterial = []
        self.SEG_Layup = []
        self.SEG_Boundary_OCC = []      #Segment Boundaries in Opencascade format (TOPO_DSwire)
        self.SEG_Boundary_DCT = []      #Segment Boundaries in discrete values (np.array) 

        if filename:
            self.read_config(filename)


    def read_config(self,filename):

        #READ FILE AND CLEAN UP COMMENTS, EMPTY LINES WITH SPACES AND NEWLINES
        a = ''
        with open(filename) as f:
            for line in f:
                line = line.partition('#')[0]
                line = line.rstrip()
                a += line 
                a += '\n'
             
        STR = ''.join([s for s in a.strip().splitlines(True) if s.strip("\r\n").strip()])
        
        
        #READ SETUP SEGMENT OF INPUT FILE
        SETUP_str,SETUP_start,SETUP_end = read_segment(STR,'Setup')  
        #print SETUP_str,SETUP_start,SETUP_end,'\n'
        self.SETUP_mat_filename = read_TXTrowSTR(SETUP_str,'mat_filename')
        self.SETUP_NbOfWebs = read_INTrowSTR(SETUP_str,'NbOfWebs')
        self.SETUP_BalanceWeight = read_BOOLrowSTR(SETUP_str,'BalanceWeight')
        self.SETUP_input_type = read_INTrowSTR(SETUP_str,'input_type')
        self.SETUP_datasource = read_TXTrowSTR(SETUP_str,'datasource')
        if self.SETUP_input_type == 3 or self.SETUP_input_type == 4:
            self.SETUP_radial_station = read_FLOATrowSTR(SETUP_str,'radial_station')
        self.SETUP_scale_factor = read_FLOATrowSTR(SETUP_str,'scale_factor')
        self.SETUP_Theta = read_FLOATrowSTR(SETUP_str,'Theta')
        self.SETUP_mesh_resolution = read_INTrowSTR(SETUP_str,'mesh_resolution')
        #self.SETUP_Airfoil = read_TXTrowSTR(SETUP_str,'Airfoil')
        #self.SETUP_chord = read_FLOATrowSTR(SETUP_str,'chord')
                                            
                     
        #======================================================================
        #CHECK and READ BALANCE WEIGHT DEFINITION!   
        if self.SETUP_BalanceWeight == True:
            #READ Balance Weight Definition
            BW_str,BW_start,BW_end = read_segment(STR,'BalanceWeight')  
            #print BW_str,BW_start,BW_end,'\n'
            self.BW_MatID   = read_INTrowSTR(BW_str,'MatID')
            self.BW_XPos    = read_FLOATrowSTR(BW_str,'XPos')
            self.BW_YPos    = read_FLOATrowSTR(BW_str,'YPos')
            self.BW_Diameter= read_FLOATrowSTR(BW_str,'Diameter')
    
    
        #======================================================================
        #CHECK and READ WEB DEFINITION!
        if self.SETUP_NbOfWebs != STR.count('&DEFN Web'):
            print 'WARNING: \t Setup variable "NbOfWebs" is not equal to the number of web definitions'
            
        if self.SETUP_NbOfWebs > 0:
            #READ Balance Weight Definition
            for j in range(0,self.SETUP_NbOfWebs):    
                WEB_str,WEB_start,WEB_end = read_segment(STR,'Web')              
        
                self.WEB_ID.append(read_INTrowSTR(WEB_str,'WebID'))
                self.WEB_Pos1.append(read_FLOATrowSTR(WEB_str,'Pos1'))
                self.WEB_Pos2.append(read_FLOATrowSTR(WEB_str,'Pos2'))

                #Replace Characters from WEB_start to WEB_end with whitspaces
                for k in range(WEB_start,WEB_end):  
                    templist = list(STR)
                    templist[k] = ' '
                    STR = ''.join(templist)
                #WEB_str,WEB_start,WEB_end = read_segment(STR,'Web')      
                
        if not allunique(self.WEB_ID):
            print 'WARNING: \t WEB IDs are not unique!'
  
  
        #======================================================================
        #CHECK and READ SEGMENT COMPOSITE LAYUP DEFINITION!
        SEG_str,SEG_start,SEG_end = read_segment(STR,'Seg')
             
                
        while SEG_str!= '':
            self.SEG_ID.append(read_INTrowSTR(SEG_str,'SegID'))
            self.SEG_CoreMaterial.append(read_INTrowSTR(SEG_str,'CoreMaterial'))
            self.SEG_Layup.append(read_layup(SEG_str))
            
             #Replace Characters from SEG_start to SEG_end with whitspaces
            for k in range(SEG_start,SEG_end):  
                templist = list(STR)
                templist[k] = ' '
                STR = ''.join(templist)
                
            SEG_str,SEG_start,SEG_end = read_segment(STR,'Seg')       
        
        #CHECK FOR SOME INPUT MISTAKES:
        if not allunique(self.SEG_ID):
            print 'WARNING: \t SEG IDs are not unique'        
        
        if not(self.SETUP_NbOfWebs == 0 and len(self.SEG_ID) == 1) and (not(self.SETUP_NbOfWebs+2 == len(self.SEG_ID))):   
            print 'WARNING: \t The number of segments does not corresponds to the number of Webs'




