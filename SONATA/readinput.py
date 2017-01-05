#==============================READ SONATA2D - INPUT FILE =====================
#This file provides the class and its function to read the input information for the structural design of the rotorblade crosssection
# - The format must be according to the sec_config.input file. 
# 
#TBD: Include more checks for input format
#TBD: Include Definition of an erosion protection   
#
#Author: Tobias Pflumm
#Date:	09/21/2016
#==============================================================================
import numpy as np 
import urllib2  

def read_segment(STR,seg2find):
    str2find = '&DEFN '+seg2find
    start1 = STR.find(str2find)
    start2 = STR[start1+len(str2find)::].find('&DEFN')
    end1  = STR[start1+len(str2find)::].find('&END' )
    end = end1
    if start2 > 0 and start2<end1:
        end  = STR[start1+len(str2find)+end1+len('&END')::].find('&END' )+start1+end1+2*len('&END')+len(str2find)
    else:
        end = start1+end1+len('&END')+len(str2find)
    temp = STR[start1:end]
    return temp,start1,end
    
def read_rowstring(STR,STR2Find):  
    Start = STR.find(STR2Find)
    End = STR[Start::].find('\n')+Start
    return Start,End    
    
def str2bool(v):
  return v.lower() in ("yes", "true", "t", "1")

def read_INTrowSTR(STR,STR2Find):
    Start,End = read_rowstring(STR,STR2Find)
    temp = STR[Start:End]
    temp = temp.split('=')[1]
    temp = temp.replace(" ", "")
    return int(temp)

def read_FLOATrowSTR(STR,STR2Find):
    Start,End = read_rowstring(STR,STR2Find)
    temp = STR[Start:End]
    temp = temp.split('=')[1]
    temp = temp.replace(" ", "")
    return float(temp)
   
def read_BOOLrowSTR(STR,STR2Find):
    Start,End = read_rowstring(STR,STR2Find)
    temp = STR[Start:End]
    temp = temp.split('=')[1]
    temp = temp.replace(" ", "")
    return str2bool(temp)    

def read_TXTrowSTR(STR,STR2Find):
    Start,End = read_rowstring(STR,STR2Find)
    temp = STR[Start:End]
    temp = temp.split('=')[1]
    temp = temp.replace(" ", "")
    return temp 

def read_layup(STR):
    str2find = '&DEFN Layup'
    start1 = STR.find(str2find)
    start2 = STR[start1+len(str2find)::].find('&DEFN')
    end1  = STR[start1+len(str2find)::].find('&END' )
    end = end1
    if start2 > 0 and start2<end1:
        end  = STR[start1+len(str2find)+end1+len('&END')::].find('&END' )+start1+end1+2*len('&END')+len(str2find)
    else:
        end = start1+end1+len('&END')+len(str2find)
    temp = STR[start1:end]
    first = temp.find('\n')
    last = temp.rfind('\n') 
    temp = temp[first+len('\n'):last]
    
    list_temp = temp.split('\n')
    NbOfLayers = temp.count('\n')+1                #The Number of Layers are the determined by the number of rows
    
    for j in range(0,NbOfLayers):  
        list_temp[j] = list_temp[j].split()
    
    #CONVERT TABLE TO np.ARRAY
    x = np.asarray(list_temp)
    nplayup = x[:,:5].astype(np.float)
    
    return nplayup
    

def UIUCAirfoil(name):
    foil_dat_url = 'http://www.ae.illinois.edu/m-selig/ads/coord_seligFmt/%s.dat' % name
    f = urllib2.urlopen(foil_dat_url)
    temp_x = []
    temp_y = []
    temp_z = []
    for line in f.readlines()[1:]:                                                                  # The first line contains info only
        line = line.lstrip().rstrip().replace('    ', ' ').replace('   ', ' ').replace('  ', ' ')   # do some cleanup on the data (mostly dealing with spaces)
        data = line.split(' ')   
        temp_y.append(float(data[0]))                                                               # data[0] = x coord.
        temp_z.append(float(data[1]))                                                               # data[1] = y coord.
        temp_x.append(float(0))                                                                     # data[2] = z coord 
    return np.array([temp_x,temp_y,temp_z])                                                         # return AirfoilCoordinate as np.arrray

def UIUCAirfoil2d(name):
    foil_dat_url = 'http://www.ae.illinois.edu/m-selig/ads/coord_seligFmt/%s.dat' % name
    f = urllib2.urlopen(foil_dat_url)
    temp_x = []
    temp_y = []
    for line in f.readlines()[1:]:                                                                  # The first line contains info only
        line = line.lstrip().rstrip().replace('    ', ' ').replace('   ', ' ').replace('  ', ' ')   # do some cleanup on the data (mostly dealing with spaces)
        data = line.split(' ')   
        temp_x.append(float(data[0]))                                                               # data[0] = x coord.
        temp_y.append(float(data[1]))                                                               # data[1] = y coord.
    return np.array([temp_x,temp_y])    

def AirfoilDat(name):
    string = "%s" %(name)
    f = open(string)
    temp_x = []
    temp_y = []
    temp_z = []
    for line in f.readlines()[1:]:                                                                  # The first line contains info only
        line = line.lstrip().rstrip().replace('    ', ' ').replace('   ', ' ').replace('  ', ' ')   # do some cleanup on the data (mostly dealing with spaces)
        data = line.split(' ')   
        temp_y.append(float(data[0]))                                                               # data[0] = x coord.
        temp_z.append(float(data[1]))                                                               # data[1] = y coord.
        temp_x.append(float(0))                                                                     # data[2] = z coord 
    return np.array([temp_x,temp_y,temp_z])                                                         # return AirfoilCoordinate as np.arrray    
    
def AirfoilDat2d(name):
    string = "%s" %(name)
    f = open(string)
    temp_x = []
    temp_y = []
    for line in f.readlines()[1:]:                                                                  # The first line contains info only
        line = line.lstrip().rstrip().replace('    ', ' ').replace('   ', ' ').replace('  ', ' ')   # do some cleanup on the data (mostly dealing with spaces)
        line = line.replace('\t', ' ')   # do some cleanup on the data (mostly dealing with spaces)
        data = line.split(' ')
        temp_x.append(float(data[0]))                                                             
        temp_y.append(float(data[1]))                                                                                                       
    return np.array([temp_x,temp_y])                                                                # return AirfoilCoordinate as np.arrray    

    

#================ SECTION-CONFIG OBJECT ========================================

class section_config(object):
    
    def __init__(self, filename):
        self.WEB_ID = []
        self.WEB_Pos1 = []
        self.WEB_Pos2 = []
        self.SEG_ID = []    
        self.SEG_CoreMaterial = []
        self.SEG_Layup = []
        self.SEG_Boundary_OCC = []      #Segment Boundaries in Opencascade format (TOPO_DSwire)
        self.SEG_Boundary_DCT = []      #Segment Boundaries in discrete values (np.array) 
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
        self.SETUP_NbOfWebs = read_INTrowSTR(SETUP_str,'NbOfWebs')
        self.SETUP_BalanceWeight = read_BOOLrowSTR(SETUP_str,'BalanceWeight')
        self.SETUP_Airfoil = read_TXTrowSTR(SETUP_str,'Airfoil')
        self.SETUP_chord = read_FLOATrowSTR(SETUP_str,'chord')             
                     
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
            print 'WARNING: Setup variable "NbOfWebs" is not equal to the number of web definitions'
            
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
                
        if (len(self.WEB_ID) > 1) and (len(set(self.WEB_ID)) == len(self.WEB_ID)):
            print 'WARNING: WEB IDs are not unique!'
  
  
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
        if (len(self.SEG_ID) > 1) and (not(len(set(self.SEG_ID)) == len(self.SEG_ID))):
            print 'WARNING: SEG IDs are not unique'        
        
        if not(self.SETUP_NbOfWebs == 0 and len(self.SEG_ID) == 1) and (not(self.SETUP_NbOfWebs+2 == len(self.SEG_ID))):   
            print 'WARNING: Make sure the number of Composite Layup Segments corresponds to your chosen number of Webs'


#======================================================
#       MAIN
#======================================================
if __name__ == '__main__':
    filename = 'sec_config.input'
    section1 = section_config()
    section1.read_config(filename)