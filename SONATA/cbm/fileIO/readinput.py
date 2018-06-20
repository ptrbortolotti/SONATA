#==============================READ SONATA2D - INPUT FILE =====================
#This file provides the class and its function to read the input information for the structural design of the rotorblade crosssection
# - The format must be according to the sec_config.input file. 
# 
#TODO: Include more checks for input format
#TODO: Include Definition of an erosion protection   
#

#Author: Tobias Pflumm
#Date:	09/21/2016
#==============================================================================
import numpy as np 
import ast
from urllib.request import urlopen

from SONATA.cbm.fileIO.material import Material


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
    if End<Start:
        End=len(STR)
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


def read_LISTrowSTR(STR,STR2Find):
    Start,End = read_rowstring(STR,STR2Find)
    temp = STR[Start:End]
    temp = temp.split('=')[1]
    temp = temp.replace(" ", "")
    temp = ast.literal_eval(temp)
    return temp 


def read_CMATRIX(STR,STR2Find):
    Start,End = read_rowstring(STR,STR2Find)
    End = STR[Start::].find(']]')+Start
    temp = STR[Start:End+len(']]')]
    temp = temp.split('=')[1]
    temp = temp.replace(" ", "")
    temp = temp.replace("\t", "")
    temp = temp.replace("\n", ",")
    temp = temp.replace("[,", "[float('nan'),")
    temp = temp.replace(",,,,,", ",float('nan'),float('nan'),float('nan'),float('nan'),")
    temp = temp.replace(",,,,", ",float('nan'),float('nan'),float('nan'),")
    temp = temp.replace(",,,", ",float('nan'),float('nan'),")
    temp = temp.replace(",,", ",float('nan'),")
    temp = temp.replace("float('nan'),", '696969,')
    a = np.asarray(ast.literal_eval(temp))
    a[a==696969] = float('nan') #replace all 696969's with nan's
    return a
    

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
    foil_dat_url = 'http://m-selig.ae.illinois.edu/ads/coord_seligFmt/%s.dat' % name
    with urlopen(foil_dat_url) as f:
        temp_x = []
        temp_y = []
        temp_z = []
        for line in f.readlines()[1:]:      
            line = line.decode('utf-8')                                                            # The first line contains info only
            line = line.lstrip().rstrip().replace('    ', ' ').replace('   ', ' ').replace('  ', ' ')   # do some cleanup on the data (mostly dealing with spaces)
            data = line.split(' ')   
            temp_y.append(float(data[0]))                                                               # data[0] = x coord.
            temp_z.append(float(data[1]))                                                               # data[1] = y coord.
            temp_x.append(float(0))                                                                     # data[2] = z coord 
    return np.array([temp_x,temp_y,temp_z])                                                         # return AirfoilCoordinate as np.arrray


def UIUCAirfoil2d(name):
    foil_dat_url = 'http://m-selig.ae.illinois.edu/ads/coord_seligFmt/%s.dat' % name
    with urlopen(foil_dat_url) as f:
        temp_x = []
        temp_y = []
        for line in f.readlines()[1:]:                                                                  # The first line contains info only
            line = line.decode('utf-8')            
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
    for line in f.readlines()[1:]:    
        line = line.decode('utf-8')                                                                # The first line contains info only
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


def read_material_input(filename):
    MaterialLst=[]
#    filename='mat_database.input'
    a = ''
    with open(filename) as f:
        for line in f:
            line = line.partition('#')[0]
            line = line.rstrip()
            a += line 
            a += '\n'
         
    STR = ''.join([s for s in a.strip().splitlines(True) if s.strip("\r\n").strip()])
    #NbMat = STR.count('&DEFN Material')

    MAT_str,MAT_start,MAT_end = read_segment(STR,'Material')
   
    while MAT_str!= '':
        MAT_str = MAT_str[len('&DEFN Material'):-len('&END')]

        #determine MatID, name and orth
        MatID = read_INTrowSTR(MAT_str,'MatID')
        name = read_TXTrowSTR(MAT_str,'name')
        orth = read_INTrowSTR(MAT_str,'orth')
        rho = read_FLOATrowSTR(MAT_str,'rho')
        
        #Remove long name string:
        Start_name,End_name = read_rowstring(MAT_str,'name')
        End_name = MAT_str[Start_name:].index('\n') + Start_name
        for k in range(Start_name,End_name):  
                templist = list(MAT_str)
                templist[k] = ' '
                MAT_str = ''.join(templist)
   
        if orth == 0:
            #print 'Isotropic material'
            E = read_FLOATrowSTR(MAT_str,'E')
            nu = read_FLOATrowSTR(MAT_str,'nu')
            alpha = read_FLOATrowSTR(MAT_str,'alpha')
            MaterialLst.append(Material(MatID,name,orth,rho,E=E,nu=nu,alpha=alpha))
            
        elif orth == 1:
            #print 'orthotropic material'
            E =  np.asarray(read_LISTrowSTR(MAT_str,'E'))
            G =  np.asarray(read_LISTrowSTR(MAT_str,'G'))
            nu =  np.asarray(read_LISTrowSTR(MAT_str,'nu'))
            alpha =  np.asarray(read_LISTrowSTR(MAT_str,'alpha'))
            MaterialLst.append(Material(MatID,name,orth,rho,E=E,G=G,nu=nu,alpha=alpha))
            
        elif orth == 2:
            #print 'general anisotropic material'
            C = read_CMATRIX(MAT_str,'C')
            alpha = np.asarray(read_LISTrowSTR(MAT_str,'alpha'))
            MaterialLst.append(Material(MatID,name,orth,rho,C=C,alpha=alpha))
        
        for k in range(MAT_start,MAT_end):  
                templist = list(STR)
                templist[k] = ' '
                STR = ''.join(templist)
                
        MAT_str,MAT_start,MAT_end = read_segment(STR,'Material') 
    
    MaterialLst.sort(key=lambda Material: (Material.id))
    return MaterialLst


#======================================================
#       MAIN
#======================================================       
if __name__ == '__main__':
    filename = 'sec_config.input'
    #section1 = Configuration(filename)
    
    filename='mat_database.input'
    MaterialLst = read_material_input(filename)
    