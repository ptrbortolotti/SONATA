


#UTILITY READ FUNCTIONS:

def read_segment(STR,seg2find):
    str2find = '&DEFN '+seg2find
    start1 = STR.find(str2find)
    start2 = STR[start1+len(str2find)::].find('&DEFN')
    end1  = STR[start1+len(str2find)::].find('&END' )
    end = end1
    if start2<end1:
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
    

#================ SECTION-CONFIG OBJECT ========================================

class section_config(object):
    def __init__(self):
        pass
    
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
        self.NbOfWebs = read_INTrowSTR(SETUP_str,'NbOfWebs')
        self.BalanceWeight = read_BOOLrowSTR(SETUP_str,'BalanceWeight')
        self.Airfoil = read_TXTrowSTR(SETUP_str,'Airfoil')
                     
                   
        if self.BalanceWeight == True:
            #READ Balance Weight Definition
            BW_str,BW_start,BW_end = read_segment(STR,'BalanceWeight')  
            #print BW_str,BW_start,BW_end,'\n'
            self.BW_MatID   = read_INTrowSTR(BW_str,'MatID')
            self.BW_XPos    = read_FLOATrowSTR(BW_str,'XPos')
            self.BW_YPos    = read_FLOATrowSTR(BW_str,'YPos')
            self.BW_Diameter= read_FLOATrowSTR(BW_str,'Diameter')
            
        if self.NbOfWebs > 0:
            #READ Balance Weight Definition
            for j in range(0,self.NbOfWebs):    
                WEB_str,WEB_start,WEB_end = read_segment(STR,'Web')
                print WEB_str
                #Replace Characters from WEB_start to WEB_end with whitspaces
                for k in range(WEB_start,WEB_end):  
                    templist = list(STR)
                    templist[k] = ' '
                    STR = ''.join(templist)
                WEB_str,WEB_start,WEB_end = read_segment(STR,'Web')      
                print WEB_str
                
                
                
#======================================================
#       MAIN
#======================================================
if __name__ == '__main__':
    
    filename = 'sec_config.input'
    section1 = section_config()
    section1.read_config(filename)










