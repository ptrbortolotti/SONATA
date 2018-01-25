"""
Created on Mon Mar 20 14:15:44 2017

@author: TPflumm
"""
import numpy as np
from StringIO import StringIO 
from SONATA.fileIO.readinput import read_rowstring, read_INTrowSTR, read_FLOATrowSTR, read_BOOLrowSTR, read_TXTrowSTR, read_LISTrowSTR

def read_VABS_Results(filename):
    a = ''
    with open(filename) as f:
        for line in f:
            line = line.partition('#')[0]
            line = line.rstrip()
            a += line 
            a += '\n'
         
    STR = ''.join([s for s in a.strip().splitlines(True) if s.strip("\r\n").strip()])
    return np.loadtxt(StringIO(STR))        


class VABS_config(object):

    def __init__(self, **kw):
        self.format_flag = 0     #0:old format, 1:new format
        self.nlayer = 0          #nlayer is not used in the old format. nlayer shoud be always given a value greater than one if format_flag=1
        self.Timoshenko_flag = 1 #VABS will construct both the classical model and the generalized Timoshenko model. If it is 0, it will only construct the classical model. 
        self.recover_flag = 0    #0: VABS will carry out the consitutive modeling. 1: non-linear beam theory, 3D stress,displacement,strain recovery, 2:linear beam theory 3D recovery
        self.thermal_flag = 0    #Either 0 or 3: 0:pure mechanical analysis. 3:one-way coupled thermoelastic analyis.
        self.curve_flag = 0      #To model initially curved or twisted beams. if 1: three real numbers for twist(k1) and curvatures (k2 and k3) shoud be provided in the very next line    
        self.trapeze_flag = 0    #to obtain the trapeze effect trapeze_flag is 1
        self.Vlasov_flag = 0     #To obtain a generalized Vlasov model, can be 1 only if Timoshenko_flag is 1. For more information see VABS USER MANUAL
        self.oblique_flag = 0    #to model oblique cross sections. See VABS USER MANUAL
        
        if 'format_flag' in kw:     self.format_flag = kw['format_flag']  #Currently no need to use this!
        if 'nlayer' in kw:          self.nlayer = kw['nlayer']  
        if 'Timoshenko_flag' in kw: self.Timoshenko_flag = kw['Timoshenko_flag']
        if 'recover_flag' in kw:    self.recover_flag = kw['recover_flag']
        if 'thermal_flag' in kw:    self.thermal_flag = kw['thermal_flag']
        if 'trapeze_flag' in kw:    self.trapeze_flag = kw.get('trapeze_flag')
        if 'Vlasov_flag'  in kw:    self.Vlasov_flag = kw['Vlasov_flag']
        
        if 'curve_flag' in kw:      
            self.curve_flag = kw['curve_flag']
            self.k1 = kw['k1']
            self.k2 = kw['k2']
            self.k3 = kw['k3']
        if 'oblique_flags' in kw:
            self.oblique_flag = kw['oblique_flag']     #to model oblique cross sections. See VABS USER MANUAL
            self.oblique_cosine1 = kw['oblique_cosine1']             #Angle between beam axis x1 and and oblique axis y1. See VABS USER MANUAL Figure 6
            self.oblique_cosine2 = kw['oblique_cosine2']             #Angle between beam axis x1 and and oblique axis y2. See VABS USER MANUAL Figure 6
            self.Timoshenko_flag = 0                             #THis feature is only enabled in the classical beam model
            
        
        self.u = [0,0,0]                                    
        self.Cij = np.array([[0,0,0],[0,0,0],[0,0,0]])
        self.F = [0,0,0]                #Timoshenko_flag == 0 -> F = [0]    
        self.M = [0,0,0]
        self.f = [0,0,0]                #Timoshenko_flag == 1:
        self.df = [0,0,0]               #Timoshenko_flag == 1:
        self.ddf =[0,0,0]               #Timoshenko_flag == 1:
        self.dddf =[0,0,0]              #Timoshenko_flag == 1:
        self.m = [0,0,0]                #Timoshenko_flag == 1:
        self.dm = [0,0,0]               #Timoshenko_flag == 1:
        self.ddm =[0,0,0]               #Timoshenko_flag == 1:
        self.dddm =[0,0,0]              #Timoshenko_flag == 1:
        self.gamma11 = [0]              #Vlasov_flag == 1:
        self.kappa = [0,0,0]            #Vlasov_flag == 1:
        self.dkappa1 = 0                #Vlasov_flag == 1:
        self.ddkapp1 = 0                #Vlasov_flag == 1:
        
    
def export_cells_for_VABS(cells,nodes,filename,VABSsetup,MaterialLst):
    f = open(filename,'w+')
    
    try:
        #PRINT HEADER!
        f.write('%i %i\n' % (VABSsetup.format_flag, VABSsetup.nlayer))
        f.write('%i %i %i\n' % (VABSsetup.Timoshenko_flag, VABSsetup.recover_flag, VABSsetup.thermal_flag))
        f.write('%i %i %i %i\n' % (VABSsetup.curve_flag,VABSsetup.oblique_flag,VABSsetup.trapeze_flag,VABSsetup.Vlasov_flag))
        if VABSsetup.curve_flag == 1:
            f.write(('%.2f %.2f %.2f\n') % (VABSsetup.k1,VABSsetup.k2,VABSsetup.k3))
        if VABSsetup.oblique_flag == 1:
            f.write('%.2f %.2f\n' % (VABSsetup.oblique_cosine1,VABSsetup.oblique_cosine2))
        f.write('\n')
        
       
        #Number of Nodes,Cells and Materials
        f.write('%i\t%i\t%i\n' % (len(nodes),len(cells),len(MaterialLst)))
        f.write('\n')
        #Node number, coordinates x_2, coordinatex x_3
        for n in nodes:
            f.write('%i\t\t%.6f\t%.5f\n' % (n.id,n.coordinates[0],n.coordinates[1]))
        f.write('\n')
        #Element number, connectivity
        for c in cells:
            f.write('%i\t\t' % (c.id))
            for i in range(0,9):
                if i<len(c.nodes):
                    f.write('%i\t' % (c.nodes[i].id))  
                else:
                    f.write('%i\t' % (0))
            f.write('\n')   
        f.write('\n')
        #Element number, Layup orientation
        for c in cells:
            f.write('%i\t\t%i\t%i\t' % (c.id,c.MatID,c.theta_3))
            for t in c.theta_1:
                f.write('%.2f\t' % (t))
            f.write('\n')  
        f.write('\n')
        #Materials 
        for m in MaterialLst:
            f.write('%i, %i\n' % (m.id,m.orth))
            if m.orth == 0:
                f.write('%.2f %.2f\n' % (m.E,m.nu))    
                f.write('%.6f\n' % (m.rho))
                if VABSsetup.thermal_flag == 1:
                    f.write('%.2f\n' % (m.alpha))
                f.write('\n')
                    
            elif m.orth == 1:
                f.write('%.2f %.2f %.2f\n' % (m.E[0],m.E[1],m.E[2]))
                f.write('%.2f %.2f %.2f\n' % (m.G[0],m.G[1],m.G[2]))
                f.write('%.2f %.2f %.2f\n' % (m.nu[0],m.nu[1],m.nu[2]))
                f.write('%.6f\n' % (m.rho))
                if VABSsetup.thermal_flag == 1:
                    f.write('%.2f %.2f %.2f\n' % (m.alpha[0],m.alpha[1],m.alpha[2]))
                f.write('\n')
                
            elif m.orth == 2:
                f.write('%.2f %.2f %.2f %.2f %.2f %.2f\n' % (m.C[0][0],m.C[0][1],m.C[0][2],m.C[0][3],m.C[0][4],m.C[0][5]))
                f.write('\t  %.2f %.2f %.2f %.2f %.2f\n' % (          m.C[1][1],m.C[1][2],m.C[1][3],m.C[1][4],m.C[1][5]))
                f.write('\t  \t   %.2f %.2f %.2f %.2f\n' % (                    m.C[2][2],m.C[2][3],m.C[2][4],m.C[2][5]))
                f.write('\t  \t   \t   %.2f %.2f %.2f\n' % (                              m.C[3][3],m.C[3][4],m.C[3][5]))
                f.write('\t  \t   \t   \t   %.2f %.2f\n' % (                                        m.C[4][4],m.C[4][5]))
                f.write('\t  \t   \t   \t   \t   %.2f\n' % (                                                  m.C[5][5]))
                f.write('%.6f\n' % (m.rho)) 
                if VABSsetup.thermal_flag == 1:
                    f.write('%.2f %.2f %.2f %.2f %.2f %.2f\n' % (m.alpha[0],m.alpha[1],m.alpha[2],m.alpha[3],m.alpha[4],m.alpha[5]))
                f.write('\n')
        f.write('\n')
        
        #LOADCASE:
        if VABSsetup.recover_flag == 1:
            
            if VABSsetup.Timoshenko_flag == 0:
                f.write('%.2f %.2f %.2f\n' % (VABSsetup.u[0],VABSsetup.u[1],VABSsetup.u[2]))
                f.write('%.2f %.2f %.2f\n' % (VABSsetup.Cij[0][0],VABSsetup.Cij[0][1],VABSsetup.Cij[0][2]))
                f.write('%.2f %.2f %.2f\n' % (VABSsetup.Cij[1][0],VABSsetup.Cij[1][1],VABSsetup.Cij[1][2]))
                f.write('%.2f %.2f %.2f\n' % (VABSsetup.Cij[2][0],VABSsetup.Cij[2][1],VABSsetup.Cij[2][2]))
                f.write('%.2f %.2f %.2f %.2f\n' % (VABSsetup.F[0],VABSsetup.M[0],VABSsetup.M[1],VABSsetup.M[2]))                
 
            elif VABSsetup.Timoshenko_flag == 1:
                f.write('%.2f %.2f %.2f\n' % (VABSsetup.u[0],VABSsetup.u[1],VABSsetup.u[2]))
                f.write('%.2f %.2f %.2f\n' % (VABSsetup.Cij[0][0],VABSsetup.Cij[0][1],VABSsetup.Cij[0][2]))
                f.write('%.2f %.2f %.2f\n' % (VABSsetup.Cij[1][0],VABSsetup.Cij[1][1],VABSsetup.Cij[1][2]))
                f.write('%.2f %.2f %.2f\n' % (VABSsetup.Cij[2][0],VABSsetup.Cij[2][1],VABSsetup.Cij[2][2]))
                f.write('%.2f %.2f %.2f %.2f\n' % (VABSsetup.F[0],VABSsetup.M[0],VABSsetup.M[1],VABSsetup.M[2])) 
                f.write('%.2f %.2f\n' % (VABSsetup.F[1],VABSsetup.F[2]))
                f.write('%.2f %.2f %.2f %.2f %.2f %.2f\n' % (VABSsetup.f[0],VABSsetup.f[1],VABSsetup.f[2],VABSsetup.m[0],VABSsetup.m[1],VABSsetup.m[2]))
                f.write('%.2f %.2f %.2f %.2f %.2f %.2f\n' % (VABSsetup.df[0],VABSsetup.df[1],VABSsetup.df[2],VABSsetup.dm[0],VABSsetup.dm[1],VABSsetup.dm[2]))
                f.write('%.2f %.2f %.2f %.2f %.2f %.2f\n' % (VABSsetup.ddf[0],VABSsetup.ddf[1],VABSsetup.ddf[2],VABSsetup.ddm[0],VABSsetup.ddm[1],VABSsetup.ddm[2]))
                f.write('%.2f %.2f %.2f %.2f %.2f %.2f\n' % (VABSsetup.dddf[0],VABSsetup.dddf[1],VABSsetup.dddf[2],VABSsetup.dddm[0],VABSsetup.dddm[1],VABSsetup.dddm[2]))
                
            elif VABSsetup.Vlasov_flag ==1:
                f.write('%.2f %.2f %.2f\n' % (VABSsetup.u[0],VABSsetup.u[1],VABSsetup.u[2]))
                f.write('%.2f %.2f %.2f\n' % (VABSsetup.Cij[0][0],VABSsetup.Cij[0][1],VABSsetup.Cij[0][2]))
                f.write('%.2f %.2f %.2f\n' % (VABSsetup.Cij[1][0],VABSsetup.Cij[1][1],VABSsetup.Cij[1][2]))
                f.write('%.2f %.2f %.2f\n' % (VABSsetup.Cij[2][0],VABSsetup.Cij[2][1],VABSsetup.Cij[2][2]))
                f.write('%.2f %.2f %.2f %.2f\n' % (VABSsetup.F[0],VABSsetup.M[0],VABSsetup.M[1],VABSsetup.M[2]))    
                
        f.write('\n')
    finally:
        f.close()
        
        

        
class XSectionalProperties(object):
    def __init__(self,filename=None):
        if filename is None:
            pass
        elif filename is not None:
            self.read_VABS_K(filename)
            self.filename = filename
 
    def read_VABS_K(self,filename):    

        a = ''
        with open(filename) as f:
            for line in f:
                line = line.partition('#')[0]
                line = line.rstrip()
                a += line 
                a += '\n'
             
        STR = ''.join([s for s in a.strip().splitlines(True) if s.strip("\r\n").strip()])
    
        i_MM = STR.find('The 6X6 Mass Matrix')
        i_MC = STR.find('The Mass Center of the Cross Section')
        i_MMatMC = STR.find('The 6X6 Mass Matrix at the Mass Center')   
        i_MP = STR.find('The Mass Properties with respect to Principal Inertial Axes')
        i_GC = STR.find('The Geometric Center of the Cross Section')
        i_CS = STR.find('Classical Stiffness Matrix (1-extension; 2-twist; 3,4-bending)')
        i_CF = STR.find('Classical Flexibility Matrix (1-extension; 2-twist; 3,4-bending)')
        i_NA = STR.find('The Neutral Axes (or Tension Center) of the Cross Section')
        i_TS = STR.find('Timoshenko Stiffness Matrix (1-extension; 2,3-shear, 4-twist; 5,6-bending)')
        i_TF =STR.find('Timoshenko Flexibility Matrix (1-extension; 2,3-shear, 4-twist; 5,6-bending)')
        i_GSC = STR.find('The Generalized Shear Center of the Cross Section in the User Coordinate System')        
        i_VS = STR.find('Vlasov Stiffness Matrix (1-extension; 2-twist; 3,4-bending; 5-twist rate)')
        i_VF = STR.find('Vlasov Flexibility Matrix (1-extension; 2-twist; 3,4-bending; 5-twist rate)'
                        )
        # Trapeze Effects Related Matrices
        i_Ag = STR.find('Ag1--Ag1--Ag1--Ag1')
        i_Bk = STR.find('Bk1--Bk1--Bk1--Bk1')
        i_Ck = STR.find('Ck2--Ck2--Ck2--Ck2')
        i_Dk = STR.find('Dk3--Dk3--Dk3--Dk3')
        
        splitpoints = [i_MM,i_MC,i_MMatMC,i_MP,i_GC,i_CS,i_CF,i_NA,i_TS,i_TF,i_GSC,i_VS,i_VF,i_Ag,i_Bk,i_Ck,i_Dk]
        splitpoints.sort()
        
        def grab_str_segment(STR,idx,splitpoints):
            k = splitpoints.index(idx)
            if k+1<len(splitpoints):
                i_end = splitpoints[k+1]
                return STR[idx:i_end]
            else:
                return STR[idx:]
            
        def read_PIA(STR,STR2Find):
            Start,End = read_rowstring(STR,STR2Find)
            temp = STR[Start:]
            temp = temp.split('\n')[1]
            temp = temp.replace(" ", "")
            return float(temp)
            
        MM = grab_str_segment(STR,i_MM,splitpoints)
        MC = grab_str_segment(STR,i_MC,splitpoints)
        MMatMC = grab_str_segment(STR,i_MMatMC,splitpoints)
        MP = grab_str_segment(STR,i_MP,splitpoints)
        GC = grab_str_segment(STR,i_GC,splitpoints)
        CS = grab_str_segment(STR,i_CS,splitpoints)
        CF = grab_str_segment(STR,i_CF,splitpoints)
        NA = grab_str_segment(STR,i_NA,splitpoints)
        TS = grab_str_segment(STR,i_TS,splitpoints)
        TF = grab_str_segment(STR,i_TF,splitpoints)
        GSC = grab_str_segment(STR,i_GSC,splitpoints)
        VS = grab_str_segment(STR,i_VS,splitpoints)
        VF = grab_str_segment(STR,i_VF,splitpoints)
        Ag = grab_str_segment(STR,i_Ag,splitpoints)
        Bk = grab_str_segment(STR,i_Bk,splitpoints)
        Ck = grab_str_segment(STR,i_Ck,splitpoints)
        Dk = grab_str_segment(STR,i_Dk,splitpoints)
        
        #The 6X6 Mass Matrix:
        self.MM = np.loadtxt(StringIO(MM),skiprows=2)
        #The Mass Center of the Cross Section:
        self.Xm2 = read_FLOATrowSTR(MC,'Xm2')
        self.Xm3 = read_FLOATrowSTR(MC,'Xm3')
        #The 6X6 Mass Matrix at the Mass Center:
        self.MMatMC = np.loadtxt(StringIO(MMatMC),skiprows=2)
        #The Mass Properties with respect to Principal Inertial Axes:
        self.MpUS = read_FLOATrowSTR(MP,'Mass Per Unit Span')
        self.I1 = read_FLOATrowSTR(MP,'Mass Moments of Intertia about x1 axis')
        self.I2 = read_FLOATrowSTR(MP,'Mass Moments of Intertia about x2 axis')
        self.I3 = read_FLOATrowSTR(MP,'Mass Moments of Intertia about x3 axis')
        
        try:
            self.PIA = read_PIA(MP, 'The Principal Inertial Axes Rotated from User Coordinate System by')
        except:
            if  'The user coordinate axes are the principal inertial axes.' in MP:
                self.PIA = float(0)
                
        self.Rg = read_FLOATrowSTR(MP,'The mass-weighted radius of gyration')
        #The Geometric Center of the Cross Section:
        self.Xg2 = read_FLOATrowSTR(GC,'Xg2')
        self.Xg3 = read_FLOATrowSTR(GC,'Xg3')
        # The Neutral Axes (or Tension Center) of the Cross Section:
        self.Xt2 = read_FLOATrowSTR(NA,'Xt2')
        self.Xt3 = read_FLOATrowSTR(NA,'Xt3')
        
        
        # Classical Stiffness and Flexibility Matrix (1-extension; 2-twist; 3,4-bending)
        try:
            self.CS = np.loadtxt(StringIO(CS),skiprows=2)
            self.CF = np.loadtxt(StringIO(CF),skiprows=2)
        except:
            self.CS =None
            self.CF =None
    
        # Timoshenko Stiffness and Flexibility Matrix (1-extension; 2,3-shear, 4-twist; 5,6-bending)
        try:
            self.TS = np.loadtxt(StringIO(TS),skiprows=2)
            self.TF = np.loadtxt(StringIO(TF),skiprows=2)
        except:
            self.TS =None
            self.TF =None
                        
        #The Generalized Shear Center of the Cross Section in the User Coordinate System    
        try:
            self.Xs2 = read_FLOATrowSTR(GSC,'Xs2')
            self.Xs3 = read_FLOATrowSTR(GSC,'Xs3')
        except:
            self.Xs2 =None
            self.Xs3 =None
        
        #Vlasov Stiffness and Flexibility Matrix (1-extension; 2-twist; 3,4-bending; 5-twist rate)
        try:
            self.VS = np.loadtxt(StringIO(VS),skiprows=2)
            self.VF = np.loadtxt(StringIO(VF),skiprows=2)
        except:
            self.VS =None
            self.VF =None
        
        #Trapeze Effects Related Matrices
        try:    
            self.Ag = np.loadtxt(StringIO(Ag),skiprows=2)
            self.Bk = np.loadtxt(StringIO(Bk),skiprows=2)
            self.Ck = np.loadtxt(StringIO(Ck),skiprows=2)
            self.Dk = np.loadtxt(StringIO(Dk),skiprows=2)
        except:
            self.Ag=None
            self.Bk=None
            self.Ck=None
            self.Dk=None
            

    def read_all_VABS_Results(self):
        self.ELE = read_VABS_Results(self.filename.replace('.K','.ELE'))
        self.U = read_VABS_Results(self.filename.replace('.K','.U'))

           
#======================================================
#       MAIN
#======================================================       
if __name__ == '__main__':
    
    test = VABS_config(recover_flag=1)
    filename = 'naca0012_cspar_mesh.vab.K'
    BeamProperties = XSectionalProperties(filename)
    filename.replace('.K','.ELE')
    BeamProperties.read_VABS_ELE(filename)
