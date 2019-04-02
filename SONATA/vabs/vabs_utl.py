# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 10:33:12 2018

@author: TPflumm
"""
from io import StringIO
import numpy as np
import pint
import re
from SONATA.cbm.fileIO.readinput import read_rowstring

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
    

def export_cells_for_VABS(cells, nodes, filename, VABSsetup, materials, rotation=0):
    """
    the export_cells for VABS function gathers the information of the 
    parameters to write the VABS input file!
    
             
    Parameters
    ----------
    cells : list
        list of cells
    nodes : list
        the corresponding list of nodes of the cells
    filename : str
        filename of the file to write
    VABSsetup : VABSconfig 
        VabsConfig object
    rotation : float
        TBD: rotation angle in rad which the mesh is rotated.
        
    """
       
    with open(filename,'w+') as f:
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
            f.write('%i\t%i\t%i\n' % (len(nodes),len(cells),len(materials)))
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
            for m in materials.values():
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
        except:
            print('ERROR in export_cells_for_VABS')
    return None


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


def assign_units_to_matrix(MMunits, dct):
    MMin =  []
    for i in MMunits:
        tmp_lst = []
        for j in i:
            s = str(j)
            for k in dct:
                s = s.replace(k, dct[k])
            tmp_lst.append(s)
        MMin.append(tmp_lst)
    return np.asarray(MMin)


def assign_pintunits_to_matrix(MMunits,dct):
    ureg = pint.UnitRegistry()
    MMin = assign_units_to_matrix(MMunits,dct)
    M_pint = []
    for i,row in enumerate(MMin):
        tmp_lst = []
        for j,item in enumerate(row):
            s = str(item)
            s = "".join(s.split())
            if s != '':
                code = '1 * '
                slst = re.split("([+-/*])", s.replace(" ", ""))
                #print slst
                for x in slst:
                    #print x
                    if re.match('^[\w-]+$', x):
                        code += 'ureg.'+x
                        
                    else:
                        code += ' '+x+' '
                #print code
                tmp_pint_from = eval(code)
                tmp_lst.append(tmp_pint_from)
        
            else: 
                tmp_lst.append(0)
        M_pint.append(tmp_lst)
    return M_pint


def transfer_matrix_unitconvertion(MMunits,in_dct,out_dct):
    M_pint_in = assign_pintunits_to_matrix(MMunits,in_dct)
    M_pint_out = assign_pintunits_to_matrix(MMunits,out_dct)
    
    TM = []
    for i,row in enumerate(M_pint_in):
        tmp_lst = []
        for j,item in enumerate(row):
            x1 = item
            x2 = M_pint_out[i][j]
            if x1 != 0:
                tmp = x1.to(x2.units)
                tmp_magn = tmp.magnitude
            else:
                tmp_magn = 0
            tmp_lst.append(tmp_magn)
        TM.append(tmp_lst)
        
    return(np.asarray(TM))  
