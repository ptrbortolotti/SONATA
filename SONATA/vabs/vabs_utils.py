# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 10:33:12 2018

@author: TPflumm
"""

import numpy as np
import pint
import re


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
