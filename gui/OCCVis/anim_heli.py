# from OCC.Display.SimpleGui import init_display
from objviewer import init_display
from OCC.TopoDS import TopoDS_Shape
from OCC.StlAPI import StlAPI_Reader
from PyQt5.QtWidgets import QSlider, QLabel
from PyQt5.QtCore import Qt
from OCC.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCC.gp import gp_Ax1, gp_Pnt, gp_Dir, gp_Vec, gp_Trsf
from OCC.TopLoc import TopLoc_Location
import os
# import sys
from numpy import deg2rad as d2r
from numpy import cos, sin, pi
import numpy as np
import pandas as pd
import pickle
import re

# set usestl to 1 to load stl file for helicopter representation,
# otherwise box will be used
usestl = 0


def rotvec(invec, ang, axis):
    if(axis == 'x'):
        rotmat = np.array([[1.0, 0.0, 0.0],
                           [0.0, cos(ang), sin(ang)],
                           [0.0, -sin(ang), cos(ang)]])
    if(axis == 'y'):
        rotmat = np.array([[cos(ang), 0.0, -sin(ang)],
                           [0.0, 1.0, 0.0],
                           [sin(ang), 0.0, cos(ang)]])
    if(axis == 'z'):
        rotmat = np.array([[cos(ang), sin(ang), 0.0],
                           [-sin(ang), cos(ang), 0.0],
                           [0.0, 0.0, 1.0]])
    return rotmat.dot(invec)


# ---------- Move object and plot indicator on slider event ----------

def move_obj():
    # Transform 3d helicopter model
    ref_pnt = np.array([1.0, 0.5, 0.5])
    idx = int(sl.value())
    trans_pnt = np.array(
        [data.loc[idx]['X_e(0)'], data.loc[idx]['X_e(1)'], data.loc[idx]['X_e(2)']])
    Psi = d2r(data.loc[idx]['alpha(2)'])
    Theta = d2r(data.loc[idx]['alpha(1)'])
    Phi = d2r(data.loc[idx]['alpha(0)'])
    x0 = np.array([1.0, 0.0, 0.0])
    y0 = np.array([0.0, 1.0, 0.0])
    z0 = np.array([0.0, 0.0, 1.0])
    # x1 = rotvec(x0, Psi, 'z')
    # y1 = rotvec(y0, Psi, 'z')
    # x2 = rotvec(x1, Theta, 'y')
    # ax1 = gp_Ax1(gp_Pnt(*(trans_pnt + ref_pnt)), gp_Dir(*z0))
    # ax2 = gp_Ax1(gp_Pnt(*(trans_pnt + ref_pnt)), gp_Dir(*y0))
    # ax3 = gp_Ax1(gp_Pnt(*(trans_pnt + ref_pnt)), gp_Dir(*x0))
    ax1 = gp_Ax1(gp_Pnt(*ref_pnt), gp_Dir(*z0))
    ax2 = gp_Ax1(gp_Pnt(*ref_pnt), gp_Dir(*y0))
    ax3 = gp_Ax1(gp_Pnt(*ref_pnt), gp_Dir(*x0))
    rot1.SetRotation(ax1, Psi)
    rot2.SetRotation(ax2, Theta)
    rot3.SetRotation(ax3, Phi)
    roty.SetRotation(ax3, - pi / 2)
    rotx.SetRotation(ax2, pi / 2)
    obj_trans.SetTranslation(gp_Vec(*trans_pnt))
    if usestl:
        objTopLoc = TopLoc_Location(
            obj_trans * rot1 * rot2 * rot3 * roty * rotx)
    else:
        objTopLoc = TopLoc_Location(obj_trans * rot1 * rot2 * rot3)
    display.Context.SetLocation(obj, objTopLoc)
    display.Context.UpdateCurrentViewer()
    #
    # update plot time indicator
    dt = data.iloc[2][0] - data.iloc[1][0]
    items = enumerate(
        zip(win.plotwidget.vlines, win.plotwidget.axs, win.plotwidget.bgs))
    for j, (line, ax, bg) in items:
        win.plotwidget.fig.canvas.restore_region(bg)
        line.set_xdata([dt * idx, dt * idx])
        ax.draw_artist(line)
        win.plotwidget.fig.canvas.blit(ax.bbox)
    # update time in toolbar
    lbl.setText("%10.2f" % (dt * idx))


# sys.stdout = os.devnull


# def setcam():
#     display.Rotation(sl2.value(), sl3.value())
#     print(sl2.value(), sl3.value())


# ---------- Load plot data ----------
def is_non_zero_file(fpath):
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0


keylist = []
fname = 'plotfile.dat'
if is_non_zero_file(fname):
    # get number of columns
    cols = pd.read_csv('plotfile.dat', delimiter=',', nrows=1,
                       header=0, skipinitialspace=True).columns
    # read all but last (empty) column from csv file
    data = pd.read_csv('plotfile.dat', delimiter=',', usecols=cols[:-1],
                       header=0, skipinitialspace=True)

    # group header names by data vectors
    vregex = '(.*)\(([0-9]*)\)'
    headers = data.columns.values
    reqset = ['X_e(0)', 'X_e(1)', 'X_e(2)', 'alpha(0)', 'alpha(1)', 'alpha(2)']
    if(not set(reqset).issubset(set(headers))):
        exit('Script requires the following variables in plotfile: ' + str(reqset))
    tmp = []
    idx = 1
    while idx < len(headers):
        prev_elem = headers[idx - 1]
        cur_elem = headers[idx]
        tmp.append(prev_elem)
        # regex match everything in front of ...(#) if present
        pm = re.search(vregex, prev_elem)
        m = re.search(vregex, cur_elem)
        # check for two consecutive matches
        if pm and m:
            # if regex group matches previous keep collecting items in tmp
            # otherwise add tmp to keylist and advance
            if pm.group(1) == m.group(1):
                next
            else:
                keylist.append(tmp)
                tmp = []
        else:
            keylist.append(tmp)
            tmp = []
        idx += 1
        # add last element to tmp and keylist
        if idx == len(headers):
            tmp.append(cur_elem)
            keylist.append(tmp)

    # remove time from keylist
    keylist = keylist[1:]
else:
    exit('No data to plot, quitting...')


# ---------- Setup window and toolbar ----------

display, start_display, add_menu, add_function_to_menu, win = init_display()
# tmin = data.iloc[0][0]
# tmax = data.iloc[-1][0]
# dt = data.iloc[1][0] - data.iloc[1][0]
sl = QSlider(Qt.Horizontal)
sl.setMinimum(0)
sl.setMaximum(data.shape[0] - 1)
sl.setValue(0)
# sl.setTickPosition(QSlider.TicksBelow)
sl.setTickInterval(1)
sl.valueChanged.connect(move_obj)
tb = win.addToolBar("File")
tb.addWidget(sl)
lbl = QLabel()
lbl.setText("0.0")
lbl.setFixedWidth(50)
lbl.setAlignment(Qt.AlignRight)
tb.addWidget(lbl)
lbl_dummy = QLabel()
tb.addWidget(lbl_dummy)
# sl2 = QSlider(Qt.Horizontal)
# sl2.setMinimum(-2000)
# sl2.setMaximum(2000)
# sl2.setValue(0)
# sl2.setTickPosition(QSlider.TicksBelow)
# sl2.setTickInterval(1)
# sl2.valueChanged.connect(setcam)
# tb.addWidget(sl2)
# sl3 = QSlider(Qt.Horizontal)
# sl3.setMinimum(-2000)
# sl3.setMaximum(2000)
# sl3.setValue(0)
# sl3.setTickPosition(QSlider.TicksBelow)
# sl3.setTickInterval(1)
# sl3.valueChanged.connect(setcam)
# tb.addWidget(sl3)


# ---------- Load helicopter 3D model ----------

pfile = 'stlpickle.dat'
if usestl:
    if not os.path.isfile(pfile):
        stl_reader = StlAPI_Reader()
        fan_shp = TopoDS_Shape()
        # stl_reader.Read(fan_shp, './Flettner.stl')
        stl_reader.Read(fan_shp, './AREA_FULL.stp')
        # stl_reader.Read(fan_shp, './a-109.stp')
        P = gp_Pnt(0, 0, 0)
        factor = 0.001
        aTrsf = gp_Trsf()
        aTrsf.SetScale(P, factor)
        aBRespTrsf = BRepBuilderAPI_Transform(fan_shp, aTrsf)
        fan_shp = aBRespTrsf.Shape()
        pickle.dump(fan_shp, open(pfile, "wb"))
    else:
        fan_shp = pickle.load(open(pfile, "rb"))
    obj = display.DisplayShape(fan_shp, update=True)
else:
    box = BRepPrimAPI_MakeBox(2.0, 1.0, 1.0).Shape()
    obj = display.DisplayShape(box, update=True)

rotx = gp_Trsf()
roty = gp_Trsf()
rot1 = gp_Trsf()
rot2 = gp_Trsf()
rot3 = gp_Trsf()
obj_trans = gp_Trsf()

# Ground plate
obj0 = display.DisplayColoredShape(BRepPrimAPI_MakeBox(
    10.0, 10.0, 0.1).Shape(), 'BLACK', update=True)
ground_trans = gp_Trsf()
ground_trans.SetTranslation(gp_Vec(-5.0, -5.0, -0.1))
objTopLoc0 = TopLoc_Location(ground_trans)
display.Context.SetLocation(obj0, objTopLoc0)


# ---------- Plot data ----------

plt = win.plotwidget.fig
cnt = 0

print(keylist)
for grp in keylist:
    win.plotwidget.axs.append(plt.add_subplot(len(keylist), 1, cnt + 1))
    for ln in grp:
        win.plotwidget.axs[cnt].plot(data['Time'].tolist(), data[ln].tolist())
    if len(grp) > 1:
        m = re.search(vregex, grp[0])
        win.plotwidget.axs[cnt].set_ylabel(m.group(1))
    else:
        win.plotwidget.axs[cnt].set_ylabel(grp[0])
    cnt += 1


win.plotwidget.axs[cnt - 1].set_xlabel('time')

plt.subplots_adjust(hspace=0.4)
plt.tight_layout()
win.plotwidget.draw()


# ---------- Start program ----------

move_obj()
display.FitAll()
display.Rotation(0, 0)
display.Rotation(-1488, -492)
display.View.SetZoom(0.8)
start_display()

# sys.stdout = sys.__stdout__
