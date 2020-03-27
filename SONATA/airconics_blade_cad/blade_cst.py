import yaml
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares
from scipy.interpolate import splrep, splev, interp1d

# import primitives
# import liftingSurface
# import AirCONICStools as act
from OCC.Display.SimpleGui import init_display
from OCC.Core.Geom import Geom_TrimmedCurve

import sys
sys.path.insert(0,'/Users/rfeil/work/6_SONATA/SONATA')  # import sys path to import 'SONATA' & 'job' modules
import SONATA.airconics_blade_cad.airconics.primitives as primitives
import SONATA.airconics_blade_cad.airconics.liftingSurface as liftingSurface
import SONATA.airconics_blade_cad.airconics.AirCONICStools as act

# class TurbineCAD:

"""
Airfoil 'class function / shape function' transformation technique (cst) used to plot 'smooth' 3D lofted shapes of 
arbitrary geometries, such as beams or blades

Inputs: blade component based on yaml file input

Outputs: Component 3D shape

"""
class blade_cst(object):

    def __init__(self, yml, flags_dict):

        self.chord_grid = yml['components']['blade']['outer_shape_bem']['chord']['grid']
        self.chord = yml['components']['blade']['outer_shape_bem']['chord']['values']
        self.pitch_axis = yml['components']['blade']['outer_shape_bem']['pitch_axis']['values']
        self.pitch_axis_grid = yml['components']['blade']['outer_shape_bem']['pitch_axis']['grid']

        if flags_dict['flag_wt_ontology']:  # WT coords
            self.twist_grid = yml['components']['blade']['outer_shape_bem']['twist']['grid']
            self.twist = yml['components']['blade']['outer_shape_bem']['twist']['values']
            self.prebend_grid = yml['components']['blade']['outer_shape_bem']['reference_axis']['x']['grid']
            self.prebend = yml['components']['blade']['outer_shape_bem']['reference_axis']['x']['values']

            self.ref_axis = {}
            self.ref_axis['x_values'] = yml['components']['blade']['outer_shape_bem']['reference_axis']['x']['values']
            self.ref_axis['x_grid'] = yml['components']['blade']['outer_shape_bem']['reference_axis']['x']['grid']
            self.ref_axis['y_values'] = yml['components']['blade']['outer_shape_bem']['reference_axis']['y']['values']
            self.ref_axis['y_grid'] = yml['components']['blade']['outer_shape_bem']['reference_axis']['y']['grid']
            self.ref_axis['z_values'] = yml['components']['blade']['outer_shape_bem']['reference_axis']['z']['values']
            self.ref_axis['z_grid'] = yml['components']['blade']['outer_shape_bem']['reference_axis']['z']['grid']

        else:   # helicopter sys HT (correct signs and coord sys to fit with WT definition from BeamDyn)
            self.twist_grid = yml['components']['blade']['outer_shape_bem']['twist']['grid']
            self.twist = list(-np.array(yml['components']['blade']['outer_shape_bem']['twist']['values']))    # sign correction
            self.prebend_grid = yml['components']['blade']['outer_shape_bem']['reference_axis']['z']['grid']
            self.prebend = yml['components']['blade']['outer_shape_bem']['reference_axis']['z']['values']

            self.ref_axis = {}
            self.ref_axis['x_grid'] = yml['components']['blade']['outer_shape_bem']['reference_axis']['z']['grid']
            self.ref_axis['x_values'] = yml['components']['blade']['outer_shape_bem']['reference_axis']['z']['values']
            self.ref_axis['y_grid'] = yml['components']['blade']['outer_shape_bem']['reference_axis']['y']['grid']
            self.ref_axis['y_values'] = list(-np.array(yml['components']['blade']['outer_shape_bem']['reference_axis']['y']['values']))
            self.ref_axis['z_grid'] = yml['components']['blade']['outer_shape_bem']['reference_axis']['x']['grid']
            self.ref_axis['z_values'] = yml['components']['blade']['outer_shape_bem']['reference_axis']['x']['values']

        self.airfoils = yml['airfoils']
        self.airfoil_labels = yml['components']['blade']['outer_shape_bem']['airfoil_position']['labels']
        self.airfoil_indices = range(np.size(self.airfoil_labels))
        self.airfoil_grid = yml['components']['blade']['outer_shape_bem']['airfoil_position']['grid']
        self.airfoil_map = np.zeros(np.size(self.airfoil_labels))  # Map to an airfoil in the airfoils table

        for i, al in enumerate(self.airfoil_labels):
            for j, af in enumerate(self.airfoils):
                try:
                    if (af['name'] == al):
                        self.airfoil_map[i] = j
                #                        print("Airfoil {} - name {} in BEM blade shape is mapped to airfoil {} - name {} in airfoils list ".format(i,al,j,af['name']))
                except:
                    pass


    def twist_func(self,eta):
        """User-defined function describing the variation of twist as a function
        of the distance along pitch axis."""
        return np.interp(eta, self.twist_grid, self.twist)

    def chord_func(self,eta):
        """User-defined function describing the variation of chord as a function of
        the distance along pitch axis"""
        return np.interp(eta, self.chord_grid, self.chord)

    def prebend_func(self,eta):
        """User-defined function describing the variation of chord as a function of
        the distance along pitch axis"""
        return np.interp(eta, self.prebend_grid, self.prebend)

    def pitch_axis_func(self, eta):
        """User-defined function describing the variation of chordwise pitch
        axis location as a function of the distance along pitch axis."""
        return np.interp(eta, self.pitch_axis_grid, self.pitch_axis)

    def cst_u_l_shape(self,sf_wt,x_ul,n1,n2,zetaT,BP,K):
        C = [(x**n1)*((1.0-x)**n2) for x in x_ul]
        S = np.array([sf_wt[0]*K[0]*((1.0-x)**BP) for i,x in enumerate(x_ul)])
        for j in range(1,BP+1):
            S += np.array([sf_wt[j]*K[j]* (x_ul[i]**j) * ((1.0-x_ul[i])**(BP-j)) for i,x in enumerate(x_ul)])
        return np.array([C[i]*S[i] + (zetaT*x) for i,x in enumerate(x_ul)])

    def af_cst_diff(self,cst_wt, x_l, x_u, y_true, BP, K, N1, N2):
        y_l = self.cst_u_l_shape(cst_wt[0:BP+1], x_l, N1, N2, cst_wt[2*(BP+1)], BP, K)
        y_u = self.cst_u_l_shape(cst_wt[BP+1:2*(BP+1)], x_u, N1, N2, cst_wt[2*(BP+1)+1], BP, K)
        return np.subtract(np.append(y_l, y_u),y_true)

    def plot_af_cst(self,cst_wt, x_l, x_u, y_true, BP, K, N1, N2):

        y_l = self.cst_u_l_shape(cst_wt[0:BP+1], x_l, N1, N2, cst_wt[2*(BP+1)], BP, K)
        y_u = self.cst_u_l_shape(cst_wt[BP+1:2*(BP+1)], x_u, N1, N2, cst_wt[2*(BP+1)+1], BP, K)
        y_pred = np.append(y_l,y_u)
        print(np.linalg.norm(y_pred-y_true))

        fig = plt.figure()
        y_l = self.cst_u_l_shape(cst_wt[0:BP+1], x_l, N1, N2, cst_wt[2*(BP+1)], BP, K)
        y_u = self.cst_u_l_shape(cst_wt[BP+1:2*(BP+1)], x_u, N1, N2, cst_wt[2*(BP+1)+1], BP, K)
        x = np.append(x_l,x_u)
        y_pred = np.append(y_l,y_u)
        plt.plot(x,y_pred,label='Pred')
        plt.plot(x,y_true,"+-",label='Orig')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.legend(loc=0)
        plt.show()
        plt.close(fig)


    def create_af_cst(self, af_name, BP=6):
        af_index = None
        # for i,af in enumerate(self.turbine['airfoils']):
        for i,af in enumerate(self.airfoils):
            if (af['name'] == af_name):
                af_index = i
        psi_true = np.array(self.airfoils[af_index]['coordinates']['x'])
        zeta_true = np.array(self.airfoils[af_index]['coordinates']['y'])
        le_loc = np.argmin(psi_true)

        psi_l = psi_true[:le_loc]
        psi_u = psi_true[le_loc:]

        K = np.array([np.math.factorial(BP)/(np.math.factorial(i)*(np.math.factorial(BP-i))) for i in range(BP+1)])
        cst_wt_0 = np.zeros( 2*(BP+1)+2)
        cst_wt_0[0:BP+1] = -1.0
        cst_wt_0[BP+1:2*(BP+1)] = 1.0

        N1 = 0.5
        N2 = 1.0
        if ('Cylinder' in af_name):
            N2 = 0.5

        af_cst = least_squares(self.af_cst_diff, cst_wt_0, jac='3-point',ftol=1e-12,args=(psi_l, psi_u, zeta_true, BP, K, N1, N2))

        #self.plot_af_cst(af_cst.x, psi_l, psi_u, zeta_true, BP, K, N1, N2)
        return np.append(af_cst.x, [N1, N2])

    def create_blade_cst(self, BPpsi=6):
        af_list = self.airfoil_labels
        # Define each airfoils with 17 values, i.e. 7 for the 6th order polynomial of each the upper and the lower surface, plus TE and beam axis location
        af_cst = np.array([ self.create_af_cst(af,BPpsi) for af in af_list ])
        # print(af_cst)
        eta_true = self.airfoil_grid
        Kpsi = np.array([np.math.factorial(BPpsi)/(np.math.factorial(i)*(np.math.factorial(BPpsi-i))) for i in range(BPpsi+1)])
        self.bld_cst = [interp1d(eta_true, af_cst[:,i],fill_value="extrapolate") for i in range( 2*(BPpsi+1)+4) ]

    def dump_all_af_cst(self):

        af_list = self.airfoil_labels[1:-1]
        for af in af_list:
            print(af)
        af_cst = np.array([ self.create_af_cst(af,8) for af in af_list])
        yaml_filename="airfoil_cst_db.yaml"
        af_cst_yaml = {
            'airfoil_cst_db': {af: np.ndarray.tolist( np.append(np.append( af_cst[i][8+1:2*(8+1)], af_cst[i][0:8+1]) ,[af_cst[i][-2], af_cst[i][-1]])  )  for i,af in enumerate(af_list)}
        }
        yaml.dump(af_cst_yaml, open(yaml_filename,'w'),default_flow_style=False)

    def airfoil_func(self,Epsilon,BPpsi=6):
        """Defines the variation of cross section as a function of Epsilon"""
        Kpsi = np.array([np.math.factorial(BPpsi)/(np.math.factorial(i)*(np.math.factorial(BPpsi-i))) for i in range(BPpsi+1)])
        N1 = 0.5
        N2 = 1.0
        af_cst = np.zeros( 2*(BPpsi+1)+2)

        # print(Epsilon)  # print to screen the individual radial stations (epsilon)
        N1 = (self.bld_cst[-2])(Epsilon)
        N2 = (self.bld_cst[-1])(Epsilon)
        for i in range( 2*(BPpsi+2) ):
            af_cst[i] = (self.bld_cst[i])(Epsilon)

        delta = 0.005
        theta = np.arange(0,np.pi+delta,np.pi*delta)

        y_u = np.cos(theta)*0.5+0.5
        y_l = -np.cos(theta)*0.5+0.5
        x_l = self.cst_u_l_shape(af_cst[0:BPpsi+1], y_l, N1, N2, -0.001, BPpsi, Kpsi)
        x_u = self.cst_u_l_shape(af_cst[BPpsi+1:2*(BPpsi+1)], y_u, N1, N2, 0.001, BPpsi, Kpsi)

        y_u = -y_u + self.pitch_axis_func(Epsilon)
        y_l = -y_l + self.pitch_axis_func(Epsilon)

        y = np.append(y_u,y_l[1:])
        x = np.append(x_u,x_l[1:])

        twist = self.twist_func(Epsilon)
        chord = self.chord_func(Epsilon)
        prebend = self.prebend_func(Epsilon)
        x_p = (x * np.cos(twist) - y * np.sin(twist)) * chord + prebend
        y_p = (x * np.sin(twist) + y * np.cos(twist)) * chord
        # import matplotlib.pyplot as plt
        # plt.plot(y_p, x_p_n)
        z_p = np.ones(np.shape(x_p)) * np.interp(Epsilon, self.ref_axis['z_grid'],self.ref_axis['z_values'])

        Af = primitives.Airfoil(np.column_stack([x_p,y_p,z_p]))
        return Af


if __name__ == "__main__":

    folder_str = '/Users/rfeil/work/6_SONATA/SONATA/jobs/RFeil/'
    job_str = 'IEA-15-240-RWT_TipShape_V1.yaml'

    yaml_file = folder_str + job_str
    with open(yaml_file,'r') as f:
        yml = yaml.load(f.read(), Loader=yaml.Loader)

    flags_dict = {}
    flag_wt_ontology = True     # define in which ontology (coord sys definition) your yaml file is
    flags_dict = {"flag_wt_ontology": flag_wt_ontology}

    CAD_shape = blade_cst(yml, flags_dict)
    CAD_shape.create_blade_cst()

    NSegments = 100
    blade = liftingSurface.LiftingSurface(CAD_shape.airfoil_func, NSegments=NSegments)
    print(["Finished creating loft using " + str(NSegments) + " segments"])
    blade.Write(job_str[:-5] + '.iges')  # writes step, stl, igs, or iges files depending on chosen extension

    # nrel5mw = TurbineCAD('example/nrel5mw_mod_update.yaml')
    # #nrel5mw.dump_all_af_cst()
    # nrel5mw.create_blade_cst()
    # # display, start_display, add_menu, add_function_to_menu = init_display()
    # # af_list = [nrel5mw.airfoil_func(eps) for eps in np.arange(0.0,1.01,0.05)]
    # # surf = act.AddSurfaceLoft(af_list)
    # # display.DisplayShape(surf, update=True)
    # # start_display()
    #
    # blade = liftingSurface.LiftingSurface(nrel5mw.airfoil_func, NSegments=100)
    # print("Finished creating loft")
    # status = blade.Write('blade.iges')
    # print(status)
    #
    # fig,ax = plt.subplots()
    # #ax = fig.gca(projection='3d')
    # for eta in np.arange(0.0,1.01,0.05):
    #     nrel5mw.blade_af(eta, ax=ax)
    # plt.show()
