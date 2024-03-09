import numpy as np 
import math
import sys
import os
sys.path.append(os.getcwd())
import method_of_characteristics.moc_mesh_engine as mesh_module

"""
Tests MOC operators for the case of a uniform flow through a duct with a single 
compression turn. Should return results which agree with analytical solution. 

Will need to comment out a lot of lines in post_process for plotting to work
"""

class Idl:
        def __init__(self, u_i, v_i, x_i, geom, nPts, GasProps):
            self.x = [x_i for _ in range (nPts)]
            self.y = list(np.linspace(geom.y_cowl(x_i), geom.y_centerbody(x_i), nPts))
            self.u = [u_i for _ in range(nPts)]
            self.v = [v_i for _ in range(nPts)]
            self.p0 = GasProps.p0
            self.cowlPoint = True 
            self.cbPoint = True
            self.p0_p0f =  1


class Geom:
        def __init__(self, deflec_deg, x_i):

            def y_x_upper(x):
                if x<1:
                    return 1
                else: 
                    return 1-math.tan(math.radians(deflec_deg))*(x-1)
            def dydx_upper(x):
                if x<1: 
                    return 0
                else:
                    return -math.tan(math.radians(deflec_deg))
            def y_x_lower(x):
                return 0.5
            def dydx_lower(x):
                return 0 
            
            self.y_cowl = y_x_upper
            self.y_centerbody = y_x_lower 
            self.dydx_cowl = dydx_upper
            self.dydx_centerbody = dydx_lower


class GasProps:
    def __init__(self):
        self.gam = 1.4
        self.R = 287.05
        self.T0 = 288.15
        self.a0 = math.sqrt(self.gam*self.R*self.T0)
        self.p0 = 101325


class Input: 

    def __init__(self, M_inf, geom, delta, pcTOL, gasProps):
        self.geom = geom 
        self.gasProps = gasProps 
        self.delta = delta
        self.pcTOL = pcTOL
        self.init_method = "STRAIGHT IDL" #hardcoded 
        self.M_inf = M_inf


def duct_2d_uniform(idl, inputObj, plotFile=None):
    """
    u_i: initial uniform x-velocity 
    v_i: initial uniform y-velocity
    geom: object containing functions for y(x) and dydx(x) for upper and 
        lower surfaces 
    nPts: number of points on initial vertical data line
    """
    geom = inputObj.geom
    mesh = mesh_module.Mesh(inputObj, idl, explicit_shocks=True)
    #mesh.generate_initial_mesh_from_straight_idl("neg")
    #try:
    #    mesh.compute_wall_to_wall_shock("neg", mesh.shockPts_frontside[-1])
    #except: 
    #    pass 
    mesh.compile_mesh_points()
    [pt.get_point_properties(mesh) for pt in mesh.meshPts]
    mesh.compile_wall_points(geom.y_cowl, geom.y_centerbody) 
    mesh.compute_local_mass_flow()

    if plotFile is not None: 
        import matplotlib.pyplot as plt 
        import post_processing.post_process as post_process
        import json
        plotDict = json.load(open(plotFile, 'r'))
        plotSettings = plotDict["global plot settings"]
        del plotDict["global plot settings"]
        
        class MainObjDummy: #dummy class for main just to get plotting working
            def __init__(self): 
                
                class obj:
                    def __init__(self): return 
                
                freestream = obj()
                freestream.mach = 2 

                inputs = obj()
                inputs.geom = geom
                inputs.geom.cowl_bounds = [0,2]
                inputs.geom.centerbody_bounds = [0,2]
                inputs.M_inf = inputObj.M_inf

                self.freestream = freestream 
                self.inputs = inputs 
                self.mesh = mesh
                self.idlObj = idl
                
        for key in plotDict.keys():
            subDict = plotDict[key] 
            #hand off subDict to the post processing module
            post_process.create_figure(subDict, MainObjDummy(), plotSettings)

        plt.figure()
        char = mesh.C_neg[-1]
        point_mach_list = [pt.mach for pt in char]
        point_thet_list = [math.degrees(math.atan(pt.v/pt.u)) for pt in char]
        #plt.plot([1.82160847,1.82160847],[0.5,1], linewidth=2, label="Analytical Solution", color="orange")
        plt.plot([-5,-5],[0.5,1], linewidth=2, label="Analytical Solution", color="orange")
        #plt.scatter(point_mach_list, [pt.y for pt in char], label="MOC")
        plt.scatter(point_thet_list, [pt.y for pt in char], label="MOC")

        plt.xlim([-5.1, -4.95])
        plt.xlabel("Downstream Flow Angle (deg)")
        plt.grid()
        #plt.legend()

        plt.show()

    return 


if __name__ == "__main__":

    gasProps = GasProps()
    x_i = 1
    u_i, v_i = 507.33, 0
    a = math.sqrt(gasProps.a0**2 - 0.5*(gasProps.gam-1)*(u_i**2 + v_i**2))
    M_inf = math.sqrt(u_i**2 + v_i**2)/a 
    print(f"freestream mach number: {M_inf}")
    nPts = 10
    delta = 0 #2D case
    pcTOL = 0.00001 #perecent change tolerance for iterative processes 
    geom = Geom(deflec_deg=5, x_i=x_i)
    inputObj = Input(M_inf, geom, delta, pcTOL, gasProps)

    idl = Idl(u_i, v_i, x_i, geom, nPts, gasProps)
    res = duct_2d_uniform(idl, inputObj, plotFile="validation_cases/mesh_plot.json")