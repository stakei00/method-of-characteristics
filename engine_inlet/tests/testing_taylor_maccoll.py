"""
This file runs unit tests on the Taylor_Maccoll.py module
"""
import sys 
import os
import unittest
import math
sys.path.append(os.getcwd()) #add path to taylor maccoll module
import aimcat.taylor_maccoll as taylor_maccoll

class Test_Taylor_Maccoll(unittest.TestCase):

    def test_TMC_Cone(self):
        """
        Expected Results from Compressible Aerodynamics Calculator 3.1 
        https://devenport.aoe.vt.edu/aoe3114/calc.html
        """
        #instantiate object: 
        cone_ang = math.radians(12.5)
        M_inf = 2.5
        
        class gasProps:
            def __init__(self, gam, R, T0): 
                self.gam, self.R, self.T0 = gam, R, T0

        gas = gasProps(1.4, 287.05, 288.15)

        cone_flow = taylor_maccoll.TaylorMaccoll_Cone(cone_ang, M_inf, gas)
        #Compare calculated shock angle with expected angle: 
        shock_ang_exp = math.radians(26.7351759) 
        self.assertAlmostEqual(cone_flow.shock_ang, shock_ang_exp, places=1)

        #compare mach number at cone surface 
        M_c_exp =  2.19924585
        self.assertAlmostEqual(cone_flow.M_c, M_c_exp, places=1)

        #compare pressure ratios
        p2_p1_exp, p02_p01_exp, pc_p1_exp, p0c_p01_exp = 1.30902568, 0.99801286, 1.59661315, 0.99801286
        self.assertAlmostEqual(cone_flow.p2_p1, p2_p1_exp,      places=1)
        self.assertAlmostEqual(cone_flow.p02_p01, p02_p01_exp,  places=1)
        self.assertAlmostEqual(cone_flow.pc_p1, pc_p1_exp,      places=1)
        self.assertAlmostEqual(cone_flow.p0c_p01, p0c_p01_exp,  places=1)

        #compare temperature ratios
        T2_T1_exp, Tc_T1_exp =  1.08058909, 1.14367828
        self.assertAlmostEqual(cone_flow.T2_T1, T2_T1_exp, places=1)
        self.assertAlmostEqual(cone_flow.Tc_T1, Tc_T1_exp, places=1)

        #compare density ratios
        rho2_rho1_exp, rhoc_rho1_exp =  1.21140005, 1.39603346
        self.assertAlmostEqual(cone_flow.rho2_rho1, rho2_rho1_exp, places=1)
        self.assertAlmostEqual(cone_flow.rhoc_rho1, rhoc_rho1_exp, places=1)

if __name__ == "__main__":
    unittest.main()