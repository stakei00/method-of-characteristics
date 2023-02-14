import unittest
import taylor_maccoll
import math

"""
This file tests the function of the Taylor_Maccoll.py module 
"""

class Test_Taylor_Maccoll(unittest.TestCase):

    def test_TMC_Cone_Shock_Angle(self):
        """
        Expected Results from Compressible Aerodynamics Calculator 3.1 
        https://devenport.aoe.vt.edu/aoe3114/calc.html
        """
        #instantiate object: 
        cone_ang = math.radians(30)
        M_inf = 3
        gam = 1.4
        cone_flow = taylor_maccoll.TaylorMaccoll_Cone(cone_ang, M_inf, gam)
        #Compare calculated shock angle with expected angle: 
        shock_ang_exp = math.radians(39.8169519) 
        self.assertAlmostEqual(cone_flow.shock_ang, shock_ang_exp, places=1)

        #compare mach number at cone surface 
        M_c_exp = 1.82959328
        self.assertAlmostEqual(cone_flow.M_c, M_c_exp, places=1)

        #compare pressure ratios
        p2_p1_exp, p02_p01_exp, pc_p1_exp, p0c_p01_exp = 4.13866355, 0.75765102, 4.62898681, 0.75765102
        self.assertAlmostEqual(cone_flow.p2_p1, p2_p1_exp,      places=1)
        self.assertAlmostEqual(cone_flow.p02_p01, p02_p01_exp,  places=1)
        self.assertAlmostEqual(cone_flow.pc_p1, pc_p1_exp,      places=1)
        self.assertAlmostEqual(cone_flow.p0c_p01, p0c_p01_exp,  places=1)

        #compare temperature ratios
        T2_T1_exp, Tc_T1_exp =  1.62436310, 1.67716661
        self.assertAlmostEqual(cone_flow.T2_T1, T2_T1_exp, places=1)
        self.assertAlmostEqual(cone_flow.Tc_T1, Tc_T1_exp, places=1)

        #compare density ratios
        rho2_rho1_exp, rhoc_rho1_exp =  2.54786848, 2.76000415
        self.assertAlmostEqual(cone_flow.rho2_rho1, rho2_rho1_exp, places=1)
        self.assertAlmostEqual(cone_flow.rhoc_rho1, rhoc_rho1_exp, places=1)

if __name__ == "__main__":
    unittest.main()