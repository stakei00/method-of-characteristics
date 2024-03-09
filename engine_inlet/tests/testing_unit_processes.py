"""
unit testing script to verify output of various method of characteristics operators is within expected bounds
TODO: add for below wall operators
"""
import unittest 
import math 
import os 
import sys 
sys.path.append(os.getcwd())
import aimcat.method_of_characteristics.unit_processes as unit_processes 

class test_operators(unittest.TestCase):
    #create input objects
    class gasProps: 
        def __init__(self, gam, R, T0):
            self.gam, self.R, self.T0 = gam, R, T0
            self.a0 = math.sqrt(gam*R*T0)
    class pointData:
        def __init__(self, u, v, x, y, thet=None):
            self.u, self.v, self.x, self.y = u, v, x, y
            if thet is not None: self.thet = thet

    gas = gasProps(1.2, 320, 3000)  
    funcs = unit_processes.operator_funcs()
    delta = 1

    #Testing Interior Operator 
    def test_interior(self): 

        pt2 = self.pointData(2473.4, 812.8, 0.13146, 0.040118)    #Point 2 (m/s, m)
        pt1 = self.pointData(2502.8, 737.6, 0.135683, 0.037123)    #point 1 (m/s, m
        
        res = unit_processes.interior_point(pt1, pt2, self.gas, self.delta, 0.0001, self.funcs) #running interior operator on single test point
        expected_res = [0.14118, 0.04056, 2510.1, 780.2] 
        self.assertAlmostEqual(res[0], expected_res[0], places=3) 
        self.assertAlmostEqual(res[1], expected_res[1], places=3) 
        self.assertAlmostEqual(res[2], expected_res[2], places=1) 
        self.assertAlmostEqual(res[3], expected_res[3], places=1) 
    
    #Testing Direct Wall Operator
    def test_dirWall_abv(self):

        y_x = lambda x : 0.0221852 + 0.71568*x - 1.0787*x**2 #wall y function 
        dydx = lambda x : 0.71568 - 2*1.0787*x #wall slope function 

        pt1 = self.pointData(1967.1, 1141.3, 0.06048, 0.059625)
        res = unit_processes.direct_wall(pt1, y_x, dydx, self.gas, self.delta, 0.0001, self.funcs, "pos")
        expected_res = [0.063485, 0.063273, 1977.4, 1144.4]

        self.assertAlmostEqual(res[0], expected_res[0], places=3)
        self.assertAlmostEqual(res[1], expected_res[1], places=3)
        self.assertAlmostEqual(res[2], expected_res[2], places=1)
        self.assertAlmostEqual(res[3], expected_res[3], places=1)
    
if __name__ == "__main__":
    unittest.main() 