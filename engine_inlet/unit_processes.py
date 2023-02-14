import math
import numpy as np
import scipy.optimize

"""
Testing out Marshall's system of equations approach to solving for the flow field properties of an interior point
Assumptions: supersonic, irrotational, isentropic flow 
"""

class operator_functions: 
    """
    hands off repeatedly-used functions object which are necessary for MOC unit processes
    """
    def __init__(self):
        self.a0 = lambda gam, R, T0: math.sqrt(gam*R*T0)
        self.a = lambda a0, gam, u, v : math.sqrt(a0**2 - 0.5*(gam-1)*(u**2 + v**2))
        self.S = lambda delta, a, v, y : delta*(a**2*v/y)
        self.Q = lambda u, a : u**2 - a**2
        self.R = lambda u, v, Q, lam: 2*u*v - Q*lam
        self.lam = lambda lam1, lam2 : 0.5*(lam1 + lam2)
        self.lam_min = lambda u, v, a : (u*v - a*math.sqrt(u**2 + v**2 - a**2))/(u**2 - a**2)
        self.lam_plus = lambda u, v, a : (u*v + a*math.sqrt(u**2 + v**2 - a**2))/(u**2 - a**2)

def interior_point(pt1, pt2, gasProps, delta, vel_TOL, funcs):
    """
    MOC interior point solution using irrotational, isentropic equations
    """
    #Unpack input objects
    u1, v1, x1, y1 = pt1.u, pt1.v, pt1.x, pt1.y #point 1 data
    u2, v2, x2, y2 = pt2.u, pt2.v, pt2.x, pt2.y #point 2 data
    gam, R_gas, T0 = gasProps.gam, gasProps.R, gasProps.T0 #gas properties

    #getting speed of sound 
    a0 = funcs.a0(gam, R_gas, T0) #stagnation speed of sound

    #initial estimates for first iteration: 
    u13, v13 = u1, v1 
    y13 = 1 if y1 == 0 else y1  
    u23, v23 = u2, v2
    y23 = 1 if y2 == 0 else y2

    def solve_interior_point(u13, v13, y13, u23, v23, y23, first_iter=None): 
        #Calculate coefficients
        a1 = funcs.a(a0, gam, u1, v1)

        if first_iter: 
            lam13 = funcs.lam_plus(u1, v1, a1)
        else: 
            a3 = funcs.a(a0, gam, u3, v3)
            lam13 = 0.5*(funcs.lam_plus(u1, v1, a1) + funcs.lam_plus(u3, v3, a3))

        a13 = funcs.a(a0, gam, u13, v13)
        S13, Q13 = funcs.S(delta, a13, v13, y13), funcs.Q(u13, a13)
        R13 = funcs.R(u13, v13, Q13, lam13)

        a2 = funcs.a(a0, gam, u2, v2)
        lam23 = funcs.lam_min(u2, v2, a2)
        a23 = funcs.a(a0, gam, u23, v23)
        S23, Q23 = funcs.S(delta, a23, v23, y23), funcs.Q(u23, a23)
        R23 = funcs.R(u23, v23, Q23, lam23)

        #Construct system and invert to find solution vector 
        coeffMat = np.array([[lam13, -1, 0, 0],[lam23, -1, 0, 0],[-S13, 0, Q13, R13],[-S23, 0, Q23, R23]])
        RHSvec = np.array([lam13*x1-y1, lam23*x2-y2, -S13*x1+Q13*u1+R13*v1, -S23*x2+Q23*u2+R23*v2])
        return np.linalg.solve(coeffMat, RHSvec) #[x3, y3, u3, v3]
    
    #first iteration
    [x3, y3, u3, v3] = solve_interior_point(u13, v13, y13, u23, v23, y23, first_iter=True)

    #Iterate above process until values converge
    delta_vel = vel_TOL
    while delta_vel >= vel_TOL:
        u3_old, v3_old = u3, v3
        u13, v13, y13 = 0.5*(u1 + u3), 0.5*(v1 + v3), 0.5*(y1 + y3)
        u23, v23, y23 = 0.5*(u2 + u3), 0.5*(v2 + v3), 0.5*(y2 + y3)
        [x3, y3, u3, v3] = solve_interior_point(u13, v13, y13, u23, v23, y23)
        delta_vel = max([abs(u3 - u3_old), abs(v3 - v3_old)])

    return [x3, y3, u3, v3]

def direct_wall(pt1, y_x, dydx, gasProps, delta, vel_TOL, funcs):
    #MOC direct wall solution using irrotational, isentropic equations 
    #currently only works for wall above 
    #TODO adjust input function for y to be in implicit form 

    #unpacking input data 
    u1, v1, x1, y1 = pt1.u, pt1.v, pt1.x, pt1.y
    gam, R_gas, T0 = gasProps.gam, gasProps.R, gasProps.T0 #gas properties 

    #getting speed of sound
    a0 = funcs.a0(gam, R_gas, T0) #stagnation speed of sound

    #initial value computations 
    u13, v13 = u1, v1 
    y13 = 1 if y1 == 0 else y1

    def solve_direct_wall(u13, v13, y13, first_iter=None):
        a1 = funcs.a(a0, gam, u1, v1)
        a13 = funcs.a(a0, gam, u13, v13)

        if first_iter: 
            lam13 = funcs.lam_plus(u1, v1, a1)
        else: 
            a3 = funcs.a(a0, gam, u3, v3)
            lam13 = 0.5*(funcs.lam_plus(u1, v1, a1) + funcs.lam_plus(u3, v3, a3))

        S13 = funcs.S(delta, a13, v13, y13)
        Q13 = funcs.Q(u13, a13)
        R13 = funcs.R(u13, v13, Q13, lam13)

        #computing the x&y location of 3' for first iteration (2 methods suggested)
        try:
            #method 2: use lam13 to project to wall to get (x,y)_3'
            x3p = float(scipy.optimize.fsolve(lambda x: lam13*(x-x1) + y1 - y_x(x), x1))
            #sol = scipy.optimize.root_scalar(lambda x: lam13*(x-x1) + y1 - y_x(x), method='bisect', bracket=(0, 0.25)) #!Not solving with number as ZH book 
            #x3p = sol.root
            y3p = y_x(x3p)  
        except: 
            #method 1: go straight up from x1 to get x_3' if no downstream solution exists given lam13 and point 1
            x3p = x1
            y3p = y_x(x3p)

        #computing partial derivatives for taylor series expansion using 2nd order centered differences
        pfpx_3p = -dydx(x3p)
        pfpy_3p = 1 #!Hard Coded
        f_w3p = 0 #!Hard Coded 

        #construct system and solve
        coeffMat = np.array([[lam13, -1, 0, 0],[pfpx_3p, pfpy_3p, 0, 0],[-S13, 0, Q13, R13],[0, 0, pfpx_3p, pfpy_3p]])
        RHSvec = np.array([lam13*x1-y1, pfpx_3p*x3p+pfpy_3p*y3p-f_w3p, -S13*x1+Q13*u1+R13*v1, 0])
        return np.linalg.solve(coeffMat, RHSvec) #[x3, y3, u3, v3]

    #first iteration
    count = 1
    [x3, y3, u3, v3] = solve_direct_wall(u13, v13, y13, first_iter=True)

    #Iterate above process until values converge
    delta_vel = vel_TOL
    while delta_vel >= vel_TOL:
        u3_old, v3_old = u3, v3
        u13, v13, y13 = 0.5*(u1 + u3), 0.5*(v1 + v3), 0.5*(y1 + y3)
        [x3, y3, u3, v3] = solve_direct_wall(u13, v13, y13)
        delta_vel = max([abs(u3 - u3_old), abs(v3 - v3_old)])
        count += 1
    
    print(f"direct wall converged after {count} iterations")

    return [x3, y3, u3, v3] 

def inverse_wall(pt1, pt2, pt3, gasProps, delta, vel_TOL, funcs):
    #!Not Working
    #unpacking input data
    u1, v1, x1, y1 = pt1.u, pt1.v, pt1.x, pt1.y
    u2, v2, x2, y2 = pt2.u, pt2.v, pt2.x, pt2.y
    x3, y3, thet3 = pt3.x, pt3.y, pt3.thet

    gam, R_gas, T0 = gasProps.gam, gasProps.R, gasProps.T0 #gas properties 

    #getting speed of sound
    a0 = funcs.a0(gam, R_gas, T0) #stagnation speed of sound

    #compute values for C- characteristic (2 -> 1)
    u12, v12, y12 = 0.5*(u1+u2), 0.5*(v1+v2), 0.5*(y1+y2)
    a1 = funcs.a(a0, gam, u1, v1)
    a2 = funcs.a(a0, gam, u2, v2)
    a12 = funcs.a(a0, gam, u12, v12)
    lam12 = 0.5*(funcs.lam_min(u1, v1, a1) + funcs.lam_min(u2, v2, a2))
    S12 = funcs.S(delta, a12, v12, y12)
    Q12 = funcs.Q(u12, a12)
    R12 = funcs.R(u12, v12, Q12, lam12)

    def solve_inverse_wall(xa, ya, ua, va, u3, v3): 
        
        ua3, va3, ya3 = 0.5*(ua+u3), 0.5*(va+v3), 0.5*(ya+y3)
        aa3 = funcs.a(a0, gam, ua3, va3)

        aa = funcs.a(a0, gam, ua, va)
        a3 = funcs.a(a0, gam, u3, v3)
        lama3 = 0.5*(funcs.lam_plus(ua, va, aa) + funcs.lam_plus(u3, v3, a3))

        Sa3 = funcs.S(delta, aa3, va3, ya3)
        Qa3 = funcs.Q(ua3, aa3)
        Ra3 = funcs.R(ua3, va3, Qa3, lama3)

        #calculate partial derivatives at point 3 
        pfpx_3 = -math.tan(thet3)
        pfpy_3 = 1

        coeffMat = np.array([[lam12,-1,0,0,0,0], [lama3,-1,0,0,0,0], [-S12,0,Q12,R12,0,0], [-Sa3,0,Qa3,Ra3,-Qa3,-Ra3], [v2-v1,0,0,-(x2-x1),0,0], [0,0,0,0,pfpx_3, pfpy_3]])
        RHSvec = np.array([lam12*x1-y1, lama3*x3-y3, -S12*x1+Q12*u1+R12*v1, -Sa3*x3, (v2-v1)*x1-(x2-x1)*v1,0])
        return np.linalg.solve(coeffMat, RHSvec) #[xa, ya ua, va, u3, v3]

    #first iteration values
    xa, ya, ua, va = 0.5*(x1+x2), 0.5*(y1+y2), 0.5*(u1+u2), 0.5*(v1+v2)
    #ua = u2 #! what zh does
    #va = v2 #! what zh does
    u3, v3 = u2, v2

    #iterate until solution converges:
    delta_vel = vel_TOL
    while delta_vel >= vel_TOL:

        u3_old, v3_old = u3, v3
        [xa, ya, ua, va, u3, v3] = solve_inverse_wall(xa, ya, ua, va, u3, v3)
        delta_vel = max([abs(u3 - u3_old), abs(v3 - v3_old)]) 

    return [xa, ya, ua, va, u3, v3]

def symmetry_boundary(pt2, gasProps, delta, vel_TOL, funcs):
    #MOC symmetry boundary solution for irrotational, isentropic equations
    #currently only works for a point above the symmetry plane

    #unpack inputs: 
    u2, v2, x2, y2 = pt2.u, pt2.v, pt2.x, pt2.y #point 2 data
    gam, R_gas, T0 = gasProps.gam, gasProps.R, gasProps.T0 #gas properties
    
    #get speed of sound
    a0 = funcs.a0(gam, R_gas, T0) 

    #make initial estimates for first iteration
    u23, v23 = u2, v2
    y23 = 1 if y2 == 0 else y2

    def solve_symmetry_boundary(u23, v23, y23, first_iter=None):
        a2 = funcs.a(a0, gam, u2, v2)
        a23 = funcs.a(a0, gam, u23, v23)
        
        if first_iter: #if running the first iteration, lam23 is calculated different from subsequent iterations 
            lam23 = funcs.lam_min(u2, v2, a2)
        else:
            a3 = funcs.a(a0, gam, u3, v3)
            lam23 = 0.5*(funcs.lam_min(u2, v2, a2) + funcs.lam_min(u3, v3, a3))

        S23 = funcs.S(delta, a23, v23, y23)
        Q23 = funcs.Q(u23, a23)
        R23 = funcs.R(u23, v23, Q23, lam23)

        #construct system and solve
        coeffMat = np.array([[lam23, -1, 0, 0],[0, 1, 0, 0],[-S23, 0, Q23, R23],[0, 0, 0, 1]])
        RHSvec = np.array([lam23*x2-y2, 0, -S23*x2+Q23*u2+R23*v2, 0])
        return np.linalg.solve(coeffMat, RHSvec) #[x3, y3, u3, v3]
    
    #first iteration 
    [x3, y3, u3, v3] = solve_symmetry_boundary(u23, v23, y23, first_iter=True)

    #iterate until convergence is attained 
    delta_vel = vel_TOL
    while delta_vel >= vel_TOL:
        u3_old, v3_old = u3, v3
        u23, v23, y23 = 0.5*(u2 + u3), 0.5*(v2 + v3), 0.5*(y2 + y3)
        [x3, y3, u3, v3] = solve_symmetry_boundary(u23, v23, y23)
        delta_vel = max([abs(u3 - u3_old), abs(v3 - v3_old)])

    return [x3, y3, u3, v3]

def obtain_convergence(X):
    """"""
    if len(X) < 2: return None
    range = max(X) - min(X)
    delta = abs(X[-1] - X[1])
    return delta/range

if __name__  == "__main__":
    #Execute if file is run directly: 
    #Initial Values: 
    class point_data: 
        def __init__(self, u, v, x, y, thet=None):
            self.u, self.v, self.x, self.y = u, v, x, y
            if thet is not None: self.thet = thet
    class gasProps:
        def __init__(self, gam, R, T0): 
            self.gam, self.R, self.T0 = gam, R, T0
    pt2 = point_data(2473.4, 812.8, 0.13146, 0.040118)    #Point 2 (m/s, m)
    pt1 = point_data(2502.8, 737.6, 0.135683, 0.037123)    #point 1 (m/s, m
    delta = 1 #axisymmetric flow
    gas = gasProps(1.2, 320, 3000)

    #Testing Interior Operator
    moc_funcs = operator_functions() #loading in heavily used functions in moc calculations
    [x3, y3, u3, v3] = interior_point(pt1, pt2, gas, delta, 0.0001, moc_funcs)

    #Testing Direct Wall Operator
    pt1 = point_data(1967.1, 1141.3, 0.06048, 0.059625)
    y_x = lambda x : 0.0221852 + 0.71568*x - 1.0787*x**2 #wall y function 
    dydx = lambda x : 0.71568 - 2*1.0787*x #wall slope function 
    #[x3, y3, u3, v3] = direct_wall(pt1, y_x, dydx, gas, delta, 0.0001, moc_funcs)

    #Testing Axis of Symmetry Operator
    pt2 = point_data(2306.1, 35.7, 0.079625, 0.001290)
    #[x3, y3, u3, v3] = symmetry_boundary(pt2, gas, delta, 0.001, moc_funcs)

    #Testing Inverse Wall Operator 
    pt1 = point_data(1578.3, 705.7, 0.005495, 0.02602) 
    pt2 = point_data(1577.5, 702.3, 0.005085, 0.02608)
    pt3 = point_data(None, None, 0.005283, 0.02617, thet=math.radians(25)) 
    [xa, ya, ua, va, u3, y3] = inverse_wall(pt1, pt2, pt3, gas, delta, 0.0001, moc_funcs)
    #expecting [u3 = 1583.4, v3 = 738.3]
    pass