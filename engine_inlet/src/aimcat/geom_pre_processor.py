from scipy.optimize import curve_fit 
#from scipy.interpolate import PchipInterpolator
from scipy import interpolate
import math 
"""
Class responsible for parametricizing discrete geometry data into polynomial 
equations which can be used with AIMCAT
"""
class Inlet_Geom:

    def __init__(self, inputDict):
        """
        creates AIMCAT geometry object from discrete data
        """
        #unpack input dictionary 
        name, init_turn_ang_deg, geom_type = None, None, None
        cowl_data_list, centerbody_data_list = None, None
        cowl_list, centerbody_list = None, None
        cowl_coord_list, centerbody_coord_list = None, None
        x_lip = 0 
        for key in inputDict.keys():
            if key == "name": name = inputDict[key]
            if key == "initial deflection angle (deg)": init_turn_ang_deg = inputDict[key]
            if key == "2D/Axi": geom_type = inputDict[key]
            if key == "cowl fit data": cowl_data_list = inputDict[key]
            if key == "centerbody fit data": centerbody_data_list = inputDict[key]
            if key == "cowl functions": cowl_list = inputDict[key]
            if key == "centerbody functions": centerbody_list = inputDict[key]
            if key == "cowl interpolant coordinates": cowl_coord_list = inputDict[key]
            if key == "centerbody interpolant coordinates": centerbody_coord_list = inputDict[key]
            if key == "cowl lip x coord": x_lip = inputDict[key]

        self.x_cowlLip = x_lip
        self.name = name
        self.init_turn_ang_deg = init_turn_ang_deg
        self.geom_type = geom_type

        #starting with centerbody:
        print("\nCenterbody Geometry Parametrization")
        if centerbody_list is not None:
            y_centerbody, dydx_centerbody, centerbody_bounds = self.unpack_surface_function_list(centerbody_list)

        elif centerbody_coord_list is not None:
            y_centerbody, dydx_centerbody, centerbody_bounds = self.interpolate_coords(centerbody_coord_list) 

        elif centerbody_data_list is not None: #if geometry includes centerbody data 
            if len(centerbody_data_list) > 1: 
                curve = Composite_Curve(centerbody_data_list)

            elif len(centerbody_data_list) == 1:
                [type_, startpoint, xdata, ydata, startpoint_dydx, x_end, \
                    endpoint] = self.read_sub_dict(centerbody_data_list[0]) 
                
                curve = Curve(type_, startpoint, xdata=xdata, ydata=ydata, \
                              startpoint_dydx=startpoint_dydx, x_end=x_end, \
                                endpoint=endpoint)

            y_centerbody = curve.y
            dydx_centerbody = curve.dydx
            centerbody_bounds = curve.x_bounds 

        if init_turn_ang_deg is None and centerbody_list is not None:
            #if not given, calculate it
            self.init_turn_ang_deg = math.degrees(math.atan(self.dydx_centerbody(curve.x_bounds[0]))) 

        #next cowl
        print("\nCowl Geometry Parameterization")
        if cowl_list is not None:
            y_cowl, dydx_cowl, cowl_bounds = self.unpack_surface_function_list(cowl_list)

        elif cowl_coord_list is not None: 
            y_cowl, dydx_cowl, cowl_bounds = self.interpolate_coords(cowl_coord_list) 

        elif cowl_data_list is not None: #if geometry includes cowl data
            
            if len(cowl_data_list) > 1: 
                curve = Composite_Curve(cowl_data_list)
            
            elif len(cowl_data_list) == 1: 
                [type_, startpoint, xdata, ydata, startpoint_dydx, startpoint_d2ydx2, x_end, \
                    endpoint] = self.read_sub_dict(cowl_data_list[0])
                
                curve = Curve(type_, startpoint, xdata=xdata, ydata=ydata,
                                                startpoint_dydx=startpoint_dydx, 
                                                startpoint_d2ydx2=startpoint_d2ydx2,
                                                x_end=x_end,
                                                endpoint=endpoint)

            y_cowl = curve.y
            dydx_cowl = curve.dydx    
            cowl_bounds = curve.x_bounds

        #set centerbody functions 
        self.y_centerbody = y_centerbody
        self.dydx_centerbody = dydx_centerbody
        self.centerbody_bounds = centerbody_bounds

        #add cowl lip offset to cowl functions and bounds
        self.y_cowl = lambda x: y_cowl(x-self.x_cowlLip)
        self.dydx_cowl = lambda x: dydx_cowl(x-self.x_cowlLip)
        self.cowl_bounds = [x+self.x_cowlLip for x in cowl_bounds]

    def unpack_surface_function_list(self, surface_list):
        """
        unpacks a subdictionary which defines the surface using functions rather 
        data points. Functions should be in lambda form
        """
        y_funcs, dydx_funcs = [], []
        endpoints = []

        for subDict in surface_list:
            endpoints.append(subDict["x bounds"])

            y_func = eval(subDict["y function"])
            dydx_func = eval(subDict["dydx function"])

            y_funcs.append(y_func), dydx_funcs.append(dydx_func)

        def y_func_combined(x):
            y = None
            for i,func in enumerate(y_funcs):
                if endpoints[i][0] <= x <= endpoints[i][-1]:
                    y = func(x)
                else:continue
            return y

        def dydx_func_combined(x):
            dydx = None
            for i,func in enumerate(dydx_funcs):
                if endpoints[i][0] <= x <= endpoints[i][-1]:
                    dydx = func(x)
                else:continue
            return dydx

        endpoints_flatten = []
        for pt in endpoints: 
            endpoints_flatten += pt

        return y_func_combined, dydx_func_combined, [min(endpoints_flatten), max(endpoints_flatten)]
       
    def interpolate_coords(self, surface_list):
        """
        uses pchip interpolation to define the surface
        """
        y_list, dydx_list = [],[]
        x_tot = []
        for curve in surface_list: 
            x_tot += curve["x coords"]
            #interpolator = PchipInterpolator(curve["x coords"], curve["y coords"], extrapolate=False)
            spl = interpolate.splrep(curve["x coords"], curve["y coords"], k=3)
            #spl = interpolate.CubicSpline(curve["x coords"], curve["y coords"])
            def y(x):
                try: 
                    #return float(interpolator(x)) #position
                    return interpolate.splev([x], spl, der=0)[0]
                    #return float(spl(x, nu=0))
                except: 
                    return None
            def dydx(x):    
                try: 
                    #return float(interpolator.derivative()(x))
                    return interpolate.splev([x], spl, der=1)[0] 
                    #return float(spl(x, nu=1))         
                except: return None

            y_list.append(y)
            dydx_list.append(dydx)

        if len(surface_list) > 1: 
            
            def y_comb(x):
                for func in y_list:
                    y = func(x)
                    if y is not None: 
                        return y
            def dydx_comb(x):
                for func in dydx_list:
                    dydx = func(x)
                    if dydx is not None: 
                        return dydx

            return y_comb, dydx_comb, [min(x_tot), max(x_tot)]

        return y_list[0], dydx_list[0], [min(x_tot), max(x_tot)] 

    def read_sub_dict(self, subDict):
        """
        takes in a curve subdictionary and unpacks it
        """
        keys = subDict.keys()
        type_ = subDict["type"]
        startpoint = subDict["startpoint"]
        startpoint_dydx, startpoint_d2ydx2, x_end,endpoint,xdata,ydata = None,None,None,None,None, None
        
        if "startpoint dydx" in keys:
            startpoint_dydx = subDict["startpoint dydx"]
        if "startpoint d2ydx2" in keys: 
            startpoint_d2ydx2 = subDict["startpoint d2ydx2"]
        if "endpoint x" in keys:
            x_end = subDict["endpoint x"]
        if "endpoint" in keys: 
            endpoint = subDict["endpoint"]
        if "xdata" in keys: 
            xdata, ydata = subDict["xdata"], subDict["ydata"]    

        return [type_, startpoint, xdata, ydata, startpoint_dydx, startpoint_d2ydx2, x_end, endpoint]


class Curve:
    """
    generates functions defining the position and first derivative of a 
    single curve from input data
    """
    def __init__(self, type_, startpoint, xdata=None, ydata=None, \
                 startpoint_dydx=None, startpoint_d2ydx2=None, x_end=None, endpoint=None):
        
        if type_=="linear":
            y, dydx, d2ydx2 = self.linear_poly(startpoint, startPoint_dydx=startpoint_dydx, \
                                  x_end=x_end, endPoint=endpoint)
            if x_end is not None: self.x_bounds = (startpoint[0], x_end)
            elif endpoint is not None: self.x_bounds = (startpoint[0], endpoint[0])

        elif type_=="polynomial":
            y,dydx, d2ydx2 = self.least_squares_poly(xdata, ydata, startpoint,\
                                                startPoint_dydx=startpoint_dydx, 
                                                startpoint_d2ydx2=startpoint_d2ydx2)
            self.x_bounds =(startpoint[0], max(xdata))

        self.y = y
        self.dydx = dydx
        self.d2ydx2 = d2ydx2

    def least_squares_poly(self, xs, ys, startPoint, startPoint_dydx=None, startpoint_d2ydx2=None):
        """
        returns polymial equation for y_x as well as dydx of a curve defined 
        """
        y_x, dydx_x = None, None
    
        x_a, y_a = startPoint    

        def func_free(x, A, B, C, D, E):
            return A*(x-x_a)**5 + B*(x-x_a)**4 + C*(x-x_a)**3 + D*(x-x_a)**2 + E*(x-x_a) + y_a

        def func_const_1st_deriv(x, A, B, C, D):
            return A*(x-x_a)**5 + B*(x-x_a)**4 + C*(x-x_a)**3 + D*(x-x_a)**2 + startPoint_dydx*(x-x_a) + y_a
        
        def func_const_1st_and_2nd_deriv(x, A, B, C, D):
            return A*(x-x_a)**6 + B*(x-x_a)**5 + C*(x-x_a)**4 + D*(x-x_a)**3 + startpoint_d2ydx2*(x-x_a)**2 + startPoint_dydx*(x-x_a) + y_a

        if startpoint_d2ydx2 is None and startPoint_dydx is None: 
            coeffs, _ = curve_fit(func_free, xs, ys)
            A,B,C,D,E = coeffs[0],coeffs[1],coeffs[2],coeffs[3],coeffs[4]

            print(f"\tleast-squares polynomial function:\n\
            {A}*(x-x_a)**5 + {B}*(x-x_a)**4 +{C}*(x-x_a)**3 + {D}*(x-x_a)**2 + {E}*(x-x_a) + y_a\n\
            x_a: {x_a}\n\
            y_a: {y_a}")

            def y_x(x):
                if x_a <= x <= max(xs):
                    return A*(x-x_a)**5 + B*(x-x_a)**4 + C*(x-x_a)**3 + D*(x-x_a)**2 + E*(x-x_a) + y_a
                else: return None
            def dydx_x(x):
                if x_a <= x <= max(xs):
                    return 5*A*(x-x_a)**4 + 4*B*(x-x_a)**3 + 3*C*(x-x_a)**2 + 2*D*(x-x_a) + E
                else: return None 
            def d2ydx2_x(x):
                if x_a <= x <= max(xs):
                    return 20*A*(x-x_a)**3 + 12*B*(x-x_a)**2 + 6*C*(x-x_a) + 2*D
                else: return None 

        if startPoint_dydx is not None and startpoint_d2ydx2 is None: 
            coeffs, _ = curve_fit(func_const_1st_deriv, xs, ys)
            A,B,C,D = coeffs[0],coeffs[1],coeffs[2],coeffs[3]
            print(f"\tleast-squares polynomial function:\n\
            {A}*(x-x_a)**5 + {B}*(x-x_a)**4 + {C}*(x-x_a)**3 + {D}*(x-x_a)**2 + {startPoint_dydx}*(x-x_a) + y_a\n\
            x_a: {x_a}\n\
            y_a: {y_a}")
            
            def y_x(x):
                if x_a <= x <= max(xs): 
                    return A*(x-x_a)**5 + B*(x-x_a)**4 + C*(x-x_a)**3 + D*(x-x_a)**2 + startPoint_dydx*(x-x_a) + y_a
                else: return None
            def dydx_x(x):
                if x_a <= x <= max(xs):
                    return 5*A*(x-x_a)**4 + 4*B*(x-x_a)**3 + 3*C*(x-x_a)**2 + 2*D*(x-x_a) + startPoint_dydx
                else: return None
            def d2ydx2_x(x):
                if x_a <= x <= max(xs):
                    return 20*A*(x-x_a)**3 + 12*B*(x-x_a)**2 + 6*C*(x-x_a) + 2*D
                else: return None

        if startpoint_d2ydx2 is not None and startPoint_dydx is not None:
            #curves enforcing a set first and 2nd derivative 
            coeffs, _ = curve_fit(func_const_1st_and_2nd_deriv, xs, ys)
            A,B,C,D = coeffs[0],coeffs[1],coeffs[2],coeffs[3]
            print(f"\tleast-squares polynomial function:\n\
            {A}*(x-x_a)**6 + {B}*(x-x_a)**5 + {C}*(x-x_a)**4 + {D}*(x-x_a)**3 + {startpoint_d2ydx2}*(x-x_a)**2 + {startPoint_dydx}*(x-x_a) + y_a\n\
            x_a: {x_a}\n\
            y_a: {y_a}")
            def y_x(x):
                if x_a <= x <= max(xs): 
                    return A*(x-x_a)**6 + B*(x-x_a)**5 + C*(x-x_a)**4 + D*(x-x_a)**3 + startpoint_d2ydx2*(x-x_a)**2 + startPoint_dydx*(x-x_a) + y_a
                else: return None
            def dydx_x(x):
                if x_a <= x <= max(xs):
                    return 6*A*(x-x_a)**5 + 5*B*(x-x_a)**4 + 4*C*(x-x_a)**3 + 3*D*(x-x_a)**2 + 2*startpoint_d2ydx2*(x-x_a) + startPoint_dydx
                else: return None 
            def d2ydx2_x(x):
                if x_a <= x <= max(xs):
                    return 30*A*(x-x_a)**4 + 20*B*(x-x_a)**3 + 12*C*(x-x_a)**2 + 6*D*(x-x_a) + 2*startpoint_d2ydx2
                else: return None 

        #compile error and print to console
        y_x_err = [abs(y_x(x) - ys[i]) for i,x in enumerate(xs)]
        print(f"\n\tmaximum least squares position error: {max(y_x_err)}\n\
              average error: {sum(y_x_err)/len(y_x_err)}")

        return y_x, dydx_x, d2ydx2_x

    def linear_poly(self, startPoint, startPoint_dydx=None, x_end=None, endPoint=None):
        """
        returns linear equation for any straight or near straight segments. 
        """
        y_x, dydx_x = None, None 
        x_a, y_a = startPoint

        if startPoint_dydx is not None and x_end is None: 
            raise ValueError("startPoint_dydx and x_end must both be specified")

        if endPoint is not None: #if an endpoint is specified 
            x_b, y_b = endPoint
            m = (y_b-y_a)/(x_b-x_a)

        elif startPoint_dydx is not None: #if start point slope is given
            x_b = x_end
            m = startPoint_dydx

        print(f"\tlinear function:\n\
            {m}*(x-x_a) + y_a\n\
            x_a: {x_a}\n\
            y_a: {y_a}")

        def y_x(x):
            if x_a <= x <= x_b:
                return m*(x-x_a) + y_a
            else: return None
        def dydx_x(x):
            if x_a <= x <= x_b:
                return m
            else: return None
        def d2ydx2_x(x):
            if x_a <= x <= x_b:
                return 0
            else: 
                return None 

        return y_x, dydx_x, d2ydx2_x


class Composite_Curve(Inlet_Geom): 
    """
    generates functions defining the position and first derivative of a curve
    made of multiple sub-curves (polynomial or linear) from input data 
    """
    def __init__(self, curves_list):
        
        y_funcs = []
        dydx_funcs = []
        endpoint_prev = None
        endpoint_dydx_prev = None
        endpoint_d2ydx2_prev = None
        x_endpoints = []

        for i,subDict in enumerate(curves_list):
            
            [type_, startpoint, xdata, ydata, startpoint_dydx, startpoint_d2ydx2, x_end, endpoint] = self.read_sub_dict(subDict)
            if i == 0: x_endpoints.append(startpoint[0])

            if startpoint == "prev endpoint":
                startpoint = endpoint_prev
            if startpoint_dydx == "prev endpoint":
                startpoint_dydx = endpoint_dydx_prev
            if startpoint_d2ydx2 == "prev endpoint":
                startpoint_d2ydx2 = endpoint_d2ydx2_prev

            curve = Curve(type_, startpoint, xdata, ydata, startpoint_dydx, startpoint_d2ydx2, x_end, endpoint)
            x_endpoints.append(curve.x_bounds[-1])
            y_funcs.append(curve.y), dydx_funcs.append(curve.dydx)
            endpoint_prev = (curve.x_bounds[-1], curve.y(curve.x_bounds[-1]))
            endpoint_dydx_prev = curve.dydx(curve.x_bounds[-1])
            endpoint_d2ydx2_prev = curve.d2ydx2(curve.x_bounds[-1])

        self.x_bounds = (min(x_endpoints), max(x_endpoints))
        self.merge_curve_functions(y_funcs, dydx_funcs)

    def merge_curve_functions(self, y_funcs, dydx_funcs):
        """
        returns functions defining y and dydx of the whole composite curve
        """
        self.y = lambda x: sum([func(x) for func in y_funcs if func(x) is not None])
        self.dydx = lambda x: sum([func(x) for func in dydx_funcs if func(x) is not None])
