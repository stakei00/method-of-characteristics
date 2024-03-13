"""
Control script for aimcat module. Specify input file locations and run AIMCAT by
calling class "Main" with the appropriate inputs filled out. 

aimcat.Main(inputFile:str, geomFile:str, plotFile:str or None, export:bool, 
    preview_geom:bool)
    
    Inputs: 
        inputFile:      file name of main input file
        geomFile:       file name of geometry file 
        plotFile:       file name of plot file 
        export:         export results to csv when solution finished 
        preview_geom:   preview of loaded geometry before running 
    Returns: 
        None
"""
#USER SETTINGS #################################################################
#Test Cases (uncomment 1) 
#test_case = "straight_Cone_M2.5"
test_case = "NASA_D6078_M3.0"
#test_case = "NASA_D6078_M3.47"
#est_case = "2D_isentropic_ramp_M2.7"

#Characterisic Mesh Type
shock_mesh = True #if true, shock waves will be computed within the mesh 

#Plotting 
display_geometry_preview = True
display_all_plots = False #if True, all plots will be displayed, otherwise just mesh and mach scalar

#Export
export_to_csv = False

################################################################################

if test_case == "straight_Cone_M2.5":
    inletFile = "data/geometry/straight_cone_M2.5_geom.json"
    if shock_mesh: inputFile = 'data/input_files/straight_cone_M2.5_shock_mesh_input.json'
    else: inputFile = 'data/input_files/straight_cone_M2.5_isentropic_mesh_input.json'

elif test_case == "NASA_D6078_M3.0":
    inletFile = "data/geometry/NASA_D6078_M3_geom.json"
    if shock_mesh: inputFile = "data/input_files/NASA_D6078_M3_shock_mesh_input.json"
    else: inputFile = "data/input_files/NASA_D6078_M3_isentropic_mesh_input.json"

elif test_case == "NASA_D6078_M3.47":
    inletFile = "data/geometry/NASA_D6078_M3.47_geom.json"
    if shock_mesh: inputFile = "data/input_files/NASA_D6078_M3.47_shock_mesh_input.json"
    else: inputFile = "data/input_files/NASA_D6078_M3.47_isentropic_mesh_input.json"

elif test_case == "2D_isentropic_ramp_M2.7":
    inletFile = "data/geometry/2D_isentropic_ramp_M2.7_geom.json"
    if shock_mesh: inputFile = "data/input_files/2D_isentropic_ramp_M2.7_shock_mesh_input.json"
    else: inputFile = "data/input_files/2D_isentropic_ramp_M2.7_isentropic_mesh_input.json"

if display_all_plots: 
    plotfile = "data/plot_settings/plot_all.json" #full plots
else: 
    plotfile = "data/plot_settings/plot_mesh.json" #only show mesh

#RUN SOLUTION###################################################################
import os, sys
sys.path.append(os.getcwd())
from aimcat import AIMCAT

AIMCAT.Main(inputFile=inputFile, 
            geomFile=inletFile, 
            plotFile=plotfile,
            export=export_to_csv, 
            preview_geom=display_geometry_preview)