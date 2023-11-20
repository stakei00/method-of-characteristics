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
#test_case = "NASA_D6078_M3.0"
#test_case = "NASA_D6078_M3.47"
test_case = "2D_isentropic_ramp_M2.7"

#Characterisic Mesh Type
shock_mesh = True #if true, shock waves will be computed within the mesh 

#Plotting 
display_geometry_preview = True
display_all_plots = False #if True, all plots will be displayed, otherwise just mesh

#Export
export_to_csv = False

################################################################################

if test_case == "straight_Cone_M2.5":
    inletFile = "straight_cone_M2.5_geom.json"
    if shock_mesh: inputFile = 'straight_cone_M2.5_shock_mesh_input.json'
    else: inputFile = 'straight_cone_M2.5_isentropic_mesh_input.json'

elif test_case == "NASA_D6078_M3.0":
    #inletFile = "NASA_D6078_M3_geom.json"
    inletFile = "NASA_D6078_M3_Interpolated.json"
    if shock_mesh: inputFile = "NASA_D6078_M3_shock_mesh_input.json"
    else: inputFile = "NASA_D6078_M3_isentropic_mesh_input.json"

elif test_case == "NASA_D6078_M3.47":
    #inletFile = "NASA_D6078_M3.47_geom.json"
    inletFile = "NASA_D6078_M3.47_Interpolated.json"
    if shock_mesh: inputFile = "NASA_D6078_M3.47_shock_mesh_input.json"
    else: inputFile = "NASA_D6078_M3.47_isentropic_mesh_input.json"

elif test_case == "2D_isentropic_ramp_M2.7":
    inletFile = "2D_isentropic_ramp_M2.7_geom.json"
    if shock_mesh: inputFile = "2D_isentropic_ramp_M2.7_shock_mesh_input.json"
    else: inputFile = "2D_isentropic_ramp_M2.7_isentropic_mesh_input.json"

if display_all_plots: 
    plotfile = "plot_all.json" #full plots
else: 
    plotfile = "plot_mesh.json" #only show mesh

#RUN SOLUTION###################################################################
import AIMCAT as aimcat
aimcat.Main(inputFile=inputFile, 
            geomFile=inletFile, 
            plotFile=plotfile,
            export=export_to_csv, 
            preview_geom=display_geometry_preview)