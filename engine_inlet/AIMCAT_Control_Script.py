import AIMCAT as aimcat
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
test_case = "straight_Cone_M2.5"
#test_case = "NASA_D6078_M3.0"
#test_case = "NASA_D6078_M3.47"
#test_case = "2D_isentropic_ramp_M2.7"

#Plotting 
display_geometry_preview = True
display_all_plots = True #if false, only the mesh will be displayed 

#Export
export_to_csv = False

################################################################################

if test_case == "straight_Cone_M2.5":
    inletFile = "straight_cone_M2.5_geom.json"
    inputFile = 'straight_cone_M2.5_input.json'

elif test_case == "NASA_D6078_M3.0":
    inletFile = "NASA_D6078_M3_geom.json"
    inputFile = "NASA_D6078_M3_input.json"

elif test_case == "NASA_D6078_M3.47":
    inletFile = "NASA_D6078_M3.47_geom.json"
    inputFile = "NASA_D6078_M3.47_input.json"

elif test_case == "2D_isentropic_ramp_M2.7":
    inletFile = "2D_isentropic_ramp_M2.7_geom.json"
    inputFile = "2D_isentropic_ramp_M2.7_input.json"

if display_all_plots: 
    plotfile = "plot_all.json" #full plots
else: 
    plotfile = "plot_mesh.json" #just the mesh

#RUN SOLUTION###################################################################
aimcat.Main(inputFile=inputFile, geomFile=inletFile, plotFile=plotfile, \
            export=export_to_csv, preview_geom=display_geometry_preview)