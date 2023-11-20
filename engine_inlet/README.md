# Air-Inlet-Method-of-Characteristics-Analysis-Tool (AIMCAT)

This python-based tool contains modules used to model axisymmetric and 
two-dimensional supersonic air inlets. The primary method used is the method of 
characteristics which models the most complex regions of the supersonic flow.

# Requirements

The following packages are required: 
numpy, scipy, matplotlib, json, pandas

I have been exclusively working with Python 3.9 for the development of this code, 
so I can't promise compatibility with previous versions.

# How to Use

The intended workflow for this tool is as follows: 

1. configure input files 
    a. run file (.json) - contains parameters including freestream properties, 
        gas constants, and mesh settings (see Input folder for examples)
    b. geometry file (.json) - class containing geometry data. Can be of multiple forms: 
        -piecewise function declaration (see 'geometry/straight_cone_M2.5_geom.json')
        -least squares polynomial fit with 1st and 2nd derivative continuity
            see ('geometry/NASA_D6078_MX.XX_geom.json)
        -spline interpolation (see 'geometry/NASA_D6078_MX.XX_Interpolated.json')
    c. plot profile (.json) - contains settings for generating figures after a 
        solution has been run (see post_processing folder)
            see 'post_processing/plot_all.json' for generating full range of plots 

2. run the main class by importing aimcat, calling aimcat.main() with all your .json
input file paths as arguments. Boolean keyword arguments include:

        export:         export results to .csv file when solution finished 
        preview_geom:  generate a preview of the loaded geometry before solving 
