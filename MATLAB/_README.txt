This directory contains MATLAB code for the modeling and simulation of wheeled mobile robots (WMRs) and tracked vehicles
Written by Neal Seegmiller in support of
Dynamic Model Formulation and Calibration for Wheeled Mobile Robots, Ph.D. thesis, CMU-RI-TR-14-

license info: see LICENSE.txt
If you use this code in your research, please cite the above thesis or related publications by Neal Seegmiller.

STYLE GUIDE
Naming:
scripts:		lower case, with underscores
functions:		Mixed case, no underscores, start with lower case letter. myFunction()
classes:		Mixed case, no underscores, start with upper case letter. MyClass
variables:		lower case, with underscores. my_variable

coords is an abbreviation for coordinates
HT_framea_to_frameb, Homogeneous Transform from frame a to frame b coords

Vector/matrix orientation:
-Cartesian vectors (position, velocity, acceleration, etc.) are column vectors, so they can be operated on by rotation matrices.
-state and its derivatives are column vectors
-in time-indexed matrices, one row per time step

Comment format for functions:
% description
% INPUT
% input_a:		size type, description
% (optional)
% input_b
% input_c
% OUTPUT
% output_a:	size type, description
% Written by Neal Seegmiller, 
% (license info)



