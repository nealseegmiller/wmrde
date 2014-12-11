WMRDE, Wheeled Mobile Robot Dynamics Engine
MATLAB code

Copyright (c) 2014, Neal Seegmiller
license info: see LICENSE.txt

Please report bugs and serious omissions in documentation

GETTING STARTED
To understand this code, start with the following demo scripts:
test/test_simulate.m
test/test_calibrate.m

Before running the demos:
*enter the command:
addpath(genpath('C:\Users\nseegmil.NREC-014635\Documents\Thesis\MATLAB'))
but replace the string with the correct path for your machine

*change the strings in the following functions to the correct paths for your machine:
CADdir.m
resourcedir.m 
datalogdir.m
Directories are included in this wmrde repo with files for the demos.

HOW TO MODEL A NEW VEHICLE
Write a function to construct the vehicle model. Follow one of the examples in /demo/model
e.g. rocky.m, zoe.m, talon.m, etc.

Several function handles must be set in the WmrModel class.
Write your own functions for the controller and (if applicable) joint constraints.
Use one of the provided functions for the wheel-ground contact and actuator models, or write your own.

You may contact me for help.

OPTIONAL MODIFICATIONS
*A method for geometric collision detection between the wheels/tracks and ground is provided in:
/collision
/surface
You may use a different method, but it must set WheelContactGeom/TrackContactGeom objects for input to the forward kinematics/dynamics functions.

*to switch between Euler angle/quaternion orientation representation, change SIZEORIENT() output
Some code assumes Euler angles.

