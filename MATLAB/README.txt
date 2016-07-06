WMRDE, Wheeled Mobile Robot Dynamics Engine
MATLAB code

Copyright (c) 2014, Neal Seegmiller
license info: see LICENSE.txt

Please report bugs and serious omissions in documentation.
You may contact me with questions.

GETTING STARTED
To understand this code, start with the following demo script:
test/test_simulate.m

Before running the demos:
*enter the commands:
addpath(genpath('/home/user/wmrde/MATLAB'))
cd '/home/user/wmrde/MATLAB'
but replace the string with the correct path for your machine

*check that the strings in the following functions are set to the correct path for your machine:
CADdir.m
resourcedir.m 
datalogdir.m (for calibration only)

HOW TO MODEL A NEW VEHICLE
Write a function to construct the vehicle model. Follow one of the examples in /demo/model
rocky.m demonstrates a rocker-bogie rover.
talon.m demonstrates a tracked vehicle

Several function handles must be set in the WmrModel class.
Write your own functions for the controller and (if applicable) joint constraints.
Use one of the provided functions for the wheel-ground contact and actuator models, or write your own.

OPTIONAL MODIFICATIONS
*A method for geometric collision detection between the wheels/tracks and ground is provided in:
/collision
/surface
You may use a different method, but it must set WheelContactGeom/TrackContactGeom objects for input to the forward kinematics/dynamics functions.

*to switch between Euler angle/quaternion orientation representation, change SIZEORIENT() output
Some code assumes Euler angles.

