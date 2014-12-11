WMRDE, Wheeled Mobile Robot Dynamics Engine
C++ code

Copyright (c) 2014, Neal Seegmiller
license info: see LICENSE.txt


BUILD INSTRUCTIONS
Dependencies:
-Eigen 3.2.1
-OGRE 1.9

later versions may work too


BUILD INSTRUCTIONS FOR VISUAL STUDIO 2013
1. Make a new Visual Studio solution. I will refer to it as WMRDE.

include code from this repository (wmrde/src) as follows:
In Solution Explorer, right click on project, select Properties
change:
Configuration: All Configurations
Platform: All Platforms

under Configuration Properties -> C/C++ -> General 
add the path to the wmrde/src directory to Additional Include Directories

2. add files to the project
In Solution Explorer, right click on Header Files (or Source Files)
Add -> Existing Item

add all files from wmrde/src, or just the ones you require

3. Download and include Eigen:
http://eigen.tuxfamily.org/index.php?title=Main_Page
add the path to the eigen directory to Additional Include Directories

4. Install OGRE
Follow the instructions here:
http://www.ogre3d.org/tikiwiki/tiki-index.php?page=Development&structure=Development
for:
"Installing the Ogre SDK"
"Building Ogre" (to use CMake if there is no SDK download for your IDE)

64 bit build is recommended
In CMake choose "Visual Studio 12 Win64" as generator
In Visual Studio choose x64 as Active solution platform

Now follow the instructions for "Setting Up An Application" in your WMRDE solution:
set the OGRE_HOME environment variable as instructed, it should point to the new /sdk directory
add the following to Additional Include Directories:
$(OGRE_HOME)\include\OGRE\Overlay
$(OGRE_HOME)\boost
$(OGRE_HOME)\Samples\Common\include
$(OGRE_HOME)\include\OGRE
$(OGRE_HOME)\include\OIS
$(OGRE_HOME)\include

under Configuration Properties -> Linker -> General
add the following to Additional Library Directories:
$(OGRE_HOME)\boost\lib
$(OGRE_HOME)\lib\$(Configuration)
$(OGRE_HOME)\lib\DebugDoubleDLL (for Debug configuration)
$(OGRE_HOME)\lib\ReleaseDoubleDLL (for Release configuration)

under Configuration Properties -> Linker -> General
add the following to Additional Dependencies:
(for Debug configuration)
OgreMain_d.lib
OIS_d.lib
OgreOverlay_d.lib
(for Release configuration)
OgreMain.lib
OIS.lib
OgreOverlay.lib

5. to enable Auto-Vectorizer (only in Release configuration):
under Configuration Properties -> C/C++ -> Command Line
add: /Qvec-report:1


BUILD INSTRUCTIONS FOR LINUX
TODO


ADDITIONAL SETUP INSTRUCTIONS

*for OGRE
edit /sdk/bin/Debug/ogre.cfg (and /sdk/bin/Release/ogre.cfg)
VSync=Yes
This limits the framerate.

edit /sdk/bin/Debug/resources_d.cfg (and /sdk/bin/Release/resources.cfg)
add the following:
FileSystem=../../Media/materials/programs/GLSL
FileSystem=../../Media/materials/scripts/MyMaterials
FileSystem=C:/Users/nseegmil.NREC-014635/Documents/Thesis/textures

but change to the correct paths for your machine.
create the /MyMaterials directory yourself, and copy MyMaterials.material from wmrde/src/animate into it.
wmrde/textures is included in this repo.

You may comment out the [Popular] FileSystem lines for faster startup


GETTING STARTED
To understand this code, start with the following demo:
test_simulate()

Before running, you will need to make the following changes:
in options.h, change the strings in the following functions to the correct paths for your machine:
ResourceDir()
CADdir()

you may change the wheel-ground contact model used by uniformWgc() in wheelgroundcontact.h

For fastest computation:
-Release build
-in collision.h, updateWheelContactGeom() use the _Root() function instead of the _Discretize() function
-#define FIXED_NROWS, FIXED_NCOLS, FIXED_N, in eigensolve.h to used fixed size matrices


HOW TO MODEL A NEW VEHICLE
Write a function to construct the vehicle model. Follow one of the examples in /demo
e.g. rockymodel, zoemodel, talonmodel

Several function handles must be set in the WmrModel class.
Write your own functions for the controller and (if applicable) joint constraints.
Use one of the provided functions for the wheel-ground contact and actuator models, or write your own.

You may contact me for help.

OPTIONAL MODIFICATIONS
*A method for geometric collision detection between the wheels/tracks and ground is provided in:
collision.h, .cpp
/surface directory
You may use a different method, but it must set ContactGeom objects for input to the forward kinematics/dynamics functions.

*To eliminate the dependency on Eigen, do not include eigensolve.h, .cpp. However, you must write your own function definitions for solve(), subset(), and chol(). Some of the functions in test.cpp also use Eigen for comparison

*To eliminate the dependency on OGRE, set #define WMRSIM_ENABLE_ANIMATION 0
do not include the /animate directory



