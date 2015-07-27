WMRDE, Wheeled Mobile Robot Dynamics Engine
C++ code

Copyright (c) 2014, Neal Seegmiller
license info: see LICENSE.txt

Please report bugs and serious omissions in documentation.
You may contact me with questions.

BUILD INSTRUCTIONS USING CMAKE IN LINUX

1. First, install the following dependencies.
-Eigen 3.2.1
http://eigen.tuxfamily.org/index.php?title=Main_Page

-OGRE 1.9 (see build instructions below, may work with other versions)
http://www.ogre3d.org/


2. Build the wmrde application with CMake
This uses the CMakeLists.txt and dist/ directory in src/

Open the CMake gui:
$ cmake-gui

Specify the following: (change to the correct paths for your machine)
Where is the source code: /home/neal/Projects/wmrde/src
Where to build the binaries: /home/neal/Projects/wmrde/build

click Configure. 
Specify the generator for this project: Unix Makefiles
Use default native compilers

in the CMake output, verify that the following are found:
OGRE, OIS, Eigen

click Generate

cd into the build directory and make the project:
$ cd /home/neal/Projects/wmrde/build
$ make
$ make install (this just copies into build directory)

execute the Ogre application:
$ cd dist/bin/
$ ./OgreApp

press p key to unpause the simulation
use mouse and w,a,s,d keys to move the camera
press escape key to exit


BUILD INSTRUCTIONS FOR OGRE 1.9

Follow the instructions here:
http://www.ogre3d.org/tikiwiki/tiki-index.php?page=Development&structure=Development
for "Building Ogre"

1. clone the OGRE repositories

make a directory /home/neal/Ogre/source and cd into it
$ hg clone https://bitbucket.org/sinbad/ogre
cd into the new /ogre directory and
$ hg pull && hg update v1-9
make a new subdirectory in /ogre called /Dependencies and cd into it
$ hg clone https://bitbucket.org/cabalistic/ogredeps

2. build the OGRE dependencies
You may not have to do this step separately, but I did to build the OIS library

open the CMake gui,
$ cmake-gui
Where is the source code: /home/neal/Ogre/source/ogre/Dependencies
Where to build the binaries: /home/neal/Ogre/build/Dependencies

click Configure.
CMAKE_BUILD_TYPE Release
click Generate.

$ cd /home/neal/Ogre/build/Dependencies
$ make
$ sudo make install
verify that the following file exists: /usr/local/lib/libOIS.so

3. build OGRE
in the CMake gui,
Where is the source code: /home/neal/Ogre/source/ogre
Where to build the binaries: /home/neal/Ogre/build/ogre-1.9

click Configure.

by default these variables should already be set correctly:
CMAKE_BUILD_TYPE RelWithDebInfo
CMAKE_INSTALL_PREFIX /usr/local
OGREDEPS_BUILD_OIS true
OGRE_BUILD_SAMPLES true
OGRE_INSTALL_SAMPLES true

in the CMake output, verify that the OIS package was found

4. install OGRE
$ cd /home/neal/Ogre/build/ogre-1.9
$ make
$ sudo make install

verify that files were added to /usr/local


GETTING STARTED
To understand this code, start with the following demo:
test_simulate()

Before running, you will need to make the following changes:
in options.h, change the strings in the following functions to the correct paths for your machine:
ResourceDir()
CADdir()

you may change the wheel-ground contact model used by uniformWgc() in wheelgroundcontact.h

For fastest computation:
-build in Release or RelWithDebugInfo configuration
-in collision.h, updateWheelContactGeom() use the _Root() function instead of the _Discretize() function
-#define FIXED_NROWS, FIXED_NCOLS, FIXED_N, in eigensolve.h to used fixed size matrices


HOW TO MODEL A NEW VEHICLE
Write a function to construct the vehicle model. Follow one of the examples in /demo
rockymodel.cpp demonstrates a rocker-bogie rover
talonmodel.cpp demonstrates a tracked vehicle

Several function handles must be set in the WmrModel class.
Write your own functions for the controller and (if applicable) joint constraints.
Use one of the provided functions for the wheel-ground contact and actuator models, or write your own.


OPTIONAL MODIFICATIONS
*A method for geometric collision detection between the wheels/tracks and ground is provided in:
collision.h and .cpp
/surface directory
You may use a different method, but it must set ContactGeom objects for input to the forward kinematics/dynamics functions.

*To eliminate the dependency on Eigen, do not include eigensolve.h However, you must write your own function definitions for solve() and chol(). Some of the functions in test.cpp also use Eigen for comparison

*To eliminate the dependency on OGRE, set #define WMRSIM_ENABLE_ANIMATION 0
do not include the /animate directory


BUILD INSTRUCTIONS FOR VISUAL STUDIO 2013
Using CMake is the preferred way to build wmrde applications, but in the past I have also built using Visual Studio following these steps.
Make a new Visual Studio solution.

1. Add files to the project:
In Solution Explorer, right click on Header Files (or Source Files)
Add -> Existing Item
add all files from wmrde/include and /src

2. Include code from this repository as follows:
In Solution Explorer, right click on project, select Properties
change:
Configuration: All Configurations
Platform: All Platforms

under Configuration Properties -> C/C++ -> General 
add the path to the wmrde/include directory to Additional Include Directories

3. Include Eigen:
add the path to eigen to Additional Include Directories

4. Include OGRE:
Follow the instructions here:
http://www.ogre3d.org/tikiwiki/tiki-index.php?page=Development&structure=Development
for "Setting Up An Application"

set the OGRE_HOME environment variable to point to your OGRE sdk directory.
.cfg and .material files in the sdk should be modified according to files in the wmrde/src/dist directory

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

