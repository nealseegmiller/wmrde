cmake_minimum_required(VERSION 1.0.0)

project(wmrde)

set(BUILD_DEPS
  cmake_modules #for FindEigen.cmake
  roscpp
  rosunit)

set(RUNTIME_DEPS
  cmake_modules
  roscpp
  rosunit)

include_directories(include)

### Catkin ###
find_package(catkin REQUIRED COMPONENTS ${BUILD_DEPS})
include_directories(${catkin_INCLUDE_DIRS})
catkin_package(CATKIN_DEPENDS ${RUNTIME_DEPS}
  INCLUDE_DIRS include
  LIBRARIES ${PROJECT_NAME})

### Eigen ###
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

#############################################################################
###                                WARNING                                ###
###  Altering CMAKE_CXX_FLAGS violates catkin CMake Style!                ###
###  http://docs.ros.org/hydro/api/catkin/html/user_guide/standards.html  ###
###  ROS Hydro and indigo are targeted at C++03.                          ###
###  http://www.ros.org/reps/rep-0003.html                                ###
#############################################################################
message(WARNING "Altering CMAKE_CXX_FLAGS violates catkin CMake Style! (http://docs.ros.org/hydro/api/catkin/html/user_guide/standards.html).
ROS Hydro and indigo are targeted at C++03. (http://www.ros.org/reps/rep-0003.html)")
set(CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}")

#enable vectorization
set(CMAKE_CXX_FLAGS "-ftree-vectorize ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "-ftree-vectorizer-verbose=1 ${CMAKE_CXX_FLAGS}") #for tests only


### Build Library ###
add_library(${PROJECT_NAME} 
  src/algebra/matrix.cpp
#  src/algebra/spatial.cpp
  src/algebra/transform.cpp
#  src/demo/rockymodel.cpp
#  src/demo/talonmodel.cpp
#  src/demo/terrains.cpp
#  src/demo/zoemodel.cpp
#  src/surface/GridSurf.cpp
#  src/surface/PlaneSurf.cpp
#  src/surface/Surface.cpp
#  src/surface/TriMeshSurf.cpp
#  src/actuator.cpp
#  src/collision.cpp
#  src/dynamics.cpp
#  src/kinematic.cpp
#  src/eigensolve.cpp
#  src/kinematics.cpp
#  src/main.cpp
#  src/state.cpp
#  src/wheelgroundcontact.cpp
#  src/WmrModel.cpp
)
target_link_libraries(${PROJECT_NAME}
  ${catkin_LIBRARIES}
)

catkin_add_gtest(algebra_test src/test/algebra_test.cpp)
target_link_libraries(algebra_test ${PROJECT_NAME})

### Build Nodes ###
#add_executable(joysim_demo_node src/nodes/joysim_demo_node.cpp)
#target_link_libraries(joysim_demo_node ${PROJECT_NAME})

### Install Libraries and Executables ###
install(TARGETS ${PROJECT_NAME}
#  joysim_demo_node
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)

### Install Header Files ###
install(DIRECTORY include/${PROJECT_NAME}/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
)

### Install Launch Files ###
#install(DIRECTORY launch/
#  DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}/launch)

### Install Python Scripts ###
#catkin_install_python(PROGRAMS nodes/TODO
#  DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION})

### Install Resource Files ###
install(DIRECTORY resource/
  DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}/resource)
