cmake_minimum_required(VERSION 1.0.0)

project(wmrde)

set(BUILD_DEPS
  cmake_modules #for FindEigen.cmake
  geometry_msgs
  roscpp
  rosunit
  sensor_msgs
  tf)

set(RUNTIME_DEPS
  cmake_modules
  geometry_msgs
  roscpp
  rosunit
  sensor_msgs
  tf)

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
#set(CMAKE_CXX_FLAGS "-ftree-vectorize ${CMAKE_CXX_FLAGS}")
#set(CMAKE_CXX_FLAGS "-ftree-vectorizer-verbose=1 ${CMAKE_CXX_FLAGS}") #for tests only


### Build Library ###
add_library(${PROJECT_NAME} 
  src/algebra/spatial.cpp
#  src/demo/rockymodel.cpp
#  src/demo/talonmodel.cpp
#  src/demo/terrains.cpp
  src/demo/zoemodel.cpp
  src/surface/grid_surface.cpp
  src/surface/plane_surface.cpp
  src/surface/surface.cpp
#  src/surface/TriMeshSurf.cpp
#  src/actuator.cpp
  src/collision.cpp
#  src/dynamics.cpp
#  src/eigensolve.cpp
#  src/kinematics.cpp
#  src/main.cpp
#  src/wheelgroundcontact.cpp
  src/wmrde_ros_interface.cpp
  src/wmrmodel.cpp
  src/wmrstate.cpp
)
target_link_libraries(${PROJECT_NAME}
  ${catkin_LIBRARIES}
)

#gtest
catkin_add_gtest(algebra_test src/test/algebra_test.cpp)
target_link_libraries(algebra_test ${PROJECT_NAME})

catkin_add_gtest(gridded_interpolant_test src/test/gridded_interpolant_test.cpp)
target_link_libraries(gridded_interpolant_test ${PROJECT_NAME})

catkin_add_gtest(wmrmodel_test src/test/wmrmodel_test.cpp)
target_link_libraries(wmrmodel_test ${PROJECT_NAME})

#benchmark
#depends on: https://github.com/google/benchmark
add_executable(algebra_benchmark src/test/algebra_benchmark.cpp)
target_link_libraries(algebra_benchmark ${PROJECT_NAME} benchmark)

### Build Nodes ###
add_executable(joy_teleop src/nodes/joy_teleop.cpp)
target_link_libraries(joy_teleop ${PROJECT_NAME})

add_executable(collision_joy_demo src/nodes/collision_joy_demo.cpp)
target_link_libraries(collision_joy_demo ${PROJECT_NAME})

add_executable(wmrmodel_joy_demo src/nodes/wmrmodel_joy_demo.cpp)
target_link_libraries(wmrmodel_joy_demo ${PROJECT_NAME})

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
