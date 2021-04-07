cmake_minimum_required(VERSION 3.10)

#Set project name
project(OrbitalDetermination)
set(CMAKE_BUILD_TYPE Release)#Debug

#set(CMAKE_CXX_FLAGS "-g -O0 -std=c++17 ${CMAKE_CXX_FLAGS} -Wno-format-security -Wno-literal-suffix")

#add the executables
add_executable(HelloWorld src/example/HelloWorld.cpp)
add_executable(BoostExample src/example/boost_example.cpp)
# add_executable(MatlabExample src/example/MatlabExample.cpp)
add_executable(HW0 src/HW0/main.cc)
# add_executable(HW1 src/HW1/main.cc
#  src/HW1/VehicleState.cc
#  src/HW1/Util.cc)
add_executable(HW2 src/HW2/main.cc
	src/HW2/VehicleState.cc
	src/HW2/Util.cc)
add_executable(HW3 src/HW3/main.cc
	src/HW3/VehicleState.cc
	src/HW3/Util.cc
	src/HW3/Estimators.cc)
add_executable(HW4 src/HW4/main.cc
	src/HW4/Util.cc)
add_executable(HW5 src/HW5/main.cc
	src/HW5/Util.cc
	src/HW5/VehicleState.cc
	src/HW5/Estimators.cc)
add_executable(project src/project/main.cc
	src/project/Util.cc
	src/project/VehicleState.cc
	src/project/Estimators.cc)

#cpp standard
set(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

message("-- PROJECT_SOURCE_DIR : " ${PROJECT_SOURCE_DIR})

# find external libraries
find_package(Boost REQUIRED)
find_package(Eigen3 3.3 REQUIRED)
# FIND_PACKAGE(Matlab REQUIRED)
find_package(Python3 COMPONENTS Development)

#link
target_link_libraries(HW0 Eigen3::Eigen)
target_include_directories(HW0 PRIVATE ${Python3_INCLUDE_DIRS})
target_include_directories(HW0 PUBLIC
  ${PROJECT_SOURCE_DIR}/external/matplotlib-cpp)
target_link_libraries(HW0 Python3::Python)

# target_link_libraries(HW1 Eigen3::Eigen)
# target_include_directories(HW1 PRIVATE ${Python3_INCLUDE_DIRS})
# target_include_directories(HW1 PUBLIC
#  ${PROJECT_SOURCE_DIR}/external/matplotlib-cpp)
# target_include_directories(HW1 PUBLIC
#  /usr/local/MATLAB/R2019b/extern/include)
# target_link_libraries(HW1 Python3::Python)
# target_include_directories(HW1 PUBLIC
#   ${Matlab_ROOT_DIR}/extern/include)
# target_link_libraries(HW1 ${Matlab_ROOT_DIR}/extern/bin/glnxa64/libMatlabEngine.so)
# target_link_libraries(HW1 ${Matlab_ROOT_DIR}/extern/bin/glnxa64/libMatlabDataArray.so)

# target_include_directories(MatlabExample PUBLIC
#   /usr/local/MATLAB/R2019b/extern/include)
# target_link_libraries(MatlabExample ${Matlab_ROOT_DIR}/extern/bin/glnxa64/libMatlabEngine.so)
# target_link_libraries(MatlabExample ${Matlab_ROOT_DIR}/extern/bin/glnxa64/libMatlabDataArray.so)

target_link_libraries(HW2 Eigen3::Eigen)

target_link_libraries(HW3 Eigen3::Eigen)

target_link_libraries(HW4 Eigen3::Eigen)

target_link_libraries(HW5 Eigen3::Eigen)

target_link_libraries(project Eigen3::Eigen)