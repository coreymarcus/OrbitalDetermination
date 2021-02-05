cmake_minimum_required(VERSION 3.10)

#Set project name
project(OrbitalDetermination)
set(CMAKE_BUILD_TYPE Release)#Debug

#set(CMAKE_CXX_FLAGS "-g -O0 -std=c++17 ${CMAKE_CXX_FLAGS} -Wno-format-security -Wno-literal-suffix")

#add the executables
add_executable(HelloWorld src/example/HelloWorld.cpp)
add_executable(BoostExample src/example/boost_example.cpp)
add_executable(HW0 src/HW0/main.cc)
add_executable(HW1 src/HW1/main.cc
  src/HW1/VehicleState.cc
  src/HW1/Util.cc)

#cpp standard
set(CMAKE_CXX_STANDARD 17)

message("-- PROJECT_SOURCE_DIR : " ${PROJECT_SOURCE_DIR})

# find external libraries
find_package(Boost REQUIRED)
find_package(Eigen3 3.3 REQUIRED)
find_package(Python3 COMPONENTS Development)

#link
target_link_libraries(HW0 Eigen3::Eigen)
target_include_directories(HW0 PRIVATE ${Python3_INCLUDE_DIRS})
target_include_directories(HW0 PUBLIC
  ${PROJECT_SOURCE_DIR}/external/matplotlib-cpp)
target_link_libraries(HW0 Python3::Python)

target_link_libraries(HW1 Eigen3::Eigen)
target_include_directories(HW1 PRIVATE ${Python3_INCLUDE_DIRS})
target_include_directories(HW1 PUBLIC
  ${PROJECT_SOURCE_DIR}/external/matplotlib-cpp)
target_link_libraries(HW1 Python3::Python)
