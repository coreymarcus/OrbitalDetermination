cmake_minimum_required(VERSION 3.10)

#Set project name
project(OrbitalDetermination)
set(CMAKE_BUILD_TYPE Release)#Debug

#set(CMAKE_CXX_FLAGS "-g -O0 -std=c++17 ${CMAKE_CXX_FLAGS} -Wno-format-security -Wno-literal-suffix")

#add the executables
add_executable(HelloWorld src/example/HelloWorld.cpp)
add_executable(BoostExample src/example/boost_example.cpp)
add_executable(HW0 src/HW0/main.cc)

# include directory for externals
include_directories(external/matplotlib-cpp)

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
target_link_libraries(HW0 Python3::Python)
