cmake_minimum_required(VERSION 3.10)

#Set project name
project(OrbitalDetermination)
set(CMAKE_BUILD_TYPE Release)#Debug

#set(CMAKE_CXX_FLAGS "-g -O0 -std=c++17 ${CMAKE_CXX_FLAGS} -Wno-format-security -Wno-literal-suffix")

#add the executable
add_executable(HelloWorld src/example/HelloWorld.cpp)
add_executable(BoostExample src/example/boost_example.cpp)

#cpp standard
set(CMAKE_CXX_STANDARD 17)

message("-- PROJECT_SOURCE_DIR : " ${PROJECT_SOURCE_DIR})

find_package(Boost)
