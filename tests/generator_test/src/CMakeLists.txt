cmake_minimum_required(VERSION 2.6)

# set the project name
project(test)

#settings
set(CMAKE_CXX_STANDARD 11)
ADD_DEFINITIONS( "-std=c++11" )

# add the executable
add_executable(test test.cpp)
