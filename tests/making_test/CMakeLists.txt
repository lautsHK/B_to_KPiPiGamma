cmake_minimum_required(VERSION 2.6)

# set the project name
project(B_to_KPiPiGamma)

# settings
set(CMAKE_CXX_STANDARD 11)
ADD_DEFINITIONS( "-std=c++11" )

# add the libraries
add_library(RandomGen STATIC include/RandomGen.hpp src/RandomGen.cpp)
target_include_directories(RandomGen PUBLIC include)
target_compile_features(RandomGen PUBLIC cxx_std_11)

# add the executable
add_executable(B_to_KPiPiGamma src/generator_test.cpp)
target_link_libraries(B_to_KPiPiGamma PUBLIC RandomGen)
