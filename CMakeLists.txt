cmake_minimum_required(VERSION 3.17)

# set the project name
project(B_to_KPiPiGamma)

# settings
set(CMAKE_CXX_STANDARD 11)
ADD_DEFINITIONS( "-std=c++11" )

# add the libraries
add_library(RandomGen STATIC src/RandomGen.cpp include/RandomGen.hpp)
target_include_directories(RandomGen PUBLIC include)
target_compile_features(RandomGen PUBLIC cxx_std_11)

add_library(Dalitz_Region STATIC src/Dalitz_Region.cpp include/Dalitz_Region.hpp)
target_include_directories(Dalitz_Region PUBLIC include)

add_library(Functions STATIC src/Functions.cpp include/Functions.hpp)
target_include_directories(Functions PUBLIC include)

add_library(Transform STATIC src/Transform.cpp include/Transform.hpp)
target_include_directories(Transform PUBLIC include)

add_library(Kinematics STATIC src/Kinematics.cpp include/Kinematics.hpp)
target_include_directories(Kinematics PUBLIC include)

add_library(ThreeBodies_Kinematics STATIC src/ThreeBodies_Kinematics.cpp include/ThreeBodies_Kinematics.hpp)
target_include_directories(ThreeBodies_Kinematics PUBLIC include)
target_link_libraries(ThreeBodies_Kinematics PUBLIC Kinematics)

add_library(Coupling_Constants STATIC src/Coupling_Constants.cpp include/Coupling_Constants.hpp)
target_include_directories(Coupling_Constants PUBLIC include)

add_library(Particle STATIC src/Particle.cpp include/Particle.hpp)
target_include_directories(Particle PUBLIC include)

add_library(QPCM_A_to_VplusP STATIC src/Quark_Pair_Creation_Model/QPCM_A_to_VplusP.cpp include/Quark_Pair_Creation_Model/QPCM_A_to_VplusP.hpp)
target_include_directories(QPCM_A_to_VplusP PUBLIC include)

add_library(Kres_to_VP STATIC src/Quark_Pair_Creation_Model/Kres_to_VP.cpp include/Quark_Pair_Creation_Model/Kres_to_VP.hpp)
target_include_directories(Kres_to_VP PUBLIC include)
target_link_libraries(Kres_to_VP PUBLIC QPCM_A_to_VplusP)

add_library(Form_Factors STATIC src/Form_Factors.cpp include/Form_Factors.hpp)
target_include_directories(Form_Factors PUBLIC include)
target_link_libraries(Form_Factors PUBLIC Kres_to_VP)

add_library(Helicity_Amplitude STATIC src/Helicity_Amplitude.cpp include/Helicity_Amplitude.hpp)
target_include_directories(Helicity_Amplitude PUBLIC include)
target_link_libraries(Helicity_Amplitude PUBLIC Functions)
target_link_libraries(Helicity_Amplitude PUBLIC ThreeBodies_Kinematics)
target_link_libraries(Helicity_Amplitude PUBLIC Coupling_Constants)
target_link_libraries(Helicity_Amplitude PUBLIC Form_Factors)

add_library(Matrix_Elements STATIC src/Matrix_Elements.cpp include/Matrix_Elements.hpp)
target_include_directories(Matrix_Elements PUBLIC include)
target_link_libraries(Matrix_Elements PUBLIC Helicity_Amplitude)

# add the executable
add_executable(run_B_to_KPiPiGamma_beta apps/Event_Generator_beta.cpp)
target_link_libraries(run_B_to_KPiPiGamma_beta PUBLIC RandomGen)
target_link_libraries(run_B_to_KPiPiGamma_beta PUBLIC Dalitz_Region)
target_link_libraries(run_B_to_KPiPiGamma_beta PUBLIC Functions)
target_link_libraries(run_B_to_KPiPiGamma_beta PUBLIC Kinematics)
target_link_libraries(run_B_to_KPiPiGamma_beta PUBLIC ThreeBodies_Kinematics)
target_link_libraries(run_B_to_KPiPiGamma_beta PUBLIC Coupling_Constants)
target_link_libraries(run_B_to_KPiPiGamma_beta PUBLIC Particle)
target_link_libraries(run_B_to_KPiPiGamma_beta PUBLIC QPCM_A_to_VplusP)
target_link_libraries(run_B_to_KPiPiGamma_beta PUBLIC Kres_to_VP)
target_link_libraries(run_B_to_KPiPiGamma_beta PUBLIC Form_Factors)
target_link_libraries(run_B_to_KPiPiGamma_beta PUBLIC Helicity_Amplitude)
target_link_libraries(run_B_to_KPiPiGamma_beta PUBLIC Matrix_Elements)
target_link_libraries(run_B_to_KPiPiGamma_beta PUBLIC Transform)

#add_executable(run_B_to_KPiPiGamma_FAST apps/Event_Generator_FAST.cpp)
#target_link_libraries(run_B_to_KPiPiGamma_FAST PUBLIC RandomGen)
#target_link_libraries(run_B_to_KPiPiGamma_FAST PUBLIC Dalitz_Region)
#target_link_libraries(run_B_to_KPiPiGamma_FAST PUBLIC Functions)
#target_link_libraries(run_B_to_KPiPiGamma_FAST PUBLIC Kinematics)
#target_link_libraries(run_B_to_KPiPiGamma_FAST PUBLIC ThreeBodies_Kinematics)
#target_link_libraries(run_B_to_KPiPiGamma_FAST PUBLIC Coupling_Constants)
#target_link_libraries(run_B_to_KPiPiGamma_FAST PUBLIC Particle)
#target_link_libraries(run_B_to_KPiPiGamma_FAST PUBLIC QPCM_A_to_VplusP)
#target_link_libraries(run_B_to_KPiPiGamma_FAST PUBLIC Kres_to_VP)
#target_link_libraries(run_B_to_KPiPiGamma_FAST PUBLIC Form_Factors)
#target_link_libraries(run_B_to_KPiPiGamma_FAST PUBLIC Helicity_Amplitude)
#target_link_libraries(run_B_to_KPiPiGamma_FAST PUBLIC Matrix_Elements)
#target_link_libraries(run_B_to_KPiPiGamma_FAST PUBLIC Transform)

#add_executable(complex_test tests/apps_test/complex_test.cpp)
#target_link_libraries(complex_test PUBLIC RandomGen)
#target_link_libraries(complex_test PUBLIC Dalitz_Region)
#target_link_libraries(complex_test PUBLIC Functions)
#target_link_libraries(complex_test PUBLIC Kinematics)
#target_link_libraries(complex_test PUBLIC ThreeBodies_Kinematics)
#target_link_libraries(complex_test PUBLIC Coupling_Constants)
#target_link_libraries(complex_test PUBLIC Particle)
#target_link_libraries(complex_test PUBLIC QPCM_A_to_VplusP)
#target_link_libraries(complex_test PUBLIC Kres_to_VP)
#target_link_libraries(complex_test PUBLIC Form_Factors)
#target_link_libraries(complex_test PUBLIC Helicity_Amplitude)
#target_link_libraries(complex_test PUBLIC Matrix_Elements)

#add_executable(parameter_test tests/apps_test/parameter_test.cpp)
#target_link_libraries(parameter_test PUBLIC RandomGen)
#target_link_libraries(parameter_test PUBLIC Dalitz_Region)
#target_link_libraries(parameter_test PUBLIC Functions)
#target_link_libraries(parameter_test PUBLIC Kinematics)
#target_link_libraries(parameter_test PUBLIC ThreeBodies_Kinematics)
#target_link_libraries(parameter_test PUBLIC Coupling_Constants)
#target_link_libraries(parameter_test PUBLIC Particle)
#target_link_libraries(parameter_test PUBLIC QPCM_A_to_VplusP)
#target_link_libraries(parameter_test PUBLIC Kres_to_VP)
#target_link_libraries(parameter_test PUBLIC Form_Factors)
#target_link_libraries(parameter_test PUBLIC Helicity_Amplitude)
#target_link_libraries(parameter_test PUBLIC Matrix_Elements)
