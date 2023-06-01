# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake3

# The command to remove a file.
RM = /usr/bin/cmake3 -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /afs/desy.de/user/t/tslau/belle2/Tak/B_to_KPiPiGamma

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /afs/desy.de/user/t/tslau/belle2/Tak/B_to_KPiPiGamma/build

# Include any dependencies generated for this target.
include CMakeFiles/Dalitz_Region.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Dalitz_Region.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Dalitz_Region.dir/flags.make

CMakeFiles/Dalitz_Region.dir/src/Dalitz_Region.cpp.o: CMakeFiles/Dalitz_Region.dir/flags.make
CMakeFiles/Dalitz_Region.dir/src/Dalitz_Region.cpp.o: ../src/Dalitz_Region.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/afs/desy.de/user/t/tslau/belle2/Tak/B_to_KPiPiGamma/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Dalitz_Region.dir/src/Dalitz_Region.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Dalitz_Region.dir/src/Dalitz_Region.cpp.o -c /afs/desy.de/user/t/tslau/belle2/Tak/B_to_KPiPiGamma/src/Dalitz_Region.cpp

CMakeFiles/Dalitz_Region.dir/src/Dalitz_Region.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Dalitz_Region.dir/src/Dalitz_Region.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /afs/desy.de/user/t/tslau/belle2/Tak/B_to_KPiPiGamma/src/Dalitz_Region.cpp > CMakeFiles/Dalitz_Region.dir/src/Dalitz_Region.cpp.i

CMakeFiles/Dalitz_Region.dir/src/Dalitz_Region.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Dalitz_Region.dir/src/Dalitz_Region.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /afs/desy.de/user/t/tslau/belle2/Tak/B_to_KPiPiGamma/src/Dalitz_Region.cpp -o CMakeFiles/Dalitz_Region.dir/src/Dalitz_Region.cpp.s

# Object files for target Dalitz_Region
Dalitz_Region_OBJECTS = \
"CMakeFiles/Dalitz_Region.dir/src/Dalitz_Region.cpp.o"

# External object files for target Dalitz_Region
Dalitz_Region_EXTERNAL_OBJECTS =

libDalitz_Region.a: CMakeFiles/Dalitz_Region.dir/src/Dalitz_Region.cpp.o
libDalitz_Region.a: CMakeFiles/Dalitz_Region.dir/build.make
libDalitz_Region.a: CMakeFiles/Dalitz_Region.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/afs/desy.de/user/t/tslau/belle2/Tak/B_to_KPiPiGamma/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libDalitz_Region.a"
	$(CMAKE_COMMAND) -P CMakeFiles/Dalitz_Region.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Dalitz_Region.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Dalitz_Region.dir/build: libDalitz_Region.a

.PHONY : CMakeFiles/Dalitz_Region.dir/build

CMakeFiles/Dalitz_Region.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Dalitz_Region.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Dalitz_Region.dir/clean

CMakeFiles/Dalitz_Region.dir/depend:
	cd /afs/desy.de/user/t/tslau/belle2/Tak/B_to_KPiPiGamma/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /afs/desy.de/user/t/tslau/belle2/Tak/B_to_KPiPiGamma /afs/desy.de/user/t/tslau/belle2/Tak/B_to_KPiPiGamma /afs/desy.de/user/t/tslau/belle2/Tak/B_to_KPiPiGamma/build /afs/desy.de/user/t/tslau/belle2/Tak/B_to_KPiPiGamma/build /afs/desy.de/user/t/tslau/belle2/Tak/B_to_KPiPiGamma/build/CMakeFiles/Dalitz_Region.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Dalitz_Region.dir/depend

