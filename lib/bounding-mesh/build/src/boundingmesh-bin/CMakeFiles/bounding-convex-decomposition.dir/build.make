# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.9

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ty/Physim/bounding-mesh-0.2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ty/Physim/bounding-mesh-0.2/build

# Include any dependencies generated for this target.
include src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/depend.make

# Include the progress variables for this target.
include src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/progress.make

# Include the compile flags for this target's objects.
include src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/flags.make

src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/bounding-convex-decomposition-cli.cpp.o: src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/flags.make
src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/bounding-convex-decomposition-cli.cpp.o: ../src/boundingmesh-bin/bounding-convex-decomposition-cli.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ty/Physim/bounding-mesh-0.2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/bounding-convex-decomposition-cli.cpp.o"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh-bin && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bounding-convex-decomposition.dir/bounding-convex-decomposition-cli.cpp.o -c /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh-bin/bounding-convex-decomposition-cli.cpp

src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/bounding-convex-decomposition-cli.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bounding-convex-decomposition.dir/bounding-convex-decomposition-cli.cpp.i"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh-bin && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh-bin/bounding-convex-decomposition-cli.cpp > CMakeFiles/bounding-convex-decomposition.dir/bounding-convex-decomposition-cli.cpp.i

src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/bounding-convex-decomposition-cli.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bounding-convex-decomposition.dir/bounding-convex-decomposition-cli.cpp.s"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh-bin && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh-bin/bounding-convex-decomposition-cli.cpp -o CMakeFiles/bounding-convex-decomposition.dir/bounding-convex-decomposition-cli.cpp.s

src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/bounding-convex-decomposition-cli.cpp.o.requires:

.PHONY : src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/bounding-convex-decomposition-cli.cpp.o.requires

src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/bounding-convex-decomposition-cli.cpp.o.provides: src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/bounding-convex-decomposition-cli.cpp.o.requires
	$(MAKE) -f src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/build.make src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/bounding-convex-decomposition-cli.cpp.o.provides.build
.PHONY : src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/bounding-convex-decomposition-cli.cpp.o.provides

src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/bounding-convex-decomposition-cli.cpp.o.provides.build: src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/bounding-convex-decomposition-cli.cpp.o


# Object files for target bounding-convex-decomposition
bounding__convex__decomposition_OBJECTS = \
"CMakeFiles/bounding-convex-decomposition.dir/bounding-convex-decomposition-cli.cpp.o"

# External object files for target bounding-convex-decomposition
bounding__convex__decomposition_EXTERNAL_OBJECTS =

bounding-convex-decomposition: src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/bounding-convex-decomposition-cli.cpp.o
bounding-convex-decomposition: src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/build.make
bounding-convex-decomposition: libboundingmesh.a
bounding-convex-decomposition: /usr/lib/x86_64-linux-gnu/libCGAL.so
bounding-convex-decomposition: /usr/lib/x86_64-linux-gnu/libboost_thread.so
bounding-convex-decomposition: /usr/lib/x86_64-linux-gnu/libboost_system.so
bounding-convex-decomposition: /usr/lib/x86_64-linux-gnu/libboost_chrono.so
bounding-convex-decomposition: /usr/lib/x86_64-linux-gnu/libboost_date_time.so
bounding-convex-decomposition: /usr/lib/x86_64-linux-gnu/libboost_atomic.so
bounding-convex-decomposition: /usr/lib/x86_64-linux-gnu/libpthread.so
bounding-convex-decomposition: /usr/lib/x86_64-linux-gnu/libgmp.so
bounding-convex-decomposition: src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ty/Physim/bounding-mesh-0.2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../bounding-convex-decomposition"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh-bin && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/bounding-convex-decomposition.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/build: bounding-convex-decomposition

.PHONY : src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/build

src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/requires: src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/bounding-convex-decomposition-cli.cpp.o.requires

.PHONY : src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/requires

src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/clean:
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh-bin && $(CMAKE_COMMAND) -P CMakeFiles/bounding-convex-decomposition.dir/cmake_clean.cmake
.PHONY : src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/clean

src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/depend:
	cd /home/ty/Physim/bounding-mesh-0.2/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ty/Physim/bounding-mesh-0.2 /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh-bin /home/ty/Physim/bounding-mesh-0.2/build /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh-bin /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/boundingmesh-bin/CMakeFiles/bounding-convex-decomposition.dir/depend
