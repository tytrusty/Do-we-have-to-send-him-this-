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
include src/boundingmesh/CMakeFiles/boundingmesh.dir/depend.make

# Include the progress variables for this target.
include src/boundingmesh/CMakeFiles/boundingmesh.dir/progress.make

# Include the compile flags for this target's objects.
include src/boundingmesh/CMakeFiles/boundingmesh.dir/flags.make

src/boundingmesh/CMakeFiles/boundingmesh.dir/ContractionUtils.cpp.o: src/boundingmesh/CMakeFiles/boundingmesh.dir/flags.make
src/boundingmesh/CMakeFiles/boundingmesh.dir/ContractionUtils.cpp.o: ../src/boundingmesh/ContractionUtils.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ty/Physim/bounding-mesh-0.2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/boundingmesh/CMakeFiles/boundingmesh.dir/ContractionUtils.cpp.o"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/boundingmesh.dir/ContractionUtils.cpp.o -c /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/ContractionUtils.cpp

src/boundingmesh/CMakeFiles/boundingmesh.dir/ContractionUtils.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/boundingmesh.dir/ContractionUtils.cpp.i"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/ContractionUtils.cpp > CMakeFiles/boundingmesh.dir/ContractionUtils.cpp.i

src/boundingmesh/CMakeFiles/boundingmesh.dir/ContractionUtils.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/boundingmesh.dir/ContractionUtils.cpp.s"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/ContractionUtils.cpp -o CMakeFiles/boundingmesh.dir/ContractionUtils.cpp.s

src/boundingmesh/CMakeFiles/boundingmesh.dir/ContractionUtils.cpp.o.requires:

.PHONY : src/boundingmesh/CMakeFiles/boundingmesh.dir/ContractionUtils.cpp.o.requires

src/boundingmesh/CMakeFiles/boundingmesh.dir/ContractionUtils.cpp.o.provides: src/boundingmesh/CMakeFiles/boundingmesh.dir/ContractionUtils.cpp.o.requires
	$(MAKE) -f src/boundingmesh/CMakeFiles/boundingmesh.dir/build.make src/boundingmesh/CMakeFiles/boundingmesh.dir/ContractionUtils.cpp.o.provides.build
.PHONY : src/boundingmesh/CMakeFiles/boundingmesh.dir/ContractionUtils.cpp.o.provides

src/boundingmesh/CMakeFiles/boundingmesh.dir/ContractionUtils.cpp.o.provides.build: src/boundingmesh/CMakeFiles/boundingmesh.dir/ContractionUtils.cpp.o


src/boundingmesh/CMakeFiles/boundingmesh.dir/Decimator.cpp.o: src/boundingmesh/CMakeFiles/boundingmesh.dir/flags.make
src/boundingmesh/CMakeFiles/boundingmesh.dir/Decimator.cpp.o: ../src/boundingmesh/Decimator.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ty/Physim/bounding-mesh-0.2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/boundingmesh/CMakeFiles/boundingmesh.dir/Decimator.cpp.o"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/boundingmesh.dir/Decimator.cpp.o -c /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/Decimator.cpp

src/boundingmesh/CMakeFiles/boundingmesh.dir/Decimator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/boundingmesh.dir/Decimator.cpp.i"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/Decimator.cpp > CMakeFiles/boundingmesh.dir/Decimator.cpp.i

src/boundingmesh/CMakeFiles/boundingmesh.dir/Decimator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/boundingmesh.dir/Decimator.cpp.s"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/Decimator.cpp -o CMakeFiles/boundingmesh.dir/Decimator.cpp.s

src/boundingmesh/CMakeFiles/boundingmesh.dir/Decimator.cpp.o.requires:

.PHONY : src/boundingmesh/CMakeFiles/boundingmesh.dir/Decimator.cpp.o.requires

src/boundingmesh/CMakeFiles/boundingmesh.dir/Decimator.cpp.o.provides: src/boundingmesh/CMakeFiles/boundingmesh.dir/Decimator.cpp.o.requires
	$(MAKE) -f src/boundingmesh/CMakeFiles/boundingmesh.dir/build.make src/boundingmesh/CMakeFiles/boundingmesh.dir/Decimator.cpp.o.provides.build
.PHONY : src/boundingmesh/CMakeFiles/boundingmesh.dir/Decimator.cpp.o.provides

src/boundingmesh/CMakeFiles/boundingmesh.dir/Decimator.cpp.o.provides.build: src/boundingmesh/CMakeFiles/boundingmesh.dir/Decimator.cpp.o


src/boundingmesh/CMakeFiles/boundingmesh.dir/FileUtils.cpp.o: src/boundingmesh/CMakeFiles/boundingmesh.dir/flags.make
src/boundingmesh/CMakeFiles/boundingmesh.dir/FileUtils.cpp.o: ../src/boundingmesh/FileUtils.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ty/Physim/bounding-mesh-0.2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object src/boundingmesh/CMakeFiles/boundingmesh.dir/FileUtils.cpp.o"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/boundingmesh.dir/FileUtils.cpp.o -c /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/FileUtils.cpp

src/boundingmesh/CMakeFiles/boundingmesh.dir/FileUtils.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/boundingmesh.dir/FileUtils.cpp.i"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/FileUtils.cpp > CMakeFiles/boundingmesh.dir/FileUtils.cpp.i

src/boundingmesh/CMakeFiles/boundingmesh.dir/FileUtils.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/boundingmesh.dir/FileUtils.cpp.s"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/FileUtils.cpp -o CMakeFiles/boundingmesh.dir/FileUtils.cpp.s

src/boundingmesh/CMakeFiles/boundingmesh.dir/FileUtils.cpp.o.requires:

.PHONY : src/boundingmesh/CMakeFiles/boundingmesh.dir/FileUtils.cpp.o.requires

src/boundingmesh/CMakeFiles/boundingmesh.dir/FileUtils.cpp.o.provides: src/boundingmesh/CMakeFiles/boundingmesh.dir/FileUtils.cpp.o.requires
	$(MAKE) -f src/boundingmesh/CMakeFiles/boundingmesh.dir/build.make src/boundingmesh/CMakeFiles/boundingmesh.dir/FileUtils.cpp.o.provides.build
.PHONY : src/boundingmesh/CMakeFiles/boundingmesh.dir/FileUtils.cpp.o.provides

src/boundingmesh/CMakeFiles/boundingmesh.dir/FileUtils.cpp.o.provides.build: src/boundingmesh/CMakeFiles/boundingmesh.dir/FileUtils.cpp.o


src/boundingmesh/CMakeFiles/boundingmesh.dir/Mesh.cpp.o: src/boundingmesh/CMakeFiles/boundingmesh.dir/flags.make
src/boundingmesh/CMakeFiles/boundingmesh.dir/Mesh.cpp.o: ../src/boundingmesh/Mesh.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ty/Physim/bounding-mesh-0.2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object src/boundingmesh/CMakeFiles/boundingmesh.dir/Mesh.cpp.o"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/boundingmesh.dir/Mesh.cpp.o -c /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/Mesh.cpp

src/boundingmesh/CMakeFiles/boundingmesh.dir/Mesh.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/boundingmesh.dir/Mesh.cpp.i"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/Mesh.cpp > CMakeFiles/boundingmesh.dir/Mesh.cpp.i

src/boundingmesh/CMakeFiles/boundingmesh.dir/Mesh.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/boundingmesh.dir/Mesh.cpp.s"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/Mesh.cpp -o CMakeFiles/boundingmesh.dir/Mesh.cpp.s

src/boundingmesh/CMakeFiles/boundingmesh.dir/Mesh.cpp.o.requires:

.PHONY : src/boundingmesh/CMakeFiles/boundingmesh.dir/Mesh.cpp.o.requires

src/boundingmesh/CMakeFiles/boundingmesh.dir/Mesh.cpp.o.provides: src/boundingmesh/CMakeFiles/boundingmesh.dir/Mesh.cpp.o.requires
	$(MAKE) -f src/boundingmesh/CMakeFiles/boundingmesh.dir/build.make src/boundingmesh/CMakeFiles/boundingmesh.dir/Mesh.cpp.o.provides.build
.PHONY : src/boundingmesh/CMakeFiles/boundingmesh.dir/Mesh.cpp.o.provides

src/boundingmesh/CMakeFiles/boundingmesh.dir/Mesh.cpp.o.provides.build: src/boundingmesh/CMakeFiles/boundingmesh.dir/Mesh.cpp.o


src/boundingmesh/CMakeFiles/boundingmesh.dir/MetricGenerator.cpp.o: src/boundingmesh/CMakeFiles/boundingmesh.dir/flags.make
src/boundingmesh/CMakeFiles/boundingmesh.dir/MetricGenerator.cpp.o: ../src/boundingmesh/MetricGenerator.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ty/Physim/bounding-mesh-0.2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object src/boundingmesh/CMakeFiles/boundingmesh.dir/MetricGenerator.cpp.o"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/boundingmesh.dir/MetricGenerator.cpp.o -c /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/MetricGenerator.cpp

src/boundingmesh/CMakeFiles/boundingmesh.dir/MetricGenerator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/boundingmesh.dir/MetricGenerator.cpp.i"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/MetricGenerator.cpp > CMakeFiles/boundingmesh.dir/MetricGenerator.cpp.i

src/boundingmesh/CMakeFiles/boundingmesh.dir/MetricGenerator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/boundingmesh.dir/MetricGenerator.cpp.s"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/MetricGenerator.cpp -o CMakeFiles/boundingmesh.dir/MetricGenerator.cpp.s

src/boundingmesh/CMakeFiles/boundingmesh.dir/MetricGenerator.cpp.o.requires:

.PHONY : src/boundingmesh/CMakeFiles/boundingmesh.dir/MetricGenerator.cpp.o.requires

src/boundingmesh/CMakeFiles/boundingmesh.dir/MetricGenerator.cpp.o.provides: src/boundingmesh/CMakeFiles/boundingmesh.dir/MetricGenerator.cpp.o.requires
	$(MAKE) -f src/boundingmesh/CMakeFiles/boundingmesh.dir/build.make src/boundingmesh/CMakeFiles/boundingmesh.dir/MetricGenerator.cpp.o.provides.build
.PHONY : src/boundingmesh/CMakeFiles/boundingmesh.dir/MetricGenerator.cpp.o.provides

src/boundingmesh/CMakeFiles/boundingmesh.dir/MetricGenerator.cpp.o.provides.build: src/boundingmesh/CMakeFiles/boundingmesh.dir/MetricGenerator.cpp.o


src/boundingmesh/CMakeFiles/boundingmesh.dir/Primitives.cpp.o: src/boundingmesh/CMakeFiles/boundingmesh.dir/flags.make
src/boundingmesh/CMakeFiles/boundingmesh.dir/Primitives.cpp.o: ../src/boundingmesh/Primitives.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ty/Physim/bounding-mesh-0.2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object src/boundingmesh/CMakeFiles/boundingmesh.dir/Primitives.cpp.o"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/boundingmesh.dir/Primitives.cpp.o -c /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/Primitives.cpp

src/boundingmesh/CMakeFiles/boundingmesh.dir/Primitives.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/boundingmesh.dir/Primitives.cpp.i"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/Primitives.cpp > CMakeFiles/boundingmesh.dir/Primitives.cpp.i

src/boundingmesh/CMakeFiles/boundingmesh.dir/Primitives.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/boundingmesh.dir/Primitives.cpp.s"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/Primitives.cpp -o CMakeFiles/boundingmesh.dir/Primitives.cpp.s

src/boundingmesh/CMakeFiles/boundingmesh.dir/Primitives.cpp.o.requires:

.PHONY : src/boundingmesh/CMakeFiles/boundingmesh.dir/Primitives.cpp.o.requires

src/boundingmesh/CMakeFiles/boundingmesh.dir/Primitives.cpp.o.provides: src/boundingmesh/CMakeFiles/boundingmesh.dir/Primitives.cpp.o.requires
	$(MAKE) -f src/boundingmesh/CMakeFiles/boundingmesh.dir/build.make src/boundingmesh/CMakeFiles/boundingmesh.dir/Primitives.cpp.o.provides.build
.PHONY : src/boundingmesh/CMakeFiles/boundingmesh.dir/Primitives.cpp.o.provides

src/boundingmesh/CMakeFiles/boundingmesh.dir/Primitives.cpp.o.provides.build: src/boundingmesh/CMakeFiles/boundingmesh.dir/Primitives.cpp.o


src/boundingmesh/CMakeFiles/boundingmesh.dir/SegmenterDownsampling.cpp.o: src/boundingmesh/CMakeFiles/boundingmesh.dir/flags.make
src/boundingmesh/CMakeFiles/boundingmesh.dir/SegmenterDownsampling.cpp.o: ../src/boundingmesh/SegmenterDownsampling.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ty/Physim/bounding-mesh-0.2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object src/boundingmesh/CMakeFiles/boundingmesh.dir/SegmenterDownsampling.cpp.o"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/boundingmesh.dir/SegmenterDownsampling.cpp.o -c /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/SegmenterDownsampling.cpp

src/boundingmesh/CMakeFiles/boundingmesh.dir/SegmenterDownsampling.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/boundingmesh.dir/SegmenterDownsampling.cpp.i"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/SegmenterDownsampling.cpp > CMakeFiles/boundingmesh.dir/SegmenterDownsampling.cpp.i

src/boundingmesh/CMakeFiles/boundingmesh.dir/SegmenterDownsampling.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/boundingmesh.dir/SegmenterDownsampling.cpp.s"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/SegmenterDownsampling.cpp -o CMakeFiles/boundingmesh.dir/SegmenterDownsampling.cpp.s

src/boundingmesh/CMakeFiles/boundingmesh.dir/SegmenterDownsampling.cpp.o.requires:

.PHONY : src/boundingmesh/CMakeFiles/boundingmesh.dir/SegmenterDownsampling.cpp.o.requires

src/boundingmesh/CMakeFiles/boundingmesh.dir/SegmenterDownsampling.cpp.o.provides: src/boundingmesh/CMakeFiles/boundingmesh.dir/SegmenterDownsampling.cpp.o.requires
	$(MAKE) -f src/boundingmesh/CMakeFiles/boundingmesh.dir/build.make src/boundingmesh/CMakeFiles/boundingmesh.dir/SegmenterDownsampling.cpp.o.provides.build
.PHONY : src/boundingmesh/CMakeFiles/boundingmesh.dir/SegmenterDownsampling.cpp.o.provides

src/boundingmesh/CMakeFiles/boundingmesh.dir/SegmenterDownsampling.cpp.o.provides.build: src/boundingmesh/CMakeFiles/boundingmesh.dir/SegmenterDownsampling.cpp.o


src/boundingmesh/CMakeFiles/boundingmesh.dir/SegmenterSimple.cpp.o: src/boundingmesh/CMakeFiles/boundingmesh.dir/flags.make
src/boundingmesh/CMakeFiles/boundingmesh.dir/SegmenterSimple.cpp.o: ../src/boundingmesh/SegmenterSimple.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ty/Physim/bounding-mesh-0.2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object src/boundingmesh/CMakeFiles/boundingmesh.dir/SegmenterSimple.cpp.o"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/boundingmesh.dir/SegmenterSimple.cpp.o -c /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/SegmenterSimple.cpp

src/boundingmesh/CMakeFiles/boundingmesh.dir/SegmenterSimple.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/boundingmesh.dir/SegmenterSimple.cpp.i"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/SegmenterSimple.cpp > CMakeFiles/boundingmesh.dir/SegmenterSimple.cpp.i

src/boundingmesh/CMakeFiles/boundingmesh.dir/SegmenterSimple.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/boundingmesh.dir/SegmenterSimple.cpp.s"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/SegmenterSimple.cpp -o CMakeFiles/boundingmesh.dir/SegmenterSimple.cpp.s

src/boundingmesh/CMakeFiles/boundingmesh.dir/SegmenterSimple.cpp.o.requires:

.PHONY : src/boundingmesh/CMakeFiles/boundingmesh.dir/SegmenterSimple.cpp.o.requires

src/boundingmesh/CMakeFiles/boundingmesh.dir/SegmenterSimple.cpp.o.provides: src/boundingmesh/CMakeFiles/boundingmesh.dir/SegmenterSimple.cpp.o.requires
	$(MAKE) -f src/boundingmesh/CMakeFiles/boundingmesh.dir/build.make src/boundingmesh/CMakeFiles/boundingmesh.dir/SegmenterSimple.cpp.o.provides.build
.PHONY : src/boundingmesh/CMakeFiles/boundingmesh.dir/SegmenterSimple.cpp.o.provides

src/boundingmesh/CMakeFiles/boundingmesh.dir/SegmenterSimple.cpp.o.provides.build: src/boundingmesh/CMakeFiles/boundingmesh.dir/SegmenterSimple.cpp.o


src/boundingmesh/CMakeFiles/boundingmesh.dir/Split.cpp.o: src/boundingmesh/CMakeFiles/boundingmesh.dir/flags.make
src/boundingmesh/CMakeFiles/boundingmesh.dir/Split.cpp.o: ../src/boundingmesh/Split.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ty/Physim/bounding-mesh-0.2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object src/boundingmesh/CMakeFiles/boundingmesh.dir/Split.cpp.o"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/boundingmesh.dir/Split.cpp.o -c /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/Split.cpp

src/boundingmesh/CMakeFiles/boundingmesh.dir/Split.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/boundingmesh.dir/Split.cpp.i"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/Split.cpp > CMakeFiles/boundingmesh.dir/Split.cpp.i

src/boundingmesh/CMakeFiles/boundingmesh.dir/Split.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/boundingmesh.dir/Split.cpp.s"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/Split.cpp -o CMakeFiles/boundingmesh.dir/Split.cpp.s

src/boundingmesh/CMakeFiles/boundingmesh.dir/Split.cpp.o.requires:

.PHONY : src/boundingmesh/CMakeFiles/boundingmesh.dir/Split.cpp.o.requires

src/boundingmesh/CMakeFiles/boundingmesh.dir/Split.cpp.o.provides: src/boundingmesh/CMakeFiles/boundingmesh.dir/Split.cpp.o.requires
	$(MAKE) -f src/boundingmesh/CMakeFiles/boundingmesh.dir/build.make src/boundingmesh/CMakeFiles/boundingmesh.dir/Split.cpp.o.provides.build
.PHONY : src/boundingmesh/CMakeFiles/boundingmesh.dir/Split.cpp.o.provides

src/boundingmesh/CMakeFiles/boundingmesh.dir/Split.cpp.o.provides.build: src/boundingmesh/CMakeFiles/boundingmesh.dir/Split.cpp.o


src/boundingmesh/CMakeFiles/boundingmesh.dir/VoxelSet.cpp.o: src/boundingmesh/CMakeFiles/boundingmesh.dir/flags.make
src/boundingmesh/CMakeFiles/boundingmesh.dir/VoxelSet.cpp.o: ../src/boundingmesh/VoxelSet.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ty/Physim/bounding-mesh-0.2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object src/boundingmesh/CMakeFiles/boundingmesh.dir/VoxelSet.cpp.o"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/boundingmesh.dir/VoxelSet.cpp.o -c /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/VoxelSet.cpp

src/boundingmesh/CMakeFiles/boundingmesh.dir/VoxelSet.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/boundingmesh.dir/VoxelSet.cpp.i"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/VoxelSet.cpp > CMakeFiles/boundingmesh.dir/VoxelSet.cpp.i

src/boundingmesh/CMakeFiles/boundingmesh.dir/VoxelSet.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/boundingmesh.dir/VoxelSet.cpp.s"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/VoxelSet.cpp -o CMakeFiles/boundingmesh.dir/VoxelSet.cpp.s

src/boundingmesh/CMakeFiles/boundingmesh.dir/VoxelSet.cpp.o.requires:

.PHONY : src/boundingmesh/CMakeFiles/boundingmesh.dir/VoxelSet.cpp.o.requires

src/boundingmesh/CMakeFiles/boundingmesh.dir/VoxelSet.cpp.o.provides: src/boundingmesh/CMakeFiles/boundingmesh.dir/VoxelSet.cpp.o.requires
	$(MAKE) -f src/boundingmesh/CMakeFiles/boundingmesh.dir/build.make src/boundingmesh/CMakeFiles/boundingmesh.dir/VoxelSet.cpp.o.provides.build
.PHONY : src/boundingmesh/CMakeFiles/boundingmesh.dir/VoxelSet.cpp.o.provides

src/boundingmesh/CMakeFiles/boundingmesh.dir/VoxelSet.cpp.o.provides.build: src/boundingmesh/CMakeFiles/boundingmesh.dir/VoxelSet.cpp.o


src/boundingmesh/CMakeFiles/boundingmesh.dir/VoxelSubset.cpp.o: src/boundingmesh/CMakeFiles/boundingmesh.dir/flags.make
src/boundingmesh/CMakeFiles/boundingmesh.dir/VoxelSubset.cpp.o: ../src/boundingmesh/VoxelSubset.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ty/Physim/bounding-mesh-0.2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object src/boundingmesh/CMakeFiles/boundingmesh.dir/VoxelSubset.cpp.o"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/boundingmesh.dir/VoxelSubset.cpp.o -c /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/VoxelSubset.cpp

src/boundingmesh/CMakeFiles/boundingmesh.dir/VoxelSubset.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/boundingmesh.dir/VoxelSubset.cpp.i"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/VoxelSubset.cpp > CMakeFiles/boundingmesh.dir/VoxelSubset.cpp.i

src/boundingmesh/CMakeFiles/boundingmesh.dir/VoxelSubset.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/boundingmesh.dir/VoxelSubset.cpp.s"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh/VoxelSubset.cpp -o CMakeFiles/boundingmesh.dir/VoxelSubset.cpp.s

src/boundingmesh/CMakeFiles/boundingmesh.dir/VoxelSubset.cpp.o.requires:

.PHONY : src/boundingmesh/CMakeFiles/boundingmesh.dir/VoxelSubset.cpp.o.requires

src/boundingmesh/CMakeFiles/boundingmesh.dir/VoxelSubset.cpp.o.provides: src/boundingmesh/CMakeFiles/boundingmesh.dir/VoxelSubset.cpp.o.requires
	$(MAKE) -f src/boundingmesh/CMakeFiles/boundingmesh.dir/build.make src/boundingmesh/CMakeFiles/boundingmesh.dir/VoxelSubset.cpp.o.provides.build
.PHONY : src/boundingmesh/CMakeFiles/boundingmesh.dir/VoxelSubset.cpp.o.provides

src/boundingmesh/CMakeFiles/boundingmesh.dir/VoxelSubset.cpp.o.provides.build: src/boundingmesh/CMakeFiles/boundingmesh.dir/VoxelSubset.cpp.o


src/boundingmesh/CMakeFiles/boundingmesh.dir/__/__/thirdparty/EigenQP.cpp.o: src/boundingmesh/CMakeFiles/boundingmesh.dir/flags.make
src/boundingmesh/CMakeFiles/boundingmesh.dir/__/__/thirdparty/EigenQP.cpp.o: ../thirdparty/EigenQP.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ty/Physim/bounding-mesh-0.2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object src/boundingmesh/CMakeFiles/boundingmesh.dir/__/__/thirdparty/EigenQP.cpp.o"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/boundingmesh.dir/__/__/thirdparty/EigenQP.cpp.o -c /home/ty/Physim/bounding-mesh-0.2/thirdparty/EigenQP.cpp

src/boundingmesh/CMakeFiles/boundingmesh.dir/__/__/thirdparty/EigenQP.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/boundingmesh.dir/__/__/thirdparty/EigenQP.cpp.i"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ty/Physim/bounding-mesh-0.2/thirdparty/EigenQP.cpp > CMakeFiles/boundingmesh.dir/__/__/thirdparty/EigenQP.cpp.i

src/boundingmesh/CMakeFiles/boundingmesh.dir/__/__/thirdparty/EigenQP.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/boundingmesh.dir/__/__/thirdparty/EigenQP.cpp.s"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ty/Physim/bounding-mesh-0.2/thirdparty/EigenQP.cpp -o CMakeFiles/boundingmesh.dir/__/__/thirdparty/EigenQP.cpp.s

src/boundingmesh/CMakeFiles/boundingmesh.dir/__/__/thirdparty/EigenQP.cpp.o.requires:

.PHONY : src/boundingmesh/CMakeFiles/boundingmesh.dir/__/__/thirdparty/EigenQP.cpp.o.requires

src/boundingmesh/CMakeFiles/boundingmesh.dir/__/__/thirdparty/EigenQP.cpp.o.provides: src/boundingmesh/CMakeFiles/boundingmesh.dir/__/__/thirdparty/EigenQP.cpp.o.requires
	$(MAKE) -f src/boundingmesh/CMakeFiles/boundingmesh.dir/build.make src/boundingmesh/CMakeFiles/boundingmesh.dir/__/__/thirdparty/EigenQP.cpp.o.provides.build
.PHONY : src/boundingmesh/CMakeFiles/boundingmesh.dir/__/__/thirdparty/EigenQP.cpp.o.provides

src/boundingmesh/CMakeFiles/boundingmesh.dir/__/__/thirdparty/EigenQP.cpp.o.provides.build: src/boundingmesh/CMakeFiles/boundingmesh.dir/__/__/thirdparty/EigenQP.cpp.o


# Object files for target boundingmesh
boundingmesh_OBJECTS = \
"CMakeFiles/boundingmesh.dir/ContractionUtils.cpp.o" \
"CMakeFiles/boundingmesh.dir/Decimator.cpp.o" \
"CMakeFiles/boundingmesh.dir/FileUtils.cpp.o" \
"CMakeFiles/boundingmesh.dir/Mesh.cpp.o" \
"CMakeFiles/boundingmesh.dir/MetricGenerator.cpp.o" \
"CMakeFiles/boundingmesh.dir/Primitives.cpp.o" \
"CMakeFiles/boundingmesh.dir/SegmenterDownsampling.cpp.o" \
"CMakeFiles/boundingmesh.dir/SegmenterSimple.cpp.o" \
"CMakeFiles/boundingmesh.dir/Split.cpp.o" \
"CMakeFiles/boundingmesh.dir/VoxelSet.cpp.o" \
"CMakeFiles/boundingmesh.dir/VoxelSubset.cpp.o" \
"CMakeFiles/boundingmesh.dir/__/__/thirdparty/EigenQP.cpp.o"

# External object files for target boundingmesh
boundingmesh_EXTERNAL_OBJECTS =

libboundingmesh.a: src/boundingmesh/CMakeFiles/boundingmesh.dir/ContractionUtils.cpp.o
libboundingmesh.a: src/boundingmesh/CMakeFiles/boundingmesh.dir/Decimator.cpp.o
libboundingmesh.a: src/boundingmesh/CMakeFiles/boundingmesh.dir/FileUtils.cpp.o
libboundingmesh.a: src/boundingmesh/CMakeFiles/boundingmesh.dir/Mesh.cpp.o
libboundingmesh.a: src/boundingmesh/CMakeFiles/boundingmesh.dir/MetricGenerator.cpp.o
libboundingmesh.a: src/boundingmesh/CMakeFiles/boundingmesh.dir/Primitives.cpp.o
libboundingmesh.a: src/boundingmesh/CMakeFiles/boundingmesh.dir/SegmenterDownsampling.cpp.o
libboundingmesh.a: src/boundingmesh/CMakeFiles/boundingmesh.dir/SegmenterSimple.cpp.o
libboundingmesh.a: src/boundingmesh/CMakeFiles/boundingmesh.dir/Split.cpp.o
libboundingmesh.a: src/boundingmesh/CMakeFiles/boundingmesh.dir/VoxelSet.cpp.o
libboundingmesh.a: src/boundingmesh/CMakeFiles/boundingmesh.dir/VoxelSubset.cpp.o
libboundingmesh.a: src/boundingmesh/CMakeFiles/boundingmesh.dir/__/__/thirdparty/EigenQP.cpp.o
libboundingmesh.a: src/boundingmesh/CMakeFiles/boundingmesh.dir/build.make
libboundingmesh.a: src/boundingmesh/CMakeFiles/boundingmesh.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ty/Physim/bounding-mesh-0.2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Linking CXX static library ../../libboundingmesh.a"
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && $(CMAKE_COMMAND) -P CMakeFiles/boundingmesh.dir/cmake_clean_target.cmake
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/boundingmesh.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/boundingmesh/CMakeFiles/boundingmesh.dir/build: libboundingmesh.a

.PHONY : src/boundingmesh/CMakeFiles/boundingmesh.dir/build

src/boundingmesh/CMakeFiles/boundingmesh.dir/requires: src/boundingmesh/CMakeFiles/boundingmesh.dir/ContractionUtils.cpp.o.requires
src/boundingmesh/CMakeFiles/boundingmesh.dir/requires: src/boundingmesh/CMakeFiles/boundingmesh.dir/Decimator.cpp.o.requires
src/boundingmesh/CMakeFiles/boundingmesh.dir/requires: src/boundingmesh/CMakeFiles/boundingmesh.dir/FileUtils.cpp.o.requires
src/boundingmesh/CMakeFiles/boundingmesh.dir/requires: src/boundingmesh/CMakeFiles/boundingmesh.dir/Mesh.cpp.o.requires
src/boundingmesh/CMakeFiles/boundingmesh.dir/requires: src/boundingmesh/CMakeFiles/boundingmesh.dir/MetricGenerator.cpp.o.requires
src/boundingmesh/CMakeFiles/boundingmesh.dir/requires: src/boundingmesh/CMakeFiles/boundingmesh.dir/Primitives.cpp.o.requires
src/boundingmesh/CMakeFiles/boundingmesh.dir/requires: src/boundingmesh/CMakeFiles/boundingmesh.dir/SegmenterDownsampling.cpp.o.requires
src/boundingmesh/CMakeFiles/boundingmesh.dir/requires: src/boundingmesh/CMakeFiles/boundingmesh.dir/SegmenterSimple.cpp.o.requires
src/boundingmesh/CMakeFiles/boundingmesh.dir/requires: src/boundingmesh/CMakeFiles/boundingmesh.dir/Split.cpp.o.requires
src/boundingmesh/CMakeFiles/boundingmesh.dir/requires: src/boundingmesh/CMakeFiles/boundingmesh.dir/VoxelSet.cpp.o.requires
src/boundingmesh/CMakeFiles/boundingmesh.dir/requires: src/boundingmesh/CMakeFiles/boundingmesh.dir/VoxelSubset.cpp.o.requires
src/boundingmesh/CMakeFiles/boundingmesh.dir/requires: src/boundingmesh/CMakeFiles/boundingmesh.dir/__/__/thirdparty/EigenQP.cpp.o.requires

.PHONY : src/boundingmesh/CMakeFiles/boundingmesh.dir/requires

src/boundingmesh/CMakeFiles/boundingmesh.dir/clean:
	cd /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh && $(CMAKE_COMMAND) -P CMakeFiles/boundingmesh.dir/cmake_clean.cmake
.PHONY : src/boundingmesh/CMakeFiles/boundingmesh.dir/clean

src/boundingmesh/CMakeFiles/boundingmesh.dir/depend:
	cd /home/ty/Physim/bounding-mesh-0.2/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ty/Physim/bounding-mesh-0.2 /home/ty/Physim/bounding-mesh-0.2/src/boundingmesh /home/ty/Physim/bounding-mesh-0.2/build /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh /home/ty/Physim/bounding-mesh-0.2/build/src/boundingmesh/CMakeFiles/boundingmesh.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/boundingmesh/CMakeFiles/boundingmesh.dir/depend
