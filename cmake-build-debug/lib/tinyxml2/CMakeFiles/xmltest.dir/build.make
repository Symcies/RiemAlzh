# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.7

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
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/clementine.fourrier/RiemAlzh

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/clementine.fourrier/RiemAlzh/cmake-build-debug

# Include any dependencies generated for this target.
include lib/tinyxml2/CMakeFiles/xmltest.dir/depend.make

# Include the progress variables for this target.
include lib/tinyxml2/CMakeFiles/xmltest.dir/progress.make

# Include the compile flags for this target's objects.
include lib/tinyxml2/CMakeFiles/xmltest.dir/flags.make

lib/tinyxml2/CMakeFiles/xmltest.dir/xmltest.cpp.o: lib/tinyxml2/CMakeFiles/xmltest.dir/flags.make
lib/tinyxml2/CMakeFiles/xmltest.dir/xmltest.cpp.o: ../lib/tinyxml2/xmltest.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/clementine.fourrier/RiemAlzh/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object lib/tinyxml2/CMakeFiles/xmltest.dir/xmltest.cpp.o"
	cd /Users/clementine.fourrier/RiemAlzh/cmake-build-debug/lib/tinyxml2 && /usr/bin/g++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/xmltest.dir/xmltest.cpp.o -c /Users/clementine.fourrier/RiemAlzh/lib/tinyxml2/xmltest.cpp

lib/tinyxml2/CMakeFiles/xmltest.dir/xmltest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/xmltest.dir/xmltest.cpp.i"
	cd /Users/clementine.fourrier/RiemAlzh/cmake-build-debug/lib/tinyxml2 && /usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/clementine.fourrier/RiemAlzh/lib/tinyxml2/xmltest.cpp > CMakeFiles/xmltest.dir/xmltest.cpp.i

lib/tinyxml2/CMakeFiles/xmltest.dir/xmltest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/xmltest.dir/xmltest.cpp.s"
	cd /Users/clementine.fourrier/RiemAlzh/cmake-build-debug/lib/tinyxml2 && /usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/clementine.fourrier/RiemAlzh/lib/tinyxml2/xmltest.cpp -o CMakeFiles/xmltest.dir/xmltest.cpp.s

lib/tinyxml2/CMakeFiles/xmltest.dir/xmltest.cpp.o.requires:

.PHONY : lib/tinyxml2/CMakeFiles/xmltest.dir/xmltest.cpp.o.requires

lib/tinyxml2/CMakeFiles/xmltest.dir/xmltest.cpp.o.provides: lib/tinyxml2/CMakeFiles/xmltest.dir/xmltest.cpp.o.requires
	$(MAKE) -f lib/tinyxml2/CMakeFiles/xmltest.dir/build.make lib/tinyxml2/CMakeFiles/xmltest.dir/xmltest.cpp.o.provides.build
.PHONY : lib/tinyxml2/CMakeFiles/xmltest.dir/xmltest.cpp.o.provides

lib/tinyxml2/CMakeFiles/xmltest.dir/xmltest.cpp.o.provides.build: lib/tinyxml2/CMakeFiles/xmltest.dir/xmltest.cpp.o


# Object files for target xmltest
xmltest_OBJECTS = \
"CMakeFiles/xmltest.dir/xmltest.cpp.o"

# External object files for target xmltest
xmltest_EXTERNAL_OBJECTS =

lib/tinyxml2/xmltest: lib/tinyxml2/CMakeFiles/xmltest.dir/xmltest.cpp.o
lib/tinyxml2/xmltest: lib/tinyxml2/CMakeFiles/xmltest.dir/build.make
lib/tinyxml2/xmltest: lib/tinyxml2/libtinyxml2.4.0.1.dylib
lib/tinyxml2/xmltest: lib/tinyxml2/CMakeFiles/xmltest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/clementine.fourrier/RiemAlzh/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable xmltest"
	cd /Users/clementine.fourrier/RiemAlzh/cmake-build-debug/lib/tinyxml2 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/xmltest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lib/tinyxml2/CMakeFiles/xmltest.dir/build: lib/tinyxml2/xmltest

.PHONY : lib/tinyxml2/CMakeFiles/xmltest.dir/build

lib/tinyxml2/CMakeFiles/xmltest.dir/requires: lib/tinyxml2/CMakeFiles/xmltest.dir/xmltest.cpp.o.requires

.PHONY : lib/tinyxml2/CMakeFiles/xmltest.dir/requires

lib/tinyxml2/CMakeFiles/xmltest.dir/clean:
	cd /Users/clementine.fourrier/RiemAlzh/cmake-build-debug/lib/tinyxml2 && $(CMAKE_COMMAND) -P CMakeFiles/xmltest.dir/cmake_clean.cmake
.PHONY : lib/tinyxml2/CMakeFiles/xmltest.dir/clean

lib/tinyxml2/CMakeFiles/xmltest.dir/depend:
	cd /Users/clementine.fourrier/RiemAlzh/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/clementine.fourrier/RiemAlzh /Users/clementine.fourrier/RiemAlzh/lib/tinyxml2 /Users/clementine.fourrier/RiemAlzh/cmake-build-debug /Users/clementine.fourrier/RiemAlzh/cmake-build-debug/lib/tinyxml2 /Users/clementine.fourrier/RiemAlzh/cmake-build-debug/lib/tinyxml2/CMakeFiles/xmltest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/tinyxml2/CMakeFiles/xmltest.dir/depend

