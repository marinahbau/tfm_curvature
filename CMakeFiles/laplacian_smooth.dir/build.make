# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /mnt/c/Users/marina/Documents/TFM_desarrollo/tfm_curvature

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/c/Users/marina/Documents/TFM_desarrollo/tfm_curvature

# Include any dependencies generated for this target.
include CMakeFiles/laplacian_smooth.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/laplacian_smooth.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/laplacian_smooth.dir/flags.make

CMakeFiles/laplacian_smooth.dir/src/laplacian.cpp.o: CMakeFiles/laplacian_smooth.dir/flags.make
CMakeFiles/laplacian_smooth.dir/src/laplacian.cpp.o: src/laplacian.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Users/marina/Documents/TFM_desarrollo/tfm_curvature/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/laplacian_smooth.dir/src/laplacian.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/laplacian_smooth.dir/src/laplacian.cpp.o -c /mnt/c/Users/marina/Documents/TFM_desarrollo/tfm_curvature/src/laplacian.cpp

CMakeFiles/laplacian_smooth.dir/src/laplacian.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/laplacian_smooth.dir/src/laplacian.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/marina/Documents/TFM_desarrollo/tfm_curvature/src/laplacian.cpp > CMakeFiles/laplacian_smooth.dir/src/laplacian.cpp.i

CMakeFiles/laplacian_smooth.dir/src/laplacian.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/laplacian_smooth.dir/src/laplacian.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/marina/Documents/TFM_desarrollo/tfm_curvature/src/laplacian.cpp -o CMakeFiles/laplacian_smooth.dir/src/laplacian.cpp.s

# Object files for target laplacian_smooth
laplacian_smooth_OBJECTS = \
"CMakeFiles/laplacian_smooth.dir/src/laplacian.cpp.o"

# External object files for target laplacian_smooth
laplacian_smooth_EXTERNAL_OBJECTS =

laplacian_smooth: CMakeFiles/laplacian_smooth.dir/src/laplacian.cpp.o
laplacian_smooth: CMakeFiles/laplacian_smooth.dir/build.make
laplacian_smooth: /usr/lib/x86_64-linux-gnu/libGLX.so
laplacian_smooth: /usr/lib/x86_64-linux-gnu/libOpenGL.so
laplacian_smooth: glad/libglad.a
laplacian_smooth: glfw/src/libglfw3.a
laplacian_smooth: /usr/lib/x86_64-linux-gnu/librt.so
laplacian_smooth: /usr/lib/x86_64-linux-gnu/libm.so
laplacian_smooth: CMakeFiles/laplacian_smooth.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/c/Users/marina/Documents/TFM_desarrollo/tfm_curvature/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable laplacian_smooth"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/laplacian_smooth.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/laplacian_smooth.dir/build: laplacian_smooth

.PHONY : CMakeFiles/laplacian_smooth.dir/build

CMakeFiles/laplacian_smooth.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/laplacian_smooth.dir/cmake_clean.cmake
.PHONY : CMakeFiles/laplacian_smooth.dir/clean

CMakeFiles/laplacian_smooth.dir/depend:
	cd /mnt/c/Users/marina/Documents/TFM_desarrollo/tfm_curvature && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/c/Users/marina/Documents/TFM_desarrollo/tfm_curvature /mnt/c/Users/marina/Documents/TFM_desarrollo/tfm_curvature /mnt/c/Users/marina/Documents/TFM_desarrollo/tfm_curvature /mnt/c/Users/marina/Documents/TFM_desarrollo/tfm_curvature /mnt/c/Users/marina/Documents/TFM_desarrollo/tfm_curvature/CMakeFiles/laplacian_smooth.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/laplacian_smooth.dir/depend

