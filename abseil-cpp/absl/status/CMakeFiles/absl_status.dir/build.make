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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.16.4/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.16.4/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/marko/repos/sima

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/marko/repos/sima

# Include any dependencies generated for this target.
include abseil-cpp/absl/status/CMakeFiles/absl_status.dir/depend.make

# Include the progress variables for this target.
include abseil-cpp/absl/status/CMakeFiles/absl_status.dir/progress.make

# Include the compile flags for this target's objects.
include abseil-cpp/absl/status/CMakeFiles/absl_status.dir/flags.make

abseil-cpp/absl/status/CMakeFiles/absl_status.dir/status.cc.o: abseil-cpp/absl/status/CMakeFiles/absl_status.dir/flags.make
abseil-cpp/absl/status/CMakeFiles/absl_status.dir/status.cc.o: abseil-cpp/absl/status/status.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/marko/repos/sima/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object abseil-cpp/absl/status/CMakeFiles/absl_status.dir/status.cc.o"
	cd /Users/marko/repos/sima/abseil-cpp/absl/status && /usr/local/opt/ccache/libexec/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/absl_status.dir/status.cc.o -c /Users/marko/repos/sima/abseil-cpp/absl/status/status.cc

abseil-cpp/absl/status/CMakeFiles/absl_status.dir/status.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/absl_status.dir/status.cc.i"
	cd /Users/marko/repos/sima/abseil-cpp/absl/status && /usr/local/opt/ccache/libexec/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/marko/repos/sima/abseil-cpp/absl/status/status.cc > CMakeFiles/absl_status.dir/status.cc.i

abseil-cpp/absl/status/CMakeFiles/absl_status.dir/status.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/absl_status.dir/status.cc.s"
	cd /Users/marko/repos/sima/abseil-cpp/absl/status && /usr/local/opt/ccache/libexec/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/marko/repos/sima/abseil-cpp/absl/status/status.cc -o CMakeFiles/absl_status.dir/status.cc.s

abseil-cpp/absl/status/CMakeFiles/absl_status.dir/status_payload_printer.cc.o: abseil-cpp/absl/status/CMakeFiles/absl_status.dir/flags.make
abseil-cpp/absl/status/CMakeFiles/absl_status.dir/status_payload_printer.cc.o: abseil-cpp/absl/status/status_payload_printer.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/marko/repos/sima/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object abseil-cpp/absl/status/CMakeFiles/absl_status.dir/status_payload_printer.cc.o"
	cd /Users/marko/repos/sima/abseil-cpp/absl/status && /usr/local/opt/ccache/libexec/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/absl_status.dir/status_payload_printer.cc.o -c /Users/marko/repos/sima/abseil-cpp/absl/status/status_payload_printer.cc

abseil-cpp/absl/status/CMakeFiles/absl_status.dir/status_payload_printer.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/absl_status.dir/status_payload_printer.cc.i"
	cd /Users/marko/repos/sima/abseil-cpp/absl/status && /usr/local/opt/ccache/libexec/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/marko/repos/sima/abseil-cpp/absl/status/status_payload_printer.cc > CMakeFiles/absl_status.dir/status_payload_printer.cc.i

abseil-cpp/absl/status/CMakeFiles/absl_status.dir/status_payload_printer.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/absl_status.dir/status_payload_printer.cc.s"
	cd /Users/marko/repos/sima/abseil-cpp/absl/status && /usr/local/opt/ccache/libexec/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/marko/repos/sima/abseil-cpp/absl/status/status_payload_printer.cc -o CMakeFiles/absl_status.dir/status_payload_printer.cc.s

# Object files for target absl_status
absl_status_OBJECTS = \
"CMakeFiles/absl_status.dir/status.cc.o" \
"CMakeFiles/absl_status.dir/status_payload_printer.cc.o"

# External object files for target absl_status
absl_status_EXTERNAL_OBJECTS =

abseil-cpp/absl/status/libabsl_status.a: abseil-cpp/absl/status/CMakeFiles/absl_status.dir/status.cc.o
abseil-cpp/absl/status/libabsl_status.a: abseil-cpp/absl/status/CMakeFiles/absl_status.dir/status_payload_printer.cc.o
abseil-cpp/absl/status/libabsl_status.a: abseil-cpp/absl/status/CMakeFiles/absl_status.dir/build.make
abseil-cpp/absl/status/libabsl_status.a: abseil-cpp/absl/status/CMakeFiles/absl_status.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/marko/repos/sima/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX static library libabsl_status.a"
	cd /Users/marko/repos/sima/abseil-cpp/absl/status && $(CMAKE_COMMAND) -P CMakeFiles/absl_status.dir/cmake_clean_target.cmake
	cd /Users/marko/repos/sima/abseil-cpp/absl/status && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/absl_status.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
abseil-cpp/absl/status/CMakeFiles/absl_status.dir/build: abseil-cpp/absl/status/libabsl_status.a

.PHONY : abseil-cpp/absl/status/CMakeFiles/absl_status.dir/build

abseil-cpp/absl/status/CMakeFiles/absl_status.dir/clean:
	cd /Users/marko/repos/sima/abseil-cpp/absl/status && $(CMAKE_COMMAND) -P CMakeFiles/absl_status.dir/cmake_clean.cmake
.PHONY : abseil-cpp/absl/status/CMakeFiles/absl_status.dir/clean

abseil-cpp/absl/status/CMakeFiles/absl_status.dir/depend:
	cd /Users/marko/repos/sima && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/marko/repos/sima /Users/marko/repos/sima/abseil-cpp/absl/status /Users/marko/repos/sima /Users/marko/repos/sima/abseil-cpp/absl/status /Users/marko/repos/sima/abseil-cpp/absl/status/CMakeFiles/absl_status.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : abseil-cpp/absl/status/CMakeFiles/absl_status.dir/depend

