# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.30

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

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/telatina/git/covtobed/test

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/telatina/git/covtobed/test/build

# Include any dependencies generated for this target.
include unit/CMakeFiles/unit_tests.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include unit/CMakeFiles/unit_tests.dir/compiler_depend.make

# Include the progress variables for this target.
include unit/CMakeFiles/unit_tests.dir/progress.make

# Include the compile flags for this target's objects.
include unit/CMakeFiles/unit_tests.dir/flags.make

unit/CMakeFiles/unit_tests.dir/test_main.cpp.o: unit/CMakeFiles/unit_tests.dir/flags.make
unit/CMakeFiles/unit_tests.dir/test_main.cpp.o: /home/telatina/git/covtobed/test/unit/test_main.cpp
unit/CMakeFiles/unit_tests.dir/test_main.cpp.o: unit/CMakeFiles/unit_tests.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/telatina/git/covtobed/test/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object unit/CMakeFiles/unit_tests.dir/test_main.cpp.o"
	cd /home/telatina/git/covtobed/test/build/unit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT unit/CMakeFiles/unit_tests.dir/test_main.cpp.o -MF CMakeFiles/unit_tests.dir/test_main.cpp.o.d -o CMakeFiles/unit_tests.dir/test_main.cpp.o -c /home/telatina/git/covtobed/test/unit/test_main.cpp

unit/CMakeFiles/unit_tests.dir/test_main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/unit_tests.dir/test_main.cpp.i"
	cd /home/telatina/git/covtobed/test/build/unit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/telatina/git/covtobed/test/unit/test_main.cpp > CMakeFiles/unit_tests.dir/test_main.cpp.i

unit/CMakeFiles/unit_tests.dir/test_main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/unit_tests.dir/test_main.cpp.s"
	cd /home/telatina/git/covtobed/test/build/unit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/telatina/git/covtobed/test/unit/test_main.cpp -o CMakeFiles/unit_tests.dir/test_main.cpp.s

unit/CMakeFiles/unit_tests.dir/test_coverage.cpp.o: unit/CMakeFiles/unit_tests.dir/flags.make
unit/CMakeFiles/unit_tests.dir/test_coverage.cpp.o: /home/telatina/git/covtobed/test/unit/test_coverage.cpp
unit/CMakeFiles/unit_tests.dir/test_coverage.cpp.o: unit/CMakeFiles/unit_tests.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/telatina/git/covtobed/test/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object unit/CMakeFiles/unit_tests.dir/test_coverage.cpp.o"
	cd /home/telatina/git/covtobed/test/build/unit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT unit/CMakeFiles/unit_tests.dir/test_coverage.cpp.o -MF CMakeFiles/unit_tests.dir/test_coverage.cpp.o.d -o CMakeFiles/unit_tests.dir/test_coverage.cpp.o -c /home/telatina/git/covtobed/test/unit/test_coverage.cpp

unit/CMakeFiles/unit_tests.dir/test_coverage.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/unit_tests.dir/test_coverage.cpp.i"
	cd /home/telatina/git/covtobed/test/build/unit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/telatina/git/covtobed/test/unit/test_coverage.cpp > CMakeFiles/unit_tests.dir/test_coverage.cpp.i

unit/CMakeFiles/unit_tests.dir/test_coverage.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/unit_tests.dir/test_coverage.cpp.s"
	cd /home/telatina/git/covtobed/test/build/unit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/telatina/git/covtobed/test/unit/test_coverage.cpp -o CMakeFiles/unit_tests.dir/test_coverage.cpp.s

unit/CMakeFiles/unit_tests.dir/test_intervals.cpp.o: unit/CMakeFiles/unit_tests.dir/flags.make
unit/CMakeFiles/unit_tests.dir/test_intervals.cpp.o: /home/telatina/git/covtobed/test/unit/test_intervals.cpp
unit/CMakeFiles/unit_tests.dir/test_intervals.cpp.o: unit/CMakeFiles/unit_tests.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/telatina/git/covtobed/test/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object unit/CMakeFiles/unit_tests.dir/test_intervals.cpp.o"
	cd /home/telatina/git/covtobed/test/build/unit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT unit/CMakeFiles/unit_tests.dir/test_intervals.cpp.o -MF CMakeFiles/unit_tests.dir/test_intervals.cpp.o.d -o CMakeFiles/unit_tests.dir/test_intervals.cpp.o -c /home/telatina/git/covtobed/test/unit/test_intervals.cpp

unit/CMakeFiles/unit_tests.dir/test_intervals.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/unit_tests.dir/test_intervals.cpp.i"
	cd /home/telatina/git/covtobed/test/build/unit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/telatina/git/covtobed/test/unit/test_intervals.cpp > CMakeFiles/unit_tests.dir/test_intervals.cpp.i

unit/CMakeFiles/unit_tests.dir/test_intervals.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/unit_tests.dir/test_intervals.cpp.s"
	cd /home/telatina/git/covtobed/test/build/unit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/telatina/git/covtobed/test/unit/test_intervals.cpp -o CMakeFiles/unit_tests.dir/test_intervals.cpp.s

unit/CMakeFiles/unit_tests.dir/test_output.cpp.o: unit/CMakeFiles/unit_tests.dir/flags.make
unit/CMakeFiles/unit_tests.dir/test_output.cpp.o: /home/telatina/git/covtobed/test/unit/test_output.cpp
unit/CMakeFiles/unit_tests.dir/test_output.cpp.o: unit/CMakeFiles/unit_tests.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/telatina/git/covtobed/test/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object unit/CMakeFiles/unit_tests.dir/test_output.cpp.o"
	cd /home/telatina/git/covtobed/test/build/unit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT unit/CMakeFiles/unit_tests.dir/test_output.cpp.o -MF CMakeFiles/unit_tests.dir/test_output.cpp.o.d -o CMakeFiles/unit_tests.dir/test_output.cpp.o -c /home/telatina/git/covtobed/test/unit/test_output.cpp

unit/CMakeFiles/unit_tests.dir/test_output.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/unit_tests.dir/test_output.cpp.i"
	cd /home/telatina/git/covtobed/test/build/unit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/telatina/git/covtobed/test/unit/test_output.cpp > CMakeFiles/unit_tests.dir/test_output.cpp.i

unit/CMakeFiles/unit_tests.dir/test_output.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/unit_tests.dir/test_output.cpp.s"
	cd /home/telatina/git/covtobed/test/build/unit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/telatina/git/covtobed/test/unit/test_output.cpp -o CMakeFiles/unit_tests.dir/test_output.cpp.s

# Object files for target unit_tests
unit_tests_OBJECTS = \
"CMakeFiles/unit_tests.dir/test_main.cpp.o" \
"CMakeFiles/unit_tests.dir/test_coverage.cpp.o" \
"CMakeFiles/unit_tests.dir/test_intervals.cpp.o" \
"CMakeFiles/unit_tests.dir/test_output.cpp.o"

# External object files for target unit_tests
unit_tests_EXTERNAL_OBJECTS =

unit/unit_tests: unit/CMakeFiles/unit_tests.dir/test_main.cpp.o
unit/unit_tests: unit/CMakeFiles/unit_tests.dir/test_coverage.cpp.o
unit/unit_tests: unit/CMakeFiles/unit_tests.dir/test_intervals.cpp.o
unit/unit_tests: unit/CMakeFiles/unit_tests.dir/test_output.cpp.o
unit/unit_tests: unit/CMakeFiles/unit_tests.dir/build.make
unit/unit_tests: /usr/lib/libCatch2Main.a
unit/unit_tests: /usr/lib/libCatch2.a
unit/unit_tests: unit/CMakeFiles/unit_tests.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/telatina/git/covtobed/test/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX executable unit_tests"
	cd /home/telatina/git/covtobed/test/build/unit && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/unit_tests.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
unit/CMakeFiles/unit_tests.dir/build: unit/unit_tests
.PHONY : unit/CMakeFiles/unit_tests.dir/build

unit/CMakeFiles/unit_tests.dir/clean:
	cd /home/telatina/git/covtobed/test/build/unit && $(CMAKE_COMMAND) -P CMakeFiles/unit_tests.dir/cmake_clean.cmake
.PHONY : unit/CMakeFiles/unit_tests.dir/clean

unit/CMakeFiles/unit_tests.dir/depend:
	cd /home/telatina/git/covtobed/test/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/telatina/git/covtobed/test /home/telatina/git/covtobed/test/unit /home/telatina/git/covtobed/test/build /home/telatina/git/covtobed/test/build/unit /home/telatina/git/covtobed/test/build/unit/CMakeFiles/unit_tests.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : unit/CMakeFiles/unit_tests.dir/depend

