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
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.30.5/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.30.5/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/alexandrshain/cpp/quasi-one-dimensional/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/alexandrshain/cpp/quasi-one-dimensional/src/build

# Include any dependencies generated for this target.
include CMakeFiles/MyExecutable.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/MyExecutable.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/MyExecutable.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/MyExecutable.dir/flags.make

CMakeFiles/MyExecutable.dir/main.cpp.o: CMakeFiles/MyExecutable.dir/flags.make
CMakeFiles/MyExecutable.dir/main.cpp.o: /Users/alexandrshain/cpp/quasi-one-dimensional/src/main.cpp
CMakeFiles/MyExecutable.dir/main.cpp.o: CMakeFiles/MyExecutable.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/alexandrshain/cpp/quasi-one-dimensional/src/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/MyExecutable.dir/main.cpp.o"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/MyExecutable.dir/main.cpp.o -MF CMakeFiles/MyExecutable.dir/main.cpp.o.d -o CMakeFiles/MyExecutable.dir/main.cpp.o -c /Users/alexandrshain/cpp/quasi-one-dimensional/src/main.cpp

CMakeFiles/MyExecutable.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/MyExecutable.dir/main.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/alexandrshain/cpp/quasi-one-dimensional/src/main.cpp > CMakeFiles/MyExecutable.dir/main.cpp.i

CMakeFiles/MyExecutable.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/MyExecutable.dir/main.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/alexandrshain/cpp/quasi-one-dimensional/src/main.cpp -o CMakeFiles/MyExecutable.dir/main.cpp.s

CMakeFiles/MyExecutable.dir/ChannelInputData.cpp.o: CMakeFiles/MyExecutable.dir/flags.make
CMakeFiles/MyExecutable.dir/ChannelInputData.cpp.o: /Users/alexandrshain/cpp/quasi-one-dimensional/src/ChannelInputData.cpp
CMakeFiles/MyExecutable.dir/ChannelInputData.cpp.o: CMakeFiles/MyExecutable.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/alexandrshain/cpp/quasi-one-dimensional/src/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/MyExecutable.dir/ChannelInputData.cpp.o"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/MyExecutable.dir/ChannelInputData.cpp.o -MF CMakeFiles/MyExecutable.dir/ChannelInputData.cpp.o.d -o CMakeFiles/MyExecutable.dir/ChannelInputData.cpp.o -c /Users/alexandrshain/cpp/quasi-one-dimensional/src/ChannelInputData.cpp

CMakeFiles/MyExecutable.dir/ChannelInputData.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/MyExecutable.dir/ChannelInputData.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/alexandrshain/cpp/quasi-one-dimensional/src/ChannelInputData.cpp > CMakeFiles/MyExecutable.dir/ChannelInputData.cpp.i

CMakeFiles/MyExecutable.dir/ChannelInputData.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/MyExecutable.dir/ChannelInputData.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/alexandrshain/cpp/quasi-one-dimensional/src/ChannelInputData.cpp -o CMakeFiles/MyExecutable.dir/ChannelInputData.cpp.s

CMakeFiles/MyExecutable.dir/ChannelOutputData.cpp.o: CMakeFiles/MyExecutable.dir/flags.make
CMakeFiles/MyExecutable.dir/ChannelOutputData.cpp.o: /Users/alexandrshain/cpp/quasi-one-dimensional/src/ChannelOutputData.cpp
CMakeFiles/MyExecutable.dir/ChannelOutputData.cpp.o: CMakeFiles/MyExecutable.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/alexandrshain/cpp/quasi-one-dimensional/src/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/MyExecutable.dir/ChannelOutputData.cpp.o"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/MyExecutable.dir/ChannelOutputData.cpp.o -MF CMakeFiles/MyExecutable.dir/ChannelOutputData.cpp.o.d -o CMakeFiles/MyExecutable.dir/ChannelOutputData.cpp.o -c /Users/alexandrshain/cpp/quasi-one-dimensional/src/ChannelOutputData.cpp

CMakeFiles/MyExecutable.dir/ChannelOutputData.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/MyExecutable.dir/ChannelOutputData.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/alexandrshain/cpp/quasi-one-dimensional/src/ChannelOutputData.cpp > CMakeFiles/MyExecutable.dir/ChannelOutputData.cpp.i

CMakeFiles/MyExecutable.dir/ChannelOutputData.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/MyExecutable.dir/ChannelOutputData.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/alexandrshain/cpp/quasi-one-dimensional/src/ChannelOutputData.cpp -o CMakeFiles/MyExecutable.dir/ChannelOutputData.cpp.s

CMakeFiles/MyExecutable.dir/ChannelParametersReader.cpp.o: CMakeFiles/MyExecutable.dir/flags.make
CMakeFiles/MyExecutable.dir/ChannelParametersReader.cpp.o: /Users/alexandrshain/cpp/quasi-one-dimensional/src/ChannelParametersReader.cpp
CMakeFiles/MyExecutable.dir/ChannelParametersReader.cpp.o: CMakeFiles/MyExecutable.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/alexandrshain/cpp/quasi-one-dimensional/src/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/MyExecutable.dir/ChannelParametersReader.cpp.o"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/MyExecutable.dir/ChannelParametersReader.cpp.o -MF CMakeFiles/MyExecutable.dir/ChannelParametersReader.cpp.o.d -o CMakeFiles/MyExecutable.dir/ChannelParametersReader.cpp.o -c /Users/alexandrshain/cpp/quasi-one-dimensional/src/ChannelParametersReader.cpp

CMakeFiles/MyExecutable.dir/ChannelParametersReader.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/MyExecutable.dir/ChannelParametersReader.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/alexandrshain/cpp/quasi-one-dimensional/src/ChannelParametersReader.cpp > CMakeFiles/MyExecutable.dir/ChannelParametersReader.cpp.i

CMakeFiles/MyExecutable.dir/ChannelParametersReader.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/MyExecutable.dir/ChannelParametersReader.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/alexandrshain/cpp/quasi-one-dimensional/src/ChannelParametersReader.cpp -o CMakeFiles/MyExecutable.dir/ChannelParametersReader.cpp.s

CMakeFiles/MyExecutable.dir/ChannelSolver.cpp.o: CMakeFiles/MyExecutable.dir/flags.make
CMakeFiles/MyExecutable.dir/ChannelSolver.cpp.o: /Users/alexandrshain/cpp/quasi-one-dimensional/src/ChannelSolver.cpp
CMakeFiles/MyExecutable.dir/ChannelSolver.cpp.o: CMakeFiles/MyExecutable.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/alexandrshain/cpp/quasi-one-dimensional/src/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/MyExecutable.dir/ChannelSolver.cpp.o"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/MyExecutable.dir/ChannelSolver.cpp.o -MF CMakeFiles/MyExecutable.dir/ChannelSolver.cpp.o.d -o CMakeFiles/MyExecutable.dir/ChannelSolver.cpp.o -c /Users/alexandrshain/cpp/quasi-one-dimensional/src/ChannelSolver.cpp

CMakeFiles/MyExecutable.dir/ChannelSolver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/MyExecutable.dir/ChannelSolver.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/alexandrshain/cpp/quasi-one-dimensional/src/ChannelSolver.cpp > CMakeFiles/MyExecutable.dir/ChannelSolver.cpp.i

CMakeFiles/MyExecutable.dir/ChannelSolver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/MyExecutable.dir/ChannelSolver.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/alexandrshain/cpp/quasi-one-dimensional/src/ChannelSolver.cpp -o CMakeFiles/MyExecutable.dir/ChannelSolver.cpp.s

# Object files for target MyExecutable
MyExecutable_OBJECTS = \
"CMakeFiles/MyExecutable.dir/main.cpp.o" \
"CMakeFiles/MyExecutable.dir/ChannelInputData.cpp.o" \
"CMakeFiles/MyExecutable.dir/ChannelOutputData.cpp.o" \
"CMakeFiles/MyExecutable.dir/ChannelParametersReader.cpp.o" \
"CMakeFiles/MyExecutable.dir/ChannelSolver.cpp.o"

# External object files for target MyExecutable
MyExecutable_EXTERNAL_OBJECTS =

MyExecutable: CMakeFiles/MyExecutable.dir/main.cpp.o
MyExecutable: CMakeFiles/MyExecutable.dir/ChannelInputData.cpp.o
MyExecutable: CMakeFiles/MyExecutable.dir/ChannelOutputData.cpp.o
MyExecutable: CMakeFiles/MyExecutable.dir/ChannelParametersReader.cpp.o
MyExecutable: CMakeFiles/MyExecutable.dir/ChannelSolver.cpp.o
MyExecutable: CMakeFiles/MyExecutable.dir/build.make
MyExecutable: CMakeFiles/MyExecutable.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/alexandrshain/cpp/quasi-one-dimensional/src/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX executable MyExecutable"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/MyExecutable.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/MyExecutable.dir/build: MyExecutable
.PHONY : CMakeFiles/MyExecutable.dir/build

CMakeFiles/MyExecutable.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/MyExecutable.dir/cmake_clean.cmake
.PHONY : CMakeFiles/MyExecutable.dir/clean

CMakeFiles/MyExecutable.dir/depend:
	cd /Users/alexandrshain/cpp/quasi-one-dimensional/src/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/alexandrshain/cpp/quasi-one-dimensional/src /Users/alexandrshain/cpp/quasi-one-dimensional/src /Users/alexandrshain/cpp/quasi-one-dimensional/src/build /Users/alexandrshain/cpp/quasi-one-dimensional/src/build /Users/alexandrshain/cpp/quasi-one-dimensional/src/build/CMakeFiles/MyExecutable.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/MyExecutable.dir/depend

