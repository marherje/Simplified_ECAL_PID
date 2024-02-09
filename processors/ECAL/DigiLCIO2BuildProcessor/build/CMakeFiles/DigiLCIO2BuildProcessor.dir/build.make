# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

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
CMAKE_COMMAND = /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/CMake/3.15.5/bin/cmake

# The command to remove a file.
RM = /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/CMake/3.15.5/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL/DigiLCIO2BuildProcessor

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL/DigiLCIO2BuildProcessor/build

# Include any dependencies generated for this target.
include CMakeFiles/DigiLCIO2BuildProcessor.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/DigiLCIO2BuildProcessor.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/DigiLCIO2BuildProcessor.dir/flags.make

CMakeFiles/DigiLCIO2BuildProcessor.dir/src/DigiLCIO2BuildProcessor.cc.o: CMakeFiles/DigiLCIO2BuildProcessor.dir/flags.make
CMakeFiles/DigiLCIO2BuildProcessor.dir/src/DigiLCIO2BuildProcessor.cc.o: ../src/DigiLCIO2BuildProcessor.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL/DigiLCIO2BuildProcessor/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/DigiLCIO2BuildProcessor.dir/src/DigiLCIO2BuildProcessor.cc.o"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/8.2.0-3fa06/x86_64-centos7/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/DigiLCIO2BuildProcessor.dir/src/DigiLCIO2BuildProcessor.cc.o -c /nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL/DigiLCIO2BuildProcessor/src/DigiLCIO2BuildProcessor.cc

CMakeFiles/DigiLCIO2BuildProcessor.dir/src/DigiLCIO2BuildProcessor.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DigiLCIO2BuildProcessor.dir/src/DigiLCIO2BuildProcessor.cc.i"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/8.2.0-3fa06/x86_64-centos7/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL/DigiLCIO2BuildProcessor/src/DigiLCIO2BuildProcessor.cc > CMakeFiles/DigiLCIO2BuildProcessor.dir/src/DigiLCIO2BuildProcessor.cc.i

CMakeFiles/DigiLCIO2BuildProcessor.dir/src/DigiLCIO2BuildProcessor.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DigiLCIO2BuildProcessor.dir/src/DigiLCIO2BuildProcessor.cc.s"
	/cvmfs/sft.cern.ch/lcg/releases/gcc/8.2.0-3fa06/x86_64-centos7/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL/DigiLCIO2BuildProcessor/src/DigiLCIO2BuildProcessor.cc -o CMakeFiles/DigiLCIO2BuildProcessor.dir/src/DigiLCIO2BuildProcessor.cc.s

# Object files for target DigiLCIO2BuildProcessor
DigiLCIO2BuildProcessor_OBJECTS = \
"CMakeFiles/DigiLCIO2BuildProcessor.dir/src/DigiLCIO2BuildProcessor.cc.o"

# External object files for target DigiLCIO2BuildProcessor
DigiLCIO2BuildProcessor_EXTERNAL_OBJECTS =

lib/libDigiLCIO2BuildProcessor.so.0.1.0: CMakeFiles/DigiLCIO2BuildProcessor.dir/src/DigiLCIO2BuildProcessor.cc.o
lib/libDigiLCIO2BuildProcessor.so.0.1.0: CMakeFiles/DigiLCIO2BuildProcessor.dir/build.make
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/Marlin/v01-17-01/lib/libMarlin.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/Marlin/v01-17-01/lib/libMarlinXML.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/lcio/v02-17/lib64/liblcio.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/gear/v01-09/lib/libgearsurf.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/gear/v01-09/lib/libgear.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/gear/v01-09/lib/libgearxml.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/CLHEP/2.3.4.3/lib/libCLHEP-2.3.4.3.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/ilcutil/v01-06-02/lib/libstreamlog.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libCore.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libImt.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libRIO.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libNet.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libHist.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libGraf.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libGraf3d.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libGpad.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libROOTDataFrame.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libTree.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libTreePlayer.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libRint.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libPostscript.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libMatrix.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libPhysics.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libMathCore.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libThread.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libMultiProc.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/lcio/v02-17/lib64/liblcio.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/gear/v01-09/lib/libgearsurf.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/gear/v01-09/lib/libgear.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/gear/v01-09/lib/libgearxml.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/CLHEP/2.3.4.3/lib/libCLHEP-2.3.4.3.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/ilcutil/v01-06-02/lib/libstreamlog.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libCore.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libImt.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libRIO.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libNet.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libHist.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libGraf.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libGraf3d.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libGpad.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libROOTDataFrame.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libTree.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libTreePlayer.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libRint.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libPostscript.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libMatrix.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libPhysics.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libMathCore.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libThread.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/root/6.18.04/lib/libMultiProc.so
lib/libDigiLCIO2BuildProcessor.so.0.1.0: CMakeFiles/DigiLCIO2BuildProcessor.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL/DigiLCIO2BuildProcessor/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library lib/libDigiLCIO2BuildProcessor.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/DigiLCIO2BuildProcessor.dir/link.txt --verbose=$(VERBOSE)
	$(CMAKE_COMMAND) -E cmake_symlink_library lib/libDigiLCIO2BuildProcessor.so.0.1.0 lib/libDigiLCIO2BuildProcessor.so.0.1 lib/libDigiLCIO2BuildProcessor.so

lib/libDigiLCIO2BuildProcessor.so.0.1: lib/libDigiLCIO2BuildProcessor.so.0.1.0
	@$(CMAKE_COMMAND) -E touch_nocreate lib/libDigiLCIO2BuildProcessor.so.0.1

lib/libDigiLCIO2BuildProcessor.so: lib/libDigiLCIO2BuildProcessor.so.0.1.0
	@$(CMAKE_COMMAND) -E touch_nocreate lib/libDigiLCIO2BuildProcessor.so

# Rule to build all files generated by this target.
CMakeFiles/DigiLCIO2BuildProcessor.dir/build: lib/libDigiLCIO2BuildProcessor.so

.PHONY : CMakeFiles/DigiLCIO2BuildProcessor.dir/build

CMakeFiles/DigiLCIO2BuildProcessor.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/DigiLCIO2BuildProcessor.dir/cmake_clean.cmake
.PHONY : CMakeFiles/DigiLCIO2BuildProcessor.dir/clean

CMakeFiles/DigiLCIO2BuildProcessor.dir/depend:
	cd /nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL/DigiLCIO2BuildProcessor/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL/DigiLCIO2BuildProcessor /nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL/DigiLCIO2BuildProcessor /nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL/DigiLCIO2BuildProcessor/build /nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL/DigiLCIO2BuildProcessor/build /nfs/dust/ilc/user/marquezh/SiWECAL_p_AHCAL_Sim/processors/ECAL/DigiLCIO2BuildProcessor/build/CMakeFiles/DigiLCIO2BuildProcessor.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/DigiLCIO2BuildProcessor.dir/depend
