# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

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
CMAKE_SOURCE_DIR = /home/rossmech/projects/Prisms-pf/phaseField/applications/Cahn_Hilliard_with_elasticity

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/rossmech/projects/Prisms-pf/phaseField/applications/Cahn_Hilliard_with_elasticity

# Include any dependencies generated for this target.
include CMakeFiles/main.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/main.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/main.dir/flags.make

CMakeFiles/main.dir/main.cc.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/main.cc.o: main.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/rossmech/projects/Prisms-pf/phaseField/applications/Cahn_Hilliard_with_elasticity/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/main.dir/main.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/main.dir/main.cc.o -c /home/rossmech/projects/Prisms-pf/phaseField/applications/Cahn_Hilliard_with_elasticity/main.cc

CMakeFiles/main.dir/main.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/main.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/rossmech/projects/Prisms-pf/phaseField/applications/Cahn_Hilliard_with_elasticity/main.cc > CMakeFiles/main.dir/main.cc.i

CMakeFiles/main.dir/main.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/main.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/rossmech/projects/Prisms-pf/phaseField/applications/Cahn_Hilliard_with_elasticity/main.cc -o CMakeFiles/main.dir/main.cc.s

CMakeFiles/main.dir/main.cc.o.requires:

.PHONY : CMakeFiles/main.dir/main.cc.o.requires

CMakeFiles/main.dir/main.cc.o.provides: CMakeFiles/main.dir/main.cc.o.requires
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/main.cc.o.provides.build
.PHONY : CMakeFiles/main.dir/main.cc.o.provides

CMakeFiles/main.dir/main.cc.o.provides.build: CMakeFiles/main.dir/main.cc.o


# Object files for target main
main_OBJECTS = \
"CMakeFiles/main.dir/main.cc.o"

# External object files for target main
main_EXTERNAL_OBJECTS =

main: CMakeFiles/main.dir/main.cc.o
main: CMakeFiles/main.dir/build.make
main: /usr/lib/x86_64-linux-gnu/libdeal.ii.g.so.8.5.1
main: ../../libprisms_pf_debug.a
main: /usr/lib/x86_64-linux-gnu/libbz2.so
main: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_usempif08.so
main: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_usempi_ignore_tkr.so
main: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_mpifh.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_pike-blackbox.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_trilinoscouplings.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_piro.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_rol.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_stokhos_muelu.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_stokhos_ifpack2.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_stokhos_amesos2.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_stokhos_tpetra.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_stokhos_sacado.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_stokhos.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_rythmos.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_muelu-adapters.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_muelu-interface.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_muelu.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_moertel.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_locathyra.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_locaepetra.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_localapack.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_loca.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_noxepetra.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_noxlapack.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_nox.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_phalanx.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_intrepid.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_teko.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_stratimikos.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_stratimikosbelos.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_stratimikosaztecoo.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_stratimikosamesos.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_stratimikosml.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_stratimikosifpack.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_ifpack2-adapters.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_ifpack2.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_anasazitpetra.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_ModeLaplace.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_anasaziepetra.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_anasazi.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_komplex.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_amesos2.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_shylu.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_belostpetra.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_belosepetra.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_belos.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_ml.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_ifpack.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_zoltan2.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_pamgen_extras.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_pamgen.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_amesos.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_galeri-xpetra.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_galeri-epetra.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_aztecoo.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_dpliris.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_isorropia.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_optipack.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_xpetra-sup.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_xpetra.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_thyratpetra.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_thyraepetraext.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_thyraepetra.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_thyracore.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_epetraext.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_trilinosss.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_tpetraext.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_tpetrainout.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_tpetra.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_kokkostsqr.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_tpetraclassiclinalg.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_tpetraclassicnodeapi.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_tpetraclassic.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_triutils.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_globipack.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_shards.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_zoltan.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_epetra.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_sacado.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_rtop.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_kokkoskernels.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_teuchoskokkoscomm.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_teuchoskokkoscompat.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_teuchosremainder.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_teuchosnumerics.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_teuchoscomm.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_teuchosparameterlist.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_teuchoscore.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_kokkosalgorithms.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_kokkoscontainers.so
main: /usr/lib/x86_64-linux-gnu/libtrilinos_kokkoscore.so
main: /usr/lib/x86_64-linux-gnu/libsmumps.so
main: /usr/lib/x86_64-linux-gnu/libdmumps.so
main: /usr/lib/x86_64-linux-gnu/libcmumps.so
main: /usr/lib/x86_64-linux-gnu/libzmumps.so
main: /usr/lib/x86_64-linux-gnu/libpord.so
main: /usr/lib/x86_64-linux-gnu/libmumps_common.so
main: /usr/lib/x86_64-linux-gnu/hdf5/openmpi/libhdf5.so
main: /usr/lib/x86_64-linux-gnu/libtbb.so
main: /usr/lib/x86_64-linux-gnu/libz.so
main: /usr/lib/x86_64-linux-gnu/libptscotch.so
main: /usr/lib/x86_64-linux-gnu/libptscotcherr.so
main: /usr/lib/x86_64-linux-gnu/libscotch.so
main: /usr/lib/x86_64-linux-gnu/libscotcherr.so
main: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_cxx.so
main: /usr/lib/x86_64-linux-gnu/libumfpack.so
main: /usr/lib/x86_64-linux-gnu/libcholmod.so
main: /usr/lib/x86_64-linux-gnu/libccolamd.so
main: /usr/lib/x86_64-linux-gnu/libcolamd.so
main: /usr/lib/x86_64-linux-gnu/libcamd.so
main: /usr/lib/x86_64-linux-gnu/libsuitesparseconfig.so
main: /usr/lib/x86_64-linux-gnu/libamd.so
main: /usr/lib/x86_64-linux-gnu/libparpack.so
main: /usr/lib/x86_64-linux-gnu/libarpack.so
main: /usr/lib/x86_64-linux-gnu/libboost_iostreams.so
main: /usr/lib/x86_64-linux-gnu/libboost_serialization.so
main: /usr/lib/x86_64-linux-gnu/libboost_system.so
main: /usr/lib/x86_64-linux-gnu/libboost_thread.so
main: /usr/lib/x86_64-linux-gnu/libboost_regex.so
main: /usr/lib/x86_64-linux-gnu/libboost_chrono.so
main: /usr/lib/x86_64-linux-gnu/libboost_date_time.so
main: /usr/lib/x86_64-linux-gnu/libboost_atomic.so
main: /usr/lib/x86_64-linux-gnu/libgsl.so
main: /usr/lib/x86_64-linux-gnu/libgslcblas.so
main: /usr/lib/x86_64-linux-gnu/hdf5/openmpi/lib/lib/libhdf5_hl.so
main: /usr/lib/x86_64-linux-gnu/hdf5/openmpi/lib/lib/libhdf5.so
main: /usr/lib/x86_64-linux-gnu/libmuparser.so
main: /usr/lib/x86_64-linux-gnu/libnetcdf_c++.so
main: /usr/lib/x86_64-linux-gnu/libnetcdf.so
main: /usr/lib/x86_64-linux-gnu/libTKBO.so
main: /usr/lib/x86_64-linux-gnu/libTKBool.so
main: /usr/lib/x86_64-linux-gnu/libTKBRep.so
main: /usr/lib/x86_64-linux-gnu/libTKernel.so
main: /usr/lib/x86_64-linux-gnu/libTKFeat.so
main: /usr/lib/x86_64-linux-gnu/libTKFillet.so
main: /usr/lib/x86_64-linux-gnu/libTKG2d.so
main: /usr/lib/x86_64-linux-gnu/libTKG3d.so
main: /usr/lib/x86_64-linux-gnu/libTKGeomAlgo.so
main: /usr/lib/x86_64-linux-gnu/libTKGeomBase.so
main: /usr/lib/x86_64-linux-gnu/libTKHLR.so
main: /usr/lib/x86_64-linux-gnu/libTKIGES.so
main: /usr/lib/x86_64-linux-gnu/libTKMath.so
main: /usr/lib/x86_64-linux-gnu/libTKMesh.so
main: /usr/lib/x86_64-linux-gnu/libTKOffset.so
main: /usr/lib/x86_64-linux-gnu/libTKPrim.so
main: /usr/lib/x86_64-linux-gnu/libTKShHealing.so
main: /usr/lib/x86_64-linux-gnu/libTKSTEP.so
main: /usr/lib/x86_64-linux-gnu/libTKSTEPAttr.so
main: /usr/lib/x86_64-linux-gnu/libTKSTEPBase.so
main: /usr/lib/x86_64-linux-gnu/libTKSTEP209.so
main: /usr/lib/x86_64-linux-gnu/libTKSTL.so
main: /usr/lib/x86_64-linux-gnu/libTKTopAlgo.so
main: /usr/lib/x86_64-linux-gnu/libTKXSBase.so
main: /usr/lib/x86_64-linux-gnu/libp4est.so
main: /usr/lib/x86_64-linux-gnu/libsc.so
main: /usr/lib/x86_64-linux-gnu/liblapack.so
main: /usr/lib/x86_64-linux-gnu/libblas.so
main: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so
main: /usr/lib/x86_64-linux-gnu/libslepc.so
main: /usr/lib/x86_64-linux-gnu/libpetsc.so
main: CMakeFiles/main.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/rossmech/projects/Prisms-pf/phaseField/applications/Cahn_Hilliard_with_elasticity/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable main"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/main.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/main.dir/build: main

.PHONY : CMakeFiles/main.dir/build

CMakeFiles/main.dir/requires: CMakeFiles/main.dir/main.cc.o.requires

.PHONY : CMakeFiles/main.dir/requires

CMakeFiles/main.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/main.dir/cmake_clean.cmake
.PHONY : CMakeFiles/main.dir/clean

CMakeFiles/main.dir/depend:
	cd /home/rossmech/projects/Prisms-pf/phaseField/applications/Cahn_Hilliard_with_elasticity && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/rossmech/projects/Prisms-pf/phaseField/applications/Cahn_Hilliard_with_elasticity /home/rossmech/projects/Prisms-pf/phaseField/applications/Cahn_Hilliard_with_elasticity /home/rossmech/projects/Prisms-pf/phaseField/applications/Cahn_Hilliard_with_elasticity /home/rossmech/projects/Prisms-pf/phaseField/applications/Cahn_Hilliard_with_elasticity /home/rossmech/projects/Prisms-pf/phaseField/applications/Cahn_Hilliard_with_elasticity/CMakeFiles/main.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/main.dir/depend
