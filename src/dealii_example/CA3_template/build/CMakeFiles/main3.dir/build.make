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
CMAKE_SOURCE_DIR = /mnt/c/Users/shenyuan/Desktop/CA3_template

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/c/Users/shenyuan/Desktop/CA3_template/build

# Include any dependencies generated for this target.
include CMakeFiles/main3.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/main3.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/main3.dir/flags.make

CMakeFiles/main3.dir/main3.cc.o: CMakeFiles/main3.dir/flags.make
CMakeFiles/main3.dir/main3.cc.o: ../main3.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Users/shenyuan/Desktop/CA3_template/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/main3.dir/main3.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/main3.dir/main3.cc.o -c /mnt/c/Users/shenyuan/Desktop/CA3_template/main3.cc

CMakeFiles/main3.dir/main3.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main3.dir/main3.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/shenyuan/Desktop/CA3_template/main3.cc > CMakeFiles/main3.dir/main3.cc.i

CMakeFiles/main3.dir/main3.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main3.dir/main3.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/shenyuan/Desktop/CA3_template/main3.cc -o CMakeFiles/main3.dir/main3.cc.s

CMakeFiles/main3.dir/main3.cc.o.requires:

.PHONY : CMakeFiles/main3.dir/main3.cc.o.requires

CMakeFiles/main3.dir/main3.cc.o.provides: CMakeFiles/main3.dir/main3.cc.o.requires
	$(MAKE) -f CMakeFiles/main3.dir/build.make CMakeFiles/main3.dir/main3.cc.o.provides.build
.PHONY : CMakeFiles/main3.dir/main3.cc.o.provides

CMakeFiles/main3.dir/main3.cc.o.provides.build: CMakeFiles/main3.dir/main3.cc.o


# Object files for target main3
main3_OBJECTS = \
"CMakeFiles/main3.dir/main3.cc.o"

# External object files for target main3
main3_EXTERNAL_OBJECTS =

main3: CMakeFiles/main3.dir/main3.cc.o
main3: CMakeFiles/main3.dir/build.make
main3: /usr/lib/x86_64-linux-gnu/libdeal.ii.g.so.8.5.1
main3: /usr/lib/x86_64-linux-gnu/libbz2.so
main3: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_usempif08.so
main3: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_usempi_ignore_tkr.so
main3: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_mpifh.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_pike-blackbox.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_trilinoscouplings.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_piro.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_rol.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_stokhos_muelu.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_stokhos_ifpack2.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_stokhos_amesos2.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_stokhos_tpetra.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_stokhos_sacado.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_stokhos.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_rythmos.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_muelu-adapters.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_muelu-interface.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_muelu.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_moertel.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_locathyra.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_locaepetra.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_localapack.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_loca.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_noxepetra.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_noxlapack.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_nox.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_phalanx.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_intrepid.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_teko.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_stratimikos.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_stratimikosbelos.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_stratimikosaztecoo.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_stratimikosamesos.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_stratimikosml.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_stratimikosifpack.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_ifpack2-adapters.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_ifpack2.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_anasazitpetra.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_ModeLaplace.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_anasaziepetra.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_anasazi.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_komplex.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_amesos2.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_shylu.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_belostpetra.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_belosepetra.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_belos.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_ml.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_ifpack.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_zoltan2.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_pamgen_extras.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_pamgen.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_amesos.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_galeri-xpetra.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_galeri-epetra.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_aztecoo.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_dpliris.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_isorropia.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_optipack.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_xpetra-sup.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_xpetra.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_thyratpetra.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_thyraepetraext.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_thyraepetra.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_thyracore.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_epetraext.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_trilinosss.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_tpetraext.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_tpetrainout.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_tpetra.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_kokkostsqr.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_tpetraclassiclinalg.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_tpetraclassicnodeapi.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_tpetraclassic.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_triutils.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_globipack.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_shards.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_zoltan.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_epetra.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_sacado.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_rtop.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_kokkoskernels.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_teuchoskokkoscomm.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_teuchoskokkoscompat.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_teuchosremainder.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_teuchosnumerics.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_teuchoscomm.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_teuchosparameterlist.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_teuchoscore.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_kokkosalgorithms.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_kokkoscontainers.so
main3: /usr/lib/x86_64-linux-gnu/libtrilinos_kokkoscore.so
main3: /usr/lib/x86_64-linux-gnu/libsmumps.so
main3: /usr/lib/x86_64-linux-gnu/libdmumps.so
main3: /usr/lib/x86_64-linux-gnu/libcmumps.so
main3: /usr/lib/x86_64-linux-gnu/libzmumps.so
main3: /usr/lib/x86_64-linux-gnu/libpord.so
main3: /usr/lib/x86_64-linux-gnu/libmumps_common.so
main3: /usr/lib/x86_64-linux-gnu/hdf5/openmpi/libhdf5.so
main3: /usr/lib/x86_64-linux-gnu/libtbb.so
main3: /usr/lib/x86_64-linux-gnu/libz.so
main3: /usr/lib/x86_64-linux-gnu/libptscotch.so
main3: /usr/lib/x86_64-linux-gnu/libptscotcherr.so
main3: /usr/lib/x86_64-linux-gnu/libscotch.so
main3: /usr/lib/x86_64-linux-gnu/libscotcherr.so
main3: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_cxx.so
main3: /usr/lib/x86_64-linux-gnu/libumfpack.so
main3: /usr/lib/x86_64-linux-gnu/libcholmod.so
main3: /usr/lib/x86_64-linux-gnu/libccolamd.so
main3: /usr/lib/x86_64-linux-gnu/libcolamd.so
main3: /usr/lib/x86_64-linux-gnu/libcamd.so
main3: /usr/lib/x86_64-linux-gnu/libsuitesparseconfig.so
main3: /usr/lib/x86_64-linux-gnu/libamd.so
main3: /usr/lib/x86_64-linux-gnu/libparpack.so
main3: /usr/lib/x86_64-linux-gnu/libarpack.so
main3: /usr/lib/x86_64-linux-gnu/libboost_iostreams.so
main3: /usr/lib/x86_64-linux-gnu/libboost_serialization.so
main3: /usr/lib/x86_64-linux-gnu/libboost_system.so
main3: /usr/lib/x86_64-linux-gnu/libboost_thread.so
main3: /usr/lib/x86_64-linux-gnu/libboost_regex.so
main3: /usr/lib/x86_64-linux-gnu/libboost_chrono.so
main3: /usr/lib/x86_64-linux-gnu/libboost_date_time.so
main3: /usr/lib/x86_64-linux-gnu/libboost_atomic.so
main3: /usr/lib/x86_64-linux-gnu/libgsl.so
main3: /usr/lib/x86_64-linux-gnu/libgslcblas.so
main3: /usr/lib/x86_64-linux-gnu/hdf5/openmpi/lib/lib/libhdf5_hl.so
main3: /usr/lib/x86_64-linux-gnu/hdf5/openmpi/lib/lib/libhdf5.so
main3: /usr/lib/x86_64-linux-gnu/libmuparser.so
main3: /usr/lib/x86_64-linux-gnu/libnetcdf_c++.so
main3: /usr/lib/x86_64-linux-gnu/libnetcdf.so
main3: /usr/lib/x86_64-linux-gnu/libTKBO.so
main3: /usr/lib/x86_64-linux-gnu/libTKBool.so
main3: /usr/lib/x86_64-linux-gnu/libTKBRep.so
main3: /usr/lib/x86_64-linux-gnu/libTKernel.so
main3: /usr/lib/x86_64-linux-gnu/libTKFeat.so
main3: /usr/lib/x86_64-linux-gnu/libTKFillet.so
main3: /usr/lib/x86_64-linux-gnu/libTKG2d.so
main3: /usr/lib/x86_64-linux-gnu/libTKG3d.so
main3: /usr/lib/x86_64-linux-gnu/libTKGeomAlgo.so
main3: /usr/lib/x86_64-linux-gnu/libTKGeomBase.so
main3: /usr/lib/x86_64-linux-gnu/libTKHLR.so
main3: /usr/lib/x86_64-linux-gnu/libTKIGES.so
main3: /usr/lib/x86_64-linux-gnu/libTKMath.so
main3: /usr/lib/x86_64-linux-gnu/libTKMesh.so
main3: /usr/lib/x86_64-linux-gnu/libTKOffset.so
main3: /usr/lib/x86_64-linux-gnu/libTKPrim.so
main3: /usr/lib/x86_64-linux-gnu/libTKShHealing.so
main3: /usr/lib/x86_64-linux-gnu/libTKSTEP.so
main3: /usr/lib/x86_64-linux-gnu/libTKSTEPAttr.so
main3: /usr/lib/x86_64-linux-gnu/libTKSTEPBase.so
main3: /usr/lib/x86_64-linux-gnu/libTKSTEP209.so
main3: /usr/lib/x86_64-linux-gnu/libTKSTL.so
main3: /usr/lib/x86_64-linux-gnu/libTKTopAlgo.so
main3: /usr/lib/x86_64-linux-gnu/libTKXSBase.so
main3: /usr/lib/x86_64-linux-gnu/libp4est.so
main3: /usr/lib/x86_64-linux-gnu/libsc.so
main3: /usr/lib/x86_64-linux-gnu/liblapack.so
main3: /usr/lib/x86_64-linux-gnu/libblas.so
main3: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so
main3: /usr/lib/x86_64-linux-gnu/libslepc.so
main3: /usr/lib/x86_64-linux-gnu/libpetsc.so
main3: CMakeFiles/main3.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/c/Users/shenyuan/Desktop/CA3_template/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable main3"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/main3.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/main3.dir/build: main3

.PHONY : CMakeFiles/main3.dir/build

CMakeFiles/main3.dir/requires: CMakeFiles/main3.dir/main3.cc.o.requires

.PHONY : CMakeFiles/main3.dir/requires

CMakeFiles/main3.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/main3.dir/cmake_clean.cmake
.PHONY : CMakeFiles/main3.dir/clean

CMakeFiles/main3.dir/depend:
	cd /mnt/c/Users/shenyuan/Desktop/CA3_template/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/c/Users/shenyuan/Desktop/CA3_template /mnt/c/Users/shenyuan/Desktop/CA3_template /mnt/c/Users/shenyuan/Desktop/CA3_template/build /mnt/c/Users/shenyuan/Desktop/CA3_template/build /mnt/c/Users/shenyuan/Desktop/CA3_template/build/CMakeFiles/main3.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/main3.dir/depend

