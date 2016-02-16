

project(EMSoftLib)

include(FortranCInterface)
FortranCInterface_HEADER(EMsoftLib_Mangling.h MACRO_NAMESPACE "EMSOFTLIB_")

configure_file(${EMsoftLib_SOURCE_DIR}/local.f90.in 
			   ${EMsoftLib_BINARY_DIR}/local.f90 @ONLY)

set(EMsoftLib_SRCS
  ${EMSoftLib_SOURCE_DIR}/apb.f90
  ${EMSoftLib_SOURCE_DIR}/constants.f90
  ${EMSoftLib_SOURCE_DIR}/crystal.f90
  ${EMSoftLib_SOURCE_DIR}/defectmodule.f90
  ${EMSoftLib_SOURCE_DIR}/diffraction.f90
  ${EMSoftLib_SOURCE_DIR}/dislocation.f90
  ${EMSoftLib_SOURCE_DIR}/eispack.f90
  ${EMSoftLib_SOURCE_DIR}/error.f90
  ${EMSoftLib_SOURCE_DIR}/files.f90
  ${EMSoftLib_SOURCE_DIR}/foilmodule.f90
  ${EMSoftLib_SOURCE_DIR}/graphics.f90
  ${EMSoftLib_SOURCE_DIR}/gvectors.f90
  #${EMSoftLib_SOURCE_DIR}/HDF5Utilities.f90 #Need HDF5 library from somewhere
  ${EMSoftLib_SOURCE_DIR}/inclusion.f90
  ${EMSoftLib_SOURCE_DIR}/io.f90
  ${EMSoftLib_SOURCE_DIR}/kvectors.f90
  ${EMSoftLib_SOURCE_DIR}/initializers.f90
  ${EMSoftLib_SOURCE_DIR}/Lambert.f90
# ${EMSoftLib_SOURCE_DIR}/local.f90
  ${EMSoftLib_BINARY_DIR}/local.f90  # NOTE THIS FILE IS GENERATED DURING CMAKE TIME
  ${EMSoftLib_SOURCE_DIR}/math.f90
  ${EMSoftLib_SOURCE_DIR}/multibeams.f90
  ${EMSoftLib_SOURCE_DIR}/EBSDmod.f90
  ${EMSoftLib_SOURCE_DIR}/EMdymod.f90
        ${EMSoftLib_SOURCE_DIR}/ECPmod.f90
  ${EMSoftLib_SOURCE_DIR}/MBmodule.f90
  ${EMSoftLib_SOURCE_DIR}/NameListTypedefs.f90
  ${EMSoftLib_SOURCE_DIR}/NameListHandlers.f90
  ${EMSoftLib_SOURCE_DIR}/NameListHDFwriters.f90
  ${EMSoftLib_SOURCE_DIR}/HDFsupport.f90
  ${EMSoftLib_SOURCE_DIR}/JSONsupport.f90
  ${EMSoftLib_SOURCE_DIR}/noise.f90
  ${EMSoftLib_SOURCE_DIR}/others.f90
  ${EMSoftLib_SOURCE_DIR}/pgm.f90
  ${EMSoftLib_SOURCE_DIR}/postscript.f90
  ${EMSoftLib_SOURCE_DIR}/quaternions.f90
# ${EMSoftLib_SOURCE_DIR}/rksuite.f
  ${EMSoftLib_SOURCE_DIR}/rng.f95
  ${EMSoftLib_SOURCE_DIR}/rotations.f90
  ${EMSoftLib_SOURCE_DIR}/so3.f90
  ${EMSoftLib_SOURCE_DIR}/dictmod.f90
  ${EMSoftLib_SOURCE_DIR}/spsort.f
  ${EMSoftLib_SOURCE_DIR}/stacking_fault.f90
  ${EMSoftLib_SOURCE_DIR}/STEMmodule.f90
  ${EMSoftLib_SOURCE_DIR}/symmetry.f90
  ${EMSoftLib_SOURCE_DIR}/tiff.f90
  ${EMSoftLib_SOURCE_DIR}/timing.f90
  ${EMSoftLib_SOURCE_DIR}/typedefs.f90
  ${EMSoftLib_SOURCE_DIR}/void.f90
  ${EMSoftLib_SOURCE_DIR}/YSHmodule.f90
)


set(install_dir "tools")
set(lib_install_dir "lib")
if(WIN32)
    set(install_dir ".")
    set(lib_install_dir "")
endif()


ADD_LIBRARY(EMsoftLib ${LIB_TYPE} ${EMsoftLib_SRCS})
TARGET_LINK_LIBRARIES(EMsoftLib ${EMsoft_BLAS_LAPACK_LIBS} jsonfortran ${FFTW3_LIBRARY})


INSTALL(TARGETS EMSoftLib
  COMPONENT Applications
  RUNTIME DESTINATION ${install_dir}
  LIBRARY DESTINATION ${lib_install_dir}
  ARCHIVE DESTINATION ${lib_install_dir}
  BUNDLE DESTINATION "."
)



#------------------------------------------------------------------------------
# Install a C/C++ header file to allow C/C++ to interface to EMsoftLib
install(FILES ${EMsoftLib_SOURCE_DIR}/EMsoftLib.h
  DESTINATION include/EMsoftLib
  COMPONENT Applications
  )


#------------------------------------------------------------------------------
# Create the EMsoftHDF5 Library which separates all the HDF5 dependent code
# into its own library.
set(EMsoftLibHDF_SRCS
	${EMsoftLib_SOURCE_DIR}/initializersHDF.f90
 	${EMsoftLib_SOURCE_DIR}/NameListHDFwriters.f90
	${EMsoftLib_SOURCE_DIR}/EBSDmod.f90
 	${EMsoftLib_SOURCE_DIR}/EMdymodHDF.f90
  ${EMsoftLib_SOURCE_DIR}/ECPmod.f90
 	${EMsoftLib_SOURCE_DIR}/HDFsupport.f90
)
message(STATUS "LIB_TYPE: ${LIB_TYPE}")
ADD_LIBRARY(EMsoftLibHDF ${LIB_TYPE} ${EMsoftLibHDF_SRCS})
TARGET_LINK_LIBRARIES(EMsoftLibHDF hdf5 hdf5_f90cstub hdf5_fortran hdf5_hl_fortran EMsoftLib)

if(BUILD_SHARED_LIBS)
	    INSTALL(TARGETS EMsoftLibHDF
        EXPORT EMsoftLibHDFTargets
        COMPONENT Applications
        RUNTIME DESTINATION ${install_dir}
        LIBRARY DESTINATION ${lib_install_dir}
        ARCHIVE DESTINATION ${lib_install_dir}
        BUNDLE DESTINATION "."
    )
endif()


add_executable(SingleEBSDPattern ${EMsoftLib_SOURCE_DIR}/SingleEBSDPattern.cpp)
TARGET_LINK_LIBRARIES(SingleEBSDPattern EMsoftLib)

# --------------------------------------------------------------------
#
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/EMsoft/EMsoftLibTargetsConfigVersion.cmake"
  VERSION ${EMsoft_VERSION}
  COMPATIBILITY AnyNewerVersion
)
if(BUILD_SHARED_LIBS)
  export(EXPORT EMsoftLibTargets
    FILE "${CMAKE_CURRENT_BINARY_DIR}/EMsoft/EMsoftLibTargets.cmake"
    NAMESPACE EMsoft::
  )
  export(EXPORT EMsoftLibHDFTargets
    FILE "${CMAKE_CURRENT_BINARY_DIR}/EMsoft/EMsoftLibHDFTargets.cmake"
    NAMESPACE EMsoft::
  )
endif()
configure_file(cmake/EMsoftConfig.cmake
  "${CMAKE_CURRENT_BINARY_DIR}/EMsoft/EMsoftConfig.cmake"
  COPYONLY
)

set(ConfigPackageLocation lib/cmake/EMsoft)

if(BUILD_SHARED_LIBS)
  install(EXPORT EMsoftLibTargets
    FILE
      EMsoftLibTargets.cmake
    NAMESPACE
      EMsoft::
    DESTINATION
      ${ConfigPackageLocation}
  )
  install(EXPORT EMsoftLibHDFTargets
    FILE
      EMsoftLibHDFTargets.cmake
    NAMESPACE
      EMsoft::
    DESTINATION
      ${ConfigPackageLocation}
  )
endif()

install(
  FILES
    cmake/EMsoftConfig.cmake
    "${CMAKE_CURRENT_BINARY_DIR}/EMsoft/EMsoftLibTargetsConfigVersion.cmake"
  DESTINATION
    ${ConfigPackageLocation}
  COMPONENT
    Devel
)
