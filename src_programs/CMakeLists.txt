
#---------------------------------------------------------------------
# Set some variables to shorten up the call to the function below
set(APP_DIR "${EMsoft_SOURCE_DIR}/src_programs")
set(TMPLT_DIR "${EMsoft_SOURCE_DIR}/templatefolder")
set(EXE_LINK_LIBRARIES EMSoftLib gomp gcc_eh jsonfortran)

#---------------------------------------------------------------------
# This function creates an executable that is to be compiled. The valid
# arguments are:
#  TARGET: the name of the executable
#  SOURCES: List of Fortran sources to be compiled
#  LINK_LIBRARIES: The list of Libraries the TARGET needs to be linked to
function(Add_EMSoft_Executable)
  set(options )
  set(oneValueArgs TARGET TEMPLATE)
  set(multiValueArgs SOURCES LINK_LIBRARIES)
  cmake_parse_arguments(Z "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  option(EMSOFT_BUILD_${Z_TARGET} "Build ${Z_TARGET}" ON)
  mark_as_advanced(EMSOFT_BUILD_${Z_TARGET})

  if(EMSOFT_BUILD_${Z_TARGET})  
    set(install_dir "bin")
    set(lib_install_dir "lib")
    if(WIN32)
      set(install_dir ".")
      set(lib_install_dir "")
    endif()

    BuildToolBundle(TARGET ${Z_TARGET}
                  DEBUG_EXTENSION ${EXE_DEBUG_EXTENSION}
                  VERSION_MAJOR ${EMsoft_VER_MAJOR}
                  VERSION_MINOR ${EMsoft_VER_MINOR}
                  VERSION_PATCH ${EMsoft_VER_PATCH}
                  BINARY_DIR ${EMsoft_BINARY_DIR}/Applications/${Z_TARGET}
                  COMPONENT Applications
                  INSTALL_DEST "${install_dir}"
                  SOLUTION_FOLDER "Applications"
                  SOURCES ${Z_SOURCES}
                  LINK_LIBRARIES ${Z_LINK_LIBRARIES}
                  LIB_SEARCH_DIRS ${EMsoft_BINARY_DIR}/Bin
    )

  endif()

  if(NOT "${Z_TEMPLATE}" STREQUAL "" AND NOT EXISTS ${Z_TEMPLATE})
    message(STATUS "Missing Template File for Exectuable ${Z_TARGET}")
    message(STATUS "  ${Z_TEMPLATE}")
  else()
    install(FILES ${Z_TEMPLATE}
      DESTINATION "templatefolder"
      COMPONENT Applications
    )
  endif()

endfunction()



#---------------------------------------------------------------------
# Add all of our executables.
#Add_EMSoft_Executable(TARGET EMECPthinfilm SOURCES ${APP_DIR}/EMECPthinfilm.f90 LINK_LIBRARIES ${EXE_LINK_LIBRARIES} TEMPLATE)
#Add_EMSoft_Executable(TARGET EMECPSMMultiLayer SOURCES ${APP_DIR}/EMECPSMMultiLayer.f90 LINK_LIBRARIES ${EXE_LINK_LIBRARIES} TEMPLATE)
#Add_EMSoft_Executable(TARGET EMECCI SOURCES ${APP_DIR}/EMECCI.f90 LINK_LIBRARIES ${EXE_LINK_LIBRARIES} TEMPLATE) #read_stacking_fault_data has wrong number of args passed to it
#Add_EMSoft_Executable(TARGET EMKossel SOURCES ${APP_DIR}/EMKossel.f90 LINK_LIBRARIES ${EXE_LINK_LIBRARIES} TEMPLATE)
#Add_EMSoft_Executable(TARGET EMEBSDoverlap SOURCES ${APP_DIR}/EMEBSDoverlap.f90 LINK_LIBRARIES ${EXE_LINK_LIBRARIES} TEMPLATE)


Add_EMSoft_Executable(TARGET EMECP SOURCES ${APP_DIR}/EMECP.f90 LINK_LIBRARIES ${EXE_LINK_LIBRARIES} TEMPLATE ${TMPLT_DIR}/EMECP.template)
Add_EMSoft_Executable(TARGET EMECPmaster SOURCES ${APP_DIR}/EMECPmaster.f90 LINK_LIBRARIES ${EXE_LINK_LIBRARIES} TEMPLATE ${TMPLT_DIR}/EMECPmaster.template)
Add_EMSoft_Executable(TARGET EMKosselMaster SOURCES ${APP_DIR}/EMKosselMaster.f90 LINK_LIBRARIES ${EXE_LINK_LIBRARIES} TEMPLATE ${TMPLT_DIR}/EMKosselMaster.template)
Add_EMSoft_Executable(TARGET EMsampleRFZ SOURCES ${APP_DIR}/EMsampleRFZ.f90 LINK_LIBRARIES ${EXE_LINK_LIBRARIES} TEMPLATE ${TMPLT_DIR}/EMsampleRFZ.template)
#Add_EMSoft_Executable(TARGET EMMC SOURCES ${APP_DIR}/EMMC.f90 LINK_LIBRARIES ${EXE_LINK_LIBRARIES} TEMPLATE ${TMPLT_DIR}/EMMC.template)
Add_EMSoft_Executable(TARGET EMEBSD SOURCES ${APP_DIR}/EMEBSD.f90 LINK_LIBRARIES ${EXE_LINK_LIBRARIES} TEMPLATE ${TMPLT_DIR}/EMEBSD.template)
Add_EMSoft_Executable(TARGET EMEBSDmaster SOURCES ${APP_DIR}/EMEBSDmaster.f90 LINK_LIBRARIES ${EXE_LINK_LIBRARIES} TEMPLATE ${TMPLT_DIR}/EMEBSDmaster.template)
#Add_EMSoft_Executable(TARGET EMPEDKIN SOURCES ${APP_DIR}/EMPEDKIN.f90 LINK_LIBRARIES ${EXE_LINK_LIBRARIES} TEMPLATE ${TMPLT_DIR}/EMPEDKIN.template)
Add_EMSoft_Executable(TARGET rotationtest SOURCES ${APP_DIR}/rotationtest.f90 LINK_LIBRARIES ${EXE_LINK_LIBRARIES} TEMPLATE ${TMPLT_DIR}/EMECP.template)
Add_EMSoft_Executable(TARGET EMmkxtal SOURCES ${APP_DIR}/EMmkxtal.f90 LINK_LIBRARIES ${EXE_LINK_LIBRARIES})
Add_EMSoft_Executable(TARGET EMlistSG SOURCES ${APP_DIR}/EMlistSG.f90 LINK_LIBRARIES ${EXE_LINK_LIBRARIES})
Add_EMSoft_Executable(TARGET EMqg SOURCES ${APP_DIR}/EMqg.f90 LINK_LIBRARIES ${EXE_LINK_LIBRARIES})
Add_EMSoft_Executable(TARGET EMfamily SOURCES ${APP_DIR}/EMfamily.f90 LINK_LIBRARIES ${EXE_LINK_LIBRARIES})
Add_EMSoft_Executable(TARGET EMstar SOURCES ${APP_DIR}/EMstar.f90 LINK_LIBRARIES ${EXE_LINK_LIBRARIES})
Add_EMSoft_Executable(TARGET EMorbit SOURCES ${APP_DIR}/EMorbit.f90 LINK_LIBRARIES ${EXE_LINK_LIBRARIES})
Add_EMSoft_Executable(TARGET EMZAgeom SOURCES ${APP_DIR}/EMZAgeom.f90 LINK_LIBRARIES ${EXE_LINK_LIBRARIES})
Add_EMSoft_Executable(TARGET EMlatgeom SOURCES ${APP_DIR}/EMlatgeom.f90 LINK_LIBRARIES ${EXE_LINK_LIBRARIES})
Add_EMSoft_Executable(TARGET EMstereo SOURCES ${APP_DIR}/EMstereo.f90 LINK_LIBRARIES ${EXE_LINK_LIBRARIES})
Add_EMSoft_Executable(TARGET EMorient SOURCES ${APP_DIR}/EMorient.f90 LINK_LIBRARIES ${EXE_LINK_LIBRARIES})
Add_EMSoft_Executable(TARGET EMxtalinfo SOURCES ${APP_DIR}/EMxtalinfo.f90 LINK_LIBRARIES ${EXE_LINK_LIBRARIES})
Add_EMSoft_Executable(TARGET EMzap SOURCES ${APP_DIR}/EMzap.f90 LINK_LIBRARIES ${EXE_LINK_LIBRARIES})
Add_EMSoft_Executable(TARGET EMdrawcell SOURCES ${APP_DIR}/EMdrawcell.f90 LINK_LIBRARIES ${EXE_LINK_LIBRARIES})


#---------------------------------------------------------------------
# Copy the Rotations.txt file from Resources into the build directory
ADD_CUSTOM_TARGET(ZZ_RotationsFile-Copy ALL
                    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${EMsoft_SOURCE_DIR}/resources/rotations.txt"
                    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${INTER_DIR}/
                    )
add_dependencies(rotationtest ZZ_RotationsFile-Copy )

#------------------------------------------------------------------------------
# Find the BLAS/LAPACK library
# On mac systems, we likely have to link against the vecLib framework
if(APPLE)
  find_library(LAPACK_LIB Accelerate)
 # include_directories(" -FvecLib")
elseif(UNIX AND NOT APPLE)
  find_package(LAPACK REQUIRED)
  if(LAPACK_FOUND)
    set(LAPACK_LIB ${LAPACK_LIBRARIES})
  else()
    message(STATUS "LAPACK NOT Found for LINUX. Install LAPACK/BLAS on your system")
  endif()
else()
	message(FATAL_ERROR "This platform needs to have CMake code inserted to find the BLAS/LAPACK Libraries.")
endif(APPLE)

#------------------------------------------------------------------------------
# Fine the OpenCL Package
find_package( OpenCL REQUIRED )
if(OpenCL_FOUND)
  message(STATUS "OpenCL_FOUND: ${OpenCL_FOUND}")
  message(STATUS "OpenCL_VERSION_STRING: ${OpenCL_VERSION_STRING}")
  message(STATUS "OpenCL_INCLUDE_DIRS: ${OpenCL_INCLUDE_DIRS}")
  message(STATUS "OpenCL_LIBRARIES: ${OpenCL_LIBRARIES}")
#  message(STATUS "OpenCL_INCLUDE_DIR: ${OpenCL_INCLUDE_DIR}")
#  message(STATUS "OpenCL_LIBRARY: ${OpenCL_LIBRARY}")
else()
  message(FATAL_ERROR "OpenCL is needed to compile some programs in EMSoft. Please install a package appropriate for your Operating System")
endif()

include_directories( ${OpenCL_INCLUDE_DIRS} )

# Figure out if the OpenCL Package has CPP bindings
if( OPENCL_HAS_CPP_BINDINGS )
  message(STATUS "OpenCL has CPP bindings. Full include is: " ${OpenCL_INCLUDE_DIRS} )
else( OPENCL_HAS_CPP_BINDINGS )
  message(STATUS "No OpenCL CPP bindings found" )
endif( OPENCL_HAS_CPP_BINDINGS )

#------------------------------------------------------------------------------
# Fine the Fortran OpenCL Bindings Package
find_package(FortranCL REQUIRED)
if( NOT FORTRANCL_FOUND)
  message(STATUS "FortranCL is REQUIRED for this project.")
  message(STATUS "FortranCL source repository is at http://code.google.com/p/fortrancl/downloads/list")
  message(FATAL_ERROR "Please Download, Build and install. After install export the environment variable FORTRANCL_INSTALL to point to the installation location.")
endif()


include_directories(${FORTRANCL_INCLUDE_DIR})
include_directories(${HDF5_INCLUDE_DIR})

# Add executables to be compiled.
Add_EMSoft_Executable(TARGET EMMCOpenCL 
                      SOURCES ${APP_DIR}/EMMCOpenCL.f90 
                      LINK_LIBRARIES ${EXE_LINK_LIBRARIES} ${LAPACK_LIB} ${FORTRANCL_LIBRARIES} ${OpenCL_LIBRARY}
                      TEMPLATE ${TMPLT_DIR}/EMMCOpenCL.template
                      )
#Add_EMSoft_Executable(TARGET EMEBSDmasterOpenCL 
#                      SOURCES ${APP_DIR}/EMEBSDmasterOpenCL.f90 
#                      LINK_LIBRARIES ${EXE_LINK_LIBRARIES} ${LAPACK_LIB} ${FORTRANCL_LIBRARIES} ${OpenCL_LIBRARY}
#                      TEMPLATE ${TMPLT_DIR}/EMEBSDmaster.template
#                      )
