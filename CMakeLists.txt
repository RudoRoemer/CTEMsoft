
if(NOT "${EMSoft_SDK}" STREQUAL "")
  include("${EMSoft_SDK}/EMSoft_SDK.cmake")
  message(STATUS "EMSoft_SDK Location: ${EMSoft_SDK}")
endif()


cmake_minimum_required (VERSION 2.8.12)
project (EMsoft C CXX Fortran)

if(${CMAKE_MAJOR_VERSION} STREQUAL "3")
  # CMake 3.0 introduces VERSION variables for each Project. We do our own version numbering so we need to have
  # cmake ignore that policy
  cmake_policy(SET CMP0048 OLD)
  # CMake 3.0 introduces some OS X RPath policies
  cmake_policy(SET CMP0042 OLD)
endif()


include(CMakeParseArguments)
  
set(EMsoft_VER_MAJOR "3")
set(EMsoft_VER_MINOR "0")
set(EMsoft_VER_PATCH "1")


# ---------- Setup output Directories -------------------------
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
  ${PROJECT_BINARY_DIR}/Bin
  CACHE PATH
  "Single Directory for all Libraries"
  )

# --------- Setup the Executable output Directory -------------
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
  ${PROJECT_BINARY_DIR}/Bin
  CACHE PATH
  "Single Directory for all Executables."
  )

# --------- Setup the Executable output Directory -------------
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
  ${PROJECT_BINARY_DIR}/Bin
  CACHE PATH
  "Single Directory for all static libraries."
  )

# Build shared libraries
OPTION(BUILD_SHARED_LIBS "Build Shared Libraries" OFF)
# changed to SHARED on 10/23/15, MDG
set(LIB_TYPE SHARED)
if(BUILD_SHARED_LIBS)
    set(LIB_TYPE SHARED)
    list(APPEND CMP_LIB_SEARCH_DIRS  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} )
    if(MSVC)
        list(APPEND CMP_LIB_SEARCH_DIRS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Debug
                                    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Release )
    endif()
endif(BUILD_SHARED_LIBS)




set(CMP_SOURCE_DIR ${EMsoft_SOURCE_DIR}/Support/cmp)
INCLUDE (${CMP_SOURCE_DIR}/cmpCMakeMacros.cmake )
INCLUDE (${CMP_SOURCE_DIR}/cmpProject.cmake)

# --------------------------------------------------------------------
# Over ride CMake's built in module directory by prepending cmp's module
# directory first
set(CMAKE_MODULE_PATH ${EMsoft_SOURCE_DIR}/Support/CMakeModules ${CMAKE_MODULE_PATH})


# set(EXE_DEBUG_EXTENSION "_debug")
set(EXE_DEBUG_EXTENSION "")
	   
# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp -fno-sign-zero")
# set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
  set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -funroll-all-loops -fno-f2c -O3 -fno-sign-zero -flat_namespace")
  set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fno-f2c -O0 -g -fno-sign-zero -flat_namespace")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3 -openmp")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g -openmp")
elseif (Fortran_COMPILER_NAME MATCHES "g77")
  # g77
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")


# ---------- Find HDF5 Headers/Libraries -----------------------
# HDF5 now comes with everything that is needed for CMake to load
# up the targets (Exported) that it needs. We just need to find where
# HDF5 is installed.
include(${CMP_SOURCE_DIR}/ExtLib/HDF5Support.cmake)

# ---------- Find Json Fortran Headers/Libraries -----------------------
# Json-Fortran comes with everything that cmake needs to determine the
# include directories, libraries and other items. One only needs to put
# the correct path onto the CMAKE_PREFIX_PATH variable.
find_package(jsonfortran-gnu REQUIRED)
if( NOT jsonfortran-gnu_FOUND)
  message(STATUS "jsonfortran is REQUIRED for this project.")
  message(STATUS "jsonfortran source repository is at https://github.com/jacobwilliams/json-fortran")
  message(FATAL_ERROR "Please Download, Build and install. After install export the environment variable JSONFORTRAN_INSTALL to point to the installation location.")
else()
  GET_FILENAME_COMPONENT (jsonfortran_LIBRARY_DIRS "${jsonfortran_INCLUDE_DIRS}" PATH)
  set(jsonfortran_LIBRARY_DIRS ${jsonfortran_LIBRARY_DIRS}/lib)
  message(STATUS "jsonfortran-gnu Location: ${JSONFORTRAN_INSTALL}")
  message(STATUS "jsonfortran-gnu Version: ${jsonfortran_VERSION}")
  message(STATUS "jsonfortran-gnu LIBRARY DIR: ${jsonfortran_LIBRARY_DIRS}")
  include_directories(${jsonfortran_INCLUDE_DIRS})
endif()

# --------------------------------------------------------------------
# Generate our ReadMe and License Files
configure_file(${PROJECT_SOURCE_DIR}/License.txt
                ${PROJECT_BINARY_DIR}/License.txt )
configure_file(${PROJECT_SOURCE_DIR}/ReadMe.md
                ${PROJECT_BINARY_DIR}/ReadMe.md )
install(FILES ${PROJECT_SOURCE_DIR}/ReadMe.md ${PROJECT_SOURCE_DIR}/License.txt
        DESTINATION .
        COMPONENT Applications)


include_directories(${EMsoft_BINARY_DIR}/src ${JSONFORTRAN_INCLUDE_DIR})
#find_package (FFTW3 COMPONENTS single double)


## THIS WILL BE NEEDED FOR WINDOWS
include(CMakeAddFortranSubdirectory)
# add the lapack subdirectory as a fortran project
# the subdir is lapack, the project is LAPACK
cmake_add_fortran_subdirectory(src
  NO_EXTERNAL_INSTALL
  PROJECT EMSoftLib  # project name in toplevel CMakeLists.txt in lapack
  #ARCHIVE_DIR lapack/lib # .lib location relative to root binary tree
  #RUNTIME_DIR lapack/bin # .dll location relative to root binary tree
  LIBRARIES EMSoftLib # target libraries created
  LINK_LIBRARIES blas lapack  # link interface libraries
  LINK_LIBS EMSoftLib blas lapack 
)

cmake_add_fortran_subdirectory(src_SampleRFZ
  NO_EXTERNAL_INSTALL
  PROJECT SampleRFZLib  # project name in toplevel CMakeLists.txt in lapack
  LIBRARIES  SampleRFZLib# target libraries created
  LINK_LIBRARIES # link interface libraries
  LINK_LIBS 
)


# --------------------------------------------------------------------
# Package the needed examples
include(${PROJECT_SOURCE_DIR}/examples/SourceList.cmake)

# --------------------------------------------------------------------
# Package the needed IDL Files
include(${PROJECT_SOURCE_DIR}/IDL/SourceList.cmake)

# --------------------------------------------------------------------
# Package the needed PDF Manuals
include(${PROJECT_SOURCE_DIR}/manuals/CMakeLists.txt)

# --------------------------------------------------------------------
# Package the OpenCL files
include(${PROJECT_SOURCE_DIR}/opencl/SourceList.cmake)

# --------------------------------------------------------------------
# Package the needed resources
include(${PROJECT_SOURCE_DIR}/resources/SourceList.cmake)


# --------------------------------------------------------------------
# Package the needed examples
include(${PROJECT_SOURCE_DIR}/templatefolder/SourceList.cmake)

# --------------------------------------------------------------------
# Add the executables
ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/src_programs ${PROJECT_BINARY_DIR}/src_programs)

# --------------------------------------------------------------------
# This should be the last line in this file:
include(${EMsoft_SOURCE_DIR}/Support/CPack/PackageProject.cmake)


