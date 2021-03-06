

cmake_minimum_required (VERSION 2.8.12)
project (CTEMsoft C CXX Fortran)
#enable_language (Fortran)

include(CMakeParseArguments)
#find_package(OpenMP_Fortran)

   
   

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
set(LIB_TYPE STATIC)
if(BUILD_SHARED_LIBS)
	message(FATAL_ERROR "Version 2 Only supports Static Libraries. Please turn OFF BUILD_SHARED_LIBS")
    set(LIB_TYPE SHARED)
    list(APPEND CMP_LIB_SEARCH_DIRS  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} )
    if(MSVC)
        list(APPEND CMP_LIB_SEARCH_DIRS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Debug
                                    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Release )
    endif()
endif(BUILD_SHARED_LIBS)




set(CMP_SOURCE_DIR ${CTEMsoft_SOURCE_DIR}/Support/cmp)
INCLUDE (${CMP_SOURCE_DIR}/cmpCMakeMacros.cmake )
INCLUDE (${CMP_SOURCE_DIR}/cmpProject.cmake)

# --------------------------------------------------------------------
# Over ride CMake's built in module directory by prepending cmp's module
# directory first
set(CMAKE_MODULE_PATH ${CTEMsoft_SOURCE_DIR}/Support/CMakeModules ${CMAKE_MODULE_PATH})


set(EXE_DEBUG_EXTENSION "_debug")

set(CTEMsoft_VER_MAJOR "2")
set(CTEMsoft_VER_MINOR "1")
set(CTEMsoft_VER_PATCH "0")


# --------------------------------------------------------------------
# Generate our ReadMe and License Files
configure_file(${PROJECT_SOURCE_DIR}/License.txt
                ${PROJECT_BINARY_DIR}/License.txt )
configure_file(${PROJECT_SOURCE_DIR}/ReadMe.md
                ${PROJECT_BINARY_DIR}/ReadMe.md )
install(FILES ${PROJECT_SOURCE_DIR}/ReadMe.md ${PROJECT_SOURCE_DIR}/License.txt
        DESTINATION .
        COMPONENT Applications)


include_directories(${CTEMsoft_BINARY_DIR}/src)
#find_package (FFTW3 COMPONENTS single double)


## THIS WILL BE NEEDED FOR WINDOWS
include(CMakeAddFortranSubdirectory)
# add the lapack subdirectory as a fortran project
# the subdir is lapack, the project is LAPACK
cmake_add_fortran_subdirectory(src
  NO_EXTERNAL_INSTALL
  PROJECT CTEMSoftLib  # project name in toplevel CMakeLists.txt in lapack
  #ARCHIVE_DIR lapack/lib # .lib location relative to root binary tree
  #RUNTIME_DIR lapack/bin # .dll location relative to root binary tree
  LIBRARIES CTEMSoftLib # target libraries created
  LINK_LIBRARIES blas lapack  # link interface libraries
  LINK_LIBS  CTEMSoftLib blas lapack 
)
	   
# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
  set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -funroll-all-loops -fno-f2c -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fno-f2c -O0 -g")
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

#ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/src ${PROJECT_BINARY_DIR}/src)
ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/src_programs ${PROJECT_BINARY_DIR}/src_programs)

# This should be the last line in this file:
include(${CTEMsoft_SOURCE_DIR}/Support/CPack/PackageProject.cmake)


