

cmake_minimum_required (VERSION 2.8.12)
project (CTEMsoft C CXX Fortran)
#enable_language (Fortran)

include(CMakeParseArguments)


   
   

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
OPTION(BUILD_SHARED_LIBS "Build Shared Libraries" ON)
set(LIB_TYPE STATIC)
if(BUILD_SHARED_LIBS)
    set(LIB_TYPE SHARED)
    list(APPEND CMP_LIB_SEARCH_DIRS  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} )
    if(MSVC)
        list(APPEND CMP_LIB_SEARCH_DIRS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Debug
                                    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Release )
    endif()
endif(BUILD_SHARED_LIBS)

include_directories(${CTEMsoft_PROJECT_DIR}/src)
include_directories(${CTEMsoft_BINARY_DIR}/src)

## THIS WILL BE NEEDED FOR WINDOWS
	#cmake_minimum_required(VERSION 2.8.7.20120206)
	#project(fortranc)
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

#ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/src ${PROJECT_BINARY_DIR}/src)
ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/src_programs ${PROJECT_BINARY_DIR}/src_programs)
