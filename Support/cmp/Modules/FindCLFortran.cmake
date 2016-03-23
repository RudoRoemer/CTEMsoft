# - Find szlib
# Find the native CLFORTRAN_LIB includes and library
#
#  CLFORTRAN_INCLUDE_DIR - where to find szlib.h, etc.
#  CLFORTRAN_LIBRARIES   - List of libraries when using szlib.
#  CLFORTRAN_FOUND       - True if szlib found.
set(CLFortran_DEBUG 1)
if(CLFortran_DEBUG)
  MESSAGE(STATUS "Finding CLFortran")
endif()

# Only set CLFORTRAN_INSTALL to the environment variable if it is blank
if("${CLFORTRAN_INSTALL}" STREQUAL "")
    set(CLFORTRAN_INSTALL  $ENV{CLFORTRAN_INSTALL})
endif()


IF (CLFORTRAN_INCLUDE_DIR)
  # Already in cache, be silent
  SET(CLFORTRAN_LIB_FIND_QUIETLY TRUE)
ENDIF (CLFORTRAN_INCLUDE_DIR)

FIND_PATH(CLFORTRAN_INCLUDE_DIR clfortran.mod
  ${CLFORTRAN_INSTALL}/include
  /usr/local/include
  /usr/include
  NO_DEFAULT_PATH
)

SET(CLFORTRAN_LIB_NAMES clfortran)
FIND_LIBRARY(CLFORTRAN_LIBRARY
  NAMES ${CLFORTRAN_LIB_NAMES}
  PATHS
  ${CLFORTRAN_INSTALL}/lib /usr/lib /usr/local/lib
  NO_DEFAULT_PATH
)

if(CLFortran_DEBUG)
  MESSAGE(STATUS "CLFORTRAN_INSTALL: ${CLFORTRAN_INSTALL}")
  MESSAGE(STATUS "CLFORTRAN_LIBRARY: ${CLFORTRAN_LIBRARY}")
  MESSAGE(STATUS "CLFORTRAN_INCLUDE_DIR: ${CLFORTRAN_INCLUDE_DIR}")
ENDIF()

IF (CLFORTRAN_INCLUDE_DIR AND CLFORTRAN_LIBRARY)
   SET(CLFORTRAN_FOUND TRUE)
   SET( CLFORTRAN_LIBRARIES ${CLFORTRAN_LIBRARY} )
ELSE (CLFORTRAN_INCLUDE_DIR AND CLFORTRAN_LIBRARY)
   SET(CLFORTRAN_FOUND FALSE)
   SET( CLFORTRAN_LIBRARIES )
ENDIF (CLFORTRAN_INCLUDE_DIR AND CLFORTRAN_LIBRARY)

IF (CLFORTRAN_FOUND)
  message(STATUS "CLFortran Location: ${CLFORTRAN_INSTALL}")
  message(STATUS "CLFortran Version: ${CLFORTRAN_VERSION}")
  message(STATUS "CLFortran LIBRARY: ${CLFORTRAN_LIBRARY}")

ELSE (CLFORTRAN_FOUND)
  IF (CLFORTRAN_LIB_FIND_REQUIRED)
    MESSAGE(STATUS "Looked for CLFortran libraries named ${CLFORTRAN_LIBS_NAMES}.")
    MESSAGE(FATAL_ERROR "Could NOT find CLFortran library")
  ENDIF (CLFORTRAN_LIB_FIND_REQUIRED)
ENDIF (CLFORTRAN_FOUND)

MARK_AS_ADVANCED(
  CLFORTRAN_LIBRARY
  CLFORTRAN_INCLUDE_DIR
  )