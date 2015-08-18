# Find the jsonfortran includes and library
#
#  JSONFORTRAN_INCLUDE_DIR 
#  JSONFORTRAN_LIBRARIES   
#  JSONFORTRAN_FOUND       
MESSAGE(STATUS "Finding jsonfortran")

#SET(JSONFORTRAN_DEBUG TRUE)


IF (JSONFORTRAN_INCLUDE_DIR)
  # Already in cache, be silent
  SET(JSONFORTRAN_LIB_FIND_QUIETLY TRUE)
ENDIF (JSONFORTRAN_INCLUDE_DIR)

FIND_PATH(JSONFORTRAN_INCLUDE_DIR json_module.mod
  /usr/local/include
  /usr/include
  $ENV{JSONFORTRAN_INSTALL}/include
)

SET(JSONFORTRAN_LIB_NAMES libjsonfortran.a)
FIND_LIBRARY(JSONFORTRAN_LIBRARY
  NAMES ${JSONFORTRAN_LIB_NAMES}
  PATHS  /usr/lib /usr/local/lib $ENV{JSONFORTRAN_INSTALL}/lib
)

if(JSONFORTRAN_DEBUG)
  MESSAGE(STATUS "JSONFORTRAN_LIBRARY: ${JSONFORTRAN_LIBRARY}")
  MESSAGE(STATUS "JSONFORTRAN_INCLUDE_DIR: ${JSONFORTRAN_INCLUDE_DIR}")
ENDIF()


IF (JSONFORTRAN_INCLUDE_DIR AND JSONFORTRAN_LIBRARY)
   SET(JSONFORTRAN_FOUND TRUE)
   SET( JSONFORTRAN_LIBRARIES ${JSONFORTRAN_LIBRARY} )
ELSE (JSONFORTRAN_INCLUDE_DIR AND JSONFORTRAN_LIBRARY)
   SET(JSONFORTRAN_FOUND FALSE)
   SET( JSONFORTRAN_LIBRARIES )
ENDIF (JSONFORTRAN_INCLUDE_DIR AND JSONFORTRAN_LIBRARY)

IF (JSONFORTRAN_FOUND)
  MESSAGE(STATUS "Found JSONFORTRAN_LIB: ${JSONFORTRAN_LIBRARY}")
  MESSAGE(STATUS "Found JSONFORTRAN_INCLUDE_DIR: ${JSONFORTRAN_INCLUDE_DIR}")
ELSE (JSONFORTRAN_FOUND)
  IF (JSONFORTRAN_LIB_FIND_REQUIRED)
    MESSAGE(STATUS "Looked for jsonfortran libraries named ${JSONFORTRAN_LIB_NAMES}.")
    MESSAGE(FATAL_ERROR "Could NOT find jsonfortran library")
  ENDIF (JSONFORTRAN_LIB_FIND_REQUIRED)
ENDIF (JSONFORTRAN_FOUND)

MARK_AS_ADVANCED(
  JSONFORTRAN_LIBRARY
  JSONFORTRAN_INCLUDE_DIR
  )
