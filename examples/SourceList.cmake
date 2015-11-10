
#---------------------------------------------------------------------
# Include all the subdirectories that should be packaged.
include(${EMsoft_SOURCE_DIR}/examples/EBSDPatterns/SourceList.cmake)
include(${EMsoft_SOURCE_DIR}/examples/KosselPatterns/SourceList.cmake)

file(COPY "${EMsoft_SOURCE_DIR}/examples" DESTINATION "${PROJECT_BINARY_DIR}/")
