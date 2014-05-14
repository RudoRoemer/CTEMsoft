project(Apps)





add_executable(CTEMEBSDmaster ${Apps_SOURCE_DIR}/CTEMEBSDmaster.f90)
TARGET_LINK_LIBRARIES(CTEMEBSDmaster CTEMSoftLib local)