project(Apps)


set(CTEMEBSDmaster_SRCS
	${Apps_SOURCE_DIR}/CTEMEBSDmaster.f90
)


add_executable(CTEMEBSDmaster ${CTEMEBSDmaster_SRCS})
TARGET_LINK_LIBRARIES(CTEMEBSDmaster CTEMSoftLib)