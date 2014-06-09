project(Apps)



include(FortranCInterface)
FortranCInterface_HEADER(SampleRFZ_Mangling.h MACRO_NAMESPACE "SAMPLERFZ_")
#This creates a "FC.h" header that defines mangling macros
# FC_GLOBAL(), FC_GLOBAL_(), FC_MODULE(), and FC_MODULE_().


INCLUDE_DIRECTORIES( ${Apps_BINARY_DIR} )

set(BridgeTest_SRCS
	${Apps_SOURCE_DIR}/BridgeTest.cpp
)


add_executable(BridgeTest ${BridgeTest_SRCS})
TARGET_LINK_LIBRARIES(BridgeTest SampleRFZLib)
