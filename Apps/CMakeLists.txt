project(Apps)


set(BridgeTest_SRCS
	${Apps_SOURCE_DIR}/BridgeTest.cpp
)


add_executable(BridgeTest ${BridgeTest_SRCS})
TARGET_LINK_LIBRARIES(BridgeTest SampleRFZLib)
