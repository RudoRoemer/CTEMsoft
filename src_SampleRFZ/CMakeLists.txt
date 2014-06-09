project(SampleRFZLib)


set(SampleRFZLib_SRCS
	${SampleRFZLib_SOURCE_DIR}/sampler.f90
)

ADD_LIBRARY(SampleRFZLib ${LIB_TYPE} ${SampleRFZLib_SRCS})
