project(EMSoftLib)

configure_file(${EMSoftLib_SOURCE_DIR}/local.f90.in 
			   ${EMSoftLib_BINARY_DIR}/local.f90 @ONLY)

set(EMSoftLib_SRCS
	${EMSoftLib_SOURCE_DIR}/apb.f90
	${EMSoftLib_SOURCE_DIR}/constants.f90
	${EMSoftLib_SOURCE_DIR}/crystal.f90
	${EMSoftLib_SOURCE_DIR}/defectmodule.f90
	${EMSoftLib_SOURCE_DIR}/diffraction.f90
	${EMSoftLib_SOURCE_DIR}/dislocation.f90
	${EMSoftLib_SOURCE_DIR}/eispack.f90
	${EMSoftLib_SOURCE_DIR}/error.f90
	${EMSoftLib_SOURCE_DIR}/files.f90
	${EMSoftLib_SOURCE_DIR}/foilmodule.f90
	${EMSoftLib_SOURCE_DIR}/graphics.f90
	${EMSoftLib_SOURCE_DIR}/gvectors.f90
	#${EMSoftLib_SOURCE_DIR}/HDF5Utilities.f90 #Need HDF5 library from somewhere
	${EMSoftLib_SOURCE_DIR}/inclusion.f90
	${EMSoftLib_SOURCE_DIR}/io.f90
	${EMSoftLib_SOURCE_DIR}/kvectors.f90
	${EMSoftLib_SOURCE_DIR}/initializers.f90
	${EMSoftLib_SOURCE_DIR}/Lambert.f90
#	${EMSoftLib_SOURCE_DIR}/local.f90
	${EMSoftLib_BINARY_DIR}/local.f90  # NOTE THIS FILE IS GENERATED DURING CMAKE TIME
	${EMSoftLib_SOURCE_DIR}/math.f90
	${EMSoftLib_SOURCE_DIR}/multibeams.f90
	${EMSoftLib_SOURCE_DIR}/EBSDmod.f90
        ${EMSoftLib_SOURCE_DIR}/ECPmod.f90
	${EMSoftLib_SOURCE_DIR}/MBmodule.f90
 	${EMSoftLib_SOURCE_DIR}/NameListTypedefs.f90
 	${EMSoftLib_SOURCE_DIR}/NameListHandlers.f90
 	${EMSoftLib_SOURCE_DIR}/NameListHDFwriters.f90
 	${EMSoftLib_SOURCE_DIR}/HDFsupport.f90
 	${EMSoftLib_SOURCE_DIR}/JSONsupport.f90
 	${EMSoftLib_SOURCE_DIR}/noise.f90
	${EMSoftLib_SOURCE_DIR}/others.f90
	${EMSoftLib_SOURCE_DIR}/pgm.f90
	${EMSoftLib_SOURCE_DIR}/postscript.f90
	${EMSoftLib_SOURCE_DIR}/quaternions.f90
	${EMSoftLib_SOURCE_DIR}/rksuite.f
	${EMSoftLib_SOURCE_DIR}/rng.f95
	${EMSoftLib_SOURCE_DIR}/rotations.f90
	${EMSoftLib_SOURCE_DIR}/so3.f90
	${EMSoftLib_SOURCE_DIR}/dictmod.f90
	${EMSoftLib_SOURCE_DIR}/spsort.f
	${EMSoftLib_SOURCE_DIR}/stacking_fault.f90
	${EMSoftLib_SOURCE_DIR}/STEMmodule.f90
	${EMSoftLib_SOURCE_DIR}/symmetry.f90
	${EMSoftLib_SOURCE_DIR}/tiff.f90
	${EMSoftLib_SOURCE_DIR}/timing.f90
	${EMSoftLib_SOURCE_DIR}/typedefs.f90
	${EMSoftLib_SOURCE_DIR}/void.f90
	${EMSoftLib_SOURCE_DIR}/YSHmodule.f90
)

# to be removed
#set(EMSoftdyLib_SRCS
#	${EMSoftLib_SOURCE_DIR}/constants.f90
#	${EMSoftLib_SOURCE_DIR}/symmetry.f90
#	${EMSoftLib_SOURCE_DIR}/Lambert.f90
#	${EMSoftLib_SOURCE_DIR}/crystal.f90
#	${EMSoftLib_SOURCE_DIR}/math.f90
#	${EMSoftLib_SOURCE_DIR}/error.f90
##	${EMSoftLib_SOURCE_DIR}/local.f90
#	${EMSoftLib_BINARY_DIR}/local.f90  # NOTE THIS FILE IS GENERATED DURING CMAKE TIME
#	${EMSoftLib_SOURCE_DIR}/EMdymod.f90
#	${EMSoftLib_SOURCE_DIR}/io.f90
# 	${EMSoftLib_SOURCE_DIR}/HDFsupport.f90
#	${EMSoftLib_SOURCE_DIR}/quaternions.f90
#	${EMSoftLib_SOURCE_DIR}/rotations.f90
#	${EMSoftLib_SOURCE_DIR}/timing.f90
#	${EMSoftLib_SOURCE_DIR}/typedefs.f90
#)

message(STATUS “JSONFORTRAN_LIBRARY = ${JSONFORTRAN_LIBRARY}”)

ADD_LIBRARY(EMSoftLib ${LIB_TYPE} ${EMSoftLib_SRCS})
TARGET_LINK_LIBRARIES(EMSoftLib blas lapack ${HDF5_LIBRARIES} ${JSONFORTRAN_LIBRARY})
#ADD_LIBRARY(EMSoftdyLib SHARED ${EMSoftdyLib_SRCS})
#ADD_LIBRARY(local ${EMSoftLib_SOURCE_DIR}/local.f90)
#TARGET_LINK_LIBRARIES(EMSoftdyLib blas lapack ${HDF5_LIBRARIES} ${JSONFORTRAN_LIBRARY})

set(install_dir "tools")
set(lib_install_dir "lib")
if(WIN32)
    set(install_dir ".")
    set(lib_install_dir "")
endif()


if(BUILD_SHARED_LIBS)
    INSTALL(TARGETS ${PROJECT_NAME}
        COMPONENT Applications
        RUNTIME DESTINATION ${install_dir}
        LIBRARY DESTINATION ${lib_install_dir}
        ARCHIVE DESTINATION ${lib_install_dir}
        BUNDLE DESTINATION "."
    )
endif()
