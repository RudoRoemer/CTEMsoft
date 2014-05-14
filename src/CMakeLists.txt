project(CTEMSoftLib)


set(CTEMSoftLib_SRCS
	${CTEMSoftLib_SOURCE_DIR}/apb.f90
	${CTEMSoftLib_SOURCE_DIR}/constants.f90
	${CTEMSoftLib_SOURCE_DIR}/crystal.f90
	${CTEMSoftLib_SOURCE_DIR}/crystalvars.f90
	${CTEMSoftLib_SOURCE_DIR}/defectmodule.f90
	${CTEMSoftLib_SOURCE_DIR}/diffraction.f90
	${CTEMSoftLib_SOURCE_DIR}/dislocation.f90
	#${CTEMSoftLib_SOURCE_DIR}/dispfield.f90
	${CTEMSoftLib_SOURCE_DIR}/doublediff.f90
	${CTEMSoftLib_SOURCE_DIR}/dynamical.f90
	#${CTEMSoftLib_SOURCE_DIR}/eispack.f90
	${CTEMSoftLib_SOURCE_DIR}/error.f90
	#${CTEMSoftLib_SOURCE_DIR}/family.f90
	${CTEMSoftLib_SOURCE_DIR}/files.f90
	${CTEMSoftLib_SOURCE_DIR}/foilmodule.f90
	${CTEMSoftLib_SOURCE_DIR}/graphics.f90
	${CTEMSoftLib_SOURCE_DIR}/gvectors.f90
	#${CTEMSoftLib_SOURCE_DIR}/HDF5Utilities.f90
	${CTEMSoftLib_SOURCE_DIR}/inclusion.f90
	${CTEMSoftLib_SOURCE_DIR}/io.f90
	${CTEMSoftLib_SOURCE_DIR}/kvectors.f90
	${CTEMSoftLib_SOURCE_DIR}/Lambert.f90
	${CTEMSoftLib_SOURCE_DIR}/local.f90
	${CTEMSoftLib_SOURCE_DIR}/math.f90
	${CTEMSoftLib_SOURCE_DIR}/multibeams.f90
	${CTEMSoftLib_SOURCE_DIR}/MBmodule.f90
	${CTEMSoftLib_SOURCE_DIR}/noise.f90
	${CTEMSoftLib_SOURCE_DIR}/others.f90
	${CTEMSoftLib_SOURCE_DIR}/pgm.f90
	${CTEMSoftLib_SOURCE_DIR}/postscript.f90
	${CTEMSoftLib_SOURCE_DIR}/quaternions.f90
	${CTEMSoftLib_SOURCE_DIR}/rksuite.f
	${CTEMSoftLib_SOURCE_DIR}/rng.f95
	${CTEMSoftLib_SOURCE_DIR}/rotations.f90
	${CTEMSoftLib_SOURCE_DIR}/spsort.f
	${CTEMSoftLib_SOURCE_DIR}/stacking_fault.f90
	${CTEMSoftLib_SOURCE_DIR}/STEMmodule.f90
	${CTEMSoftLib_SOURCE_DIR}/symmetry.f90
	${CTEMSoftLib_SOURCE_DIR}/symmetryvars.f90
	${CTEMSoftLib_SOURCE_DIR}/tiff.f90
	${CTEMSoftLib_SOURCE_DIR}/timing.f90
	${CTEMSoftLib_SOURCE_DIR}/void.f90
	${CTEMSoftLib_SOURCE_DIR}/YSHmodule.f90
)

ADD_LIBRARY(CTEMSoftLib ${LIB_TYPE} ${CTEMSoftLib_SRCS})
TARGET_LINK_LIBRARIES(CTEMSoftLib blas lapack)

