project(Apps)


macro (Create_exe_and_tll INPUT)
	ADD_EXECUTABLE(${INPUT} ${Apps_SOURCE_DIR}/${INPUT}.f90)
 	TARGET_LINK_LIBRARIES(${INPUT} EMSoftLib gomp gcc_eh)
	
	
	set(install_dir "bin")
	set(lib_install_dir "lib")
	if(WIN32)
		set(install_dir ".")
		set(lib_install_dir "")
	endif()

	INSTALL(TARGETS ${INPUT}
		COMPONENT Applications
		RUNTIME DESTINATION ${install_dir}
		LIBRARY DESTINATION ${lib_install_dir}
		ARCHIVE DESTINATION ${lib_install_dir}
		BUNDLE DESTINATION "."
	)


endmacro (Create_exe_and_tll)


Create_exe_and_tll(EMlacbed)
Create_exe_and_tll(EMmbcbed)
Create_exe_and_tll(EMSRdefect) #MP 
Create_exe_and_tll(EMZAdefect) 
Create_exe_and_tll(EMECP)
Create_exe_and_tll(EMECPmaster)
Create_exe_and_tll(EMECPthinfilm)
Create_exe_and_tll(EMECPSMMultiLayer)
Create_exe_and_tll(EMECCI) #read_stacking_fault_data has wrong number of args passed to it
Create_exe_and_tll(EMKossel)
Create_exe_and_tll(EMKosselMaster)
Create_exe_and_tll(EMsampleRFZ)
Create_exe_and_tll(EMMC)  #MP
Create_exe_and_tll(EMEBSD)
Create_exe_and_tll(EMEBSDoverlap)
Create_exe_and_tll(EMquatSP)
Create_exe_and_tll(EMEBSDmaster)
Create_exe_and_tll(EMEBSDmasterOpenCL)
Create_exe_and_tll(EMmkxtal)
Create_exe_and_tll(EMlistSG)
Create_exe_and_tll(EMqg)
Create_exe_and_tll(EMped)
Create_exe_and_tll(EMfamily)
Create_exe_and_tll(EMstar)
Create_exe_and_tll(EMorbit)
Create_exe_and_tll(EMZAgeom)
Create_exe_and_tll(EMlatgeom)
Create_exe_and_tll(EMstereo)
Create_exe_and_tll(EMorient)
Create_exe_and_tll(EMxtalinfo)
Create_exe_and_tll(EMzap)
Create_exe_and_tll(EMdrawcell)
Create_exe_and_tll(EMBethe)
Create_exe_and_tll(EMPED)
Create_exe_and_tll(EMPEDKIN)

Create_exe_and_tll(rotationtest)
Create_exe_and_tll(hdftest)



 # On mac systems, we likely have to link against the vecLib framework
if(APPLE)
  find_library(LAPACK_LIB Accelerate)
 # include_directories(" -FvecLib")
else()
	message(FATAL_ERROR "This platform needs to have CMake code inserted to find the BLAS/LAPACK Libraries.")
endif(APPLE)


find_package( OpenCL REQUIRED )

include_directories( ${OPENCL_INCLUDE_DIRS} )

if( OPENCL_HAS_CPP_BINDINGS )
  message( "OpenCL has CPP bindings. Full include is: " ${OPENCL_INCLUDE_DIRS} )
else( OPENCL_HAS_CPP_BINDINGS )
  message( "No OpenCL CPP bindings found" )
endif( OPENCL_HAS_CPP_BINDINGS )

find_package(FortranCL REQUIRED)
if( NOT FORTRANCL_FOUND)
  message(STATUS "FortranCL is REQUIRED for this project.")
  message(STATUS "FortranCL source repository is at http://code.google.com/p/fortrancl/downloads/list")
  message(FATAL_ERROR "Please Download, Build and install. After install export the environment variable FORTRANCL_INSTALL to point to the installation location.")
endif()

include_directories(${FORTRANCL_INCLUDE_DIR})

message(STATUS "FORTRANCL_LIBRARIES: ${FORTRANCL_LIBRARIES}")

ADD_EXECUTABLE(EMMCOpenCL EMMCOpenCL.f90)
TARGET_LINK_LIBRARIES(EMMCOpenCL ${OPENCL_LIBRARIES} ${LAPACK_LIB} ${FORTRANCL_LIBRARIES} EMSoftLib)


#Create_exe_and_tll(EMBWshow)
#Create_exe_and_tll(EMECPSM)
#Create_exe_and_tll(EMMCTKD) 
##Create_exe_and_tll(EMSRBW) #Explicit interface required for 'extractbwdata' 
#Create_exe_and_tll(EMSRCBED)
#Create_exe_and_tll(EMTBBFDF)
#Create_exe_and_tll(EMTBBW)
#Create_exe_and_tll(EMTBSM)
#Create_exe_and_tll(EMexpeditions)
#Create_exe_and_tll(EMgamma)
#Create_exe_and_tll(EMgammadefect)
#Create_exe_and_tll(EMholz)
##Create_exe_and_tll(EMlens) #can't open rksuite
##Create_exe_and_tll(EMoverlap) #doesn't look like dbdiffB has been initialized anywhere
#Create_exe_and_tll(EMvr)
##Create_exe_and_tll(EMvrfft)  #need fft? 
