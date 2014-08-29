project(Apps)


macro (Create_exe_and_tll INPUT)
	ADD_EXECUTABLE(${INPUT} ${Apps_SOURCE_DIR}/${INPUT}.f90)
 TARGET_LINK_LIBRARIES(${INPUT} CTEMSoftLib gomp)
	
	
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
#add_executable(CTEMEBSDmaster ${Apps_SOURCE_DIR}/CTEMEBSDmaster.f90)
#TARGET_LINK_LIBRARIES(CTEMEBSDmaster CTEMSoftLib)


Create_exe_and_tll(CTEMlacbed)
Create_exe_and_tll(CTEMmbcbed)
Create_exe_and_tll(CTEMSRdefect) #MP 
Create_exe_and_tll(CTEMZAdefect) 
Create_exe_and_tll(CTEMECP)
Create_exe_and_tll(CTEMECPpattern)
Create_exe_and_tll(CTEMECCI) #read_stacking_fault_data has wrong number of args passed to it
Create_exe_and_tll(CTEMKossel)
Create_exe_and_tll(CTEMMC)  #MP
Create_exe_and_tll(CTEMEBSD)
Create_exe_and_tll(CTEMEBSDmaster)
Create_exe_and_tll(CTEMmkxtal)
Create_exe_and_tll(CTEMlistSG)
Create_exe_and_tll(CTEMqg)
Create_exe_and_tll(CTEMfamily)
Create_exe_and_tll(CTEMstar)
Create_exe_and_tll(CTEMorbit)
Create_exe_and_tll(CTEMZAgeom)
Create_exe_and_tll(CTEMlatgeom)
Create_exe_and_tll(CTEMstereo)
Create_exe_and_tll(CTEMorient)
Create_exe_and_tll(CTEMxtalinfo)
Create_exe_and_tll(CTEMzap)
Create_exe_and_tll(CTEMdrawcell)

Create_exe_and_tll(structest)


#Create_exe_and_tll(CTEMBWshow)
#Create_exe_and_tll(CTEMECPSM)
#Create_exe_and_tll(CTEMMCTKD) 
##Create_exe_and_tll(CTEMSRBW) #Explicit interface required for 'extractbwdata' 
#Create_exe_and_tll(CTEMSRCBED)
#Create_exe_and_tll(CTEMTBBFDF)
#Create_exe_and_tll(CTEMTBBW)
#Create_exe_and_tll(CTEMTBSM)
#Create_exe_and_tll(CTEMexpeditions)
#Create_exe_and_tll(CTEMgamma)
#Create_exe_and_tll(CTEMgammadefect)
#Create_exe_and_tll(CTEMholz)
##Create_exe_and_tll(CTEMlens) #can't open rksuite
##Create_exe_and_tll(CTEMoverlap) #doesn't look like dbdiffB has been initialized anywhere
#Create_exe_and_tll(CTEMped)
#Create_exe_and_tll(CTEMvr)
##Create_exe_and_tll(CTEMvrfft)  #need fft? 
