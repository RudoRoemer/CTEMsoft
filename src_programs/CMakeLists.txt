project(Apps)

#-------------------------------------------------------------------------------
# Function COMPILE_TOOL to help alleviate lots of extra code below for adding
# simple command line tools that just need one or two source files
#
function(COMPILE_TOOL)
    set(options)
    set(oneValueArgs TARGET DEBUG_EXTENSION BINARY_DIR COMPONENT INSTALL_DEST DEFINITION)
    set(multiValueArgs SOURCES LINK_LIBRARIES)
    cmake_parse_arguments(FTOOL "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    if( ${FTOOL_DEFINITION} )
    add_definitions(-D${DEFINITION})
    endif()

    BuildToolBundle(
        TARGET ${FTOOL_TARGET}
        SOURCES ${FTOOL_SOURCES}
        DEBUG_EXTENSION ${FTOOL_DEBUG_EXTENSION}
        VERSION_MAJOR ${FTOOL_CTEMsoft_VER_MAJOR}
        VERSION_MINOR ${FTOOL_CTEMsoft_VER_MINOR}
        VERSION_PATCH ${FTOOL_CTEMsoft_VER_PATCH}
        BINARY_DIR    ${FTOOL_BINARY_DIR}
        LINK_LIBRARIES "${FTOOL_LINK_LIBRARIES}"
        LIB_SEARCH_DIRS ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
        COMPONENT     Applications
        INSTALL_DEST  "${FTOOL_INSTALL_DEST}"
        SOLUTION_FOLDER "Tools"
    )

endfunction()


#-----------------------------------------------------------------------------------------
# Set the install directory for all the tools apps that are getting created
set(install_dir "bin")
	set(lib_install_dir "lib")
	if(WIN32)
		set(install_dir ".")
		set(lib_install_dir "")
	endif()



#-----------------------------------------------------------------------------------------
# List all the executables here that only have a SINGLE source file
set(TOOLS
  CTEMmbcbed
  CTEMSRdefect #MP 
  CTEMZAdefect 
  CTEMECP
  CTEMECCI #read_stacking_fault_data has wrong number of args passed to it
  CTEMKossel
  CTEMMC  #MP
  CTEMEBSD
  CTEMEBSDmaster
  CTEMmkxtal
  CTEMlistSG
  CTEMqg
  CTEMfamily
  CTEMstar
  CTEMorbit
  CTEMZAgeom
  CTEMlatgeom
  CTEMstereo
  CTEMorient
  CTEMxtalinfo
  CTEMzap
  CTEMdrawcell
)



#-----------------------------------------------------------------------------------------
# Now loop over each of the tool executables and create targets
FOREACH(TOOL ${TOOLS} )

  COMPILE_TOOL(
	TARGET ${TOOL}
	SOURCES ${Apps_SOURCE_DIR}/${TOOL}.f90
	DEBUG_EXTENSION ${EXE_DEBUG_EXTENSION}
    VERSION_MAJOR ${CTEMsoft_VER_MAJOR}
    VERSION_MINOR ${CTEMsoft_VER_MINOR}
    VERSION_PATCH ${CTEMsoft_VER_PATCH}
    BINARY_DIR    ${${PROJECT_NAME}_BINARY_DIR}
    COMPONENT     Applications
    INSTALL_DEST  "${install_dir}"
    LINK_LIBRARIES CTEMSoftLib
	)


ENDFOREACH()



#  CTEMBWshow)
#  CTEMECPSM)
#  CTEMMCTKD) 
##  CTEMSRBW) #Explicit interface required for 'extractbwdata' 
#  CTEMSRCBED)
#  CTEMTBBFDF)
#  CTEMTBBW)
#  CTEMTBSM)
#  CTEMexpeditions)
#  CTEMgamma)
#  CTEMgammadefect)
#  CTEMholz)
##  CTEMlens) #can't open rksuite
##  CTEMoverlap) #doesn't look like dbdiffB has been initialized anywhere
#  CTEMped)
#  CTEMvr)
##  CTEMvrfft)  #need fft? 
