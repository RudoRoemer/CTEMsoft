
#---------------------------------------------------------------------
# Set some variables to shorten up the call to the function below
set(APP_DIR ${EMsoft_SOURCE_DIR}/templatefolder)

#---------------------------------------------------------------------
# Aggregate all the OpenCL files that are needed
set(EMSoft_RESOURCE_FILES
  ${APP_DIR}/BetheParameters.template
  ${APP_DIR}/EMEBSD.template
  ${APP_DIR}/EMEBSDmaster.template
  ${APP_DIR}/EMECP.template
  ${APP_DIR}/EMECPmaster.template
  ${APP_DIR}/EMKosselmaster.template
  ${APP_DIR}/EMMCOpenCL.template
  ${APP_DIR}/EMsampleRFZ.template
)

file(COPY "${EMsoft_SOURCE_DIR}/templatefolder" DESTINATION "${PROJECT_BINARY_DIR}/")


#---------------------------------------------------------------------
# Create the Installation Rules
INSTALL(FILES ${EMSoft_RESOURCE_FILES}
  COMPONENT Applications
  DESTINATION "templatefolder"
)
