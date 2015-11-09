
#---------------------------------------------------------------------
# Set some variables to shorten up the call to the function below
set(APP_DIR ${EMsoft_SOURCE_DIR}/templatefolder)

#---------------------------------------------------------------------
# Aggregate all the OpenCL files that are needed
set(EMSoft_RESOURCE_FILES
  ${APP_DIR}/BetheParameters.template
)


#---------------------------------------------------------------------
# Create the Installation Rules
INSTALL(FILES ${EMSoft_RESOURCE_FILES}
  COMPONENT Applications
  DESTINATION "templatefolder"
)
