
#---------------------------------------------------------------------
# Set some variables to shorten up the call to the function below
set(APP_DIR ${EMsoft_SOURCE_DIR}/resources)

#---------------------------------------------------------------------
# Aggregate all the OpenCL files that are needed
set(EMSoft_CL_SRCS
  ${APP_DIR}/rotations.txt
  ${APP_DIR}/templatecodes.txt
#  ${APP_DIR}/RandomSeed.data
)

#---------------------------------------------------------------------
# Create the Installation Rules
INSTALL(FILES ${EMSoft_CL_SRCS}
  COMPONENT Applications
  DESTINATION "resources"
)
