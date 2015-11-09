
#---------------------------------------------------------------------
# Set some variables to shorten up the call to the function below
set(APP_DIR ${EMsoft_SOURCE_DIR}/examples_todo/DefectSimulations)

#---------------------------------------------------------------------
# Aggregate all the OpenCL files that are needed
set(EMSoft_DefectSim_FILES
  ${APP_DIR}/CTEMZAdefect.template
  ${APP_DIR}/dislocation.template
  ${APP_DIR}/FOIL_rundata.template
  ${APP_DIR}/inclusion.template
  ${APP_DIR}/stackingfault.template
  ${APP_DIR}/STEM_rundata.template
  ${APP_DIR}/void.template
)

#---------------------------------------------------------------------
# Create the Installation Rules
INSTALL(FILES ${EMSoft_DefectSim_FILES}
  COMPONENT Applications
  DESTINATION "examples/DefectSimulations"
)
