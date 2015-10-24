#!/bin/bash


HOST_SYSTEM=`uname`
echo "Host System: $HOST_SYSTEM"
PARALLEL_BUILD=8
SCRIPT_DIR=`pwd`


SDK_PARENT_DIR=/opt
SDK_INSTALL=${SDK_PARENT_DIR}/EMSoft_SDK

# Create the actual SDK directory and set the ownership
if [ ! -e "$SDK_INSTALL" ];
then
  sudo mkdir -p ${SDK_INSTALL}
  sudo chmod ugo+rwx ${SDK_INSTALL}
fi

WGET=`type -P wget`
CURL=`type -P curl`

if [[ "$WGET" == "" ]];
   then
  if [[ "$CURL" == "" ]];
     then
    echo "wget and curl are NOT present on your machine. One of them is needed to download sources from the internet."
    exit 1
  fi
fi

DOWNLOAD_PROG=""
DOWNLOAD_ARGS=""

if [[ "$WGET" != "" ]];
then
  DOWNLOAD_PROG=$WGET
fi

if [[ "$CURL" != "" ]];
then
  DOWNLOAD_PROG=$CURL
  DOWNLOAD_ARGS=""
fi

# Only look for Doxygen.app on OS X systems.
if [[ "$HOST_SYSTEM" = "Darwin" ]];
then
  if [ ! -e "/Applications/Doxygen.app" ];
  then
    echo "--------------------------------------------"
    echo "Doxygen is missing from your system."
    echo "Downloading Doxygen 1.8.10 for you."
    $DOWNLOAD_PROG  "http://ftp.stack.nl/pub/users/dimitri/Doxygen-1.8.10.dmg" -o "${EMSoft_SDK}/Doxygen-1.8.10.dmg"
    open "${EMSoft_SDK}/Doxygen-1.8.10.dmg"
    echo "Please Copy the Doxygen.app from the mounted disk image into the /Applications directory. CMake can most"
    echo "easily find it in this location."
  fi
fi

archiveName="EMSoft_SDK_Linux.tar.gz"
cmakename="Linux"
if [[ "$HOST_SYSTEM" = "Darwin" ]];
then
  archiveName="EMSoft_SDK_OSX.tar.gz"
  cmakename="Darwin"
fi

# If we are missing the actual source archives then download from the web site
if [ ! -e "${SDK_PARENT_DIR}/EMSoft_SDK_OSX.tar.gz" ];
  then
  echo "-----------------------------------------------------------"
  echo "An archive named ${archiveName} should be located at "
  echo "${SDK_PARENT_DIR}/${archiveName} but was not found. This "
  echo "archive contains all the dependent library codes that will be"
  echo "compiled for the EMSoft SDK. The SDK Archive will be downloaded"
  echo "from http://dream3d.bluequartz.net"
  $DOWNLOAD_PROG  "http://dream3d.bluequartz.net/binaries/EMSoft_SDK/${archiveName}" -O "${SDK_PARENT_DIR}/${archiveName}"

fi

#-------------------------------------------------
# Move one Directory Above the SDK Folder and untar the
if [ -e "$SDK_PARENT_DIR/${archiveName}" ];
  then
    cd "$SDK_PARENT_DIR"
    tar -xvzf ${archiveName}
fi

#-------------------------------------------------
# Copy our scripts over to the SDK directory
cp ${SCRIPT_DIR}/Build_HDF5.sh ${SDK_INSTALL}/.
cp ${SCRIPT_DIR}/Build_JsonFortran.sh ${SDK_INSTALL}/.
cp ${SCRIPT_DIR}/Build_FortranCL.sh ${SDK_INSTALL}/.

#-------------------------------------------------
# Move into the SDK directory
cd ${SDK_INSTALL}

#-------------------------------------------------
# Unpack CMake
tar -xvzf ${SDK_INSTALL}/cmake-3.3.1-${cmakename}-x86_64.tar.gz

#-------------------------------------------------
# Get CMake on our path
if [[ "$HOST_SYSTEM" = "Darwin" ]];
then
  export PATH=$PATH:${SDK_INSTALL}/cmake-3.3.1-${cmakename}-x86_64/CMake.app/Contents/bin
else
  export PATH=$PATH:${SDK_INSTALL}/cmake-3.3.1-${cmakename}-x86_64/bin
fi

#-------------------------------------------------
# Create the EMSoft_SKD.cmake file, but back up any existing one first
if [ -e "$SDK_INSTALL/EMSoft_SDK.cmake" ]
  then
  mv "$SDK_INSTALL/EMSoft_SDK.cmake" "$SDK_INSTALL/EMSoft_SDK.cmake.bak"
fi
echo "# This is the EMSoft_SDK File. This file contains all the paths to the dependent libraries." > "$SDK_INSTALL/EMSoft_SDK.cmake"

echo "if(NOT DEFINED EMSoft_FIRST_CONFIGURE)" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "  message(STATUS \"*******************************************************\")" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "  message(STATUS \"* EMSoft First Configuration Run                    *\")" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "  message(STATUS \"* EMSoft_SDK Loading from \${CMAKE_CURRENT_LIST_DIR}  *\")" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "  message(STATUS \"*******************************************************\")" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "endif()" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "set(CMAKE_CXX_FLAGS \"-Wmost -Wno-four-char-constants -Wno-unknown-pragmas -mfpmath=sse\" CACHE STRING \"\" FORCE)" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "#--------------------------------------------------------------------------------------------------" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "# These settings are specific to DREAM3D. DREAM3D needs these variables to" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "# configure properly." >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "set(BUILD_TYPE \${CMAKE_BUILD_TYPE})" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "if(\"\${BUILD_TYPE}\" STREQUAL \"\")" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "    set(BUILD_TYPE \"Release\" CACHE STRING \"\" FORCE)" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "endif()" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "message(STATUS \"The Current Build type being used is \${BUILD_TYPE}\")" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "set(BUILD_SHARED_LIBS ON CACHE BOOL \"\")" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "set(EMSoft_SDK_ROOT \"$SDK_INSTALL\")" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "set(EMSoft_DATA_DIR \${EMSoft_SDK_ROOT}/EMSoft_Data CACHE PATH \"\")" >> "$SDK_INSTALL/EMSoft_SDK.cmake"

# Write out the Qt5 directory/installation
echo "#--------------------------------------------------------------------------------------------------" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "# Currently EMSoft does not Depend on Qt5, but if it did, this line is needed." >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "# Qt 5.4.x Library" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "# set(Qt5_DIR \"\${EMSoft_SDK_ROOT}/Qt5.4.2/5.4/clang_64/lib/cmake/Qt5\" CACHE PATH \"\")" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo ""  >> "$SDK_INSTALL/EMSoft_SDK.cmake"

#-------------------------------------------------
# Start building all the packages

${SDK_INSTALL}/Build_JsonFortran.sh "${SDK_INSTALL}" ${PARALLEL_BUILD}
rm ${SDK_INSTALL}/Build_JsonFortran.sh


${SDK_INSTALL}/Build_HDF5.sh "${SDK_INSTALL}" ${PARALLEL_BUILD}
rm ${SDK_INSTALL}/Build_HDF5.sh


${SDK_INSTALL}/Build_FortranCL.sh "${SDK_INSTALL}" ${PARALLEL_BUILD}
rm ${SDK_INSTALL}/Build_FortranCL.sh

# Continue writing the EMSoft_SDK.cmake file after all those libraries were compiled
echo "" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "#--------------------------------------------------------------------------------------------------" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "# Update CMake Module Path with additional paths in order to better find the libraries." >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "set(CMAKE_MODULE_PATH \${CMAKE_MODULE_PATH} \${HDF5_DIR} \${JSONFORTRAN_DIR})" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "#--------------------------------------------------------------------------------------------------" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "# Only Run this the first time when configuring DREAM.3D. After that the values" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "# are cached properly and the user can add additional plugins through the normal" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "# CMake GUI or CCMake programs." >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "if(NOT DEFINED EMSoft_FIRST_CONFIGURE)" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "  set(EMSoft_FIRST_CONFIGURE \"ON\" CACHE STRING \"Determines if DREAM3D has already been configured\")" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "endif()" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "" >> "$SDK_INSTALL/EMSoft_SDK.cmake"


sudo chmod -R ugo+rw *
