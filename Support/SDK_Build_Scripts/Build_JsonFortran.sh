#!/bin/bash
# This script requires 2 arguments. The root of the EMSoft_SDK (/Users/Shared/EMSoft_SDK
# or /opt/EMSoft_SDK) and the number of parallel processes to use to compile. This
# is typically 2x the number of physical cores in the machine.

SDK_INSTALL=$1

cd $SDK_INSTALL

PARALLEL_BUILD=$2


HOST_SYSTEM=`uname`
echo "Host System: $HOST_SYSTEM"

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


CMAKE=`type -P cmake`
if [[ $CMAKE == "" ]];
  then
  echo "CMake is needed for this script. Please install it on your system and be sure it is on your path."
  exit 1
fi

version="4.2.0"
# Build the jsonfortran libraries we need and set our Environment Variable.
jsonfortranArchiveName="json-fortran"

if [ ! -e "$SDK_INSTALL/${jsonfortranArchiveName}.tar.gz" ];
then
  echo "-------------------------------------------"
  echo " Downloading jsonfortran Version ${version}"
  echo "-------------------------------------------"
  $DOWNLOAD_PROG  "" -o ${jsonfortranArchiveName}.tar.gz
fi

if [ ! -e "$SDK_INSTALL/${jsonfortranArchiveName}" ];
then
  tar -xvzf ${jsonfortranArchiveName}.tar.gz
# mv jsonfortran-1.8.15 jsonfortran-1.8.15_source
fi

# We assume we already have downloaded the source for json-fortran and have it in a folder
# called json-fortran
cd ${jsonfortranArchiveName}
mkdir Build
cd Build
cmake -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Release -DSKIP_DOC_GEN=TRUE -DCMAKE_INSTALL_PREFIX=${SDK_INSTALL} ../
make -j${PARALLEL_BUILD}
make install
cd ../

#------------------------------------------------------------------------------
# This next bit of code sets the install name of the dylib to the full absolute
# path of the library. This will come in handy when packagin EMSoft with CMake
# by allowing CMake to more easily find the library and adjust its internal paths
cd ${SDK_INSTALL}/jsonfortran-gnu-${version}/lib
install_name_tool -id ${SDK_INSTALL}/jsonfortran-gnu-${version}/lib/libjsonfortran.4.2.dylib ${SDK_INSTALL}/jsonfortran-gnu-${version}/lib/libjsonfortran.4.2.dylib 

echo "#--------------------------------------------------------------------------------------------------" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "# jsonfortran Library" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "set(JSONFORTRAN_INSTALL \"\${EMSoft_SDK_ROOT}/jsonfortran-gnu-${version}\")" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "set(JSONFORTRAN_DIR \"\${JSONFORTRAN_INSTALL}/cmake\")" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "set(jsonfortran-gnu_DIR \"\${JSONFORTRAN_DIR}\")" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "" >> "$SDK_INSTALL/EMSoft_SDK.cmake"

