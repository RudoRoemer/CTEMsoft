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

version="0.1alpha4"
# Build the HDF5 libraries we need and set our Environment Variable.
fortranClArchiveName="fortrancl-${version}"

if [ ! -e "$SDK_INSTALL/${fortranClArchiveName}.tar.gz" ];
then
  echo "-------------------------------------------"
  echo " Downloading HDF5 Version ${version}"
  echo "-------------------------------------------"
  #$DOWNLOAD_PROG  "http://www.hdfgroup.org/ftp/HDF5/current/src/${fortranClArchiveName}.tar.gz" -o ${fortranClArchiveName}.tar.gz
fi

if [ ! -e "$SDK_INSTALL/${fortranClArchiveName}" ];
then
  tar -xvzf ${fortranClArchiveName}.tar.gz
fi
# We assume we already have downloaded the source for FortranCL and have it extracted
cd $SDK_INSTALL/${fortranClArchiveName}

./configure --prefix=$SDK_INSTALL/fortrancl
make
make install


echo "#--------------------------------------------------------------------------------------------------" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "# FORTRAN CL Library" >> "$SDK_INSTALL/EMSoft_SDK.cmake"
echo "set(FORTRANCL_INSTALL \"\${EMSoft_SDK_ROOT}/fortrancl\")" >> "$SDK_INSTALL/EMSoft_SDK.cmake"


