# EMSoft SDK Setup #

These instructions have been tested on OS X  10.10.x (Yosemite) but should also be generally applicable to LINUX systems as well. Separate instructions will be made available for Windows installations.

## Preliminaries ##

EMSoft relies on the following component libraries:


| Package | Version | Notes |
|---------|---------|-------|
| GFortran | 4.9.0 or Greater | [https://gcc.gnu.org/wiki/GFortranBinaries] (https://gcc.gnu.org/wiki/GFortranBinaries) |
| JSON Fortran | 4.2.0 | [https://github.com/jacobwilliams/json-fortran] (https://github.com/jacobwilliams/json-fortran)|
| HDF5 1.8.15 | 1.8.15 | Must build the Fortran libraries as well as C/C++ libraries |
| FortranCL | fortrancl-0.1alpha4 | [https://code.google.com/p/fortrancl/downloads/detail?name=fortrancl-0.1alpha4.tar.gz&can=2&q=](https://code.google.com/p/fortrancl/downloads/detail?name=fortrancl-0.1alpha4.tar.gz&can=2&q=) |

## Install FORTRAN Compiler ##

Install gfortran or Intel Fortran (untested as of this writing) on your system. Go to the web site and download/install the appropriate version of gfortran for your OS version.  Make sure that you have at least version 4.9; EMSoft will not compile with earlier gfortran versions.

## EMSoft_SDK Location ##

On OS X in particular we like to *sandbox* our installation of the above libraries in order to not affect any other software packages on the system. This is optional and can be changed if the user has the knowledge of how to do so. For this example we will locate our **sandbox** at */Users/Shared/EMSoft_SDK*.

## Build SDK ##

There is a shell script in EMSoft/Support/OSX_Build_Scripts/Build_SDK.sh that can be run by the user to download and build the dependent libraries in such a way that EMSoft will compile. For OS X systems the EMSoft_SDK is coded to be in **/Users/Shared/EMSoft_SDK**. If you want this in a different location then the script file can be adjusted as needed. Simply run the script as **sudo** so that the EMSoft_SDK can be created. Note any errors that occur during the process.

## Build EMSoft ##

We will use CMake to configure a build system for EMSoft. CMake has been included in the EMSoft_SDK for use with EMSoft development. The easiest way to proceed is with a terminal/command prompt on your system of choice. First Navigate to the location of EMSoft.

	[user@system] $ export PATH=$PATH:/Users/Shared/EMSoft_SDK/cmake-3.3.1-Darwin-x86_64/CMake.app/Contents/bin/
	[user@system] $ cd /Path/to/EMSoft
	[user@system] $ mkdir Build
	[user@system] $ cd Build
	[user@system] $ cmake -DEMSoft_SDK=/Users/Shared/EMSoft_SDK -DCMAKE_BUILD_TYPE=Debug ../
	[user@system] $ make -j
	
After compilation the various programs will be available to execute.

# Linux Installation #

These instructions were developed on Ubuntu 14.04 with an nVidea 650 GTX video card.

## GFortran 4.9.2 ##

Ubuntu 14.04 comes with GFortran 4.8 which is not new enough for EMSoft to use. We need to install 
a newer version (Preferably 5.x, but 4.9.2 should work.). The following commands will get you
the newer GFortran library:

  [user] $ sudo add-apt-repository ppa:ubuntu-toolchain-r/test
  [user] $ sudo apt-get update && sudo apt-get upgrade && sudo apt-get install gfortran-4.9

## OpenCL ##
Download from (https://software.intel.com/en-us/articles/opencl-drivers)[https://software.intel.com/en-us/articles/opencl-drivers]
The download we are looking for is:
  Intel OpenCL™ Code Builder for Ubuntu (157MB)



