# CTEMsoft Version 2.0 #

This package contains a series of programs along with a library, mostly written in fortran-90, for the computation and visualization of electron microscopy images and diffraction patterns. The programs can be used as standalone command-line programs, and produce binary output files that can then be visualized using IDL routines (Interactive Data Language).

## IMPORTANT NOTICE ##
As of Friday November 13, 2015, CTEMsoft 2.0 is no longer supported.  The SEM portion of this package has been completely rewritten and is made available as the new EMsoft 3.0 package here on GitHub.  The TEM portion of CTEMsoft 2.0 will be rewritten during the coming months and will become part of EMsoft 3.2 in the Fall of 2016.

## Financial Support ##
Current development of this package is supported by an AFOSR/MURI grant, FA9550-12-1-0458; the original EBSD code was developed with support from an ONR grant, N00014-12-1-0075

## Source ##
Source codes are here on GitHub. 
## Installation/Compilation ##
The code uses the CMake approach for compilation, and requires the gfortran compiler [tested on Mac and PC; Intel compiler on PC currently being tested].  A CMake run will produce a library and a series of executables.  The package requires the BLAS and LAPACK libraries, so you should install those before compilation.

## Feedback ##
Feedback on the code and visualization routines is always welcome!

## Precompiled Binaries ##
Precompiled binaries, manuals, and IDL Virtual Machine apps, are located at [Here](http://muri.materials.cmu.edu/binaries/CTEMsoft-2.0-OSX.zip); they are known to work on OS X 10.8.9 and 10.9.3

## History ##
This code has a long history. The very first version was written in the mid-1990s in fortran-77 and contained only a few simple crystallography and diffraction routines.  Around 2000, all the code was ported to fortran-90, and served as the basis for a TEM textbook (Introduction to Conventional Transmission Electron Microscopy, M. De Graef, 2003 Cambridge University Press).  Since 2003, many additions have been made.  During 2013, I was on sabbatical at the Ohio State University (Center for Electron Microscopy and MicroAnalysis), and the entire package was rewritten, visualization routines were created, all of which resulted in the present Release; most importantly, we started to add developers (it was pretty much a one-man-party until then), and the source code is now managed by a private git-repository, from which we will periodically push to github.com. 

## License ##

!###################################################################
! Copyright (c) 2013, Marc De Graef/Carnegie Mellon University
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are 
! permitted provided that the following conditions are met:
!
!     - Redistributions of source code must retain the above copyright notice, this list 
!        of conditions and the following disclaimer.
!     - Redistributions in binary form must reproduce the above copyright notice, this 
!        list of conditions and the following disclaimer in the documentation and/or 
!        other materials provided with the distribution.
!     - Neither the names of Marc De Graef, Carnegie Mellon University nor the names 
!        of its contributors may be used to endorse or promote products derived from 
!        this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE 
! USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! ###################################################################

