# EMsoft Version 3.0 #

This package contains a series of programs along with a library, mostly written in fortran-90 with some OpenCL bits, for the computation and visualization of scanning electron microscopy diffraction patterns, notable EBSD, ECP, and EKP. The programs can be used as standalone command-line programs, and produce HDF5 output files that can then be visualized using an IDL virtual app (Interactive Data Language) or read and processed from any other package with HDF5 capability.

## Financial Support ##
Current development of this package is supported by an AFOSR/MURI grant, FA9550-12-1-0458; the original EBSD code from CTEMsoft 2.0 was developed with support from an ONR grant, N00014-12-1-0075.

## Source ##
Source codes are here on GitHub. Precompiled binaries, manuals, and IDL Virtual Machine apps, as well as additional information about the package are located at [Here](http://muri.materials.cmu.edu/?p=858). 

## Installation/Compilation ##
See [here](http://muri.materials.cmu.edu/wp-content/uploads/2015/11/EMsoft.pdf).

## Future Releases ##
Our plan is to provide Release 3.1 in the Spring of 2016, and 3.2 in the Fall.  3.1 will include ECCI defect image simulations as well as dictionary-based indexing codes for EBSD and ECP.  In 3.2 we will re-introduce all of the TEM simulations that were available in the original [CTEMsoft 2.0](https://github.com/marcdegraef/CTEMsoft) release. 
 
## Feedback ##
Feedback on the code and visualization routines is always welcome!

## History ##
This code has a long history. The very first version was written in the mid-1990s in fortran-77 and contained only a few simple crystallography and diffraction routines.  Around 2000, all the code was ported to fortran-90, and served as the basis for a TEM textbook (Introduction to Conventional Transmission Electron Microscopy, M. De Graef, 2003 Cambridge University Press).  Since 2003, many additions have been made.  During 2013, I was on sabbatical at the Ohio State University (Center for Electron Microscopy and MicroAnalysis), and the entire package was rewritten, visualization routines were created, all of which resulted in the present Release; most importantly, we started to add developers (it was pretty much a one-man-party until then), and the source code is now managed by a private git-repository, from which we will periodically push to GitHub.com. 

## License ##

	!###################################################################
	! Copyright (c) 2013-2015, Marc De Graef/Carnegie Mellon University
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

