! ###################################################################
! Copyright (c) 2013-2014, Marc De Graef/Carnegie Mellon University
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

!--------------------------------------------------------------------------
! CTEMsoft2013:NameListTypedefs.f90
!--------------------------------------------------------------------------
!
! PROGRAM: NameListTypedefs
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief collection of namelist type declarations
!
!> @date 06/13/14 MDG 1.0 initial version
!--------------------------------------------------------------------------
module NameListTypedefs

use local

IMPLICIT NONE

! namelist for the CTEMKossel program
type KosselNameListType
        integer(kind=irg)       :: stdout
        integer(kind=irg)       :: numthick
        integer(kind=irg)       :: npix
        integer(kind=irg)       :: maxHOLZ
        integer(kind=irg)       :: nthreads
        integer(kind=irg)       :: k(3)
        integer(kind=irg)       :: fn(3)
        real(kind=sgl)          :: voltage
        real(kind=sgl)          :: dmin
        real(kind=sgl)          :: convergence
        real(kind=sgl)          :: startthick
        real(kind=sgl)          :: thickinc
        real(kind=sgl)          :: minten
        character(fnlen)        :: xtalname
        character(fnlen)        :: outname
end type KosselNameListType

! namelist for the CTEMMC program
type MCNameListType
        integer(kind=irg)       :: stdout
        integer(kind=irg)       :: numsx
        integer(kind=irg)       :: primeseed
        integer(kind=irg)       :: num_el
        integer(kind=irg)       :: nthreads
        real(kind=dbl)          :: sig
        real(kind=dbl)          :: omega
        real(kind=dbl)          :: EkeV
        real(kind=dbl)          :: Ehistmin
        real(kind=dbl)          :: Ebinsize
        real(kind=dbl)          :: depthmax
        real(kind=dbl)          :: depthstep
        character(4)            :: MCmode
        character(fnlen)        :: xtalname
        character(fnlen)        :: dataname
end type MCNameListType

! namelist for the CTEMEBSDmaster program
type EBSDMasterNameListType
        integer(kind=irg)       :: stdout
        integer(kind=irg)       :: npx
        integer(kind=irg)       :: Esel
        integer(kind=irg)       :: nthreads
        real(kind=sgl)          :: dmin
        character(fnlen)        :: energyfile
        character(fnlen)        :: outname
end type EBSDMasterNameListType

! namelist for the CTEMEBSD program
! note that not all of these are actually entered via a namelist file
! some of them are used to facilitate passing of subroutine arguments in EBSDmod.f90
type EBSDNameListType
        integer(kind=irg)       :: stdout
        integer(kind=irg)       :: numsx
        integer(kind=irg)       :: numsy
        integer(kind=irg)       :: binning
        integer(kind=irg)       :: nthreads
	real(kind=sgl)          :: L
	real(kind=sgl)          :: thetac
	real(kind=sgl)          :: delta
	real(kind=sgl)          :: xpc
	real(kind=sgl)          :: ypc
	real(kind=sgl)          :: energymin
	real(kind=sgl)          :: energymax
	real(kind=sgl)          :: gammavalue
	real(kind=sgl)          :: axisangle(4)
	real(kind=dbl)          :: beamcurrent
	real(kind=dbl)          :: dwelltime
        character(3)            :: scalingmode
        character(3)            :: eulerconvention
        character(fnlen)        :: anglefile
        character(fnlen)        :: masterfile
        character(fnlen)        :: energyfile 
        character(fnlen)        :: datafile
! everything below here is not part of the namelist input structure, but is used to pass arguments to subroutines
        integer(kind=irg)       :: numangles
        integer(kind=irg)       :: numEbins
        integer(kind=irg)       :: numzbins 
        integer(kind=irg)       :: nsx
        integer(kind=irg)       :: nsy
        integer(kind=irg)       :: num_el
        integer(kind=irg)       :: MCnthreads
        integer(kind=irg)       :: npx
        integer(kind=irg)       :: npy
        integer(kind=irg)       :: nE
        integer(kind=irg)       :: numset
        real(kind=dbl)          :: EkeV
        real(kind=dbl)          :: Ehistmin 
        real(kind=dbl)          :: Ebinsize 
        real(kind=dbl)          :: depthmax
        real(kind=dbl)          :: depthstep
        real(kind=dbl)          :: MCsig
        real(kind=dbl)          :: MComega
        character(4)            :: MCmode       ! Monte Carlo mode
        character(5)            :: anglemode    ! 'quats' or 'euler' for angular input
        character(6)            :: sqorhe       ! from Master file, square or hexagonal Lambert projection
        character(8)            :: MCscversion
        character(8)            :: Masterscversion
        character(fnlen)        :: Masterprogname
        character(fnlen)        :: Masterxtalname 
        character(fnlen)        :: Masterenergyfile
        character(fnlen)        :: MCprogname 
        character(fnlen)        :: MCxtalname
  
end type EBSDNameListType


end module NameListTypedefs