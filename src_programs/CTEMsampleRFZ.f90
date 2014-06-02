! ###################################################################
! Copyright (c) 2014, Marc De Graef/Carnegie Mellon University
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
!
! PROGRAM: CTEMsampleRFZ
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Basic program to generate a uniform sampling of Rodrigues Fundamental Zone
!
!> @details This program calls the SampleRFZ routine of the so3 module to generate
!> an angle file of euler angles for points that uniformly sample an RFZ for a given
!> crystal symmetry.  
!
!> @date 5/12/14   MDG 1.0 original
!> @date 5/29/14   MDG 1.1 integrated with CTEMsoft package (started from standalone program)
!--------------------------------------------------------------------------
program CTEMsampleRFZ

use local
use files
use io

IMPLICIT NONE

character(fnlen)	:: nmldeffile

! deal with the command line arguments, if any
nmldeffile = 'CTEMsampleRFZ.nml'
progname = 'CTEMsampleRFZ.f90'
call Interpret_Program_Arguments(nmldeffile,1,(/ 60 /) )

! perform the zone axis computations
call CreateSampling(nmldeffile)

end program CTEMsampleRFZ

!--------------------------------------------------------------------------
!
! SUBROUTINE:CreateSampling
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Generate a sampling of the Rodrigues Fundamental Zone for a given xtal symmetry
!
!> @param nmlfile namelist file name
!
!> @date 05/29/14  MDG 1.0 original
!--------------------------------------------------------------------------
subroutine CreateSampling(nmlfile)

use local
use constants
use rotations
use so3

IMPLICIT NONE

character(fnlen),INTENT(IN)::nmlfile

integer(kind=irg)	    :: i, pgnum, nsteps
real(kind=dbl)             :: eud(3), rtod
character(fnlen)           :: outname

namelist /RFZlist/ pgnum, nsteps, outname

rtod = 180.D0/cPi

pgnum = 1
nsteps = 10
outname = 'anglefile.txt'

write (*,*) ro2eu( (/ 0.0, 0.0, 0.0 /) )

! read the namelist file
open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
read(UNIT=dataunit,NML=RFZlist)
close(UNIT=dataunit,STATUS='keep')

write (*,*) pgnum, nsteps, outname

! get the linked list for the FZ for point group symmetry pgnum for nsteps along the cubic semi-edge
call SampleRFZ(nsteps,pgnum)

! now we have the linked list so we can do anything we want with it.
! in this test program, we create a VTK file so that we can visualize the RFZ with ParaView
open (UNIT=22,FILE=trim(outname),FORM='formatted',STATUS='unknown')
write (22,"(A)") 'eu'
write (22,"(I8)") FZcnt

! scan through the list
FZtmp => FZlist
do i = 1, FZcnt
! convert each rodrigues triplet into euler angles and write them to the file
  eud = ro2eu(FZtmp%rod) * rtod
  write (22,"(3F16.4)") eud
  FZtmp => FZtmp%next
end do

close(UNIT=22,STATUS='keep')

end subroutine CreateSampling
