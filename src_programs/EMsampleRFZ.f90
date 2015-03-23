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
! PROGRAM: EMsampleRFZ
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
!> @date 5/29/14   MDG 1.1 integrated with EMsoft package (started from standalone program)
!--------------------------------------------------------------------------
program EMsampleRFZ

use local
use typedefs
use NameListTypedefs
use NameListHandlers
use files
use io

IMPLICIT NONE

character(fnlen)                :: nmldeffile, progname, progdesc
type(RFZNameListType)           :: rfznl

! deal with the command line arguments, if any
nmldeffile = 'EMsampleRFZ.nml'
progname = 'EMsampleRFZ.f90'
progdesc = 'Create a uniform sampling of Rodrigues space and output Euler angles list'
call Interpret_Program_Arguments(nmldeffile,1,(/ 60 /), progname )

! deal with the namelist stuff
call GetRFZNameList(nmldeffile,rfznl)

! print some information
call EMsoft(progname, progdesc)

! perform the zone axis computations
call CreateSampling(rfznl,progname)

end program EMsampleRFZ

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
!> @date 05/29/14 MDG 1.0 original
!> @date 12/09/14 MDG 2.0 changed rfznl handling
!--------------------------------------------------------------------------
subroutine CreateSampling(rfznl, progname)

use local
use typedefs
use NameListTypedefs
use constants
use rotations
use io
use so3

IMPLICIT NONE

type(RFZNameListType),INTENT(IN)        :: rfznl
character(fnlen),INTENT(IN)             :: progname

integer(kind=irg)                       :: i, FZcnt, io_int(1)
real(kind=dbl)                          :: eud(3), rtod
type(FZpointd),pointer                  :: FZlist, FZtmp

rtod = 180.D0/cPi

! a bit of output
call Message('Starting computation for point group '//PGTHD(rfznl%pgnum))

! get the linked list for the FZ for point group symmetry pgnum for nsteps along the cubic semi-edge
nullify(FZlist)
FZcnt = 0
call SampleRFZ(rfznl%nsteps,rfznl%pgnum,FZcnt,FZlist)

io_int(1) = FZcnt
call WriteValue('Total number of unique orientations generated = ',io_int,1,"(I10)")

! generate a list of all orientations in Euler angle format
open (UNIT=22,FILE=trim(rfznl%outname),FORM='formatted',STATUS='unknown')
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

call Message('Euler angles stored in file '//rfznl%outname)


end subroutine CreateSampling
