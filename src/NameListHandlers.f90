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
! CTEMsoft2013:NameListHandlers.f90
!--------------------------------------------------------------------------
!
! PROGRAM: NameListHandlers
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief routines for reading and returning name list type structures
!
!> @date 06/13/14 MDG 1.0 original
!--------------------------------------------------------------------------
module NameListHandlers

use local
use NameListTypedefs

contains

!--------------------------------------------------------------------------
!
! SUBROUTINE:GetKosselNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read namelist file and fill knl structure (used by CTEMKossel.f90)
!
!> @param nmlfile namelist file name
!> @param knl Kossel name list structure
!
!> @date 06/13/14  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine GetKosselNameList(nmlfile, knl)

use error

IMPLICIT NONE

character(fnlen),INTENT(IN)             :: nmlfile
type(KosselNameListType),INTENT(INOUT)  :: knl

integer(kind=irg)                       :: stdout, numthick, npix, maxHOLZ, nthreads, k(3), fn(3)
real(kind=sgl)                          :: voltage, dmin, convergence, startthick, thickinc, minten
character(fnlen)                        :: xtalname, outname

namelist /Kossellist/ stdout, xtalname, voltage, k, fn, dmin, convergence, minten, nthreads, &
                              startthick, thickinc, numthick, outname, npix, maxHOLZ

! set the input parameters to default values (except for xtalname, which must be present)
stdout = 6                      ! standard output
numthick = 10                   ! number of increments
npix = 256                      ! output arrays will have size npix x npix
maxHOLZ = 3                     ! output arrays will have size npix x npix
nthreads = 4                    ! default number of threads for OpenMP
k = (/ 0, 0, 1 /)               ! beam direction [direction indices]
fn = (/ 0, 0, 1 /)              ! foil normal [direction indices]
voltage = 200000.0              ! acceleration voltage [V]
dmin = 0.025                    ! smallest d-spacing to include in dynamical matrix [nm]
convergence = 25.0              ! beam convergence angle [mrad]
startthick = 10.0               ! starting thickness [nm]
thickinc = 10.0                 ! thickness increment
minten = 1.0E-5                 ! minimum intensity in diffraction disk to make it into the output file
xtalname = 'undefined'          ! initial value to check that the keyword is present in the nml file
outname = 'Kosselout.data'      ! output filename

! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=Kossellist)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(xtalname).eq.'undefined') then
  call FatalError('CTEMKossel:',' structure file name is undefined in '//nmlfile)
 end if

! if we get here, then all appears to be ok, and we need to fill in the knl fields
knl%stdout = stdout
knl%numthick = numthick
knl%npix = npix
knl%maxHOLZ = maxHOLZ
knl%nthreads = nthreads
knl%k = k
knl%fn = fn
knl%voltage = voltage
knl%dmin = dmin
knl%convergence = convergence
knl%startthick = startthick
knl%thickinc = thickinc
knl%minten = minten
knl%xtalname = xtalname
knl%outname = outname

end subroutine GetKosselNameList

end module NameListHandlers
