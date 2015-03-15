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
! CTEMsoft2013:void.f90
!--------------------------------------------------------------------------
!
! MODULE: void
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Provides routines to deal with spherical voids
! 
!> @date   04/29/11 MDG 1.0 original
!> @date   06/04/13 MDG 2.0 rewrite + quaternions instead of rotations
!--------------------------------------------------------------------------
module void

use local
use typedefs

contains

!--------------------------------------------------------------------------
!
! SUBROUTINE: read_void_data
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  read void parameters from file
! 
!> @param defects defects structure
!> @param foil foil structure
!> @param DF_L column edge length 
!> @param DF_npix number of x-pixels
!> @param DF_npiy number of y-pixels
!> @param dinfo logical to trigger verbose output
! 
!> @date  01/05/99 MDG 1.0 original
!> @date  05/19/01 MDG 2.0 f90 version
!> @date  11/27/01 MDG 2.1 added kind support
!> @date  03/25/13 MDG 3.0 updated IO
!> @date  06/09/14 MDG 4.0 added defects argument
!> @date  06/10/14 MDG 4.1 added foil argument
!--------------------------------------------------------------------------
subroutine read_void_data(defects,foil,DF_L,DF_npix,DF_npiy,dinfo)

use io
use files
use foilmodule
use quaternions

IMPLICIT NONE

type(defecttype),INTENT(INOUT)      :: defects
type(foiltype),INTENT(INOUT)        :: foil
integer(kind=irg),INTENT(IN)	      :: dinfo,DF_npix,DF_npiy
real(kind=sgl),INTENT(IN)	      :: DF_L

integer(kind=irg) 			:: i, io_int(1)
real(kind=sgl)      			:: Vx,Vy,Vz,Vrad,tmp(3)

! open the void data file
call Message('Opening '//trim(defects%voidname), frm = "(A)")
open(unit=dataunit,file=trim(defects%voidname),form='formatted')
read(unit=dataunit) defects%numvoids

allocate(defects%voids(defects%numvoids))

if (dinfo.eq.1) then
  io_int(1) = defects%numvoids
  call WriteValue(' Number of voids',io_int, 1, "(I)")
end if

! read each subsequent line 
do i=1,defects%numvoids
  read(unit=dataunit) Vx,Vy,Vz,Vrad
  defects%voids(i)%xpos = Vx * 0.5 * float(DF_npix) * DF_L
  defects%voids(i)%ypos = Vy * 0.5 * float(DF_npiy) * DF_L
  defects%voids(i)%zpos = Vz * foil%z0
  defects%voids(i)%radius = Vrad    ! radius in nanometers
! transform to the foil reference frame  
  tmp = quat_Lp( conjg(foil%a_fc), dble((/ defects%voids(i)%xpos, defects%voids(i)%ypos, defects%voids(i)%zpos /)) )  
  defects%voids(i)%xpos = tmp(1)
  defects%voids(i)%ypos = tmp(2)
  defects%voids(i)%zpos = tmp(3)
  if (dinfo.eq.1) write (*,*) i,Vx,Vy,Vz,Vrad
end do

! close datafile
close(unit=dataunit,status='keep')
end subroutine read_void_data


end module void
