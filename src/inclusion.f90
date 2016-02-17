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
! EMsoft:inclusion.f90
!--------------------------------------------------------------------------
!
! MODULE: inclusion
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Provides routines to deal with spherical inclusions
! 
!> @date   04/29/11 MDG 1.0 original
!> @date   06/04/13 MDG 2.0 rewrite + quaternions instead of rotations
!--------------------------------------------------------------------------
module inclusion

use local
use typedefs

contains

!--------------------------------------------------------------------------
!
! SUBROUTINE: read_inclusion_data
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  read inclusion parameters from file
! 
!> @param defects defect structure
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
subroutine read_inclusion_data(defects,foil,DF_L,DF_npix,DF_npiy,dinfo)

use io
use files
use foilmodule
use quaternions

IMPLICIT NONE

type(defecttype),INTENT(INOUT) :: defects
type(foiltype),INTENT(INOUT)   :: foil
integer(kind=irg),INTENT(IN)	:: dinfo,DF_npix,DF_npiy
real(kind=sgl),INTENT(IN)	:: DF_L

integer(kind=irg) :: i, io_int(1)
real(kind=sgl)      :: Vx,Vy,Vz,Vrad,C,tmp(3)

! open the inclusion data file
call Message('Opening '//trim(defects%incname), frm = "(A)")
open(unit=dataunit,file=trim(EMsoft_toNativePath(defects%incname)),form='formatted')
read(unit=dataunit) defects%numinc

allocate(defects%inclusions(defects%numinc))

if (dinfo.eq.1) then
  io_int(1) = defects%numinc
  call WriteValue(' Number of inclusions',io_int, 1, "(I)")
end if


! read each subsequent line 
do i=1,defects%numinc
  read(unit=dataunit) Vx,Vy,Vz,Vrad,C
  defects%inclusions(i)%xpos = Vx * 0.5 * float(DF_npix)*DF_L
  defects%inclusions(i)%ypos = Vy * 0.5 * float(DF_npiy)*DF_L
  defects%inclusions(i)%zpos = Vz * foil%z0         ! vertical fractional location in interval [-1,1]
  defects%inclusions(i)%radius = Vrad    ! radius in nanometers
  defects%inclusions(i)%C = C                 ! this is the parameter defined in equation (8.36) of the EM book.
  tmp = quat_Lp( conjg(foil%a_fc), dble((/ defects%inclusions(i)%xpos, defects%inclusions(i)%ypos, &
        defects%inclusions(i)%zpos /)) )  
  defects%inclusions(i)%xpos = tmp(1)
  defects%inclusions(i)%ypos = tmp(2)
  defects%inclusions(i)%zpos = tmp(3)
  if (dinfo.eq.1) write (*,*) i,Vx,Vy,Vz,Vrad,C
end do

! close datafile
close(unit=dataunit,status='keep')
end subroutine read_inclusion_data

end module inclusion
