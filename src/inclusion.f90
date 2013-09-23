! ###################################################################
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

!--------------------------------------------------------------------------
! CTEMsoft2013:inclusion.f90
!--------------------------------------------------------------------------
!
! MODULE: inclusion
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief Provides routines to deal with spherical inclusions
! 
!> @date   04/29/11 MDG 1.0 original
!> @date   06/04/13 MDG 2.0 rewrite + quaternions instead of rotations
!--------------------------------------------------------------------------
module inclusion

use local

type inclusiontype
	real(kind=sgl)       ::  xpos,ypos,zpos,radius,C
end type inclusiontype

type (inclusiontype), allocatable  :: inclusions(:)

contains

!--------------------------------------------------------------------------
!
! SUBROUTINE: read_inclusion_data
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief  read inclusion parameters from file
! 
!> @param numinc number of inclusions
!> @param incname name of inclusion file
!> @param DF_L column edge length 
!> @param DF_npix number of x-pixels
!> @param DF_npiy number of y-pixels
!> @param dinfo logical to trigger verbose output
! 
!> @date 1/5/99   MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   03/25/13 MDG 3.0 updated IO
!--------------------------------------------------------------------------
subroutine read_inclusion_data(numinc,incname,DF_L,DF_npix,DF_npiy,dinfo)

use local
use io
use files
use foilmodule

IMPLICIT NONE

character(50),INTENT(IN)       	:: incname
integer(kind=irg),INTENT(OUT):: numinc
integer(kind=irg),INTENT(IN)	:: dinfo,DF_npix,DF_npiy
real(kind=sgl),INTENT(IN)	:: DF_L

integer(kind=irg) :: i, io_int(1)
real(kind=sgl)      :: Vx,Vy,Vz,Vrad,C,tmp(3)

! open the inclusion data file
mess = 'Opening '//incname; call Message("(A)")
open(unit=dataunit,file=incname,form='formatted')
read(unit=dataunit,*) numinc
allocate(inclusions(numinc))
if (dinfo.eq.1) then
  io_int(1) = numinc
  call WriteValue(' Number of inclusions',io_int, 1, "(I)")
end if


! read each subsequent line 
do i=1,numinc
  read(unit=dataunit,*) Vx,Vy,Vz,Vrad,C
  inclusions(i)%xpos = Vx * 0.5 * float(DF_npix)*DF_L
  inclusions(i)%ypos = Vy * 0.5 * float(DF_npiy)*DF_L
  inclusions(i)%zpos = Vz * foil%z0         ! vertical fractional location in interval [-1,1]
  inclusions(i)%radius = Vrad    ! radius in nanometers
  inclusions(i)%C = C                 ! this is the parameter defined in equation (8.36) of the CTEM book.
  tmp = quat_rotate_vector( foil%a_fc, dble((/ inclusions(i)%xpos, inclusions(i)%ypos, inclusions(i)%zpos /)) )  
  inclusions(i)%xpos = tmp(1)
  inclusions(i)%ypos = tmp(2)
  inclusions(i)%zpos = tmp(3)
  if (dinfo.eq.1) write (*,*) i,Vx,Vy,Vz,Vrad,C
end do

! close datafile
close(unit=dataunit,status='keep')
end subroutine read_inclusion_data

end module inclusion
