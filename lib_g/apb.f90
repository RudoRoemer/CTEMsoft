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
! CTEMsoft2013:apb.f90
!--------------------------------------------------------------------------
!
! MODULE: inclusion
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Provides routines to deal with cylindrical APBs
! 
!> @date   02/10/14 MDG 1.0 original
!--------------------------------------------------------------------------
module apb

use local

type apbtype
	real(kind=sgl)       ::  xpos,ypos,zpos,radius,w,Rdisp(3)
end type apbtype

type (apbtype), allocatable  :: apbs(:)

contains

!--------------------------------------------------------------------------
!
! SUBROUTINE: read_apb_data
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  read apb parameters from file
! 
!> @param numapb number of inclusions
!> @param apbname name of inclusion file
!> @param DF_L column edge length 
!> @param DF_npix number of x-pixels
!> @param DF_npiy number of y-pixels
!> @param dinfo logical to trigger verbose output
! 
!> @date   02/10/14 MDG 1.0 new code
!--------------------------------------------------------------------------
subroutine read_apb_data(numapb,apbname,DF_L,DF_npix,DF_npiy,dinfo)

use local
use io
use files
use crystal
use foilmodule

IMPLICIT NONE

character(fnlen),INTENT(IN)    :: apbname
integer(kind=irg),INTENT(OUT)  :: numapb
integer(kind=irg),INTENT(IN)   :: dinfo,DF_npix,DF_npiy
real(kind=sgl),INTENT(IN)      :: DF_L

integer(kind=irg)              :: i, io_int(1)
real(kind=sgl)                 :: Vx,Vy,Vz,Vrad,w,tmp(3), tmp2(3), Rx, Ry, Rz

! open the apbdata file
mess = 'Opening '//apbname; call Message("(A)")
open(unit=dataunit,file=apbname,form='formatted')
read(dataunit,*) numapb ! PGC unit=dataunit -> dataunit
allocate(apbs(numapb))
if (dinfo.eq.1) then
  io_int(1) = numapb
  call WriteValue(' Number of APBs ',io_int, 1, "(I)")
end if


! read each subsequent line 
do i=1,numapb
  read(dataunit,*) Vx,Vy,Vz,Vrad,w,Rx,Ry,Rz ! PGC unit=dataunit -> dataunit
  apbs(i)%xpos = Vx * 0.5 * float(DF_npix)*DF_L
  apbs(i)%ypos = Vy * 0.5 * float(DF_npiy)*DF_L
  apbs(i)%zpos = Vz * foil%z0         ! vertical fractional location in interval [-1,1]
  apbs(i)%radius = Vrad               ! radius in nanometers
  apbs(i)%w = w   
  tmp = (/ Rx, Ry, Rz /)
  call TransSpace(tmp,tmp2,'d','c') 
  apbs(i)%Rdisp = tmp2 

  tmp = quat_rotate_vector( foil%a_fc, dble((/ apbs(i)%xpos, apbs(i)%ypos, apbs(i)%zpos /)) )  
  apbs(i)%xpos = tmp(1)
  apbs(i)%ypos = tmp(2)
  apbs(i)%zpos = tmp(3)
  write (*,*) apbs(i)
end do

! close datafile
close(unit=dataunit,status='keep')
end subroutine read_apb_data

end module apb
