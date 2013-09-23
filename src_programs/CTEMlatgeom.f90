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
! CTEMsoft2013:CTEMlatgeom.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMlatgeom 
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief general lattice geometry computations
!
!> @details  not meant to be an all encompassing program; simply an illustration of how
!> to do crystallographic computations.
! 
!> @date 1/5/99   MDG 1.0 original
!> @date  5/21/01  MDG 2.0 f90
!> @date 4/16/13 MDG 3.0 rewrite
!--------------------------------------------------------------------------
program CTEMlatgeom

use local
use io
use files
use crystalvars
use crystal
use constants

integer(kind=irg)        	:: isel,another, oi_int(3)
real(kind=sgl)           	:: v1(3),v2(3),vc(3),p,q,r, oi_real(1)
character(1)             	:: sp,sp2

 progname = 'CTEMlatgeom.f90'
 progdesc = 'Simple lattice geometry program'
 call CTEMsoft
 
! load crystal structure data
 call CrystalData

! set up loop
 another=1
 do while (another.eq.1) 
  call PrintMenu(isel)
  sp='d'
  if (mod(isel,2).eq.0) sp='r'
  if (isel.lt.3) then
   call ReadValue('    Enter vector components : ', v1, 3)
   oi_real(1) = CalcLength(v1,sp)
   if (sp.eq.'d') then 
    call WriteValue(' -> Length [nm] = ', oi_real, 1, "(2x,F10.6)")
   else
    call WriteValue(' -> Length [nm-1] = ', oi_real, 1, "(2x,F10.6)")
   end if
  else
   call ReadValue('    Enter first vector components : ', v1, 3)
   call ReadValue('    Enter second vector components : ', v2, 3)
   if (isel.lt.5) then 
    p = CalcLength(v1,sp)
    q = CalcLength(v2,sp)
    r = CalcDot(v1,v2,sp)
    oi_real(1) = CalcAngle(v1,v2,sp)*180.0/cPi
    call WriteValue(' -> Angle [deg] = ', oi_real, 1, "(2x,F8.4)")
   else
    if (sp.eq.'d') sp2='r'
    if (sp.eq.'r') sp2='d'
    call CalcCross(v1,v2,vc,sp,sp2,0)
    oi_int(1:3)=int(vc(1:3))
    if (sp.eq.'d') then
     call WriteValue(' ', oi_int, 3, "(' -> p x q = (',2(i3,1x),i3,')')")
    else
     call WriteValue(' ', oi_int, 3, "(' -> p x q = [',2(i3,1x),i3,']')")
    end if
   end if
  end if
  call ReadValue(' Another computation? (1/0) : ', oi_int, 1)
  another = oi_int(1)
 end do 
 
end program CTEMlatgeom

!--------------------------------------------------------------------------
!
! SUBROUTINE: PrintMenu 
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief Print a menu and return a selection
!
!> @param isel seletion number
! 
!> @date 1/5/99   MDG 1.0 original
!> @date  5/21/01  MDG 2.0 f90
!> @date 4/16/13 MDG 3.0 rewrite
!--------------------------------------------------------------------------
subroutine PrintMenu(isel)

use local
use io
        
integer(kind=irg),INTENT(OUT)  	:: isel	!< selection 
integer(kind=irg)				:: io_int(1)

 mess = ' Select from the following options: '; call Message("(A/)")
 mess = ' [1] length of direct space vector'; call Message("(A)")
 mess = ' [2] length of reciprocal space vector'; call Message("(A)")
 mess = ' [3] angle between direct space vectors'; call Message("(A)")
 mess = ' [4] angle between reciprocal space vectors'; call Message("(A)")
 mess = ' [5] cross product, direct space vectors'; call Message("(A)")
 mess = ' [6] cross product, reciprocal space vectors'; call Message("(A/)")
 call ReadValue('    Enter selection: ', io_int, 1)
 isel = io_int(1)

end subroutine PrintMenu
