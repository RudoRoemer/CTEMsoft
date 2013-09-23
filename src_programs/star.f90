!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CTEMsoft: star.f90                                                            !
! Copyright (c) 2001, 2002  Marc De Graef/Carnegie Mellon University            !
!                                                                               !
!     This program is free software; you can redistribute it and/or modify      !
!     it under the terms of the GNU General Public License as published by      !
!     the Free Software Foundation; either version 2 of the License, or         !
!     (at your option) any later version.                                       !
!                                                                               !
!     This program is distributed in the hope that it will be useful,           !
!     but WITHOUT ANY WARRANTY; without even the implied warranty of            !
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             !
!     GNU General Public License for more details.                              !
!                                                                               !
!     You should have received a copy of the GNU General Public License         !
!     along with this program; if not, write to the Free Software               !
!     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA !
!                                                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! -Fortran-90 Source Code-
! ###################################################################
!  Introduction to Conventional Transmission Electron Microscopy
!
!  by Marc De Graef
!
!  Cambridge University Press
!  SBN 0521629950 (Paperback)
!  SBN 0521620066 (Hardback)
! 
!  FILE: "star.f90"
!                                    created: 1/5/99 {11:26:55 AM}
!                                last update: 6/5/2001 {7:24:17 PM}
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description: print star of a reciprocal lattice position for a given space group;
!               also list the Fourier coefficients of the lattice potential.
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  1/5/99  MDG 1.0 original
!  5/18/01 MDG 2.0 f90
!  3/07/02 MDG 2.1 added CalcUcg support
! ###################################################################
program star

use local
use io
use files
use symmetryvars
use symmetry
use diffraction
use dynamical

IMPLICIT NONE

integer(kind=irg) :: g(3),gg(3),ans,n,i,j
real(kind=dbl)    :: kk(3),stmp(0:47,3)
logical           :: first
character(1)      :: space

 progname = 'star.f90'
 progdesc = 'Compute the star of a reciprocal lattice vector'
 call CTEMsoft

 SG % SYM_reduce=.FALSE.
 space = 'r'
 call CrystalData
 call GetVoltage
 call CalcPositions('v')

 if (SG % SYM_centrosym) then 
  mess = 'structure is centrosymmetric'; call Message("(A/)") 
 else 
  mess = 'structure is non-centrosymmetric'; call Message("(A/)") 
 end if
 ans = 1
 first = .TRUE.

 do while (ans.eq.1)
  mess = ''; call Message("(/A)")
  mess = 'Enter reciprocal lattice coordinates [I] : '; call GetInt(3)
  g(1:3) = io_int(1:3)
  kk = dble(g)
  call CalcStar(kk,n,stmp,space)
  mess = 'number of equivalent planes in star = '; oi_int(1)=n; call WriteInt(1,"(I3)")
! compute and display the structure factor
  do i=0,n-1
   do j=1,3 
    gg(j) = nint(stmp(i,j))
   end do
   call CalcUcg(gg)
   if (i.eq.0) then 
    call Printrlp(first)
   else
    call Printrlp
   endif
  end do
  mess = 'Another star ? (1/0) '; call GetInt(1)
  ans = io_int(1)
 end do
end program
