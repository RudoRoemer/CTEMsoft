!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CTEMsoft:  family.f90                                                         !
! Copyright (c) 2001, 2002  Marc De Graef/Carnegie Mellon University (CMU)      !
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
! ###################################################################
!  Introduction to Conventional Transmission Electron Microscopy
!
!  by Marc De Graef
!
!  Cambridge University Press
!  SBN 0521629950 (Paperback)
!  SBN 0521620066 (Hardback)
! 
!  FILE: "family.f90"
!                                    created: 10/16/98 {9:29:46 AM} 
!                                last update: 7/16/99 {10:37:32 AM} 
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description:  This program draws a stereographic projection of 
!                either real or reciprocal space for an arbitrary
!                crystal system and viewing direction.  The program
!                is different from stereo.f in that it only draws
!                the requested families <uvw> or {hkl}.  Multiple
!                output files can be generated.
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/16/98 MDG 1.0 original
!   5/21/01 MDG 2.0 f90
! ###################################################################
program family

use local
use crystalvars
use crystal
use symmetryvars
use symmetry
use postscript
use io
use math
use graphics
use files

IMPLICIT NONE

character(1)         :: sp
logical              :: nn,topbot
real(kind=sgl)       :: rr(3),g(3),r(3),M(3,3), CX, CY, CRad, negthresh,xst,yst
integer(kind=irg)    :: h,k,l,hkl(3),iview(3),cr,ans,sgn,i,num

 progname = 'family.f90'
 progdesc = 'Stereographic projection of family'
 call CTEMsoft
 
 SG % SYM_reduce=.TRUE.
 topbot=.TRUE.
! 20cm radius projection circle [inches]
 CRad = 3.937
 CX = 3.25
 CY = 3.5
 negthresh=-0.0001
! read crystal information
 call CrystalData
 sgn = 1 
! main loop 
 do while (sgn.eq.1) 
! real space or reciprocal space?
  call GetDrawingSpace(sp)
! viewing direction (watch for hexagonal indices !)
  call GetViewingDirection(iview)
! create transformation matrix
  call ProjectionMatrix(iview,M)
! open PostScript file
  call PS_openfile
! write text and draw projection circle
  call DrawSPFrame(CX, CY, CRad, iview, sp)
  ans = 1
! loop over families
  do while (ans.eq.1)
! loop over all points and draw projection+label
   call GetIndex(hkl,sp)
   call CalcFamily(hkl,num,sp)
   do i=1,num
    h=itmp(i,1)
    k=itmp(i,2)
    l=itmp(i,3)
    hkl(1)=h
    hkl(2)=k
    hkl(3)=l
! reduce to smallest integers to avoid overlap
! of indices, such as (111) and (222)
    call IndexReduce(hkl)
    g(1)=float(h)
    g(2)=float(k)
    g(3)=float(l)
    h=hkl(1)
    k=hkl(2)
    l=hkl(3)
    call TransSpace(g,r,sp,'c')
    call NormVec(r,'c')
! apply viewing tansformation
    rr = matmul(M,r)
! compute stereographic projection coordinates
    xst=CX+CRad*rr(1)/(1.0+abs(rr(3)))
    yst=CY+CRad*rr(2)/(1.0+abs(rr(3)))
! spit them out
    write (*,*) h,k,l,xst,yst
! and draw the projection point along with its label
    cr=1
    if (rr(3).gt.negthresh) then
     call PS_filledcircle(xst,yst,0.015/PS % psscale,0.0)
     nn = .TRUE.
     call DumpIndices(sp,h,k,l,cr,xst,yst,nn)
    else if (topbot) then
     call PS_circle(xst,yst,0.035/PS % psscale)
     nn = .FALSE.
     call DumpIndices(sp,h,k,l,cr,xst,yst,nn)
    endif
   end do
! another family on the same drawing ?
   mess = ' Another family (1/0) ? '; call GetInt(1)
   ans = io_int(1)
  end do
! close Postscript file
  call PS_closefile
! loop for another drawing ?
   mess = ' Another pattern (1/0) ? '; call GetInt(1)
   sgn = io_int(1)
 end do
end program
