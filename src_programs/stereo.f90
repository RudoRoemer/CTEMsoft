!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CTEMsoft:   stereo.f90                                                        !
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
!  FILE: "stereo.f90"
!                                    created: 10/13/98 {9:29:46 AM} 
!                                
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description: standard stereographic projection (direct & reciprocal)
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  10/13/98 MDG 1.0 original
!   5/21/01 MDG 2.0 f90
! ###################################################################
program stereo

use local
use crystalvars
use crystal
use symmetryvars
use graphics
use files
use postscript
use io

IMPLICIT NONE

character(1)   :: sp
logical        :: topbot
integer        :: hm,km,lm,i,iview(3)

 progname = 'stereo.f90'
 progdesc = 'Stereographic projections'
 call CTEMsoft
 
 SG % SYM_reduce=.TRUE.
 topbot=.FALSE.
! read crystal information
 call CrystalData
! real space or reciprocal space
 call GetDrawingSpace(sp)
! viewing direction
 call GetViewingDirection(iview)
! open PostScript file
 call PS_openfile
! get index ranges
 mess = 'Enter the maximum index for h,k and l, or for '
 mess = 'u,v, and w. For a hexagonal system, please use'
 mess = '4-index notation [uv.w] or (hk.l) to determine'
 mess = 'the largest index.'
 mess = 'Enter maximum indices (h,k,l) : '; call GetInt(3)
 hm = io_int(1)
 km = io_int(2)
 lm = io_int(3)
! call the drawing routine
 call StereoProj(sp,iview,hm,km,lm,topbot)
! close Postscript file
 call PS_closefile
end program

