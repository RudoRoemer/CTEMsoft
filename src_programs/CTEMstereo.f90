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
! CTEMsoft2013:CTEMstereo.f90
!--------------------------------------------------------------------------
!
! PROGRAM:CTEMstereo 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief STandard stereographic projections
!
!> @todo fix a bug that causes a segmentation fault when different range parameters are
!> entered.  The program does produce a correct ps file, but does not end gracefully....
!
!> @date   10/13/98 MDG 1.0 original
!> @date    5/21/01 MDG 2.0 f90
!> @date  4/16/13 MDG 3.0 rewrite
!--------------------------------------------------------------------------
program CTEMstereo

use local
use crystalvars
use crystal
use symmetryvars
use graphics
use files
use postscript
use io

IMPLICIT NONE

character(1)   			:: sp
logical        			:: topbot
integer(kind=irg)        	:: hm,km,lm,i,iview(3), io_int(3)

 progname = 'CTEMstereo.f90'
 progdesc = 'Stereographic projections (direct/ reciprocal space)'
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
 mess = '  Enter the maximum index for h,k and l, or for '; call Message("(/A)")
 mess = '  u,v, and w. For a hexagonal system, please use'; call Message("(A)")
 mess = '  4-index notation [uv.w] or (hk.l) to determine'; call Message("(A)")
 mess = '  the largest index.'; call Message("(A/)")
 call ReadValue(' Enter maximum indices (h,k,l) : ', io_int, 3)

! call the drawing routine
 call StereoProj(sp,iview,io_int(1),io_int(2),io_int(3),topbot)

! close Postscript file
 call PS_closefile

end program CTEMstereo

