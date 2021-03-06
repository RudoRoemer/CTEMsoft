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
! CTEMsoft2013:CTEMzap.f90
!--------------------------------------------------------------------------
!
! PROGRAM:CTEMzap 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief PostScript output of kinematical zone axis diffraction patterns
!
! 
!> @date 12/11/98   MDG 1.0 original
!> @date 04/08/13 MDG 2.0 rewrite
!--------------------------------------------------------------------------
program CTEMzap

use local
use crystalvars
use crystal
use symmetryvars
use symmetry
use graphics
use files
use postscript
use io
use diffraction

IMPLICIT NONE

real(kind=sgl)			:: io_real(1)


 progname = 'CTEMzap.f90'
 progdesc = 'Kinematical Zone Axis Diffraction Patterns'
 call CTEMsoft

 SG % SYM_reduce=.TRUE.

! read crystal information, microscope voltage, and camera length
 call CrystalData
 call GetVoltage
 call ReadValue(' Camera length L  [mm, real] ', io_real, 1)
 camlen = io_real(1)

! generate all atom positions in the fundamental unit cell
 call CalcPositions('v')

! open PostScript file
 call PS_openfile

! generate a set of zone axis patterns
 call DiffPage

! close Postscript file
 call PS_closefile

end program CTEMzap
