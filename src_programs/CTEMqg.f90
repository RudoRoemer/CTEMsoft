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
! CTEMsoft2013:CTEMqg.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMqg 
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief compute the complex extinction distance q_g
!
!> @date   10/13/98 MDG 1.0 original
!> @date    5/27/01 MDG 2.0 f90
!> @date  4/16/13 MDG 3.0 rewrite
!--------------------------------------------------------------------------
program CTEMqg

use local
use io
use crystal
use files
use diffraction
use symmetry
use constants
use postscript
use dynamical

IMPLICIT NONE

integer(kind=irg)        	:: ind(3),ans, oi_int(3)
real(kind=sgl)			:: oi_real(7)
complex(kind=sgl)		:: oi_cmplx(1)
real(kind=sgl)           	:: preg

 progname = 'CTEMqg.f90'
 progdesc = 'Display potential coefficient values'
 call CTEMsoft
 
 call CrystalData
 call GetVoltage
 call CalcPositions('v')
 preg = 2.0 * sngl(cRestmass*cCharge/cPlanck**2)*1.0E-18

 ans = 1
 do while (ans.eq.1)
  mess = 'Enter Miller indices :'; call Message("(/A)")
  call GetIndex(ind,'r')
  call CalcUcg(ind)
  
  mess = '   h  k  l    |g|    Ucg_r     Ucg_i      |Ug|      phase     |Ugp|     phase     xi_g   xi_gp    ratio    Re-1/q_g-Im'
  call Message("(A)") 
 oi_int(1:3) = rlp%hkl(1:3)
 call WriteValue('',oi_int, 3, "(1x,3I3,1x,$)")
 oi_real(1) = rlp%g
 call WriteValue('',oi_real, 1, "(F9.4,$)")
 oi_cmplx(1) = rlp%Ucg
 call WriteValue('',oi_cmplx, 1, "(2F10.6,1x,$)")
 oi_real(1:7)  = (/ rlp%Umod,rlp%Vphase*180.0/sngl(cPi),rlp%Upmod,rlp%Vpphase*180.0/sngl(cPi),rlp%xg,rlp%xgp,rlp%ar /)
 call WriteValue('',oi_real, 7, "(4F10.5,3F8.1,$)")
 oi_cmplx(1) = rlp%qg
 call WriteValue('',oi_cmplx, 1, "(2F8.5)")



  call ReadValue(' Another one ? (1/0) : ', oi_int, 1)
  ans = oi_int(1)
 end do 

end  program CTEMqg
       



