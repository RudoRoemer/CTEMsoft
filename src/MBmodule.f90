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
! CTEMsoft2013:MBmodule.f90
!--------------------------------------------------------------------------
!
! MODULE: MBmodule
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief routine for multi-beam diffraction calculations, used by more than one program
!
!> @date 7/10/13   MDG 1.0 original
!> @date 7/12/13   MDG 1.1 added forward cube to ball to quaternion mappings
!> @date 8/01/13   MDG 1.2 added standard Lambert projection
!> @date 8/12/13   MDG 1.3 added inverse Lambert projections for Ball to Cube
!> @date 9/20/13   MDG 1.4 added ApplyLaueSymmetry
!--------------------------------------------------------------------------
module MBmodule

use local

IMPLICIT NONE

contains


!--------------------------------------------------------------------------
!
! SUBROUTINE: CalcBWint
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute the scattered intensities for a range of thicknesses
!
!> @param nn number of strong beams
!> @param nw number of weak beams
!> @param nt number of thickness values
!> @param thick thickness array
!> @param inten output intensity list, including weak beams
!
!> @date  10/13/98 MDG 1.0 original
!> @date   7/04/01 MDG 2.0 f90
!> @date  04/29/13 MDG 3.0 inclusion of Bethe weak beams
!--------------------------------------------------------------------------
subroutine CalcBWint(nn,nw,nt,thick,inten)

use local
use io
use diffraction
use kvectors
use gvectors
use dynamical
use constants

IMPLICIT NONE

integer(kind=irg),INTENT(IN)	:: nn			!< number of strong beams
integer(kind=irg),INTENT(IN)	:: nw			!< number of weak beams
integer(kind=irg),INTENT(IN)	:: nt			!< number of thickness values
real(kind=sgl),INTENT(IN)	:: thick(nt)		!< thickness array
real(kind=sgl),INTENT(INOUT)	:: inten(nt,nn+nw)	!< output intensities (both strong and weak)

integer(kind=irg)		:: i,j,IPIV(nn), ll(3), jp
complex(kind=dbl)		:: CGinv(nn,nn), Minp(nn,nn),diag(nn),Wloc(nn), lCG(nn,nn), lW(nn), &
				lalpha(nn), delta(nn,nn), weak(nw,nn), Ucross(nw,nn), tmp(nw,nn), c
real(kind=sgl) 			:: th

! compute the eigenvalues and eigenvectors
 Minp = DynMat
 IPIV = 0
 call BWsolve(Minp,Wloc,lCG,CGinv,nn,IPIV)

! the alpha coefficients are in the first column of the inverse matrix
! the minus sign in W(i) stems from the fact that k_n is in the direction
! opposite to the foil normal
 lW = cPi*Wloc/cmplx(ktmp%kn,0.0)
 do i=1,nn
  lalpha(i) = CGinv(i,1)
 end do

! in preparation for the intensity computation, we need the prefactor array for the
! weak beam amplitude computation...
! we will also need a potential coefficient array for the cross coefficients, which we 
! compute in the same loop
 do j=1,nn   ! strong beam loop
   do jp=1,nw  ! weak beam loop
! prefactor value
     c = cmplx(2.D0*BetheParameter%weaksg(jp)/mLambda) - 2.D0*ktmp%kn*Wloc(j)
     weak(jp,j) = cmplx(-1.D0,0.D0)/c
! cross potential coefficient
     ll(1:3) = BetheParameter%weakhkl(1:3,jp) - BetheParameter%stronghkl(1:3,j)
     Ucross(jp,j) = LUT( ll(1),ll(2),ll(3) )
   end do
 end do

! compute the strong beam intensities, stored in the first nn slots of inten 
! we can also compute the weak beams, since they make use of the same diag(1:nn) expression
! as the strong beams, plus a few other factors (excitation error, wave length, Fourier coefficients)
 do i=1,nt
  th = thick(i)
  diag(1:nn)=exp(-th*imag(lW(1:nn)))*cmplx(cos(th*real(lW(1:nn))),sin(th*real(lW(1:nn))))*lalpha(1:nn)
! the delta array is common to the strong and weak beam intensity computation, so we compute it first
  do j=1,nn
   delta(j,1:nn) = lCG(j,1:nn)*diag(1:nn)
  end do
! strong beams
  do j=1,nn
   inten(i,j) = cabs(sum(delta(j,1:nn)))**2
  end do 
! weak beams
  tmp = matmul(Ucross,delta)
  do jp=1,nw
   inten(i,nn+jp) = cabs( sum(weak(jp,1:nn)*tmp(jp,1:nn)) )**2
  end do  
 end do
   
end subroutine CalcBWint



end module MBmodule
