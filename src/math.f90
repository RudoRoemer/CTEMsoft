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
! CTEMsoft2013:math.f90
!--------------------------------------------------------------------------
!
! MODULE: math
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief collection of mathematical/numerical routines that don't fit anywhere else
! 
!> @date   10/13/98 MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   03/19/13 MDG 3.0 updated all routines
!> @date   11/13/13 MDG 4.0 added MatrixExponential routine
!--------------------------------------------------------------------------
! ###################################################################
!  

module math

use local

contains

!--------------------------------------------------------------------------
!
! SUBROUTINE:mInvert
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Invert a 3x3 matrix
!
!> @details  Invert a 3x3 matrix; if unitary, simply transpose
!
!> @param a input matrix
!> @param b output matrix
!> @param uni .TRUE. if unitary matrix, .FALSE. otherwise
!
!> @todo this should really be replaced by a BLAS call
! 
!> @date   10/13/98 MDG 1.0 original
!> @date    4/ 5/00 MDG 1.1 added inverse of unitary matrix
!> @date    5/19/01 MDG 2.0 f90
!> @date   11/27/01 MDG 2.1 added kind support
!
!--------------------------------------------------------------------------
subroutine mInvert(a,b,uni)

use local
use error

IMPLICIT NONE

real(kind=dbl),INTENT(IN)		:: a(3,3) 		!< input matrix
real(kind=dbl),INTENT(OUT)		:: b(3,3)		!< output matrix
logical,INTENT(IN)				:: uni		!< unitary logical
real(kind=dbl)					:: d			!< auxiliary variable

! it is a regular (non-unitary) matrix
 if (.not.uni) then 
  d = a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+ &
         a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)- &
         a(1,2)*a(2,1)*a(3,3)-a(1,1)*a(2,3)*a(3,2)
  if (d.ne.0.0) then
   b(1,1)=a(2,2)*a(3,3)-a(2,3)*a(3,2)
   b(1,2)=a(1,3)*a(3,2)-a(1,2)*a(3,3)
   b(1,3)=a(1,2)*a(2,3)-a(1,3)*a(2,2)
   b(2,1)=a(2,3)*a(3,1)-a(2,1)*a(3,3)
   b(2,2)=a(1,1)*a(3,3)-a(1,3)*a(3,1)
   b(2,3)=a(1,3)*a(2,1)-a(1,1)*a(2,3)
   b(3,1)=a(2,1)*a(3,2)-a(2,2)*a(3,1)
   b(3,2)=a(1,2)*a(3,1)-a(1,1)*a(3,2)
   b(3,3)=a(1,1)*a(2,2)-a(1,2)*a(2,1)
   b = b/d
  else
   call FatalError('mInvert','matrix has zero determinant')
  end if
 else
! it is a unitary matrix, so simply get the transpose
  b = transpose(a)
 endif

end subroutine mInvert
      
!--------------------------------------------------------------------------
!
! SUBROUTINE:cInvert
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Invert a 3x3 complex matrix
!
!> @param a input matrix
!> @param b output matrix
!
!> @todo this should really be replaced by a BLAS call
! 
!> @date   10/13/98 MDG 1.0 original
!> @date    4/ 5/00 MDG 1.1 added inverse of unitary matrix
!> @date    5/19/01 MDG 2.0 f90
!> @date   11/27/01 MDG 2.1 added kind support
!
!--------------------------------------------------------------------------
subroutine cInvert(a,b)

use local
use error

IMPLICIT NONE

complex(kind=dbl),INTENT(IN)		:: a(3,3) 		!< input matrix
complex(kind=dbl),INTENT(OUT)		:: b(3,3)		!< output matrix
complex(kind=dbl)					:: d			!< auxiliary variable

  d = a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+ &
      a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)- &
      a(1,2)*a(2,1)*a(3,3)-a(1,1)*a(2,3)*a(3,2)
  if (abs(d).ne.0.D0) then
   b(1,1)=a(2,2)*a(3,3)-a(2,3)*a(3,2)
   b(1,2)=a(1,3)*a(3,2)-a(1,2)*a(3,3)
   b(1,3)=a(1,2)*a(2,3)-a(1,3)*a(2,2)
   b(2,1)=a(2,3)*a(3,1)-a(2,1)*a(3,3)
   b(2,2)=a(1,1)*a(3,3)-a(1,3)*a(3,1)
   b(2,3)=a(1,3)*a(2,1)-a(1,1)*a(2,3)
   b(3,1)=a(2,1)*a(3,2)-a(2,2)*a(3,1)
   b(3,2)=a(1,2)*a(3,1)-a(1,1)*a(3,2)
   b(3,3)=a(1,1)*a(2,2)-a(1,2)*a(2,1)
   b = b/d
  else
   call FatalError('cInvert','Matrix has complex zero determinant')
  end if

end subroutine cInvert


!--------------------------------------------------------------------------
!
! SUBROUTINE: MatrixExponential
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute the exponential of a dynamical matrix
!
!> @details This routine uses two different methods, one based on the Taylor
!> expansion, the other on Pade approximants, both with "scaling & squaring".
!> Currently, the routine targets an accuracy level of 10^{-9}.
!> This routine uses table 1 in "Nineteen dubious ways to compute the exponential of a matrix,
!> twenty-five years later", C. Moler, C. Van Loan, SIAM Review, 45, 1 (2003)
!
!> @param A input matrix
!> @param E output matrix
!> @param z0 slice thickness in [nm]
!> @param TP 'Tayl' or 'Pade', to select method
!> @param nn number of row/column entries in A
!
!> @date   09/16/13 MDG 1.0 original, tested against analytical version for small array
!--------------------------------------------------------------------------
recursive subroutine MatrixExponential(A,E,z0,TP,nn)

use local
use io
use error

IMPLICIT NONE

integer(kind=irg),INTENT(IN)		:: nn
complex(kind=dbl),INTENT(IN)		:: A(nn,nn)
complex(kind=dbl),INTENT(OUT)		:: E(nn,nn)
real(kind=dbl),INTENT(IN)		:: z0
character(4),INTENT(IN)		:: TP

real(kind=dbl)				:: modA, pref, sgn
complex(kind=dbl),allocatable		:: B(:,:), add(:,:), Nqq(:,:), Dqq(:,:), C(:,:)

integer(kind=irg)			:: i, k, j, icnt, q, istat, ilev
integer(kind=irg)    			:: INFO, LDA, MILWORK
integer(kind=irg),allocatable		:: JPIV(:)
complex(kind=dbl),allocatable 		:: MIWORK(:)
integer(kind=irg),parameter		:: kTaylor(6) = (/ 3, 4, 6, 8, 7, 6 /)		! from table 1 in reference above
integer(kind=irg),parameter		:: jTaylor(6) = (/ 1, 2, 3, 5, 9, 13 /)
integer(kind=irg),parameter		:: qPade(6) = (/ 2, 3, 4, 4, 4, 4 /)
integer(kind=irg),parameter		:: jPade(6) = (/ 0, 0, 1, 5, 8, 11 /)

! set output array to zero
E = dcmplx(0.D0,0.D0)

!get the max( ||A|| ) value and determine index into (k,j) or (q,j) pairs
modA = maxval(cabs(A)*z0)
ilev = nint(alog10(modA))+3
if (ilev.le.0) ilev = 1		! can not be smaller than 1

! if modA gets to be too large, abort with a message
if (modA.gt.1000.D0) then
  mess = 'MatrixExponential routine can not deal with ||A|| > 1000.0'
  call Message("(/A/)")
  stop 'Program aborted'
end if


if (TP.eq.'Tayl') then ! use scaling and squaring for the Taylor expansion
	k = kTaylor(ilev)
	j = jTaylor(ilev)

	! allocate an auxiliary array
	allocate( B(nn,nn), add(nn,nn), stat=istat )
	if (istat.ne.0) then 
	  call FatalError('MatrixExponential','Error allocating arrays for Taylor approximation')
	  stop
	end if
	
	! perform the scaling step
	B = (A * z0) / 2.0D0**j ! dcmplx(2.0**j,0.0)
	
	! initialize the diagonal of E
	forall (i=1:nn) E(i,i) = dcmplx(1.0D0,0.D0)
	
	! loop over the Taylor series
	add = B
	E = E + add
	do icnt=2,k
	  add = matmul( add, B/dcmplx(icnt) )
	  E = E + add
	end do
	
	! and deallocate the auxiliary arrays
	deallocate(add, B)

else ! Pade approximation for target accuracy 10^(-9)
	q = qPade(ilev)
	j = jPade(ilev)

	! allocate auxiliary arrays
	allocate(B(nn,nn),C(nn,nn), Nqq(nn,nn), Dqq(nn,nn), stat=istat )
	if (istat.ne.0) then 
	  call FatalError('MatrixExponential','Error allocating arrays for Pade approximation')
	  stop
	end if
	
	! perform the scaling step
	B = (A * z0) / 2.D0**j  ! dcmplx(2.0**j,0.0)
	C = B
		
	! initialize the diagonal of both arrays
	Nqq = dcmplx(0.D0,0.D0)
	forall (i=1:nn) Nqq(i,i) = dcmplx(1.0D0,0.D0)
	Dqq = Nqq

	! init some constants
	pref = 1.D0
	sgn = -1.D0
	
	! and loop
	do icnt=1,q
	  pref = pref * dble(q-icnt+1) / dble(icnt) / dble(2*q-icnt+1)
	  Nqq = Nqq + pref * C
	  Dqq = Dqq + sgn * pref * C
	  sgn = -sgn
	  C = matmul( C, B )
	end do
	
	! get the inverse of Dqq using the LAPACK routines zgetrf and zgetri
	LDA = nn
	allocate( JPIV(nn) )
 	call zgetrf(nn,nn,Dqq,LDA,JPIV,INFO)
	if (INFO.ne.0) call FatalError('Error in MatrixExponential: ','ZGETRF return not zero')

	MILWORK = 64*nn 
 	allocate(MIWORK(MILWORK))

	MIWORK = dcmplx(0.0_dbl,0.0_dbl)
 	call zgetri(nn,Dqq,LDA,JPIV,MIWORK,MILWORK,INFO)
 	if (INFO.ne.0) call FatalError('Error in MatrixExponential: ','ZGETRI return not zero')

	! and compute E
	E = matmul( Dqq, Nqq )
	
	! clean up
	deallocate(Nqq, Dqq, C, B, JPIV, MIWORK)
end if

! and finally compute the power 2^j of the matrix E (i.e. the squaring step)
do icnt = 1,j
  E = matmul( E, E )
end do

end subroutine MatrixExponential









end module math
