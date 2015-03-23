! ###################################################################
! Copyright (c) 2013-2014, Marc De Graef/Carnegie Mellon University
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
! EMsoft:math.f90
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

use error

IMPLICIT NONE

real(kind=dbl),INTENT(IN)               :: a(3,3)               !< input matrix
real(kind=dbl),INTENT(OUT)              :: b(3,3)               !< output matrix
logical,INTENT(IN)                              :: uni          !< unitary logical
real(kind=dbl)                                  :: d                    !< auxiliary variable

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

use error

IMPLICIT NONE

complex(kind=dbl),INTENT(IN)            :: a(3,3)               !< input matrix
complex(kind=dbl),INTENT(OUT)           :: b(3,3)               !< output matrix
complex(kind=dbl)                                       :: d                    !< auxiliary variable

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
!> @date 09/16/13 MDG 1.0 original, tested against analytical version for small array
!> @date 06/05/14 MDG 1.1 updated IO
!--------------------------------------------------------------------------
recursive subroutine MatrixExponential(A,E,z0,TP,nn)

use io
use error

IMPLICIT NONE

integer(kind=irg),INTENT(IN)            :: nn
complex(kind=dbl),INTENT(IN)            :: A(nn,nn)
complex(kind=dbl),INTENT(OUT)           :: E(nn,nn)
real(kind=dbl),INTENT(IN)               :: z0
character(4),INTENT(IN)         :: TP

real(kind=dbl)                          :: modA, pref, sgn
complex(kind=dbl),allocatable           :: B(:,:), add(:,:), Nqq(:,:), Dqq(:,:), C(:,:)

integer(kind=irg)                       :: i, k, j, icnt, q, istat, ilev

integer(kind=irg)                       :: INFO, LDA, MILWORK
integer(kind=irg),allocatable           :: JPIV(:)
complex(kind=dbl),allocatable           :: MIWORK(:)

integer(kind=irg),parameter             :: kTaylor(6) = (/ 3, 4, 6, 8, 7, 6 /)          ! from table 1 in reference above
integer(kind=irg),parameter             :: jTaylor(6) = (/ 1, 2, 3, 5, 9, 13 /)
integer(kind=irg),parameter             :: qPade(6) = (/ 2, 3, 4, 4, 4, 4 /)
integer(kind=irg),parameter             :: jPade(6) = (/ 0, 0, 1, 5, 8, 11 /)

! set output array to zero
E = dcmplx(0.D0,0.D0)

!get the max( ||A|| ) value and determine index into (k,j) or (q,j) pairs
modA = maxval(cdabs(A)*z0)
ilev = nint(log10(modA))+3
if (ilev.le.0) ilev = 1         ! can not be smaller than 1

! if modA gets to be too large, abort with a message
if (modA.gt.10000.D0) call FatalError('MatrixExponential','  routine can not deal with ||A|| > 10000.0')

if (TP.eq.'Tayl') then ! use scaling and squaring for the Taylor expansion
        k = kTaylor(ilev)
        j = jTaylor(ilev)

        ! allocate an auxiliary array
        allocate( B(nn,nn), add(nn,nn), stat=istat )
        if (istat.ne.0) call FatalError('MatrixExponential',' Error allocating arrays for Taylor approximation')
        
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
        if (istat.ne.0) call FatalError('MatrixExponential',' Error allocating arrays for Pade approximation')
        
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
        if (INFO.ne.0) call FatalError('Error in MatrixExponential: ',' ZGETRF return not zero')

        MILWORK = 64*nn 
        allocate(MIWORK(MILWORK))

        MIWORK = dcmplx(0.0_dbl,0.0_dbl)
        call zgetri(nn,Dqq,LDA,JPIV,MIWORK,MILWORK,INFO)
        if (INFO.ne.0) call FatalError('Error in MatrixExponential: ',' ZGETRI return not zero')

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


!--------------------------------------------------------------------------
!
! FUNCTION BesselIn
!
!> @author Marc De Graef, Carnegie Mellon University/ J-P Moreau (www.jpmoreau.fr)
!
!> @brief compute the first order modified Bessel function I_n(x)
!
!> @details original source: http://jean-pierre.moreau.pagesperso-orange.fr/Fortran/tbessi_f90.txt (12/31/14)
!>* -------------------------------------------------------------------- *
!>*   Reference: From Numath Library By Tuan Dang Trong in Fortran 77.   *
!>*                                                                      *
!>*                               F90 Release 1.2 By J-P Moreau, Paris.  *
!>*                                        (www.jpmoreau.fr)             *
!> ----------------------------------------------------------------------
!> routine tested and compared with both Matlab and IDL output.
!
!> @param X input argument
!> @param N Bessel order
!
!> @date 12/31/14 MDG 1.0 original, rewritten with EMsoft module calls and renamed
!--------------------------------------------------------------------------
recursive function BesselIn(X,N) result(BESSI)
! original comment:
!     This subroutine calculates the first kind modified Bessel function
!     of integer order N, for any REAL X. We use here the classical
!     recursion formula, when X > N. For X < N, the Miller's algorithm
!     is used to avoid overflows. 
!     REFERENCE:
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.
!
use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN)       :: X
integer(kind=irg),INTENT(IN)    :: N

integer(kind=irg),parameter     :: IACC = 40
real(kind=dbl),parameter        :: BIGNO = 1.D10, BIGNI = 1.D-10
integer(kind=irg)               :: M, J
real(kind=dbl)                  :: BESSI,BESSI0,BESSI1,TOX,BIM,BI,BIP

! special cases
if (N.EQ.0) then
  BESSI = BesselI0(X)
  return
end if

if (N.EQ.1) then
  BESSI = BesselI1(X)
  return
endif

if(X.EQ.0.D0) then
  BESSI=0.D0
  return
endif

! set up the recursion 
TOX = 2.D0/X
BIP = 0.D0
BI  = 1.D0
BESSI = 0.D0
M = 2*((N+int(sqrt(float(IACC*N)))))

DO J = M,1,-1
   BIM = BIP+dble(J)*TOX*BI
   BIP = BI
   BI  = BIM
   if (dabs(BI).GT.BIGNO) then
     BI  = BI*BIGNI
     BIP = BIP*BIGNI
     BESSI = BESSI*BIGNI
   end if
   if (J.EQ.N) BESSI = BIP
end do
BESSI = BESSI*BesselI0(X)/BI

end function BesselIn

!--------------------------------------------------------------------------
!
! FUNCTION BesselI0
!
!> @author Marc De Graef, Carnegie Mellon University/ J-P Moreau (www.jpmoreau.fr)
!
!> @brief compute the first order modified Bessel function I_0(x)
!
!> @details original source: http://jean-pierre.moreau.pagesperso-orange.fr/Fortran/tbessi_f90.txt (12/31/14)
!>* -------------------------------------------------------------------- *
!>*   Reference: From Numath Library By Tuan Dang Trong in Fortran 77.   *
!>*                                                                      *
!>*                               F90 Release 1.2 By J-P Moreau, Paris.  *
!>*                                        (www.jpmoreau.fr)             *
!> ----------------------------------------------------------------------
!> routine tested and compared with both Matlab and IDL output.
!
!> @param X input argument
!
!> @date 12/31/14 MDG 1.0 original, rewritten with EMsoft module calls and renamed
!--------------------------------------------------------------------------
recursive function BesselI0(X) result(BESSI0)

use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN)       :: X
real(kind=dbl)                  :: BESSI0

real(kind=dbl)                  :: Y,Ax,BX
real(kind=dbl),parameter        :: P1=1.D0, P2=3.5156229D0, P3=3.0899424D0, P4=1.2067492D0,  &
                                   P5=0.2659732D0, P6=0.360768D-1, P7=0.45813D-2, Q1=0.39894228D0, &
                                   Q2=0.1328592D-1, Q3=0.225319D-2, Q4=-0.157565D-2, Q5=0.916281D-2, &
                                   Q6=-0.2057706D-1,  Q7=0.2635537D-1, Q8=-0.1647633D-1, Q9=0.392377D-2

if(dabs(X).lt.3.75D0) then
  Y=(X/3.75D0)**2
  BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
else
  AX=dabs(X)
  Y=3.75D0/AX
  BX=dexp(AX)/dsqrt(AX)
  AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
  BESSI0=AX*BX
ENDIF

end function BesselI0


!--------------------------------------------------------------------------
!
! FUNCTION BesselI1
!
!> @author Marc De Graef, Carnegie Mellon University/ J-P Moreau (www.jpmoreau.fr)
!
!> @brief compute the first order modified Bessel function I_1(x)
!
!> @details original source: http://jean-pierre.moreau.pagesperso-orange.fr/Fortran/tbessi_f90.txt (12/31/14)
!>* -------------------------------------------------------------------- *
!>*   Reference: From Numath Library By Tuan Dang Trong in Fortran 77.   *
!>*                                                                      *
!>*                               F90 Release 1.2 By J-P Moreau, Paris.  *
!>*                                        (www.jpmoreau.fr)             *
!> ----------------------------------------------------------------------
!> routine tested and compared with both Matlab and IDL output.
!
!> @param X input argument
!
!> @date 12/31/14 MDG 1.0 original, rewritten with EMsoft module calls and renamed
!--------------------------------------------------------------------------
recursive function BesselI1(X) result(BESSI1)

use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN)       :: X
real(kind=dbl)                  :: BESSI1

real(kind=dbl)                  :: Y,AX,BX
real(kind=dbl),parameter        :: P1=0.5D0, P2=0.87890594D0, P3=0.51498869D0, P4=0.15084934D0, &
                                   P5=0.2658733D-1, P6=0.301532D-2, P7=0.32411D-3, Q1=0.39894228D0, &
                                   Q2=-0.3988024D-1, Q3=-0.362018D-2, Q4=0.163801D-2, Q5=-0.1031555D-1, &
                                   Q6=0.2282967D-1, Q7=-0.2895312D-1, Q8=0.1787654D-1, Q9=-0.420059D-2


if(dabs(X).lt.3.75D0) then
      Y=(X/3.75D0)**2
      BESSI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
else
      AX=dabs(X)
      Y=3.75D0/AX
      BX=dexp(AX)/dsqrt(AX)
      AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
      BESSI1=AX*BX
end if

end function BesselI1

!*****************************************************************************80
!*****************************************************************************80
!*****************************************************************************80
! the functions below are taken from the normal.f90 file posted on 
! http://people.sc.fsu.edu/~jburkardt/f_src/normal/normal.html
!*****************************************************************************80
!*****************************************************************************80
!*****************************************************************************80


function c4_normal_01 ( seed )

!*****************************************************************************80
!
!! C4_NORMAL_01 returns a unit pseudonormal C4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, complex ( kind = 4 ) C4_NORMAL_01, a unit pseudonormal value.
!
  implicit none

  complex ( kind = 4 ) c4_normal_01
  real ( kind = 4 ), parameter :: r4_pi = 3.141592653589793E+00
! real ( kind = 4 ) r4_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 4 ) v1
  real ( kind = 4 ) v2
  real ( kind = 4 ) x_c
  real ( kind = 4 ) x_r

  v1 = r4_uniform_01 ( seed )
  v2 = r4_uniform_01 ( seed )

  x_r = sqrt ( - 2.0E+00 * log ( v1 ) ) * cos ( 2.0E+00 * r4_pi * v2 )
  x_c = sqrt ( - 2.0E+00 * log ( v1 ) ) * sin ( 2.0E+00 * r4_pi * v2 )

  c4_normal_01 = cmplx ( x_r, x_c )

  return
end
function c8_normal_01 ( seed )

!*****************************************************************************80
!
!! C8_NORMAL_01 returns a unit pseudonormal C8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, complex ( kind = 8 ) C8_NORMAL_01, a sample of the PDF.
!
  implicit none

  complex ( kind = 8 ) c8_normal_01
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
! real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) v1
  real ( kind = 8 ) v2
  real ( kind = 8 ) x_c
  real ( kind = 8 ) x_r

  v1 = r8_uniform_01 ( seed )
  v2 = r8_uniform_01 ( seed )

  x_r = sqrt ( - 2.0D+00 * log ( v1 ) ) * cos ( 2.0D+00 * r8_pi * v2 )
  x_c = sqrt ( - 2.0D+00 * log ( v1 ) ) * sin ( 2.0D+00 * r8_pi * v2 )

  c8_normal_01 = cmplx ( x_r, x_c, kind = 8 )

  return
end
function i4_huge ( )

!*****************************************************************************80
!
!! I4_HUGE returns a "huge" I4.
!
!  Discussion:
!
!    On an IEEE 32 bit machine, I4_HUGE should be 2^31 - 1, and its
!    bit pattern should be
!
!     01111111111111111111111111111111
!
!    In this case, its numerical value is 2147483647.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) I4_HUGE, a "huge" I4.
!
  implicit none

  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_huge

  i4_huge = 2147483647

  return
end
function i4_normal_ab ( a, b, seed )

!*****************************************************************************80
!
!! I4_NORMAL_AB returns a scaled pseudonormal I4.
!
!  Discussion:
!
!    The normal probability distribution function (PDF) is sampled,
!    with mean A and standard deviation B.
!
!    The result is then rounded to the nearest integer.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A, the mean of the PDF.
!
!    Input, real ( kind = 4 ) B, the standard deviation of the PDF.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the
!    random number generator.
!
!    Output, integer ( kind = 4 ) I4_NORMAL_AB, a sample of the normal PDF.
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) b
  integer ( kind = 4 ) i4_normal_ab
  real ( kind = 4 ) r1
  real ( kind = 4 ) r2
  real ( kind = 4 ), parameter :: r4_pi = 3.141592653589793E+00
! real ( kind = 4 ) r4_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 4 ) x

  r1 = r4_uniform_01 ( seed )
  r2 = r4_uniform_01 ( seed )
  x = sqrt ( - 2.0E+00 * log ( r1 ) ) * cos ( 2.0E+00 * r4_pi * r2 )

  i4_normal_ab = nint ( a + b * x )

  return
end
function i8_normal_ab ( a, b, seed )

!*****************************************************************************80
!
!! I8_NORMAL_AB returns a scaled pseudonormal I8.
!
!  Discussion:
!
!    The normal probability distribution function (PDF) is sampled,
!    with mean A and standard deviation B.
!
!    The result is then rounded to the nearest integer.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the mean of the PDF.
!
!    Input, real ( kind = 8 ) B, the standard deviation of the PDF.
!
!    Input/output, integer ( kind = 8 ) SEED, a seed for the
!    random number generator.
!
!    Output, integer ( kind = 8 ) I8_NORMAL_AB, a sample of the normal PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 8 ) i8_normal_ab
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
! real ( kind = 8 ) r8_uniform_01
! integer ( kind = 8 ) seed     ! compiler error, type mismatch; corrected MDG 01/02/15
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  r1 = r8_uniform_01 ( seed )
  r2 = r8_uniform_01 ( seed )
  x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * r8_pi * r2 )

  i8_normal_ab = nint ( a + b * x )

  return
end
function r4_normal_01 ( seed )

!*****************************************************************************80
!
!! R4_NORMAL_01 returns a unit pseudonormal R4.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 4 ) R4_NORMAL_01, a sample of the standard
!    normal PDF.
!
  implicit none

  real ( kind = 4 ) r1
  real ( kind = 4 ) r2
  real ( kind = 4 ) r4_normal_01
  real ( kind = 4 ), parameter :: r4_pi = 3.141592653589793E+00
! real ( kind = 4 ) r4_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 4 ) x

  r1 = r4_uniform_01 ( seed )
  r2 = r4_uniform_01 ( seed )
  x = sqrt ( - 2.0E+00 * log ( r1 ) ) * cos ( 2.0E+00 * r4_pi * r2 )

  r4_normal_01 = x

  return
end
function r4_normal_ab ( a, b, seed )

!*****************************************************************************80
!
!! R4_NORMAL_AB returns a scaled pseudonormal R4.
!
!  Discussion:
!
!    The normal probability distribution function (PDF) is sampled,
!    with mean A and standard deviation B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A, the mean of the PDF.
!
!    Input, real ( kind = 4 ) B, the standard deviation of the PDF.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 4 ) R4_NORMAL_AB, a sample of the normal PDF.
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) b
  real ( kind = 4 ) r1
  real ( kind = 4 ) r2
  real ( kind = 4 ) r4_normal_ab
  real ( kind = 4 ), parameter :: r4_pi = 3.141592653589793E+00
! real ( kind = 4 ) r4_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 4 ) x

  r1 = r4_uniform_01 ( seed )
  r2 = r4_uniform_01 ( seed )
  x = sqrt ( - 2.0E+00 * log ( r1 ) ) * cos ( 2.0E+00 * r4_pi * r2 )

  r4_normal_ab = a + b * x

  return
end

function r4_uniform_01 ( seed )

!*****************************************************************************80
!
!! R4_UNIFORM_01 returns a unit pseudorandom R4.
!
!  Discussion:
!
!    An R4 is a real ( kind = 4 ) value.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r4_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R4_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 4 ) R4_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 4 ) r4_uniform_01

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r4_uniform_01 = real ( seed, kind = 4 ) * 4.656612875E-10

  return
end
subroutine r4vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R4VEC_UNIFORM_01 returns a unit pseudorandom R4VEC.
!
!  Discussion:
!
!    An R4VEC is an array of real ( kind = 4 ) values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value,
!    which should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, real ( kind = 4 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 4 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r(i) = real ( seed, kind = 4 ) * 4.656612875E-10

  end do

  return
end
subroutine r4vec_normal_ab ( n, a, b, seed, x )

!*****************************************************************************80
!
!! R4VEC_NORMAL_AB returns a scaled pseudonormal R4VEC.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    An R4VEC is a vector of R4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values desired.
!
!    Input, real ( kind = 4 ) A, B, the mean and standard deviation.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 4 ) X(N), a sample of the standard normal PDF.
!
!  Local parameters:
!
!    Local, real ( kind = 4 ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    Local, integer ( kind = 4 ) X_LO_INDEX, X_HI_INDEX, records the range
!    of entries of X that we need to compute.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 4 ) a
  real ( kind = 4 ) b
  integer ( kind = 4 ) m
  real ( kind = 4 ) r(n+1)
  real ( kind = 4 ), parameter :: r4_pi = 3.141592653589793E+00
! real ( kind = 4 ) r4_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 4 ) x(n)
  integer ( kind = 4 ) x_hi_index
  integer ( kind = 4 ) x_lo_index
!
!  Record the range of X we need to fill in.
!
  x_lo_index = 1
  x_hi_index = n
!
!  If we need just one new value, do that here to avoid null arrays.
!
  if ( x_hi_index - x_lo_index + 1 == 1 ) then

    r(1) = r4_uniform_01 ( seed )

    if ( r(1) == 0.0E+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R4VEC_NORMAL_AB - Fatal error!'
      write ( *, '(a)' ) '  R4_UNIFORM_01 returned a value of 0.'
      stop 1
    end if

    r(2) = r4_uniform_01 ( seed )

    x(x_hi_index) = &
      sqrt ( - 2.0E+00 * log ( r(1) ) ) * cos ( 2.0E+00 * r4_pi * r(2) )
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index, 2 ) == 1 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    call r4vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( - 2.0E+00 * log ( r(1:2*m-1:2) ) ) &
      * cos ( 2.0E+00 * r4_pi * r(2:2*m:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( - 2.0E+00 * log ( r(1:2*m-1:2) ) ) &
      * sin ( 2.0E+00 * r4_pi * r(2:2*m:2) )
!
!  If we require an odd number of values, we generate an even number,
!  and handle the last pair specially, storing one in X(N), and
!  saving the other for later.
!
  else

    x_hi_index = x_hi_index - 1

    m = ( x_hi_index - x_lo_index + 1 ) / 2 + 1

    call r4vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( - 2.0E+00 * log ( r(1:2*m-3:2) ) ) &
      * cos ( 2.0E+00 * r4_pi * r(2:2*m-2:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( - 2.0E+00 * log ( r(1:2*m-3:2) ) ) &
      * sin ( 2.0E+00 * r4_pi * r(2:2*m-2:2) )

    x(n) = sqrt ( - 2.0E+00 * log ( r(2*m-1) ) ) &
      * cos ( 2.0E+00 * r4_pi * r(2*m) )

  end if

  x(1:n) = a + b * x(1:n)

  return
end
function r8_normal_01 ( seed )

!*****************************************************************************80
!
!! R8_NORMAL_01 returns a unit pseudonormal R8.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) R8_NORMAL_01, a normally distributed
!    random value.
!
  implicit none

  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r8_normal_01
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
! real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  r1 = r8_uniform_01 ( seed )
  r2 = r8_uniform_01 ( seed )
  x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * r8_pi * r2 )

  r8_normal_01 = x

  return
end
function r8_normal_ab ( a, b, seed )

!*****************************************************************************80
!
!! R8_NORMAL_AB returns a scaled pseudonormal R8.
!
!  Discussion:
!
!    The normal probability distribution function (PDF) is sampled,
!    with mean A and standard deviation B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the mean of the PDF.
!
!    Input, real ( kind = 8 ) B, the standard deviation of the PDF.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) R8_NORMAL_AB, a sample of the normal PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r8_normal_ab
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
! real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  r1 = r8_uniform_01 ( seed )
  r2 = r8_uniform_01 ( seed )
  x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * r8_pi * r2 )

  r8_normal_ab = a + b * x

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, 1969, pages 136-143.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine r8mat_normal_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R8MAT_NORMAL_01 returns a unit pseudonormal R8MAT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the array.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(M,N), the array of pseudonormal values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(m,n)

  call r8vec_normal_01 ( m * n, seed, r )

  return
end
subroutine r8mat_normal_ab ( m, n, a, b, seed, r )

!*****************************************************************************80
!
!! R8MAT_NORMAL_AB returns a scaled pseudonormal R8MAT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the array.
!
!    Input, real ( kind = 8 ) A, B, the mean and standard deviation.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(M,N), the array of pseudonormal values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(m,n)

  call r8vec_normal_ab ( m * n, a, b, seed, r )

  return
end
subroutine r8vec_normal_01 ( n, seed, x )

!*****************************************************************************80
!
!! R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values desired.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) X(N), a sample of the standard normal PDF.
!
!  Local parameters:
!
!    Local, real ( kind = 8 ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    Local, integer ( kind = 4 ) X_LO_INDEX, X_HI_INDEX, records the range
!    of entries of X that we need to compute
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) m
  real ( kind = 8 ) r(n+1)
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
! real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)
  integer ( kind = 4 ) x_hi_index
  integer ( kind = 4 ) x_lo_index
!
!  Record the range of X we need to fill in.
!
  x_lo_index = 1
  x_hi_index = n
!
!  If we need just one new value, do that here to avoid null arrays.
!
  if ( x_hi_index - x_lo_index + 1 == 1 ) then

    r(1) = r8_uniform_01 ( seed )

    if ( r(1) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC_NORMAL_01 - Fatal error!'
      write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
      stop 1
    end if

    r(2) = r8_uniform_01 ( seed )

    x(x_hi_index) = &
      sqrt ( - 2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * r8_pi * r(2) )
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index, 2 ) == 1 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    call r8vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2:2*m:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * sin ( 2.0D+00 * r8_pi * r(2:2*m:2) )
!
!  If we require an odd number of values, we generate an even number,
!  and handle the last pair specially, storing one in X(N), and
!  saving the other for later.
!
  else

    x_hi_index = x_hi_index - 1

    m = ( x_hi_index - x_lo_index + 1 ) / 2 + 1

    call r8vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2:2*m-2:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * sin ( 2.0D+00 * r8_pi * r(2:2*m-2:2) )

    x(n) = sqrt ( - 2.0D+00 * log ( r(2*m-1) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2*m) )

  end if

  return
end
subroutine r8vec_normal_ab ( n, a, b, seed, x )

!*****************************************************************************80
!
!! R8VEC_NORMAL_AB returns a scaled pseudonormal R8VEC.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values desired.
!
!    Input, real ( kind = 8 ) A, B, the mean and standard deviation.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) X(N), a sample of the standard normal PDF.
!
!  Local parameters:
!
!    Local, real ( kind = 8 ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    Local, integer ( kind = 4 ) X_LO_INDEX, X_HI_INDEX, records the range
!    of entries of X that we need to compute. 
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) m
  real ( kind = 8 ) r(n+1)
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
! real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)
  integer ( kind = 4 ) x_hi_index
  integer ( kind = 4 ) x_lo_index
!
!  Record the range of X we need to fill in.
!
  x_lo_index = 1
  x_hi_index = n
!
!  If we need just one new value, do that here to avoid null arrays.
!
  if ( x_hi_index - x_lo_index + 1 == 1 ) then

    r(1) = r8_uniform_01 ( seed )

    if ( r(1) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC_NORMAL - Fatal error!'
      write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
      stop 1
    end if

    r(2) = r8_uniform_01 ( seed )

    x(x_hi_index) = &
      sqrt ( - 2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * r8_pi * r(2) )
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index, 2 ) == 1 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    call r8vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2:2*m:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * sin ( 2.0D+00 * r8_pi * r(2:2*m:2) )
!
!  If we require an odd number of values, we generate an even number,
!  and handle the last pair specially, storing one in X(N), and
!  saving the other for later.
!
  else

    x_hi_index = x_hi_index - 1

    m = ( x_hi_index - x_lo_index + 1 ) / 2 + 1

    call r8vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2:2*m-2:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * sin ( 2.0D+00 * r8_pi * r(2:2*m-2:2) )

    x(n) = sqrt ( - 2.0D+00 * log ( r(2*m-1) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2*m) )

  end if

  x(1:n) = a + b * x(1:n)

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end







end module math
