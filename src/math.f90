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
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief collection of mathematical/numerical routines that don't fit anywhere else
! 
!> @date   10/13/98 MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   03/19/13 MDG 3.0 updated all routines
!
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
!> @author Marc De Graef, Carnegie Melon University
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
!> @author Marc De Graef, Carnegie Melon University
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

end module math
