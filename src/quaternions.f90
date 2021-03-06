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
! CTEMsoft2013:quaternions.f90
!--------------------------------------------------------------------------
!
! MODULE: quaternions
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief module with basic quaternion functions (some overloaded operators)
!                    
!> @details   [verified against Mathematica Quaternion package on 3/15/12]
!>
!> REMEMBER THAT QUATERNION MULTIPLICATION IS NON-COMMUTATIVE !!!!!
!>
!> quaternions are defined as arrays of 4 single or double precision reals;
!> the first entry is the scalar part, the remaining three form the vector part.
!>
!> If you want to try it out, here is an example test program\n
!>\n
!> program qtest\n
!>\n
!>use local\n
!>use quaternions\n
!>\n
!>IMPLICIT NONE\n
!>\n
!>! for single precision, use the following lines\n
!>!real(kind=sgl)  ::  u(4), v(4), w(4)\n
!>!real(kind=sgl) :: x, a=2\n
!>! define two quaternions (single)\n
!>!u = (/1.0,2.0,3.0,4.0/)\n
!>!v = (/5.0,6.0,7.0,8.0/)\n
!>\n
!>! for double precision, uncomment the next set and comment the previous set lines\n
!>real(kind=dbl)  ::  u(4), v(4), w(4)\n
!>real(kind=dbl) :: x, a=2.D0\n
!>! double\n
!>u = (/1.D0,2.D0,3.D0,4.D0/)\n
!>v = (/5.D0,6.D0,7.D0,8.D0/)\n
!>\n
!>\n
!>write (stdout,*) ' quaternion u '\n
!>call quaternion_print(u)\n
!>write (stdout,*) ' quaternion v '\n
!>call quaternion_print(v)\n
!>\n
!>! next, do all the operations to make sure that they are correct\n
!>\n
!>write (stdout,*) '   addition u+v '\n
!>w = u+v\n
!>call quaternion_print(w)\n
!>write (stdout,*) ' '\n
!>\n
!>write (stdout,*) '   subtraction u-v '\n
!>w = u-v\n
!>call quaternion_print(w)\n
!>write (stdout,*) ' '\n
!>\n
!>write (stdout,*) '   scalar multiplication (both orderings)  '\n
!>w = a*u\n
!>call quaternion_print(w)\n
!>w = u*a\n
!>call quaternion_print(w)\n
!>write (stdout,*) ' '\n
!>\n
!>write (stdout,*) '   multiplication uv '\n
!>w = quat_mult(u,v)\n
!>call quaternion_print(w)\n
!>write (stdout,*) ' '\n
!>\n
!write (stdout,*) '   conjugate u '\n
!>w = conjg(u)\n
!>call quaternion_print(w)\n
!>write (stdout,*) ' '\n
!>\n
!>write (stdout,*) '   norm(u) '\n
!>x = cabs(u)\n
!>write (stdout,*) x\n
!>write (stdout,*) ' '\n
!>\n
!>write (stdout,*) '   division u/v '\n
!>w = quat_div(u,v)\n
!>call quaternion_print(w)\n
!>write (stdout,*) ' '\n
!>\n
!>end program qtest\n
!>
!
!> @note Quaternions are defined with the scalar part in position 1, and the 
!> vector part in positions 2:4.
! 
!> @date 3/15/12   MDG 1.0 original
!> @date 8/04/13   MDG 1.1 moved rotation conversion functions to rotations.f90
!> @date 8/12/13   MDG 2.0 re-defined quaternions to be arrays of 4 reals rather than by letter ... 
!--------------------------------------------------------------------------
module quaternions

use local

IMPLICIT NONE

! only this routine is public; all the others are done via overloaded operators
public :: quaternion_print
interface quaternion_print
	module procedure quaternion_print
	module procedure quaternion_print_d
 end interface

! quaternion multiplication (single and double precision)
public :: quat_mult
interface quat_mult
     module procedure quat_mult
     module procedure quat_mult_d
  end interface

! complex conjugation (single and double precision)
intrinsic :: conjg
public :: conjg
interface conjg
     module procedure quat_conjg
     module procedure quat_conjg_d
  end interface

! quaternion norm (single and double precision)
intrinsic :: cabs
public :: cabs
interface cabs
     module procedure quat_norm
     module procedure quat_norm_d
  end interface

! quaternion division (single and double precision)
public :: quat_div
interface quat_div
     module procedure quat_div
     module procedure quat_div_d
  end interface

! quaternion inner product (single and double precision)
public :: quat_innerproduct
interface quat_innerproduct
     module procedure quat_innerproduct
     module procedure quat_innerproduct_d
  end interface

! interquaternion angle (single and double precision)
public :: quat_angle
interface quat_angle
     module procedure quat_angle
     module procedure quat_angle_d
  end interface

! quaternion rotation of a unit vector  q v q-1
public :: quat_rotate_vector
interface quat_rotate_vector
	module procedure quat_rotate_vector
	module procedure quat_rotate_vector_d
end interface




contains

!--------------------------------------------------------------------------
!
! SUBROUTINE: quaternion_print
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief print a quaternion (for debugging purposes mostly)
!
!> @param q quaternion to be printed (single precision)  
! 
!> @date 3/15/12   MDG 1.0 original
!> @date 8/12/13   MDG 2.0 rewrite
!--------------------------------------------------------------------------
subroutine quaternion_print(q)

use local
use io

    real(kind=sgl), intent(in) 	:: q(4)		!< input quaternion (single precision)

    call WriteValue('', q, 4, "('(',4f12.6,')')")

end subroutine quaternion_print

!--------------------------------------------------------------------------
!
! SUBROUTINE:quaternion_print_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief print a quaternion (for debugging purposes mostly)
!
!> @param q quaternion to be printed (double precision)  
! 
!> @date 3/15/12   MDG 1.0 original
!> @date 8/12/13   MDG 2.0 rewrite
!--------------------------------------------------------------------------
subroutine quaternion_print_d(q)

use local
use io

    real(kind=dbl), intent(in) 	:: q(4)		!< input quaternion (double precision)

    call WriteValue('', q, 4, "('(',4f15.9,')')")

end subroutine quaternion_print_d

!--------------------------------------------------------------------------
!
! FUNCTION: quat_mult
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief quaternion multiplication   (single precision)
!
!> @param x first quaternion 
!> @param y second quaternion 
! 
!> @date 3/15/12   MDG 1.0 original
!> @date 8/12/13   MDG 2.0 rewrite
!--------------------------------------------------------------------------
function quat_mult(x,y) result (res)

    real(kind=sgl), intent(in) 		:: x(4), y(4)		!< input quaternions
    real(kind=sgl) 				:: res(4)

    res = (/ x(1)*y(1) - x(2)*y(2) - x(3)*y(3) - x(4)*y(4), &
             x(1)*y(2) + x(2)*y(1) + x(3)*y(4) - x(4)*y(3), &
             x(1)*y(3) - x(2)*y(4) + x(3)*y(1) + x(4)*y(2), &
             x(1)*y(4) + x(2)*y(3) - x(3)*y(2) + x(4)*y(1) /)
    
end function quat_mult

!--------------------------------------------------------------------------
!
! FUNCTION: quat_mult_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief quaternion multiplication   (double precision)
!
!> @param x first quaternion 
!> @param y second quaternion 
! 
!> @date 3/15/12   MDG 1.0 original
!> @date 8/12/13   MDG 2.0 rewrite
!--------------------------------------------------------------------------
function quat_mult_d(x,y) result (res)

    real(kind=dbl), intent(in) 		:: x(4), y(4)		!< input quaternions
    real(kind=dbl) 				:: res(4)

    res = (/ x(1)*y(1) - x(2)*y(2) - x(3)*y(3) - x(4)*y(4), &
             x(1)*y(2) + x(2)*y(1) + x(3)*y(4) - x(4)*y(3), &
             x(1)*y(3) - x(2)*y(4) + x(3)*y(1) + x(4)*y(2), &
             x(1)*y(4) + x(2)*y(3) - x(3)*y(2) + x(4)*y(1) /)

end function quat_mult_d

!--------------------------------------------------------------------------
!
! FUNCTION: quat_conjg
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief quaternion complex conjugation (extends intrinsic routine conjg)
!
!> @param x quaternion to be conjugated (single precision)
! 
!> @date 3/15/12   MDG 1.0 original
!> @date 8/12/13   MDG 2.0 rewrite
!--------------------------------------------------------------------------
function quat_conjg(x) result (res)

    real(kind=sgl), intent(in) 		:: x(4)		!< input quaternion
    real(kind=sgl) 				:: res(4)

    res = (/ x(1), -x(2), -x(3), -x(4) /)

end function quat_conjg

!--------------------------------------------------------------------------
!
! FUNCTION: quat_conjg_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief quaternion complex conjugation (extends intrinsic routine conjg)
!
!> @param x quaternion to be conjugated (double precision)
! 
!> @date 3/15/12   MDG 1.0 original
!> @date 8/12/13   MDG 2.0 rewrite
!--------------------------------------------------------------------------
function quat_conjg_d(x) result (res)

    real(kind=dbl), intent(in) 		:: x(4)		!< input quaternion
    real(kind=dbl) 				:: res(4)

    res = (/ x(1), -x(2), -x(3), -x(4) /)

end function quat_conjg_d


!--------------------------------------------------------------------------
!
! FUNCTION: quat_norm
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief quaternion norm (extends intrinsic routine cabs)
!
!> @param x quaternion to be normed (single precision)
! 
!> @date 3/15/12   MDG 1.0 original
!> @date 8/12/13   MDG 2.0 rewrite
!--------------------------------------------------------------------------
function quat_norm(x) result (res)

use local

   real(kind=sgl), intent(in) 		:: x(4)		!< input quaternion
   real(kind=sgl) 			:: res

    res =  sqrt( sum(x*x) )

end function quat_norm

!--------------------------------------------------------------------------
!
! FUNCTION: quat_norm_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief quaternion norm (extends intrinsic routine cabs)
!
!> @param x quaternion to be normed (double precision)
! 
!> @date 3/15/12   MDG 1.0 original
!> @date 8/12/13   MDG 2.0 rewrite
!--------------------------------------------------------------------------
function quat_norm_d(x) result (res)

use local

    real(kind=dbl), intent(in) 	:: x(4)		!< input quaternion
    real(kind=dbl) 			:: res

    res =  dsqrt( sum(x*x) )

end function quat_norm_d

!--------------------------------------------------------------------------
!
! FUNCTION: quat_div
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief quaternion division (single precision)
!
!> @param x nominator quaternion 
!> @param y denominator quaternion 
! 
!> @date 3/15/12   MDG 1.0 original
!> @date 8/12/13   MDG 2.0 rewrite
!--------------------------------------------------------------------------
function quat_div(x,y) result (res)

    real(kind=sgl), intent(in) 		:: x(4),y(4)		!< input quaternions
    real(kind=sgl) 				:: res(4), p(4), q

    q = quat_norm(y)
    p = quat_conjg(y)/(q*q)
    res =  quat_mult(x,p)

end function quat_div

!--------------------------------------------------------------------------
!
! FUNCTION: quat_div_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief quaternion division (double precision)
!
!> @param x nominator quaternion 
!> @param y denominator quaternion 
! 
!> @date 3/15/12   MDG 1.0 original
!> @date 8/12/13   MDG 2.0 rewrite
!--------------------------------------------------------------------------
function quat_div_d(x,y) result (res)

    real(kind=dbl), intent(in) 		:: x(4),y(4)		!< input quaternions
    real(kind=dbl) 				:: res(4), p(4), q

    q = quat_norm_d(y)
    p = quat_conjg_d(y)/(q*q)
    res =  quat_mult_d(x,p)

end function quat_div_d

!--------------------------------------------------------------------------
!
! FUNCTION: quat_innerproduct
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  quaternion inner product  (single precision)
!
!> @param x nominator quaternion 
!> @param y denominator quaternion 
! 
!> @date 3/15/12   MDG 1.0 original
!> @date 8/12/13   MDG 2.0 rewrite
!--------------------------------------------------------------------------
function quat_innerproduct(x,y) result (res)

use local

    real(kind=sgl), intent(in) 		:: x(4),y(4)		!< input quaternions
    real(kind=sgl) 				:: res

    res = sum(x * y)

end function quat_innerproduct

!--------------------------------------------------------------------------
!
! FUNCTION: quat_innerproduct_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  quaternion inner product  (double precision)
!
!> @param x nominator quaternion 
!> @param y denominator quaternion 
! 
!> @date 3/15/12   MDG 1.0 original
!> @date 8/12/13   MDG 2.0 rewrite
!--------------------------------------------------------------------------
function quat_innerproduct_d(x,y) result (res)

use local

    real(kind=dbl), intent(in) 		:: x(4),y(4)		!< input quaternions
    real(kind=dbl) 				:: res

    res = sum(x * y) 

end function quat_innerproduct_d

!--------------------------------------------------------------------------
!
! FUNCTION: quat_angle
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief   interquaternion angle   (single precision)
!
!> @param x first quaternion 
!> @param y second quaternion 
! 
!> @date 3/15/12   MDG 1.0 original
!> @date 8/12/13   MDG 2.0 rewrite
!--------------------------------------------------------------------------!
function quat_angle(x,y) result (res)

use local

    real(kind=sgl), intent(in) 			:: x(4),y(4)		!< input quaternions
    real(kind=sgl) 					:: res, q

    q = quat_innerproduct(x,y)
    res = acos( 2.0_sgl*q*q - 1.0_sgl )

end function quat_angle

!--------------------------------------------------------------------------
!
! FUNCTION: quat_angle_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief   interquaternion angle   (double precision)
!
!> @param x first quaternion 
!> @param y second quaternion 
! 
!> @date 3/15/12   MDG 1.0 original
!> @date 8/12/13   MDG 2.0 rewrite
!--------------------------------------------------------------------------!
function quat_angle_d(x,y) result (res)

use local

    real(kind=dbl), intent(in) 		:: x(4),y(4)		!< input quaternions
    real(kind=dbl) 				:: res, q

    q = quat_innerproduct_d(x,y)
    res = dacos( 2.0_dbl*q*q - 1.0_dbl )

end function quat_angle_d

!--------------------------------------------------------------------------
!
! FUNCTION: quat_rotate_vector
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief   rotate a unit vector by a unit quaternion
!
!> @param q quaternion 
!> @param v vector to be rotated 
! 
!> @date 3/15/12   MDG 1.0 original
!> @date 8/12/13   MDG 2.0 rewrite
!--------------------------------------------------------------------------!
recursive function quat_rotate_vector(q,v) result (res)

use local

    real(kind=sgl),intent(in) 			:: q(4)		!< input quaternion
    real(kind=sgl),intent(in)			:: v(3)		!< input vector (must be normalized)
    real(kind=sgl)				:: qv(4), rqv(4)
    real(kind=sgl) 				:: res(3)

    qv = (/ 0.0, v(1), v(2), v(3) /)   
    rqv = quat_mult(q,quat_mult(qv,quat_conjg(q) ) )
    res = rqv(2:4)

end function quat_rotate_vector

!--------------------------------------------------------------------------
!
! FUNCTION: quat_rotate_vector_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief   rotate a unit vector by a unit quaternion (double precision)
!
!> @param q quaternion 
!> @param v vector to be rotated 
! 
!> @date 3/15/12   MDG 1.0 original
!> @date 8/12/13   MDG 2.0 rewrite
!--------------------------------------------------------------------------!
recursive function quat_rotate_vector_d(q,v) result (res)

use local

    real(kind=dbl), intent(in) 		:: q(4)		!< input quaternion
    real(kind=dbl),intent(in)			:: v(3)		!< input vector (must be normalized)
    real(kind=dbl)				:: qv(4), rqv(4)
    real(kind=dbl) 				:: res(3)

    qv = (/ 0.D0, v(1), v(2), v(3) /)   
    rqv = quat_mult_d(q,quat_mult_d(qv,quat_conjg_d(q) ) )
    res = rqv(2:4)

end function quat_rotate_vector_d




end module quaternions
