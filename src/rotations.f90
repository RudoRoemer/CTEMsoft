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
! CTEMsoft2013:rotations.f90
!--------------------------------------------------------------------------
!
! MODULE: rotations
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief everything that has to do with rotations and conversions between rotations
!
!> @details This file relies a lot on the relations listed in the book "Orientations
!> and Rotations" by Adam Morawiec [Springer 2004].  I've tried to implement every
!> available representation for rotations in a way that makes it easy to convert 
!> between any pair.  Needless to say, this needs extensive testing and debugging...
!>
!> Instead of converting all the time between representations, I've opted to 
!> "waste" a little more memory and time and precompute all the representations.
!> This way all representations are available via a single data structure.
!>
!> Obviously, the individual conversion routines also exist and can be called either in
!> single or in double precision (using a function interface for each call, so that only
!> one function name is used).  The conversion routines use the following format for their
!> call name:  ab2cd, where (ab and cd are two-characters strings selected from the following
!> possiblities: [the number in parenthesis lists the number of entries that need to be provided] 
!>
!> eu : euler angle representation (3)
!> om : orientation matrix representation (3x3)
!> ax : axis angle representation (4)
!> ro : Rodrigues vector representation (3)
!> qu : unit quaternion representation (4)
!> ho : homochoric representation (3)
!> cu : cubochoric representation (3).
!>
!> hence, conversion from homochoric to euler angle is called as ho2eu(); the argument of 
!> each routine must have the correct number of dimensions and entries.
!> All 42 conversion routines exist.

!> @date 8/04/13   MDG 1.0 original
!--------------------------------------------------------------------------
module rotations

use local
use quaternions

! the "orientation" type contains entries for all rotation and orientation representations
type orientationtype
  real(kind=sgl)	:: eulang(3)		! Bunge Euler angles in radians
  real(kind=sgl)	:: om(3,3)		! 3x3 matrix
  real(kind=sgl)	:: axang(4)		! axis-angle pair (angle in rad, component 4; axis in direction cosines)
  real(kind=sgl)	:: rodrigues(3)		! Rodrigues vector
  real(kind=sgl)	:: quat(4)		! quaternion representation (q(1) is scalar part, q(2:4) vector part)
  real(kind=sgl)	:: homochoric(3)	! homochoric representation according to Frank's paper  
  real(kind=sgl)	:: cubochoric(3)	! cubic grid representation (derived from homochoric)
end type orientationtype


! double precision version
type orientationtyped
  real(kind=dbl)	:: eulang(3)		! Bunge Euler angles in radians
  real(kind=dbl)	:: om(3,3)		! 3x3 matrix
  real(kind=dbl)	:: axang(4)		! axis-angle pair (angle in rad, component 4; axis in direction cosines)
  real(kind=dbl)	:: rodrigues(3)		! Rodrigues vector
  real(kind=dbl)	:: quat(4)		! quaternion representation (q(1) is scalar part, q(2:4) vector part)
  real(kind=dbl)	:: homochoric(3)	! homochoric representation according to Frank's paper  
  real(kind=dbl)	:: cubochoric(3)	! cubic grid representation (derived from homochoric)
end type orientationtyped



! general interface routine to populate the orientation type
public:: init_orientation
interface init_orientation
	module procedure init_orientation
	module procedure init_orientation_om
	module procedure init_orientation_d
	module procedure init_orientation_om_d
end interface


! convert Euler angles to 3x3 orientation matrix
public :: eu2om
interface eu2om
	module procedure eu2om
	module procedure eu2om_d
end interface

! convert Euler angles to axis angle
public :: eu2ax
interface eu2ax
	module procedure eu2ax
	module procedure eu2ax_d
end interface

! convert Euler angles to Rodrigues vector
public :: eu2ro
interface eu2ro
	module procedure eu2ro
	module procedure eu2ro_d
end interface

! convert Euler angles to quaternion
public :: eu2qu
interface eu2qu
	module procedure eu2qu
	module procedure eu2qu_d
end interface

! convert Euler angles to homochoric
public :: eu2ho
interface eu2ho
	module procedure eu2ho
	module procedure eu2ho_d
end interface

! convert Euler angles to cubochoric
public :: eu2cu
interface eu2cu
	module procedure eu2cu
	module procedure eu2cu_d
end interface

!--------------------------------
! convert 3x3 orientation matrix to Euler angles
public :: om2eu
interface om2eu
	module procedure om2eu
	module procedure om2eu_d
end interface

! convert 3x3 orientation matrix to axis angle
public :: om2ax
interface om2ax
	module procedure om2ax
	module procedure om2ax_d
end interface

! convert 3x3 orientation matrix to Rodrigues
public :: om2ro
interface om2ro
	module procedure om2ro
	module procedure om2ro_d
end interface

! convert 3x3 rotation matrix to quaternion
public :: om2qu
interface om2qu
	module procedure om2qu
	module procedure om2qu_d
end interface

! convert 3x3 rotation matrix to homochoric
public :: om2ho
interface om2ho
	module procedure om2ho
	module procedure om2ho_d
end interface

! convert 3x3 rotation matrix to cubochoric
public :: om2cu
interface om2cu
	module procedure om2cu
	module procedure om2cu_d
end interface

!--------------------------------
! convert axis angle pair to euler
public :: ax2eu
interface ax2eu
	module procedure ax2eu
	module procedure ax2eu_d
end interface

! convert axis angle pair to orientation matrix
public :: ax2om
interface ax2om
	module procedure ax2om
	module procedure ax2om_d
end interface

! convert axis angle pair to Rodrigues
public :: ax2ro
interface ax2ro
	module procedure ax2ro
	module procedure ax2ro_d
end interface

! convert axis angle pair to quaternion
public :: ax2qu
interface ax2qu
	module procedure ax2qu
	module procedure ax2qu_d
end interface

! convert axis angle pair to homochoric representation
public :: ax2ho
interface ax2ho
	module procedure ax2ho
	module procedure ax2ho_d
end interface

! convert axis angle pair to cubochoric
public :: ax2cu
interface ax2cu
	module procedure ax2cu
	module procedure ax2cu_d
end interface

!--------------------------------
! convert Rodrigues vector to Euler angles
public :: ro2eu
interface ro2eu
	module procedure ro2eu
	module procedure ro2eu_d
end interface

! convert Rodrigues vector to orientation matrix
public :: ro2om
interface ro2om
	module procedure ro2om
	module procedure ro2om_d
end interface

! convert Rodrigues vector to axis angle pair
public :: ro2ax
interface ro2ax
	module procedure ro2ax
	module procedure ro2ax_d
end interface

! convert Rodrigues vector to quaternion
public :: ro2qu
interface ro2qu
	module procedure ro2qu
	module procedure ro2qu_d
end interface

! convert Rodrigues vector to homochoric
public :: ro2ho
interface ro2ho
	module procedure ro2ho
	module procedure ro2ho_d
end interface

! convert Rodrigues vector to cubochoric
public :: ro2cu
interface ro2cu
	module procedure ro2cu
	module procedure ro2cu_d
end interface

!--------------------------------
! convert quaternion to Euler angles
public :: qu2eu
interface qu2eu
	module procedure qu2eu
	module procedure qu2eu_d
end interface

! convert quaternion to orientation matrix 
public :: qu2om
interface qu2om
	module procedure qu2om
	module procedure qu2om_d
end interface

! convert quaternion to axis angle
public :: qu2ax
interface qu2ax
	module procedure qu2ax
	module procedure qu2ax_d
end interface

! convert quaternion to Rodrigues
public :: qu2ro
interface qu2ro
	module procedure qu2ro
	module procedure qu2ro_d
end interface

! convert quaternion to homochoric
public :: qu2ho
interface qu2ho
	module procedure qu2ho
	module procedure qu2ho_d
end interface

! convert quaternion to cubochoric
public :: qu2cu
interface qu2cu
	module procedure qu2cu
	module procedure qu2cu_d
end interface

!--------------------------------
! convert homochoric to euler
public :: ho2eu
interface ho2eu
	module procedure ho2eu
	module procedure ho2eu_d
end interface

! convert homochoric to orientation matrix 
public :: ho2om
interface ho2om
	module procedure ho2om
	module procedure ho2om_d
end interface

! convert homochoric to axis angle pair 
public :: ho2ax
interface ho2ax
	module procedure ho2ax
	module procedure ho2ax_d
end interface

! convert homochoric to Rodrigues 
public :: ho2ro
interface ho2ro
	module procedure ho2ro
	module procedure ho2ro_d
end interface

! convert homochoric to quaternion
public :: ho2qu
interface ho2qu
	module procedure ho2qu
	module procedure ho2qu_d
end interface

! convert homochoric to cubochoric
public :: ho2cu
interface ho2cu
	module procedure ho2cu
	module procedure ho2cu_d
end interface

!--------------------------------
! convert cubochoric to euler
public :: cu2eu
interface cu2eu
	module procedure cu2eu
	module procedure cu2eu_d
end interface

! convert cubochoric to orientation matrix
public :: cu2om
interface cu2om
	module procedure cu2om
	module procedure cu2om_d
end interface

! convert cubochoric to axis angle
public :: cu2ax
interface cu2ax
	module procedure cu2ax
	module procedure cu2ax_d
end interface

! convert cubochoric to Rodrigues
public :: cu2ro
interface cu2ro
	module procedure cu2ro
	module procedure cu2ro_d
end interface

! convert cubochoric to quaternion
public :: cu2qu
interface cu2qu
	module procedure cu2qu
	module procedure cu2qu_d
end interface

! convert cubochoric to homochoric
public :: cu2ho
interface cu2ho
	module procedure cu2ho
	module procedure cu2ho_d
end interface




!--------------------------------
! print quaternion and equivalent 3x3 rotation matrix
public :: print_orientation
interface print_orientation
	module procedure print_orientation
	module procedure print_orientation_d
end interface



contains


!--------------------------------------------------------------------------
!
! Function: init_orientation
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief take an orientation representation with 3 component and init all others (single precision)
!
!> @param orient 3-component vector (single precision)  
!> @param intype input type ['eu', 'ro', 'ho', 'cu']
! 
!> @date 8/04/13   MDG 1.0 original
!--------------------------------------------------------------------------
function init_orientation(orient,intype) result(res)

use local

IMPLICIT NONE

real(kind=sgl),INTENT(IN)	:: orient(*)		! 3 or 4 component orientation descriptor
character(2),INTENT(IN)	:: intype		! describes input type

type(orientationtype)		:: res

real(kind=sgl)			:: r,neworient(4)

select case (intype)
	case ('eu')	! Euler angles
		res%eulang = orient(1:3)
		res%om = eu2om(orient(1:3))
		res%quat = eu2qu(orient(1:3))
		res%rodrigues = eu2ro(orient(1:3))
		res%axang = ro2ax(res%rodrigues)
		res%homochoric = ax2ho(res%axang)
		res%cubochoric = ho2cu(res%homochoric)
	case ('ro') 	! Rodrigues vector
		res%rodrigues = orient(1:3)
		res%eulang = ro2eu(orient(1:3))
		res%om = eu2om(res%eulang)
		res%quat = eu2qu(res%eulang)
		res%axang = ro2ax(res%rodrigues)
		res%homochoric = ax2ho(res%axang)
		res%cubochoric = ho2cu(res%homochoric)
	case ('ho')	! homochoric
		res%homochoric = orient(1:3)
		res%axang = ho2ax(orient(1:3))
		res%om = ax2om(res%axang)
		res%eulang = om2eu(res%om)
		res%rodrigues = eu2ro(res%eulang)
		res%quat = eu2qu(res%eulang)
		res%cubochoric = ho2cu(res%homochoric)
	case ('cu')	! cubochoric
		res%cubochoric = orient(1:3)
		res%homochoric = cu2ho(res%cubochoric)
		res%eulang = cu2eu(res%cubochoric)
		res%om = cu2om(res%cubochoric)
		res%quat = cu2qu(res%cubochoric)
		res%axang = cu2ax(res%cubochoric)
		res%rodrigues = cu2ro(res%cubochoric)
	case ('qu')	! quaternion
		res%quat = orient(1:4)
		res%eulang = qu2eu(res%quat)
		res%om = eu2om(res%eulang)
		res%rodrigues = eu2ro(res%eulang)
		res%axang = ro2ax(res%rodrigues)
		res%homochoric = ax2ho(res%axang)
		res%cubochoric = ho2cu(res%homochoric)
	case ('ax') 	! axis angle pair
		r = sqrt(sum(orient(1:3)**2))				! normalize, just in case...
		neworient(1:3) = orient(1:3)/r
		neworient(4) = orient(4)
		res%axang = neworient
		res%om = ax2om(neworient)
		res%eulang = om2eu(res%om)
		res%rodrigues = eu2ro(res%eulang)
		res%quat = eu2qu(res%eulang)
		res%homochoric = ax2ho(res%axang)
		res%cubochoric = ho2cu(res%homochoric)
end select 

end function init_orientation

!--------------------------------------------------------------------------
!
! Function: init_orientation_3_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief take an orientation representation with 3 component and init all others (double precision)
!
!> @param orient 3-component vector (double precision)  
!> @param intype input type ['eu', 'ro', 'ho', 'cu']
! 
!> @date 8/04/13   MDG 1.0 original
!--------------------------------------------------------------------------
function init_orientation_d(orient,intype) result(res)

use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN)	:: orient(*)		! 3-component orientation descriptor
character(2),INTENT(IN)	:: intype		! describes input type

type(orientationtyped)		:: res

real(kind=dbl)			:: r, neworient(4)

select case (intype)
	case ('eu')	! Euler angles
		res%eulang = orient(1:3)
		res%om = eu2om_d(orient(1:3))
		res%quat = eu2qu_d(orient(1:3))
		res%rodrigues = eu2ro_d(orient(1:3))
		res%axang = ro2ax_d(res%rodrigues)
		res%homochoric = ax2ho_d(res%axang)
		res%cubochoric = ho2cu_d(res%homochoric)
	case ('ro') 	! Rodrigues vector
		res%rodrigues = orient(1:3)
		res%eulang = ro2eu_d(orient(1:3))
		res%om = eu2om_d(res%eulang)
		res%quat = eu2qu_d(res%eulang)
		res%axang = ro2ax_d(res%rodrigues)
		res%homochoric = ax2ho_d(res%axang)
		res%cubochoric = ho2cu_d(res%homochoric)
	case ('ho')	! homochoric
		res%homochoric = orient(1:3)
		res%axang = ho2ax_d(orient(1:3))
		res%om = ax2om_d(res%axang)
		res%eulang = om2eu_d(res%om)
		res%rodrigues = eu2ro_d(res%eulang)
		res%quat = eu2qu_d(res%eulang)
		res%cubochoric = ho2cu_d(res%homochoric)
	case ('cu')	! cubochoric
		res%cubochoric = orient(1:3)
		res%homochoric = cu2ho_d(res%cubochoric)
		res%eulang = cu2eu_d(res%cubochoric)
		res%om = cu2om_d(res%cubochoric)
		res%quat = cu2qu_d(res%cubochoric)
		res%axang = cu2ax_d(res%cubochoric)
		res%rodrigues = cu2ro_d(res%cubochoric)
	case ('qu')	! quaternion
		res%quat = orient(1:4)
		res%eulang = qu2eu_d(res%quat)
		res%om = eu2om_d(res%eulang)
		res%rodrigues = eu2ro_d(res%eulang)
		res%axang = ro2ax_d(res%rodrigues)
		res%homochoric = ax2ho_d(res%axang)
		res%cubochoric = ho2cu_d(res%homochoric)
	case ('ax') 	! axis angle pair
		r = dsqrt(sum(orient(1:3)**2))				! normalize, just in case...
		neworient(1:3) = orient(1:3)/r
		neworient(4) = orient(4)
		res%axang = neworient
		res%om = ax2om_d(neworient)
		res%eulang = om2eu_d(res%om)
		res%rodrigues = eu2ro_d(res%eulang)
		res%quat = eu2qu_d(res%eulang)
		res%homochoric = ax2ho_d(res%axang)
		res%cubochoric = ho2cu_d(res%homochoric)
end select 

end function init_orientation_d

!--------------------------------------------------------------------------
!
! Function: init_orientation_om
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief take an orientation representation with 3x3 components and init all others (single precision)
!
!> @param orient r-component vector (single precision)  
!> @param intype input type ['om']
! 
!> @date 8/04/13   MDG 1.0 original
!--------------------------------------------------------------------------
function init_orientation_om(orient,intype) result(res)

use local

IMPLICIT NONE

real(kind=sgl),INTENT(IN)	:: orient(3,3)		! 3x3-component orientation descriptor
character(2),INTENT(IN)	:: intype		! describes input type

type(orientationtype)		:: res

select case (intype)
	case ('om')	! orientation matrix
		res%om = orient
		res%quat = om2qu(orient)
		res%eulang = qu2eu(res%quat)
		res%rodrigues = eu2ro(res%eulang)
		res%axang = ro2ax(res%rodrigues)
		res%homochoric = ax2ho(res%axang)
		res%cubochoric = ho2cu(res%homochoric)
end select 

end function init_orientation_om


!--------------------------------------------------------------------------
!
! Function: init_orientation_om_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief take an orientation representation with 3x3 components and init all others (double precision)
!
!> @param orient r-component vector (double precision)  
!> @param intype input type ['om']
! 
!> @date 8/04/13   MDG 1.0 original
!--------------------------------------------------------------------------
function init_orientation_om_d(orient,intype) result(res)

use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN)	:: orient(3,3)		! 3x3-component orientation descriptor
character(2),INTENT(IN)	:: intype		! describes input type

type(orientationtyped)		:: res

select case (intype)
	case ('om')	! orientation matrix
		res%om = orient
		res%quat = om2qu_d(orient)
		res%eulang = qu2eu_d(res%quat)
		res%rodrigues = eu2ro_d(res%eulang)
		res%axang = ro2ax_d(res%rodrigues)
		res%homochoric = ax2ho_d(res%axang)
		res%cubochoric = ho2cu_d(res%homochoric)
end select 

end function init_orientation_om_d





!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! here we start with a series of conversion routines between representations
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------


!--------------------------------------------------------------------------
!
! Function: om2eu
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief orientation matrix to euler angles (single precision)
!
!> @note verified 8/5/13.
!
!> @param o orientation matrix (single precision)  
!
!> @date 8/04/13   MDG 1.0 original
!--------------------------------------------------------------------------
function om2eu(o) result(res)

use local

IMPLICIT NONE

real(kind=sgl),INTENT(IN)	:: o(3,3)		!< orientation matrix
real(kind=sgl)			:: res(3)

res = qu2eu(om2qu(o))

end function om2eu

!--------------------------------------------------------------------------
!
! Function: om2eu_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief orientation matrix to euler angles (double precision)
!
!> @note verified 8/5/13.
!
!> @param o orientation matrix (double precision)  
!
!> @date 8/04/13   MDG 1.0 original
!--------------------------------------------------------------------------
function om2eu_d(o) result(res)

use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN)	:: o(3,3)		!< orientation matrix
real(kind=dbl)			:: res(3)

res = qu2eu_d(om2qu_d(o))

end function om2eu_d

!--------------------------------------------------------------------------
!
! Function: ax2om
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Axis angle pair to orientation matrix (single precision)
!
!> @note verified 8/5/13.
!
!> @param a axis angle pair (single precision)  
!
!> @date 8/04/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ax2om(a) result(res)

use local

IMPLICIT NONE

real(kind=sgl),INTENT(IN)	:: a(4)		!< axis angle pair
real(kind=sgl)			:: res(3,3)

real(kind=sgl)			:: q, c, s, omc
integer(kind=irg)		:: i

c = cos(a(4))
s = sin(a(4))
omc = 1.0-c

do i=1,3
  res(i,i) = a(i)**2*omc + c
end do

q = omc*a(1)*a(2)
res(1,2) = q + s*a(3)
res(2,1) = q - s*a(3)

q = omc*a(2)*a(3)
res(2,3) = q + s*a(1)
res(3,2) = q - s*a(1)

q = omc*a(3)*a(1)
res(3,1) = q + s*a(2)
res(1,3) = q - s*a(2)

end function ax2om

!--------------------------------------------------------------------------
!
! Function: ax2om_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Axis angle pair to orientation matrix (double precision)
!
!> @note verified 8/5/13.
!
!> @param a axis angle pair (double precision)  
!
!> @date 8/04/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ax2om_d(a) result(res)

use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN)	:: a(4)		!< axis angle pair
real(kind=dbl)			:: res(3,3)

real(kind=dbl)			:: q, c, s, omc
integer(kind=irg)		:: i

c = dcos(a(4))
s = dsin(a(4))
omc = 1.D0-c

do i=1,3
  res(i,i) = a(i)**2*omc + c
end do


q = omc*a(1)*a(2)
res(1,2) = q + s*a(3)
res(2,1) = q - s*a(3)

q = omc*a(2)*a(3)
res(2,3) = q + s*a(1)
res(3,2) = q - s*a(1)

q = omc*a(3)*a(1)
res(3,1) = q + s*a(2)
res(1,3) = q - s*a(2)

end function ax2om_d

!--------------------------------------------------------------------------
!
! Function: qu2eu
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Quaternion to Euler angles (single precision) [Morawiec page 40, with errata !!!! ]
!
!> @param q quaternion (single precision)  
!
!> @date 8/04/13   MDG 1.0 original
!--------------------------------------------------------------------------
function qu2eu(q) result(res)

use local

IMPLICIT NONE

real(kind=sgl),INTENT(IN)	:: q(4)		!< quaternion
real(kind=sgl)			:: res(3)

real(kind=sgl)			:: q12, q03, chi, Phi, phi1, phi2
real(kind=sgl),parameter	:: pi = 3.14159265

q03 = q(1)**2+q(4)**2
q12 = q(2)**2+q(3)**2
chi = sqrt(q03*q12)

if (chi.eq.0.0) then
  if (q12.eq.0.0) then 
    Phi = 0.0
    phi2 = 0.0   		! arbitrarily due to degeneracy
    phi1 = atan2(-2.0*q(1)*q(4),q(1)**2-q(4)**2)
  else
    Phi = pi
    phi2 = 0.0   		! arbitrarily due to degeneracy
    phi1 = atan2(2.0*q(2)*q(3),q(2)**2-q(3)**2)
  end if
else  		! this is not a special degenerate case
  Phi = atan2( 2.0*chi, q03-q12 )
  chi = 1.0/chi
  phi1 = atan2( (-q(1)*q(3)+q(2)*q(4))*chi, (-q(1)*q(2)-q(3)*q(4))*chi )
  phi2 = atan2( (q(1)*q(3)+q(2)*q(4))*chi, (-q(1)*q(2)+q(3)*q(4))*chi )
end if

res = (/ phi1, Phi, phi2 /)

end function qu2eu

!--------------------------------------------------------------------------
!
! Function: qu2eu_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Quaternion to Euler angles (double precision) [Morawiec page 40, with errata !!!! ]
!
!> @param q quaternion (double precision)  
!
!> @date 8/04/13   MDG 1.0 original
!--------------------------------------------------------------------------
function qu2eu_d(q) result(res)

use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN)	:: q(4)		!< quaternion
real(kind=dbl)			:: res(3)

real(kind=dbl)			:: q12, q03, chi, Phi, phi1, phi2
real(kind=dbl),parameter	:: pi = 3.14159262217387D0

q03 = q(1)**2+q(4)**2
q12 = q(2)**2+q(3)**2
chi = dsqrt(q03*q12)

if (chi.eq.0.D0) then
  if (q12.eq.0.D0) then 
    Phi = 0.D0
    phi2 = 0.D0   		! arbitrarily due to degeneracy
    phi1 = datan2(-2.D0*q(1)*q(4),q(1)**2-q(4)**2)
  else
    Phi = pi
    phi2 = 0.D0   		! arbitrarily due to degeneracy
    phi1 = datan2(2.D0*q(2)*q(3),q(2)**2-q(3)**2)
  end if
else  		! this is not a special degenerate case
  Phi = datan2( 2.D0*chi, q03-q12 )
  chi = 1.D0/chi
  phi1 = datan2( (-q(1)*q(3)+q(2)*q(4))*chi, (-q(1)*q(2)-q(3)*q(4))*chi )
  phi2 = datan2( (q(1)*q(3)+q(2)*q(4))*chi, (-q(1)*q(2)+q(3)*q(4))*chi )
end if

res = (/ phi1, Phi, phi2 /)

end function qu2eu_d

!--------------------------------------------------------------------------
!
! Function: ro2eu
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Rodrigues vector to Euler angles (single precision)
!
!> @param r Rodrigues vector (single precision)  
! 
!> @todo what to do if r(1) = 0
!
!> @date 8/04/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ro2eu(r) result(res)

use local

IMPLICIT NONE

real(kind=sgl),INTENT(IN)	:: r(3)		!< Rodrigues vector
real(kind=sgl)			:: res(3)

real(kind=sgl)			:: s, d

s = atan(r(3))
d = atan(r(2)/r(1))

res = (/ s+d, 2.0*atan(r(1)*cos(s)/cos(d)), s-d /)

end function ro2eu

!--------------------------------------------------------------------------
!
! Function: ro2eu_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Rodrigues vector to Euler angles (double precision)
!
!> @param r Rodrigues vector (double precision)  
! 
!> @todo what to do if r(1) = 0
!
!> @date 8/04/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ro2eu_d(r) result(res)

use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN)	:: r(3)		!< Rodrigues vector
real(kind=dbl)			:: res(3)

real(kind=dbl)			:: s, d

s = datan(r(3))
d = datan(r(2)/r(1))

res = (/ s+d, 2.D0*datan(r(1)*dcos(s)/dcos(d)), s-d /)

end function ro2eu_d




!--------------------------------------------------------------------------
!
! Function: ax2ho
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Axis angle pair to homochoric (single precision)
!
!> @param a axis-angle pair (single precision)  
! !
!> @date 8/04/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ax2ho(a) result(res)

use local

IMPLICIT NONE

real(kind=sgl),INTENT(IN)	:: a(4)		!< axis angle pair
real(kind=sgl)			:: res(3)

real(kind=sgl)			:: f

f = 0.75 * ( a(4) - sin(a(4)) )
f = f**(1.0/3.0)

res = a(1:3) * f

end function ax2ho

!--------------------------------------------------------------------------
!
! Function: ax2ho_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Axis angle pair to homochoric (double precision)
!
!> @param a axis-angle pair (double precision)  
! !
!> @date 8/04/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ax2ho_d(a) result(res)

use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN)	:: a(4)		!< axis angle pair
real(kind=dbl)			:: res(3)

real(kind=dbl)			:: f

f = 0.75D0 * ( a(4) - dsin(a(4)) )
f = f**(1.D0/3.D0)

res = a(1:3) * f

end function ax2ho_d

!--------------------------------------------------------------------------
!
! Function: ho2ax
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Homochoric to axis angle pair (single precision)
!
!> @param h homochoric coordinates (single precision)  
! !
!> @date 8/04/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ho2ax(h) result(res)

use local

IMPLICIT NONE

real(kind=sgl),INTENT(IN)	:: h(3)		!< homochoric coordinates
real(kind=sgl)			:: res(4)

integer(kind=irg)		:: i
real(kind=sgl)			:: hn(3), hmag, s, hm

! fit parameters determined with Mathematica
real(kind=dbl),parameter	:: c(7) = (/ -0.5000096149170321D0, -0.02486606148871731D0, &
              				-0.004549381779362819D0, 0.0005118668366387526D0, &
              				-0.0016500827333575548D0, 0.0007593352203388718D0, &
              				-0.0002040422502566876D0 /)

! normalize h and store the magnitude
hmag = sum(h*h)
if (hmag.eq.0.0) then
  res = (/ 0.0, 0.0, 0.0, 0.0 /)
else
  hm = hmag
  hn = h/sqrt(hmag)

! convert the magnitude to the rotation angle
  s = c(1) * hmag
  do i=2,7
    hm = hm*hmag
    s = s + c(i) * hm
  end do

  res = (/ hn(1), hn(2), hn(3), 2.0*acos(1.0+s) /)
end if

end function ho2ax

!--------------------------------------------------------------------------
!
! Function: ho2ax_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Homochoric to axis angle pair (double precision)
!
!> @param h homochoric coordinates (double precision)  
! !
!> @date 8/04/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ho2ax_d(h) result(res)

use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN)	:: h(3)		!< homochoric coordinates
real(kind=dbl)			:: res(4)

integer(kind=irg)		:: i
real(kind=dbl)			:: hn(3), hmag, s, hm

! fit parameters determined with Mathematica
real(kind=dbl),parameter	:: c(7) = (/ -0.5000096149170321D0, -0.02486606148871731D0, &
              				-0.004549381779362819D0, 0.0005118668366387526D0, &
              				-0.0016500827333575548D0, 0.0007593352203388718D0, &
              				-0.0002040422502566876D0 /)

! normalize h and store the magnitude
hmag = sum(h*h)
if (hmag.eq.0.D0) then
  res = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
else
  hm = hmag
  hn = h/dsqrt(hmag)

! convert the magnitude to the rotation angle
  s = c(1) * hmag
  do i=2,7
    hm = hm*hmag
    s = s + c(i) * hm
  end do

  res = (/ hn(1), hn(2), hn(3), 2.D0*dacos(1.D0+s) /)
end if

end function ho2ax_d

!--------------------------------------------------------------------------
!
! Function: ro2ax
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Rodrigues vector to axis angle pair (single precision)
!
!> @param r Rodrigues vector (single precision)  
! 
!> @todo figure out what to do when ta equals zero
!
!> @date 8/04/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ro2ax(r) result(res)

use local

IMPLICIT NONE

real(kind=sgl),INTENT(IN)	:: r(3)		!< input Rodrigues vector
real(kind=sgl)			:: res(4)	!< output axis-angle pair

real(kind=sgl)			:: ta, angle

ta = sqrt(sum(r**2))		! tan(omega/2)
angle = 2.0*atan(ta)

res = (/ r(1)/ta, r(2)/ta, r(3)/ta, angle /)

! normalize the direction cosines
ta = sqrt(sum(res(1:3)**2))
res(1:3) = res(1:3)/ta

end function ro2ax

!--------------------------------------------------------------------------
!
! Function: ro2ax_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Rodrigues vector to axis angle pair (double precision)
!
!> @param r Rodrigues vector (double precision)  
! 
!> @todo figure out what to do when ta equals zero
!
!> @date 8/04/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ro2ax_d(r) result(res)

use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN)	:: r(3)		!< input Rodrigues vector
real(kind=dbl)			:: res(4)	!< output axis-angle pair

real(kind=dbl)			:: ta, angle

ta = dsqrt(sum(r**2))		! tan(omega/2)
angle = 2.D0*datan(ta)

res = (/ r(1)/ta, r(2)/ta, r(3)/ta, angle /)

! normalize the direction cosines
ta = dsqrt(sum(res(1:3)**2))
res(1:3) = res(1:3)/ta

end function ro2ax_d



!--------------------------------------------------------------------------
!
! Function: eu2ro
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Euler angles to Rodrigues vector (single precision) [Morawiec, page 40]
!
!> @param e 3 Euler angles in radians (single precision)  
! 
!> @todo figure out what to do when cs equals zero
!
!> @date 8/04/13   MDG 1.0 original
!--------------------------------------------------------------------------
function eu2ro(e) result(res)

use local

IMPLICIT NONE

real(kind=sgl),INTENT(IN)	:: e(3)		!< input Euler angles (radians)
real(kind=sgl)			:: res(3)	!< output Rodrigues vector

real(kind=sgl)			:: ee(3), sm, df, cs, t

ee = e*0.5

sm = ee(1)+ee(3)
df = ee(1)-ee(3)

cs = cos(sm)
t = tan(ee(2))

res = (/ t * cos(df)/cs, t * sin(df)/cs, tan(sm) /)

end function eu2ro

!--------------------------------------------------------------------------
!
! Function: eu2ro_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Euler angles to Rodrigues vector (double precision) [Morawiec, page 40]
!
!> @param e 3 Euler angles in radians (double precision)  
! 
!> @todo figure out what to do when cs equals zero
!
!> @date 8/04/13   MDG 1.0 original
!--------------------------------------------------------------------------
function eu2ro_d(e) result(res)

use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN)	:: e(3)		!< input Euler angles (radians)
real(kind=dbl)			:: res(3)	!< output Rodrigues vector

real(kind=dbl)			:: ee(3), sm, df, cs, t

ee = e*0.5D0

sm = ee(1)+ee(3)
df = ee(1)-ee(3)

cs = dcos(sm)
t = dtan(ee(2))

res = (/ t * dcos(df)/cs, t * dsin(df)/cs, tan(sm) /)

end function eu2ro_d


!--------------------------------------------------------------------------
!
! Function: eu2qu
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Euler angles to quaternion (single precision) [Morawiec, page 40]
!
!> @note verified 8/5/13
!
!> @param e 3 Euler angles in radians (single precision)  
! 
!> @date 8/04/13   MDG 1.0 original
!--------------------------------------------------------------------------
function eu2qu(e) result(res)

use local

IMPLICIT NONE

real(kind=sgl),INTENT(IN)	:: e(3)		!< input Euler angles in radians
real(kind=sgl)			:: res(4)	!< output quaternion

real(kind=sgl)			:: ee(3), phip, phim, cPhi, cp, cm, sPhi, sp, sm, q

ee = 0.5*e

cPhi = cos(ee(2))
sPhi = sin(ee(2))
cm = cos(ee(1)-ee(3))
sm = sin(ee(1)-ee(3))
cp = cos(ee(1)+ee(3))
sp = sin(ee(1)+ee(3))

q = cPhi*cp

! get the sign correct (Northern hemisphere)
if (q.lt.0.0) then
  res = (/ -q, sPhi*cm, sPhi*sm, cPhi*sp /)
else
  res = (/ q, -sPhi*cm, -sPhi*sm, -cPhi*sp /)
end if

end function eu2qu

!--------------------------------------------------------------------------
!
! Function: eu2qu_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Euler angles to quaternion (double precision) [Morawiec, page 40]
!
!> @note verified 8/5/13
!
!> @param e 3 Euler angles in radians (double precision)  
! 
!> @date 8/04/13   MDG 1.0 original
!--------------------------------------------------------------------------
function eu2qu_d(e) result(res)

use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN)	:: e(3)		!< input Euler angles in radians
real(kind=dbl)			:: res(4)	!< output quaternion

real(kind=dbl)			:: ee(3), phip, phim, cPhi, cp, cm, sPhi, sp, sm, q

ee = 0.5D0*e

cPhi = dcos(ee(2))
sPhi = dsin(ee(2))
cm = dcos(ee(1)-ee(3))
sm = dsin(ee(1)-ee(3))
cp = dcos(ee(1)+ee(3))
sp = dsin(ee(1)+ee(3))

q = cPhi*cp

! get the sign correct (Northern hemisphere)
if (q.lt.0.D0) then
  res = (/ -q, sPhi*cm, sPhi*sm, cPhi*sp /)
else
  res = (/ q, -sPhi*cm, -sPhi*sm, -cPhi*sp /)
end if

end function eu2qu_d



!--------------------------------------------------------------------------
!
! Function: eu2om
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Euler angles to orientation matrix (single precision) [Morawiec, page 28]
!
!> @param e 3 Euler angles in radians (single precision)  
! 
!> @date 8/04/13   MDG 1.0 original
!--------------------------------------------------------------------------
function eu2om(e) result(res)

use local

IMPLICIT NONE

real(kind=sgl),INTENT(IN)	:: e(3)		!< Euler angles in radians
real(kind=sgl)			:: res(3,3)	!< output orientation matrix

real(kind=sgl)			:: c1, c2, c3, s1, s2, s3

c1 = cos(e(1))
c2 = cos(e(2))
c3 = cos(e(3))
s1 = sin(e(1))
s2 = sin(e(2))
s3 = sin(e(3))

res(1,1) = c1*c3-s1*s3*c2
res(1,2) = s1*c3+c1*s3*c2
res(1,3) = s3*s2
res(2,1) = -c1*s3-s1*c3*c2
res(2,2) = -s1*s3+c1*c3*c2
res(2,3) = c3*s2
res(3,1) = s1*s2
res(3,2) = -c1*s2
res(3,3) = c2

end function eu2om

!--------------------------------------------------------------------------
!
! Function: eu2om_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Euler angles to orientation matrix (double precision) [Morawiec, page 28]
!
!> @param e 3 Euler angles in radians (double precision)  
! 
!> @date 8/04/13   MDG 1.0 original
!--------------------------------------------------------------------------
function eu2om_d(e) result(res)

use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN)	:: e(3)		!< Euler angles in radians
real(kind=dbl)			:: res(3,3)	!< output orientation matrix

real(kind=dbl)			:: c1, c2, c3, s1, s2, s3

c1 = dcos(e(1))
c2 = dcos(e(2))
c3 = dcos(e(3))
s1 = dsin(e(1))
s2 = dsin(e(2))
s3 = dsin(e(3))

res(1,1) = c1*c3-s1*s3*c2
res(1,2) = s1*c3+c1*s3*c2
res(1,3) = s3*s2
res(2,1) = -c1*s3-s1*c3*c2
res(2,2) = -s1*s3+c1*c3*c2
res(2,3) = c3*s2
res(3,1) = s1*s2
res(3,2) = -c1*s2
res(3,3) = c2

end function eu2om_d


!--------------------------------------------------------------------------
!
! FUNCTION: qu2om
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert a quaternion to a 3x3 matrix
!
!> @param q quaternion (single precision)  
! 
!> @note verified 8/5/13
!
!> @date 6/03/13   MDG 1.0 original
!--------------------------------------------------------------------------
function qu2om(q) result(res)

use local

real(kind=sgl),INTENT(IN)		:: q(4)
real(kind=sgl)				:: res(3,3)

real(kind=sgl)				:: qq

qq=q(1)*q(1)-(q(2)*q(2)+q(3)*q(3)+q(4)*q(4))

res(1,1) = qq+2.0*q(2)*q(2)
res(2,2) = qq+2.0*q(3)*q(3)
res(3,3) = qq+2.0*q(4)*q(4)

res(1,2) = 2.0*(q(2)*q(3)-q(1)*q(4))
res(2,3) = 2.0*(q(3)*q(4)-q(1)*q(2))
res(3,1) = 2.0*(q(4)*q(2)-q(1)*q(3))
res(2,1) = 2.0*(q(3)*q(2)+q(1)*q(4))
res(3,2) = 2.0*(q(4)*q(3)+q(1)*q(2))
res(1,3) = 2.0*(q(2)*q(4)+q(1)*q(3))

end function qu2om


!--------------------------------------------------------------------------
!
! FUNCTION: qu2om_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert a quaternion to a 3x3 matrix (double precision)
!
!> @param q quaternion (double precision)  
! 
!> @note verified 8/5/13
!
!> @date 6/03/13   MDG 1.0 original
!--------------------------------------------------------------------------
function qu2om_d(q) result(res)

use local

real(kind=dbl),INTENT(IN)		:: q(4)
real(kind=dbl)				:: res(3,3)

real(kind=dbl)				:: qq

qq=q(1)*q(1)-(q(2)*q(2)+q(3)*q(3)+q(4)*q(4))

res(1,1) = qq+2.D0*q(2)*q(2)
res(2,2) = qq+2.D0*q(3)*q(3)
res(3,3) = qq+2.D0*q(4)*q(4)

res(1,2) = 2.D0*(q(2)*q(3)-q(1)*q(4))
res(2,3) = 2.D0*(q(3)*q(4)-q(1)*q(2))
res(3,1) = 2.D0*(q(4)*q(2)-q(1)*q(3))
res(2,1) = 2.D0*(q(3)*q(2)+q(1)*q(4))
res(3,2) = 2.D0*(q(4)*q(3)+q(1)*q(2))
res(1,3) = 2.D0*(q(2)*q(4)+q(1)*q(3))

end function qu2om_d


!--------------------------------------------------------------------------
!
! FUNCTION: om2qu
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert a 3x3 rotation matrix to a unit quaternion (see Morawiec, page 37)
!
!> @param x 3x3 matrix to be converted (single precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function om2qu(x) result (res)

use local 

real(kind=sgl), intent(in) 		:: x(3,3)		!< input matrix
real(kind=sgl) 				:: res(4)

real(kind=sgl)				:: s
    
s = x(1,1) + x(2,2) + x(3,3) + 1.0
res(1) = sqrt(s)/2.0				! use the positive sign here

s = abs(x(1,1) - x(2,2) - x(3,3) + 1.0)
res(2) = sqrt(s)/2.0

s = abs(-x(1,1) + x(2,2) - x(3,3) + 1.0)
res(3) = sqrt(s)/2.0

s = abs(-x(1,1) - x(2,2) + x(3,3) + 1.0)
res(4) = sqrt(s)/2.0

if (x(3,2)-x(2,3) .lt. 0.0)  res(2) = -res(2)
if (x(1,3)-x(3,1) .lt. 0.0)  res(3) = -res(3)
if (x(2,1)-x(1,2) .lt. 0.0)  res(4) = -res(4)

s = 1.0/cabs(res)

res = s * res

end function om2qu

!--------------------------------------------------------------------------
!
! FUNCTION: om2qu_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert a 3x3 rotation matrix to a unit quaternion (see Morawiec, page 37)
!
!> @param x 3x3 matrix to be converted (double precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function om2qu_d(x) result (res)

use local 

real(kind=dbl), intent(in) 		:: x(3,3)		!< input matrix
real(kind=dbl) 				:: res(4)

real(kind=dbl)				:: s
    
s = dabs(x(1,1) + x(2,2) + x(3,3) + 1.D0)
res(1) = dsqrt(s)/2.D0				! use the positive sign here

s = dabs(x(1,1) - x(2,2) - x(3,3) + 1.D0)
res(2) = dsqrt(s)/2.D0

s = dabs(-x(1,1) + x(2,2) - x(3,3) + 1.D0)
res(3) = dsqrt(s)/2.D0

s = dabs(-x(1,1) - x(2,2) + x(3,3) + 1.D0)
res(4) = dsqrt(s)/2.D0

if (x(3,2)-x(2,3) .lt. 0.D0)  res(2) = -res(2)
if (x(1,3)-x(3,1) .lt. 0.D0)  res(3) = -res(3)
if (x(2,1)-x(1,2) .lt. 0.D0)  res(4) = -res(4)

s = 1.D0/cabs(res)

res = s * res

end function om2qu_d


!--------------------------------------------------------------------------
!
! FUNCTION: ho2cu
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert homochoric to cubochoric
!
!> @param h homochoric coordinates (single precision)
! 
!> @note calling program MUST have initialized the Lambert parameters first!!!
!
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ho2cu(h) result (res)

use local 
use Lambert, only: LambertBallToCube

real(kind=sgl), intent(in) 		:: h(3)		!< input coordinates
real(kind=sgl) 				:: res(3)

integer(kind=irg)			:: ierr

! calling program must have initialized the Lambert parameters!!!!
res = LambertBallToCube(h,ierr)

end function ho2cu

!--------------------------------------------------------------------------
!
! FUNCTION: ho2cu_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert homochoric to cubochoric
!
!> @param h homochoric coordinates (double precision)
! 
!> @note calling program MUST have initialized the Lambert parameters first!!!
!
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ho2cu_d(h) result (res)

use local 
use Lambert, only: LambertBallToCube

real(kind=dbl), intent(in) 		:: h(3)		!< input coordinates
real(kind=dbl) 				:: res(3)

integer(kind=irg)			:: ierr

! calling program must have initialized the Lambert parameters!!!!
res = LambertBallToCube(h,ierr)

end function ho2cu_d

!--------------------------------------------------------------------------
!
! FUNCTION: cu2ho
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert cubochoric to homochoric
!
!> @param c cubochoric coordinates (single precision)
! 
!> @note calling program MUST have initialized the Lambert parameters first!!!
!
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function cu2ho(c) result (res)

use local 
use Lambert, only: LambertCubeToBall

real(kind=sgl), intent(in) 		:: c(3)		!< input coordinates
real(kind=sgl) 				:: res(3)

integer(kind=irg)			:: ierr

! calling program must have initialized the Lambert parameters!!!!
res = LambertCubeToBall(c,ierr)

end function cu2ho

!--------------------------------------------------------------------------
!
! FUNCTION: cu2ho_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert cubochoric to homochoric
!
!> @param c cubochoric coordinates (double precision)
! 
!> @note calling program MUST have initialized the Lambert parameters first!!!
!
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function cu2ho_d(c) result (res)

use local 
use Lambert, only: LambertCubeToBall

real(kind=dbl), intent(in) 		:: c(3)		!< input coordinates
real(kind=dbl) 				:: res(3)

integer(kind=irg)			:: ierr

! calling program must have initialized the Lambert parameters!!!!
res = LambertCubeToBall(c,ierr)

end function cu2ho_d


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! and here are a bunch of transformation routines that are derived from the others
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!
! FUNCTION: eu2ax
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert euler to axis angle
!
!> @param e 3 euler angles (single precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function eu2ax(e) result (res)

use local 

real(kind=sgl), intent(in) 		:: e(3)
real(kind=sgl)				:: res(4)

res = ro2ax(eu2ro(e))

end function eu2ax

!--------------------------------------------------------------------------
!
! FUNCTION: eu2ax_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert euler to axis angle
!
!> @param e 3 euler angles (double precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function eu2ax_d(e) result (res)

use local 

real(kind=dbl), intent(in) 		:: e(3)
real(kind=dbl)				:: res(4)

res = ro2ax_d(eu2ro_d(e))

end function eu2ax_d

!--------------------------------------------------------------------------
!
! FUNCTION: eu2ho
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert euler to homochoric
!
!> @param e 3 euler angles (single precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function eu2ho(e) result (res)

use local 

real(kind=sgl), intent(in) 		:: e(3)
real(kind=sgl)				:: res(3)

res = ax2ho(eu2ax(e))

end function eu2ho

!--------------------------------------------------------------------------
!
! FUNCTION: eu2ho_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert euler to homochoric
!
!> @param e 3 euler angles (double precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function eu2ho_d(e) result (res)

use local 

real(kind=dbl), intent(in) 		:: e(3)
real(kind=dbl)				:: res(3)

res = ax2ho_d(eu2ax_d(e))

end function eu2ho_d


!--------------------------------------------------------------------------
!
! FUNCTION: om2ax
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert orientation matrix to axis angle
!
!> @param om 3x3 orientation matrix (single precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function om2ax(om) result (res)

use local 

real(kind=sgl), intent(in) 		:: om(3,3)
real(kind=sgl)				:: res(4)

res = eu2ax(om2eu(om))

end function om2ax

!--------------------------------------------------------------------------
!
! FUNCTION: om2ax_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert orientation matrix to axis angle
!
!> @param om 3x3 orientation matrix (double precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function om2ax_d(om) result (res)

use local 

real(kind=dbl), intent(in) 		:: om(3,3)
real(kind=dbl)				:: res(4)

res = eu2ax_d(om2eu_d(om))

end function om2ax_d

!--------------------------------------------------------------------------
!
! FUNCTION: om2ro
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert orientation matrix to Rodrigues
!
!> @param om 3x3 orientation matrix (single precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function om2ro(om) result (res)

use local 

real(kind=sgl), intent(in) 		:: om(3,3)
real(kind=sgl)				:: res(3)

res = eu2ro(om2eu(om))

end function om2ro

!--------------------------------------------------------------------------
!
! FUNCTION: om2ro_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert orientation matrix to Rodrigues
!
!> @param om 3x3 orientation matrix (double precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function om2ro_d(om) result (res)

use local 

real(kind=dbl), intent(in) 		:: om(3,3)
real(kind=dbl)				:: res(3)

res = eu2ro_d(om2eu_d(om))

end function om2ro_d

!--------------------------------------------------------------------------
!
! FUNCTION: om2ho
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert orientation matrix to homochoric
!
!> @param om 3x3 orientation matrix (single precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function om2ho(om) result (res)

use local 

real(kind=sgl), intent(in) 		:: om(3,3)
real(kind=sgl)				:: res(3)

res = eu2ho(om2eu(om))

end function om2ho

!--------------------------------------------------------------------------
!
! FUNCTION: om2ho_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert orientation matrix to homochoric
!
!> @param om 3x3 orientation matrix (double precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function om2ho_d(om) result (res)

use local 

real(kind=dbl), intent(in) 		:: om(3,3)
real(kind=dbl)				:: res(3)

res = eu2ho(om2eu(om))

end function om2ho_d

!--------------------------------------------------------------------------
!
! FUNCTION: ax2eu
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert axis angle to euler
!
!> @param a axis angle pair (single precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ax2eu(a) result (res)

use local 

real(kind=sgl), intent(in) 		:: a(4)
real(kind=sgl)				:: res(3)

res = om2eu(ax2om(a))

end function ax2eu

!--------------------------------------------------------------------------
!
! FUNCTION: ax2eu_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert axis angle to euler
!
!> @param a axis angle pair (double precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ax2eu_d(a) result (res)

use local 

real(kind=dbl), intent(in) 		:: a(4)
real(kind=dbl)				:: res(3)

res = om2eu_d(ax2om_d(a))

end function ax2eu_d

!--------------------------------------------------------------------------
!
! FUNCTION: ax2ro
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert axis angle to Rodrigues
!
!> @param a axis angle pair (single precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ax2ro(a) result (res)

use local 

real(kind=sgl), intent(in) 		:: a(4)
real(kind=sgl)				:: res(3)

res = om2ro(ax2om(a))

end function ax2ro

!--------------------------------------------------------------------------
!
! FUNCTION: ax2ro_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert axis angle to Rodrigues
!
!> @param a axis angle pair (double precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ax2ro_d(a) result (res)

use local 

real(kind=dbl), intent(in) 		:: a(4)
real(kind=dbl)				:: res(3)

res = om2ro_d(ax2om_d(a))

end function ax2ro_d

!--------------------------------------------------------------------------
!
! FUNCTION: ax2qu
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert axis angle to quaternion
!
!> @param a axis angle pair (single precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ax2qu(a) result (res)

use local 

real(kind=sgl), intent(in) 		:: a(4)
real(kind=sgl)				:: res(4)

res = om2qu(ax2om(a))

end function ax2qu

!--------------------------------------------------------------------------
!
! FUNCTION: ax2qu_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert axis angle to quaternion
!
!> @param a axis angle pair (double precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ax2qu_d(a) result (res)

use local 

real(kind=dbl), intent(in) 		:: a(4)
real(kind=dbl)				:: res(4)

res = om2qu_d(ax2om_d(a))

end function ax2qu_d

!--------------------------------------------------------------------------
!
! FUNCTION: ro2om
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert rodrigues to orientation matrix
!
!> @param r Rodrigues vector (single precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ro2om(r) result (res)

use local 

real(kind=sgl), intent(in) 		:: r(3)
real(kind=sgl)				:: res(3,3)

res = ax2om(ro2ax(r))

end function ro2om

!--------------------------------------------------------------------------
!
! FUNCTION: ro2om_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert rodrigues to orientation matrix
!
!> @param r Rodrigues vector (double precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ro2om_d(r) result (res)

use local 

real(kind=dbl), intent(in) 		:: r(3)
real(kind=dbl)				:: res(3,3)

res = ax2om_d(ro2ax_d(r))

end function ro2om_d


!--------------------------------------------------------------------------
!
! FUNCTION: ro2qu
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert rodrigues to quaternion
!
!> @param r Rodrigues vector (single precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ro2qu(r) result (res)

use local 

real(kind=sgl), intent(in) 		:: r(3)
real(kind=sgl)				:: res(4)

res = ax2qu(ro2ax(r))

end function ro2qu


!--------------------------------------------------------------------------
!
! FUNCTION: ro2qu_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert rodrigues to quaternion
!
!> @param r Rodrigues vector (double precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ro2qu_d(r) result (res)

use local 

real(kind=dbl), intent(in) 		:: r(3)
real(kind=dbl)				:: res(4)

res = ax2qu_d(ro2ax_d(r))

end function ro2qu_d

!--------------------------------------------------------------------------
!
! FUNCTION: ro2ho
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert rodrigues to homochoric
!
!> @param r Rodrigues vector (single precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ro2ho(r) result (res)

use local 

real(kind=sgl), intent(in) 		:: r(3)
real(kind=sgl)				:: res(3)

res = ax2ho(ro2ax(r))

end function ro2ho

!--------------------------------------------------------------------------
!
! FUNCTION: ro2ho_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert rodrigues to homochoric
!
!> @param r Rodrigues vector (double precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ro2ho_d(r) result (res)

use local 

real(kind=dbl), intent(in) 		:: r(3)
real(kind=dbl)				:: res(3)

res = ax2ho_d(ro2ax_d(r))

end function ro2ho_d


!--------------------------------------------------------------------------
!
! FUNCTION: qu2ax
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert quaternion to axis angle
!
!> @param q quaternion (single precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function qu2ax(q) result (res)

use local 

real(kind=sgl), intent(in) 		:: q(4)
real(kind=sgl)				:: res(4)

res = eu2ax(qu2eu(q))

end function qu2ax

!--------------------------------------------------------------------------
!
! FUNCTION: qu2ax_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert quaternion to axis angle
!
!> @param q quaternion (double precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function qu2ax_d(q) result (res)

use local 

real(kind=dbl), intent(in) 		:: q(4)
real(kind=dbl)				:: res(4)

res = eu2ax_d(qu2eu_d(q))

end function qu2ax_d


!--------------------------------------------------------------------------
!
! FUNCTION: qu2ro
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert quaternion to Rodrigues
!
!> @param q quaternion (single precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function qu2ro(q) result (res)

use local 

real(kind=sgl), intent(in) 		:: q(4)
real(kind=sgl)				:: res(3)

res = eu2ro(qu2eu(q))

end function qu2ro

!--------------------------------------------------------------------------
!
! FUNCTION: qu2ro_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert quaternion to Rodrigues
!
!> @param q quaternion (double precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function qu2ro_d(q) result (res)

use local 

real(kind=dbl), intent(in) 		:: q(4)
real(kind=dbl)				:: res(3)

res = eu2ro_d(qu2eu_d(q))

end function qu2ro_d

!--------------------------------------------------------------------------
!
! FUNCTION: qu2ho
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert quaternion to homochoric
!
!> @param q quaternion (single precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function qu2ho(q) result (res)

use local 

real(kind=sgl), intent(in) 		:: q(4)
real(kind=sgl)				:: res(3)

res = eu2ho(qu2eu(q))

end function qu2ho

!--------------------------------------------------------------------------
!
! FUNCTION: qu2ho_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert quaternion to homochoric
!
!> @param q quaternion (double precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function qu2ho_d(q) result (res)

use local 

real(kind=dbl), intent(in) 		:: q(4)
real(kind=dbl)				:: res(3)

res = eu2ho_d(qu2eu_d(q))

end function qu2ho_d

!--------------------------------------------------------------------------
!
! FUNCTION: ho2eu
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert homochoric to euler
!
!> @param h homochoric coordinates (single precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ho2eu(h) result (res)

use local 

real(kind=sgl), intent(in) 		:: h(3)
real(kind=sgl)				:: res(3)

res = ax2eu(ho2ax(h))

end function ho2eu

!--------------------------------------------------------------------------
!
! FUNCTION: ho2eu_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert homochoric to euler
!
!> @param h homochoric coordinates (double precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ho2eu_d(h) result (res)

use local 

real(kind=dbl), intent(in) 		:: h(3)
real(kind=dbl)				:: res(3)

res = ax2eu_d(ho2ax_d(h))

end function ho2eu_d

!--------------------------------------------------------------------------
!
! FUNCTION: ho2om
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert homochoric to orientation matrix
!
!> @param h homochoric coordinates (single precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ho2om(h) result (res)

use local 

real(kind=sgl), intent(in) 		:: h(3)
real(kind=sgl)				:: res(3,3)

res = ax2om(ho2ax(h))

end function ho2om

!--------------------------------------------------------------------------
!
! FUNCTION: ho2om_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert homochoric to orientation matrix
!
!> @param h homochoric coordinates (double precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ho2om_d(h) result (res)

use local 

real(kind=dbl), intent(in) 		:: h(3)
real(kind=dbl)				:: res(3,3)

res = ax2om_d(ho2ax_d(h))

end function ho2om_d

!--------------------------------------------------------------------------
!
! FUNCTION: ho2ro
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert homochoric to Rodrigues
!
!> @param h homochoric coordinates (single precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ho2ro(h) result (res)

use local 

real(kind=sgl), intent(in) 		:: h(3)
real(kind=sgl)				:: res(3)

res = ax2ro(ho2ax(h))

end function ho2ro

!--------------------------------------------------------------------------
!
! FUNCTION: ho2ro_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert homochoric to Rodrigues
!
!> @param h homochoric coordinates (double precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ho2ro_d(h) result (res)

use local 

real(kind=dbl), intent(in) 		:: h(3)
real(kind=dbl)				:: res(3)

res = ax2ro_d(ho2ax_d(h))

end function ho2ro_d

!--------------------------------------------------------------------------
!
! FUNCTION: ho2qu
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert homochoric to quaternion
!
!> @param h homochoric coordinates (single precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ho2qu(h) result (res)

use local 

real(kind=sgl), intent(in) 		:: h(3)
real(kind=sgl)				:: res(4)

res = ax2qu(ho2ax(h))

end function ho2qu

!--------------------------------------------------------------------------
!
! FUNCTION: ho2qu_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert homochoric to quaternion
!
!> @param h homochoric coordinates (double precision)
! 
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ho2qu_d(h) result (res)

use local 

real(kind=dbl), intent(in) 		:: h(3)
real(kind=dbl)				:: res(4)

res = ax2qu_d(ho2ax_d(h))

end function ho2qu_d

!--------------------------------------------------------------------------
!
! FUNCTION: eu2cu
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert euler angles to cubochoric
!
!> @param e euler angles (single precision)
! 
!> @note calling program MUST have initialized the Lambert parameters first!!!
!
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function eu2cu(e) result (res)

use local 

real(kind=sgl), intent(in) 		:: e(3)		!< input coordinates
real(kind=sgl) 				:: res(3)

res = ho2cu(eu2ho(e))

end function eu2cu

!--------------------------------------------------------------------------
!
! FUNCTION: eu2cu_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert euler angles to cubochoric
!
!> @param e euler angles (double precision)
! 
!> @note calling program MUST have initialized the Lambert parameters first!!!
!
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function eu2cu_d(e) result (res)

use local 

real(kind=dbl), intent(in) 		:: e(3)		!< input coordinates
real(kind=dbl) 				:: res(3)

res = ho2cu_d(eu2ho_d(e))

end function eu2cu_d

!--------------------------------------------------------------------------
!
! FUNCTION: om2cu
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert orientation matrix to cubochoric
!
!> @param o orientation matrix (single precision)
! 
!> @note calling program MUST have initialized the Lambert parameters first!!!
!
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function om2cu(o) result (res)

use local 

real(kind=sgl), intent(in) 		:: o(3,3)		!< input coordinates
real(kind=sgl) 				:: res(3)

res = ho2cu(om2ho(o))

end function om2cu

!--------------------------------------------------------------------------
!
! FUNCTION: om2cu_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert orientation matrix to cubochoric
!
!> @param o orientation matrix (double precision)
! 
!> @note calling program MUST have initialized the Lambert parameters first!!!
!
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function om2cu_d(o) result (res)

use local 

real(kind=dbl), intent(in) 		:: o(3,3)		!< input coordinates
real(kind=dbl) 				:: res(3)

res = ho2cu_d(om2ho_d(o))

end function om2cu_d

!--------------------------------------------------------------------------
!
! FUNCTION: ax2cu
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert axis angle to cubochoric
!
!> @param a axis angle (single precision)
! 
!> @note calling program MUST have initialized the Lambert parameters first!!!
!
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ax2cu(a) result (res)

use local 

real(kind=sgl), intent(in) 		:: a(4)		!< input coordinates
real(kind=sgl) 				:: res(3)

res = ho2cu(ax2ho(a))

end function ax2cu

!--------------------------------------------------------------------------
!
! FUNCTION: ax2cu_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert axis angle to cubochoric
!
!> @param a axis angle (double precision)
! 
!> @note calling program MUST have initialized the Lambert parameters first!!!
!
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ax2cu_d(a) result (res)

use local 

real(kind=dbl), intent(in) 		:: a(4)		!< input coordinates
real(kind=dbl) 				:: res(3)

res = ho2cu_d(ax2ho_d(a))

end function ax2cu_d

!--------------------------------------------------------------------------
!
! FUNCTION: ro2cu
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert Rodrigues to cubochoric
!
!> @param r Rodrigues (single precision)
! 
!> @note calling program MUST have initialized the Lambert parameters first!!!
!
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ro2cu(r) result (res)

use local 

real(kind=sgl), intent(in) 		:: r(3)		!< input coordinates
real(kind=sgl) 				:: res(3)

res = ho2cu(ro2ho(r))

end function ro2cu

!--------------------------------------------------------------------------
!
! FUNCTION: ro2cu_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert Rodrigues to cubochoric
!
!> @param r Rodrigues (double precision)
! 
!> @note calling program MUST have initialized the Lambert parameters first!!!
!
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ro2cu_d(r) result (res)

use local 

real(kind=dbl), intent(in) 		:: r(3)		!< input coordinates
real(kind=dbl) 				:: res(3)

res = ho2cu_d(ro2ho_d(r))

end function ro2cu_d

!--------------------------------------------------------------------------
!
! FUNCTION: qu2cu
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert quaternion to cubochoric
!
!> @param q quaternion (single precision)
! 
!> @note calling program MUST have initialized the Lambert parameters first!!!
!
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function qu2cu(q) result (res)

use local 

real(kind=sgl), intent(in) 		:: q(4)		!< input coordinates
real(kind=sgl) 				:: res(3)

res = ho2cu(qu2ho(q))

end function qu2cu

!--------------------------------------------------------------------------
!
! FUNCTION: qu2cu_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert quaternion to cubochoric
!
!> @param q quaternion (double precision)
! 
!> @note calling program MUST have initialized the Lambert parameters first!!!
!
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function qu2cu_d(q) result (res)

use local 

real(kind=dbl), intent(in) 		:: q(4)		!< input coordinates
real(kind=dbl) 				:: res(3)

res = ho2cu_d(qu2ho_d(q))

end function qu2cu_d


!--------------------------------------------------------------------------
!
! FUNCTION: cu2eu
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert cubochoric to euler angles
!
!> @param c cubochoric coordinates (single precision)
! 
!> @note calling program MUST have initialized the Lambert parameters first!!!
!
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function cu2eu(c) result (res)

use local 

real(kind=sgl), intent(in) 		:: c(3)		!< input coordinates
real(kind=sgl) 				:: res(3)

res = ho2eu(cu2ho(c))

end function cu2eu

!--------------------------------------------------------------------------
!
! FUNCTION: cu2eu_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert cubochoric to euler angles
!
!> @param c cubochoric coordinates (double precision)
! 
!> @note calling program MUST have initialized the Lambert parameters first!!!
!
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function cu2eu_d(c) result (res)

use local 

real(kind=dbl), intent(in) 		:: c(3)		!< input coordinates
real(kind=dbl) 				:: res(3)

res = ho2eu_d(cu2ho_d(c))

end function cu2eu_d

!--------------------------------------------------------------------------
!
! FUNCTION: cu2om
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert cubochoric to orientation matrix
!
!> @param c cubochoric coordinates (single precision)
! 
!> @note calling program MUST have initialized the Lambert parameters first!!!
!
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function cu2om(c) result (res)

use local 

real(kind=sgl), intent(in) 		:: c(3)		!< input coordinates
real(kind=sgl) 				:: res(3,3)

res = ho2om(cu2ho(c))

end function cu2om

!--------------------------------------------------------------------------
!
! FUNCTION: cu2om_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert cubochoric to orientation matrix
!
!> @param c cubochoric coordinates  (double precision)
! 
!> @note calling program MUST have initialized the Lambert parameters first!!!
!
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function cu2om_d(c) result (res)

use local 

real(kind=dbl), intent(in) 		:: c(3)		!< input coordinates
real(kind=dbl) 				:: res(3,3)

res = ho2om_d(cu2ho_d(c))

end function cu2om_d

!--------------------------------------------------------------------------
!
! FUNCTION: cu2ax
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert cubochoric to axis angle
!
!> @param c cubochoric coordinates (single precision)
! 
!> @note calling program MUST have initialized the Lambert parameters first!!!
!
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function cu2ax(c) result (res)

use local 

real(kind=sgl), intent(in) 		:: c(3)		!< input coordinates
real(kind=sgl) 				:: res(4)

res = ho2ax(cu2ho(c))

end function cu2ax

!--------------------------------------------------------------------------
!
! FUNCTION: cu2ax_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert cubochoric to axis angle
!
!> @param c cubochoric coordinates  (double precision)
! 
!> @note calling program MUST have initialized the Lambert parameters first!!!
!
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function cu2ax_d(c) result (res)

use local 

real(kind=dbl), intent(in) 		:: c(3)		!< input coordinates
real(kind=dbl) 				:: res(4)

res = ho2ax_d(cu2ho_d(c))

end function cu2ax_d

!--------------------------------------------------------------------------
!
! FUNCTION: cu2ro
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert cubochoric to Rodrigues
!
!> @param c cubochoric coordinates (single precision)
! 
!> @note calling program MUST have initialized the Lambert parameters first!!!
!
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function cu2ro(c) result (res)

use local 

real(kind=sgl), intent(in) 		:: c(3)		!< input coordinates
real(kind=sgl) 				:: res(3)

res = ho2ro(cu2ho(c))

end function cu2ro

!--------------------------------------------------------------------------
!
! FUNCTION: cu2ro_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert cubochoric to Rodrigues
!
!> @param c cubochoric coordinates  (double precision)
! 
!> @note calling program MUST have initialized the Lambert parameters first!!!
!
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function cu2ro_d(c) result (res)

use local 

real(kind=dbl), intent(in) 		:: c(3)		!< input coordinates
real(kind=dbl) 				:: res(3)

res = ho2ro_d(cu2ho_d(c))

end function cu2ro_d

!--------------------------------------------------------------------------
!
! FUNCTION: cu2qu
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert cubochoric to quaternion
!
!> @param c cubochoric coordinates (single precision)
! 
!> @note calling program MUST have initialized the Lambert parameters first!!!
!
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function cu2qu(c) result (res)

use local 

real(kind=sgl), intent(in) 		:: c(3)		!< input coordinates
real(kind=sgl) 				:: res(4)

res = ho2qu(cu2ho(c))

end function cu2qu

!--------------------------------------------------------------------------
!
! FUNCTION: cu2qu_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert cubochoric to quaternion
!
!> @param c cubochoric coordinates  (double precision)
! 
!> @note calling program MUST have initialized the Lambert parameters first!!!
!
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function cu2qu_d(c) result (res)

use local 

real(kind=dbl), intent(in) 		:: c(3)		!< input coordinates
real(kind=dbl) 				:: res(4)

res = ho2qu_d(cu2ho_d(c))

end function cu2qu_d



!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! and finally some printing routines, mostly used for debugging
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!
! SUBROUTINE: print_orientation
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  prints a complete orientationtype record or a single entry
! 
!> @param o orientationtype record
!> @param outtype (optional) indicates which representation to print
!> @param pretext (optional) up to 10 characters that will precede each line
! 
!> @date  8/4/13   MDG 1.0 original
!--------------------------------------------------------------------------
subroutine print_orientation(o,outtype,pretext)

use local
use io
use constants

IMPLICIT NONE

type(orientationtype),INTENT(IN)	:: o
character(2),INTENT(IN),OPTIONAL	:: outtype
character(10),INTENT(IN),OPTIONAL	:: pretext

real(kind=sgl)				:: ioreal(4)
character(10)				:: pret

pret = ''
if (present(pretext)) pret=trim(pretext)

if (present(outtype)) then
  select case (outtype)
  	case ('eu')
	  ioreal(1:3) = o%eulang(1:3)*180.0/sngl(cPi)
  	  call WriteValue(trim(pret)//'Euler angles			: ', ioreal, 3, "(3(F8.4,' '))")

  	case ('ax')
	  ioreal(1:4) = o%axang(1:4)
 	  ioreal(4) = ioreal(4)*180.0/sngl(cPi)
	  call WriteValue(trim(pret)//'Axis angle pair [n; angle]	: ', ioreal, 4, "(3(F8.4,' '),'; ',F8.4)")

  	case ('ro')
	  ioreal(1:3) = o%rodrigues(1:3)
  	  call WriteValue(trim(pret)//'Rodigues vector		: ', ioreal, 3, "(3(F8.4,' '))")

  	case ('ho')
	  ioreal(1:3) = o%homochoric(1:3)
  	  call WriteValue(trim(pret)//'Homochoric representation	: ', ioreal, 3, "(3(F8.4,' '))")

  	case ('cu')
  	  ioreal(1:3) = o%cubochoric(1:3)
  	  call WriteValue(trim(pret)//'Cubochoric representation	: ', ioreal, 3, "(3(F8.4,' '))")

  	case ('qu')
	  ioreal(1:4) = o%quat
	  call WriteValue(trim(pret)//'Quaternion			: ', ioreal, 4, "(4(F8.4,' '))")

  	case ('om')
	  ioreal(1:3) = o%om(1,1:3)
	  call WriteValue('					  /', ioreal, 3, "(2(F8.4,' '),F8.4,' \')")
	  ioreal(1:3) = o%om(2,1:3)
	  call WriteValue(trim(pret)//'Orientation Matrix		: |', ioreal, 3, "(2(F8.4,' '),F8.4,' |')")
	  ioreal(1:3) = o%om(3,1:3)
	  call WriteValue('					  \', ioreal, 3, "(2(F8.4,' '),F8.4,' /')")
  		  		
  end select
else
! print the entire record with all representations
  ioreal(1:3) = o%eulang(1:3)*180.0/sngl(cPi)
  call WriteValue(trim(pret)//'Euler angles			: ', ioreal, 3, "(3(F8.4,' '))")
  ioreal(1:4) = o%axang(1:4)
  ioreal(4) = ioreal(4)*180.0/sngl(cPi)
  call WriteValue(trim(pret)//'Axis angle pair [n; angle]	: ', ioreal, 4, "(3(F8.4,' '),'; ',F8.4)")
  ioreal(1:3) = o%rodrigues(1:3)
  call WriteValue(trim(pret)//'Rodigues vector		: ', ioreal, 3, "(3(F8.4,' '))")
  ioreal(1:3) = o%homochoric(1:3)
  call WriteValue(trim(pret)//'Homochoric representation	: ', ioreal, 3, "(3(F8.4,' '))")
  ioreal(1:3) = o%cubochoric(1:3)
  call WriteValue(trim(pret)//'Cubochoric representation	: ', ioreal, 3, "(3(F8.4,' '))")
  ioreal(1:4) = o%quat
  call WriteValue(trim(pret)//'Quaternion			: ', ioreal, 4, "(4(F8.4,' '))")
  ioreal(1:3) = o%om(1,1:3)
  call WriteValue('					  /', ioreal, 3, "(2(F8.4,' '),F8.4,' \')")
  ioreal(1:3) = o%om(2,1:3)
  call WriteValue(trim(pret)//'Orientation Matrix		: |', ioreal, 3, "(2(F8.4,' '),F8.4,' |')")
  ioreal(1:3) = o%om(3,1:3)
  call WriteValue('					  \', ioreal, 3, "(2(F8.4,' '),F8.4,' /')")
end if

end subroutine print_orientation

!--------------------------------------------------------------------------
!
! SUBROUTINE: print_orientation_d
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  prints a complete orientationtype record or a single entry (double precision)
! 
!> @param o orientationtype record
!> @param outtype (optional) indicates which representation to print
!> @param pretext (optional) up to 10 characters that will precede each line
!
!> @date  8/4/13   MDG 1.0 original
!--------------------------------------------------------------------------
subroutine print_orientation_d(o,outtype,pretext)

use local
use io
use constants

IMPLICIT NONE

type(orientationtyped),INTENT(IN)	:: o
character(2),INTENT(IN),OPTIONAL	:: outtype
character(10),INTENT(IN),OPTIONAL	:: pretext

real(kind=dbl)				:: ioreal(4)
character(10)				:: pret

pret = ''
if (present(pretext)) pret=trim(pretext)

if (present(outtype)) then
  select case (outtype)
  	case ('eu')
	  ioreal(1:3) = o%eulang(1:3)*180.D0/cPi
  	  call WriteValue(trim(pret)//'Euler angles			: ', ioreal, 3, "(3(F12.7,' '))")

  	case ('ax')
	  ioreal(1:4) = o%axang(1:4)
	  ioreal(4) = ioreal(4)*180.D0/cPi
	  call WriteValue(trim(pret)//'Axis angle pair [n; angle]	: ', ioreal, 4, "(3(F12.7,' '),'; ',F12.7)")

  	case ('ro')
	  ioreal(1:3) = o%rodrigues(1:3)
  	  call WriteValue(trim(pret)//'Rodigues vector		: ', ioreal, 3, "(3(F12.7,' '))")

  	case ('ho')
	  ioreal(1:3) = o%homochoric(1:3)
  	  call WriteValue(trim(pret)//'Homochoric representation	: ', ioreal, 3, "(3(F12.7,' '))")

  	case ('cu')
  	  ioreal(1:3) = o%cubochoric(1:3)
  	  call WriteValue(trim(pret)//'Cubochoric representation	: ', ioreal, 3, "(3(F12.7,' '))")

  	case ('qu')
	  ioreal(1:4) = o%quat
	  call WriteValue(trim(pret)//'Quaternion			: ', ioreal, 4, "(4(F12.7,' '))")

  	case ('om')
	  ioreal(1:3) = o%om(1,1:3)
	  call WriteValue('					  /', ioreal, 3, "(2(F8.4,' '),F8.4,' \')")
	  ioreal(1:3) = o%om(2,1:3)
	  call WriteValue(trim(pret)//'Orientation Matrix		: |', ioreal, 3, "(2(F8.4,' '),F8.4,' |')")
	  ioreal(1:3) = o%om(3,1:3)
	  call WriteValue('					  \', ioreal, 3, "(2(F8.4,' '),F8.4,' /')")
  		  		
  end select
else
! print the entire record with all representations
  ioreal(1:3) = o%eulang(1:3)*180.D0/cPi
  call WriteValue(trim(pret)//'Euler angles			: ', ioreal, 3, "(3(F12.7,' '))")
  ioreal(1:4) = o%axang(1:4)
  ioreal(4) = ioreal(4)*180.D0/cPi
  call WriteValue(trim(pret)//'Axis angle pair [n; angle]	: ', ioreal, 4, "(3(F12.7,' '),'; ',F12.7)")
  ioreal(1:3) = o%rodrigues(1:3)
  call WriteValue(trim(pret)//'Rodigues vector		: ', ioreal, 3, "(3(F12.7,' '))")
  ioreal(1:3) = o%homochoric(1:3)
  call WriteValue(trim(pret)//'Homochoric representation	: ', ioreal, 3, "(3(F12.7,' '))")
  ioreal(1:3) = o%cubochoric(1:3)
  call WriteValue(trim(pret)//'Cubochoric representation	: ', ioreal, 3, "(3(F12.7,' '))")
  ioreal(1:4) = o%quat
  call WriteValue(trim(pret)//'Quaternion			: ', ioreal, 4, "(4(F12.7,' '))")
  ioreal(1:3) = o%om(1,1:3)
  call WriteValue('					  /', ioreal, 3, "(2(F8.4,' '),F8.4,' \')")
  ioreal(1:3) = o%om(2,1:3)
  call WriteValue(trim(pret)//'Orientation Matrix		: |', ioreal, 3, "(2(F8.4,' '),F8.4,' |')")
  ioreal(1:3) = o%om(3,1:3)
  call WriteValue('					  \', ioreal, 3, "(2(F8.4,' '),F8.4,' /')")
end if

end subroutine print_orientation_d



end module rotations
