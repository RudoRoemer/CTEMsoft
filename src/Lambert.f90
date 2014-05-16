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
! CTEMsoft2013:Lambert.f90
!--------------------------------------------------------------------------
!
! MODULE: Lambert
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief everything that has to do with the modified Lambert projections
!
!> @details This module contains a number of projection functions for the modified 
!> Lambert projection between square lattice and 2D hemisphere, hexagonal lattice
!> and 2D hemisphere, as well as the more complex mapping between a 3D cubic grid
!> and the unit quaternion hemisphere with positive scalar comoonent.  In addition, there 
!> are some other projections, such as the stereographic one.  Each function is named
!> by the projection, the dimensionality of the starting grid, and the forward or inverse
!> character.  For each function, there is also a single precision and a double precision
!> version, but we use the interface formalism to have only a single call.  The Forward
!> mapping is taken to be the one from the simple grid to the curved grid.  Since the module
!> deals with various grids, we also add a few functions/subroutines that apply symmetry
!> operations on those grids.
!>
!> To use any of these routines, one must first call the InitLambertParameters routine
!> to initialize a bunch of constants.
!
!> @date 7/10/13   MDG 1.0 original
!> @date 7/12/13   MDG 1.1 added forward cube to ball to quaternion mappings
!> @date 8/01/13   MDG 1.2 added standard Lambert projection
!> @date 8/12/13   MDG 1.3 added inverse Lambert projections for Ball to Cube
!> @date 9/20/13   MDG 1.4 added ApplyLaueSymmetry
!--------------------------------------------------------------------------
module Lambert

use local
use quaternions

IMPLICIT NONE

! these are a bunch of constants used for the projections; they are all in double precision
type LambertParameters
	real(kind=dbl)		:: Pi    	!  pi
	real(kind=dbl)		:: iPi    	!  1/pi
	real(kind=dbl)		:: sPi		!  sqrt(pi)
	real(kind=dbl)		:: sPio2	!  sqrt(pi/2)
	real(kind=dbl)		:: sPi2		!  sqrt(pi)/2
	real(kind=dbl)		:: srt    	!  sqrt(3)/2
	real(kind=dbl)		:: isrt    	!  1/sqrt(3)
	real(kind=dbl)		:: alpha    	!  sqrt(pi)/3^(1/4)
	real(kind=dbl)		:: rtt  	!  sqrt(3)
	real(kind=dbl)		:: prea    	!  3^(1/4)/sqrt(2pi)
	real(kind=dbl)		:: preb   	!  3^(1/4)sqrt(2/pi)
	real(kind=dbl)		:: prec    	!  pi/2sqrt(3)
	real(kind=dbl)		:: pred   	!  2pi/3
	real(kind=dbl)		:: pree   	!  3^(-1/4)
	real(kind=dbl)		:: pref		!  sqrt(6/pi)
! the following constants are used for the cube to quaternion hemisphere mapping
	real(kind=dbl)		:: a		! pi^(5/6)/6^(1/6)
	real(kind=dbl)		:: ap		! pi^(2/3)
	real(kind=dbl)		:: sc		! a/ap
	real(kind=dbl)		:: beta		! pi^(5/6)/6^(1/6)/2
	real(kind=dbl)		:: R1		! (3pi/4)^(1/3)
	real(kind=dbl)		:: r2		! sqrt(2)
	real(kind=dbl)		:: pi12		! pi/12
	real(kind=dbl)		:: prek		! R1 2^(1/4)/beta
	real(kind=dbl)		:: r24		! sqrt(24)
	real(kind=dbl)		:: tfit(7)	! fit parameters
end type LambertParameters

type(LambertParameters)	:: LPs

!------------
! public functions and subroutines
!------------

! this routine should be called to init all the Lambert parameters
public::InitLambertParameters

! mappings from 2D square grid to the Northern hemisphere of a 2D sphere
public :: LambertSquareToSphere
interface LambertSquareToSphere
	module procedure Lambert2DSquareForwardSingle
	module procedure Lambert2DSquareForwardDouble
end interface

public :: LambertSphereToSquare
interface LambertSphereToSquare
	module procedure Lambert2DSquareInverseSingle
	module procedure Lambert2DSquareInverseDouble
end interface


! mappings from 2D hexagonal grid to the Northern hemisphere of a 2D sphere
public :: LambertHexToSphere
interface LambertHexToSphere
	module procedure Lambert2DHexForwardSingle
	module procedure Lambert2DHexForwardDouble
end interface

public :: LambertSphereToHex
interface LambertSphereToHex
	module procedure Lambert2DHexInverseSingle
	module procedure Lambert2DHexInverseDouble
end interface


! auxiliary private functions for the hexagonal mappings
private :: GetSextantSingle
private :: GetSextantDouble


! mappings from the 3D cubic grid to the 3D spherical grid
public :: LambertCubeToBall
interface LambertCubeToBall
	module procedure Lambert3DCubeForwardSingle
	module procedure Lambert3DCubeForwardDouble
end interface

public :: LambertBallToCube
interface LambertBallToCube
	module procedure Lambert3DCubeInverseSingle
	module procedure Lambert3DCubeInverseDouble
end interface


! mappings from the 3D spherical grid to the unit quaternion sphere
public :: LambertBallToQuaternion
interface LambertBallToQuaternion
	module procedure Lambert3DBallToQuaternionSingle
	module procedure Lambert3DBallToQuaternionDouble
end interface

!public :: LambertQuaternionToBall
!interface LambertQuaternionToBall
!	module procedure Lambert3DBallToQuaternionSingleInverse
!	module procedure Lambert3DBallToQuaternionDoubleInverse
!end interface


! mappings from the 3D cube grid to the unit quaternion sphere
public :: LambertCubeToQuaternion
interface LambertCubeToQuaternion
	module procedure Lambert3DCubeToQuaternionSingle
	module procedure Lambert3DCubeToQuaternionDouble
end interface

!public :: LambertQuaternionToCube
!interface LambertQuaternionToCube
!	module procedure Lambert3DCubeToQuaternionSingleInverse
!	module procedure Lambert3DCubeToQuaternionDoubleInverse
!end interface
!

! auxiliary private functions for the cube to sphere mappings
private :: GetPyramidSingle
private :: GetPyramidDouble


! here we also add the simple stereographic projection
public :: StereoGraphicForward
interface StereoGraphicForward
	module procedure StereoGraphicForwardSingle
	module procedure StereoGraphicForwardDouble
end interface

public :: StereoGraphicInverse
interface StereoGraphicInverse
	module procedure StereoGraphicInverseSingle
	module procedure StereoGraphicInverseDouble
end interface

! as well as the simple Lambert projection
public :: LambertForward
interface LambertForward
	module procedure LambertForwardSingle
	module procedure LambertForwardDouble
end interface

public :: LambertInverse
interface LambertInverse
	module procedure LambertInverseSingle
	module procedure LambertInverseDouble
end interface

contains

!--------------------------------------------------------------------------
!
! SUBROUTINE: InitLambertParameters
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief initialize the various constants used for Lambert projections
!
!> @note We could define these parameters explicitly, but this way they are
!> computed according to machine precision.
!
!> @date 7/10/13   MDG 1.0 original
!--------------------------------------------------------------------------
subroutine InitLambertParameters

use local

IMPLICIT NONE

real(kind=dbl)		:: dpi

dpi = 4.D0*datan(1.D0)

LPs%Pi = dpi
LPs%iPi = 1.D0/dpi
LPs%sPi = dsqrt(dpi)
LPs%sPio2 = dsqrt(dpi*0.5D0)
LPs%sPi2 = LPs%sPi*0.5D0
LPs%srt = dsqrt(3.D0)/2.D0
LPs%isrt = 1.D0/dsqrt(3.D0)
LPs%alpha = dsqrt(dpi)/3.D0**(0.25D0)
LPs%rtt = dsqrt(3.D0)
LPs%prea = 3.D0**(0.25D0)/dsqrt(2.D0*dpi)
LPs%preb = 3.D0**(0.25D0)*dsqrt(2.D0/dpi)
LPs%prec = dpi*0.5D0/dsqrt(3.D0)
LPs%pred = 2.D0*dpi/3.D0
LPs%pree = 3.D0**(-0.25D0)
LPs%pref = dsqrt(6.D0/dpi)

LPs%a = dpi**(5.D0/6.D0)/6.D0**(1.D0/6.D0)
LPs%ap = dpi**(2.D0/3.D0)
LPs%sc = LPs%a/LPs%ap
LPs%beta = 0.5D0*LPs%a
LPs%R1 = (3.D0*dpi/4.D0)**(1.D0/3.D0)
LPs%r2 = dsqrt(2.D0)
LPs%pi12 = dpi/12.D0
LPs%prek =  LPs%R1*2.D0**(0.25D0)/LPs%beta
LPs%r24 = dsqrt(24.D0)

! fit parameters determined with Mathematica
LPs%tfit = (/ -0.5000096149170321D0, -0.02486606148871731D0, &
              -0.004549381779362819D0, 0.0005118668366387526D0, &
              -0.0016500827333575548D0, 0.0007593352203388718D0, &
              -0.0002040422502566876D0 /)

end subroutine InitLambertParameters

!--------------------------------------------------------------------------
!
! FUNCTION: Lambert2DSquareForwardSingle
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief forward projection from 2D square to 3D hemisphere, single precision
!
!> @param xy 2D coordinates to be transformed (single precision)  
!> @param ierr error status: 0=OK, 1=input point lies outside square grid bounds
!
!> @date 7/10/13   MDG 1.0 original
!--------------------------------------------------------------------------
recursive function Lambert2DSquareForwardSingle(xy,ierr) result(res)

use local

IMPLICIT NONE

real(kind=sgl),INTENT(IN)		:: xy(2)
integer(kind=irg),INTENT(INOUT)	:: ierr
real(kind=sgl)				:: res(3), q, qq
real(kind=sgl),parameter		:: eps = 1.0E-6

ierr = 0
! check to make sure that the input point lies inside the square of edge length 2 sqrt(pi/2)
if ((maxval(abs(xy))-LPs%sPio2).gt.eps) then
  res = (/ 0.0, 0.0, 0.0 /)		
  ierr = 1
else
! Forward projection from square grid to Northern hemisphere.
! Equations (8) and (9) from D. Rosca, "New uniform grids on the sphere," 
! Astronomy & Astrophysics, 520, A63 (2010)

! deal with the origin first:
 if (maxval(abs(xy)).eq.0.0) then 
   res = (/ 0.0, 0.0, 1.0 /)
 else
  if (abs(xy(1)).le.abs(xy(2))) then
   q = 2.0*xy(2)*LPs%iPi*sqrt(LPs%Pi-xy(2)*xy(2))
   qq = xy(1)*LPs%Pi*0.25/xy(2)
   res = (/ q*sin(qq), q*cos(qq), 1.0-2.0*xy(2)*xy(2)*sngl(LPs%iPi) /)  
  else
   q = 2.0*xy(1)*LPs%iPi*sqrt(LPs%Pi-xy(1)*xy(1))
   qq = xy(2)*LPs%Pi*0.25/xy(1)
   res = (/ q*cos(qq), q*sin(qq), 1.0-2.0*xy(1)*xy(1)*sngl(LPs%iPi) /)  
  end if
 end if
end if

end function Lambert2DSquareForwardSingle


!--------------------------------------------------------------------------
!
! FUNCTION: Lambert2DSquareForwardDouble
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief forward projection from 2D square to 3D hemisphere, double precision
!
!> @param xy 2D coordinates to be transformed (double precision)  
!> @param ierr error status: 0=OK, 1=input point lies outside square grid bounds
! 
!> @date 7/10/13   MDG 1.0 original
!--------------------------------------------------------------------------
recursive function Lambert2DSquareForwardDouble(xy,ierr) result(res)

use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN)		:: xy(2)
integer(kind=irg),INTENT(INOUT)	:: ierr
real(kind=dbl)				:: res(3), q, qq
real(kind=dbl),parameter		:: eps = 1.0D-12

ierr = 0
! check to make sure that the input point lies inside the square of edge length 2 sqrt(pi)
if ((maxval(dabs(xy))-LPs%sPio2).gt.eps) then
  res = (/ 0.D0, 0.D0, 0.D0 /)
  ierr = 1   ! input point does not lie inside square with edge length 2 sqrt(pi/2)
else
! Forward projection from square grid to Northern hemisphere.
! Equations (8) and (9) from D. Rosca, "New uniform grids on the sphere," 
! Astronomy & Astrophysics, 520, A63 (2010)

! deal with the origin first:
 if (maxval(abs(xy)).eq.0.0) then 
   res = (/ 0.D0, 0.D0, 1.D0 /)
 else 
  if (dabs(xy(1)).le.dabs(xy(2))) then
   q = 2.D0*xy(2)*LPs%iPi*dsqrt(LPs%Pi-xy(2)*xy(2))
   qq = xy(1)*LPs%Pi*0.25D0/xy(2)
   res = (/ q*dsin(qq), q*dcos(qq), 1.D0-2.D0*xy(2)*xy(2)*LPs%iPi /)  
  else
   q = 2.D0*xy(1)*LPs%iPi*dsqrt(LPs%Pi-xy(1)*xy(1))
   qq = xy(2)*LPs%Pi*0.25D0/xy(1)
   res = (/ q*dcos(qq), q*dsin(qq), 1.D0-2.D0*xy(1)*xy(1)*LPs%iPi /)  
  end if
 end if
end if

end function Lambert2DSquareForwardDouble


!--------------------------------------------------------------------------
!
! FUNCTION: Lambert2DSquareInverseSingle
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief inverse projection from 3D hemisphere to 2D square, single precision
!
!> @param xyz 3D coordinates to be transformed (single precision)  
!> @param ierr error status: 0=OK, 1=input point has norm different from 1
!
!> @date 7/10/13   MDG 1.0 original
!--------------------------------------------------------------------------
recursive function Lambert2DSquareInverseSingle(xyz,ierr) result(res)

use local

IMPLICIT NONE

real(kind=sgl),INTENT(IN)		:: xyz(3)
integer(kind=irg),INTENT(INOUT)	:: ierr
real(kind=sgl)				:: res(2), q
real(kind=sgl),parameter		:: eps = 1.0E-6

ierr = 0
! check to make sure that the input point lies on the unit sphere
if (abs(1.0-sum(xyz**2)).gt.eps) then
  res = (/ 0.0, 0.0 /)
  ierr = 1
else
! intercept the points (0,0,+-1)
  if (abs(xyz(3)).eq.1.0) then 
    res = (/ 0.0, 0.0 /)
  else
    if (abs(xyz(2)).le.abs(xyz(1))) then
      q = abs(xyz(1))/xyz(1) * sqrt(2.0*(1.0-xyz(3)))
      res = (/ q * sngl(LPs%sPi2), q * atan(xyz(2)/xyz(1))/sngl(LPs%sPi2) /)
    else
      q = abs(xyz(2))/xyz(2) * sqrt(2.0*(1.0-xyz(3)))
      res = (/  q * atan(xyz(1)/xyz(2))/sngl(LPs%sPi2), q * sngl(LPs%sPi2) /)
    end if
  end if
end if

end function Lambert2DSquareInverseSingle


!--------------------------------------------------------------------------
!
! FUNCTION: Lambert2DSquareInverseDouble
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief inverse projection from 3D hemisphere to 2D square, double precision
!
!> @param xyz 3D coordinates to be transformed (double precision)  
!> @param ierr error status: 0=OK, 1=input point has norm different from 1
!
!> @date 7/10/13   MDG 1.0 original
!--------------------------------------------------------------------------
recursive function Lambert2DSquareInverseDouble(xyz,ierr) result(res)

use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN)		:: xyz(3)
integer(kind=irg),INTENT(INOUT)	:: ierr
real(kind=dbl)				:: res(2), q
real(kind=dbl),parameter		:: eps = 1.0D-12

ierr = 0
! check to make sure that the input point lies on the unit sphere
if (dabs(1.D0-sum(xyz**2)).gt.eps) then
  res = (/ 0.D0, 0.D0 /)
  ierr = 1
else
! intercept the points (0,0,+-1)
  if (dabs(xyz(3)).eq.1.D0) then 
    res = (/ 0.D0, 0.D0 /)
  else
    if (dabs(xyz(2)).le.dabs(xyz(1))) then
      q = dabs(xyz(1))/xyz(1) * dsqrt(2.D0*(1.D0-xyz(3)))
      res = (/ q * LPs%sPi2, q * datan(xyz(2)/xyz(1))/LPs%sPi2 /)
    else
      q = dabs(xyz(2))/xyz(2) * dsqrt(2.D0*(1.D0-xyz(3)))
      res = (/  q * datan(xyz(1)/xyz(2))/LPs%sPi2, q * LPs%sPi2 /)
    end if
  end if
end if

end function Lambert2DSquareInverseDouble

!--------------------------------------------------------------------------
! the functions below deal with the hexagonal to 2D hemisphere projection
!
! all derivations and equations can be found in 
!
! D. Rosca and M. De Graef, "Area-preserving projections from hexagonal and triangular
! domains to the sphere and applications to electron back-scatter diffraction pattern simulations,"
! Modelling Simul. Mater. Sci. Eng. 21 (2013) 055021.
!--------------------------------------------------------------------------


!--------------------------------------------------------------------------
!
! FUNCTION: Lambert2DHexForwardSingle
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief forward projection from 2D hexagon to 3D hemisphere, single precision
!
!> @param xy 2D coordinates to be transformed (single precision)  
!> @param ierr error status: 0=OK, 1=input point outside hexagon
!
!> @date 7/10/13   MDG 1.0 original
!--------------------------------------------------------------------------
recursive function Lambert2DHexForwardSingle(xy,ierr) result(res)

use local

IMPLICIT NONE

real(kind=sgl),INTENT(IN)	:: xy(2) 
integer(kind=irg),INTENT(INOUT):: ierr

real(kind=sgl)			:: res(3), q, XX, YY, xp, yp
integer(kind=irg)		:: ks

 ierr = 0
  
 if (maxval(abs(xy)).eq.0.0) then
  res = (/ 0.0, 0.0, 1.0 /)
 else
! determine in which sextant this point lies
  ks = GetSextantSingle(xy)
  
  select case (ks)
    case (0,3)
  	XX = LPs%preb*xy(1)*cos(xy(2)*LPs%prec/xy(1))
  	YY = LPs%preb*xy(1)*sin(xy(2)*LPs%prec/xy(1))
    case (1,4)
  	xp = xy(1)+LPs%rtt*xy(2)
  	yp = xy(1)*LPs%pred/xp
  	XX = LPs%prea*xp*sin(yp)
  	YY = LPs%prea*xp*cos(yp)
    case (2,5)
  	xp = xy(1)-LPs%rtt*xy(2)
  	yp = xy(1)*LPs%pred/xp
  	XX = LPs%prea*xp*sin(yp)
  	YY = -LPs%prea*xp*cos(yp)	  
  end select
  q = XX**2+YY**2
! does the point lie outside the hexagon ?
  if (q.gt.4.0) then
    res = (/ 0.0, 0.0, 0.0 /)
    ierr = 1
  else
    res = (/ 0.50*XX*sqrt(4.0-q), 0.50*YY*sqrt(4.0-q), 1.0-0.5*q /)
  end if
 end if
 
end function Lambert2DHexForwardSingle

!--------------------------------------------------------------------------
!
! FUNCTION: Lambert2DHexForwardDouble
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief forward projection from 2D hexagon to 3D hemisphere, double precision
!
!> @param xy 2D coordinates to be transformed (double precision)  
!> @param ierr error status: 0=OK, 1=input point outside hexagon
!
!> @date 7/10/13   MDG 1.0 original
!--------------------------------------------------------------------------
recursive function Lambert2DHexForwardDouble(xy,ierr) result(res)

use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN)	:: xy(2) 
integer(kind=irg),INTENT(INOUT):: ierr

real(kind=dbl)			:: res(3), q, XX, YY, xp, yp
integer(kind=irg)		:: ks

  ierr = 0

 if (maxval(dabs(xy)).eq.0.D0) then
  res = (/ 0.D0, 0.D0, 1.D0 /)
 else
! determine in which sextant this point lies
  ks = GetSextantDouble(xy)
  
  select case (ks)
    case (0,3)
  	XX = LPs%preb*xy(1)*dcos(xy(2)*LPs%prec/xy(1))
  	YY = LPs%preb*xy(1)*dsin(xy(2)*LPs%prec/xy(1))
    case (1,4)
  	xp = xy(1)+LPs%rtt*xy(2)
  	yp = xy(1)*LPs%pred/xp
  	XX = LPs%prea*xp*dsin(yp)
  	YY = LPs%prea*xp*dcos(yp)
    case (2,5)
  	xp = xy(1)-LPs%rtt*xy(2)
  	yp = xy(1)*LPs%pred/xp
  	XX = LPs%prea*xp*dsin(yp)
  	YY = -LPs%prea*xp*dcos(yp)	  
  end select
  q = XX**2+YY**2
! does the point lie outside the hexagon ?
  if (q.gt.4.D0) then
    res = (/ 0.D0, 0.D0, 0.D0 /)
    ierr = 1
  else
    res = (/ 0.5D0*XX*dsqrt(4.D0-q), 0.5D0*YY*dsqrt(4.D0-q), 1.D0-0.5D0*q /)
  end if
 end if
 
end function Lambert2DHexForwardDouble


!--------------------------------------------------------------------------
!
! FUNCTION: Lambert2DHexInverseSingle
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief inverse projection from 3D hemisphere to 2D hexagon, single precision
!
!> @param xyz 3D coordinates to be transformed (single precision)  
!> @param ierr error status: 0=OK, 1=input point not normalized
!
!> @date 7/10/13   MDG 1.0 original
!--------------------------------------------------------------------------
recursive function Lambert2DHexInverseSingle(xyz,ierr) result(res)

use local

IMPLICIT NONE

real(kind=sgl),INTENT(IN)	:: xyz(3) 
integer(kind=irg),INTENT(INOUT):: ierr

real(kind=sgl)			:: res(2), q, qq, XX, YY, lxx, lyy
integer(kind=irg)		:: ks
real(kind=sgl),parameter	:: eps = 1.0E-7

ierr = 0
! check to make sure that the input point lies on the unit sphere
if (abs(1.0-sum(xyz**2)).gt.eps) then
  res = (/ 0.0, 0.0 /)
  ierr = 1
else
 if (abs(xyz(3)).eq.1.0) then
  res = (/ 0.0, 0.0 /)
 else
! first do the Lambert projection
  q = sqrt(2.0/(1.0+xyz(3)))
  XX = q * xyz(1)
  YY = q * xyz(2)

! determine in which sextant this point lies
  ks = GetSextantSingle( (/ XX, YY /) )

! then perform the inverse to the hexagonal grid
  select case (ks)
    case (0,3)
        q = LPs%pree * (abs(XX)/XX) * sqrt(XX**2+YY**2)
  	lxx = q * LPs%sPio2 
  	lyy = q * LPs%pref * atan(YY/XX)
    case (1,4)
    	q = LPs%prea * (abs(XX)/XX) * sqrt(XX**2+YY**2)
    	qq= atan((YY-LPs%rtt*XX)/(XX+LPs%rtt*YY))
  	lxx = q * LPs%rtt *( LPs%Pi/6.0 - qq )
  	lyy = q * ( 0.5*LPs%Pi + qq )
    case (2,5)
    	q = LPs%prea * (abs(XX)/XX) * sqrt(XX**2+YY**2)
    	qq= atan((YY+LPs%rtt*XX)/(XX-LPs%rtt*YY))
  	lxx = q * LPs%rtt *( LPs%Pi/6.0 + qq )
  	lyy = q * ( -0.5*LPs%Pi + qq )
  end select
  res = (/ lxx, lyy /)
 end if
end if

end function Lambert2DHexInverseSingle


!--------------------------------------------------------------------------
!
! FUNCTION: Lambert2DHexInverseDouble
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief inverse projection from 3D hemisphere to 2D hexagon, double precision
!
!> @param xyz 3D coordinates to be transformed (double precision)  
!> @param ierr error status: 0=OK, 1=input point not normalized
!
!> @date 7/10/13   MDG 1.0 original
!--------------------------------------------------------------------------
recursive function Lambert2DHexInverseDouble(xyz,ierr) result(res)

use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN)	:: xyz(3) 
integer(kind=irg),INTENT(INOUT):: ierr

real(kind=dbl)			:: res(2), q, qq, XX, YY, lxx, lyy
integer(kind=irg)		:: ks
real(kind=dbl),parameter	:: eps = 1.0E-12

ierr = 0
! check to make sure that the input point lies on the unit sphere
if (dabs(1.0-sum(xyz**2)).gt.eps) then
  res = (/ 0.D0, 0.D0 /)
  ierr = 1
else
 if (dabs(xyz(3)).eq.1.D0) then
  res = (/ 0.D0, 0.D0 /)
 else
! first do the Lambert projection
  q = dsqrt(2.D0/(1.D0+xyz(3)))
  XX = q * xyz(1)
  YY = q * xyz(2)

! determine in which sextant this point lies
  ks = GetSextantDouble( (/ XX, YY /) )

! then perform the inverse to the hexagonal grid
  select case (ks)
    case (0,3)
        q = LPs%pree * (dabs(XX)/XX) * dsqrt(XX**2+YY**2)
  	lxx = q * LPs%sPio2 
  	lyy = q * LPs%pref * datan(YY/XX)
    case (1,4)
    	q = LPs%prea * (dabs(XX)/XX) * dsqrt(XX**2+YY**2)
    	qq= datan((YY-LPs%rtt*XX)/(XX+LPs%rtt*YY))
  	lxx = q * LPs%rtt *( LPs%Pi/6.D0 - qq )
  	lyy = q * ( 0.5D0*LPs%Pi + qq )
    case (2,5)
    	q = LPs%prea * (dabs(XX)/XX) * dsqrt(XX**2+YY**2)
    	qq= datan((YY+LPs%rtt*XX)/(XX-LPs%rtt*YY))
  	lxx = q * LPs%rtt *( LPs%Pi/6.D0 + qq )
  	lyy = q * ( -0.5D0*LPs%Pi + qq )
  end select
    res = (/ lxx, lyy /)
end if
end if

end function Lambert2DHexInverseDouble

!--------------------------------------------------------------------------
!
! FUNCTION: GetSextantSingle
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief determine to which sextant a point in hexagonal coordinates belongs
!
!> @param xy 2D coordinates to be considered (single precision)  
! 
!> @date 11/21/12    MDG 1.0 original
!--------------------------------------------------------------------------
recursive function GetSextantSingle(xy) result(res)

use local

IMPLICIT NONE

real(kind=sgl),INTENT(IN)	:: xy(2) 
integer(kind=irg)		:: res

real(kind=sgl)			:: xx

xx = abs(xy(1)*LPs%isrt)    	! |x| / sqrt(3)

if (xy(1).gt.0.0) then 
  if (abs(xy(2)).le.xx) then
    res = 0
  else 
    if (xy(2).gt.xx) then
      res = 1
    else
      res = 5
    end if
  end if
else
  if (abs(xy(2)).le.xx) then
    res = 3
  else 
    if (xy(2).gt.xx) then
      res = 2
    else
      res = 4
    end if
  end if
end if

end function GetSextantSingle

!--------------------------------------------------------------------------
!
! FUNCTION: GetSextantDouble
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief determine to which sextant a point in hexagonal coordinates belongs
!
!> @param xy 2D coordinates to be considered (double precision)  
! 
!> @date 11/21/12    MDG 1.0 original
!--------------------------------------------------------------------------
recursive function GetSextantDouble(xy) result(res)

use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN)	:: xy(2) 
integer(kind=irg)		:: res

real(kind=dbl)			:: xx

xx = dabs(xy(1)*LPs%isrt)    	! |x| / sqrt(3)

if (xy(1).gt.0.D0) then 
  if (dabs(xy(2)).le.xx) then
    res = 0
  else 
    if (xy(2).gt.xx) then
      res = 1
    else
      res = 5
    end if
  end if
else
  if (dabs(xy(2)).le.xx) then
    res = 3
  else 
    if (xy(2).gt.xx) then
      res = 2
    else
      res = 4
    end if
  end if
end if

end function GetSextantDouble


!--------------------------------------------------------------------------
! the functions below deal with the cubic grid to the 3D ball, and then from
! the 3D ball to the unit quaternion hemisphere projection
!
! all derivations and equations can be found in 
!
! D. Rosca and M. De Graef, ... in preparation
!--------------------------------------------------------------------------


!--------------------------------------------------------------------------
!
! FUNCTION: Lambert3DCubeForwardSingle
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief map from 3D cubic grid to 3D ball 
!
!> @param lxyz 3D coordinates to be considered (single precision)  
!> @param ierr error flag 0 = OK, 1 = 
! 
!> @date 7/12/13    MDG 1.0 original
!--------------------------------------------------------------------------
recursive function Lambert3DCubeForwardSingle(lxyz,ierr) result(res)

use local
use constants

IMPLICIT NONE

real(kind=sgl),INTENT(IN)	:: lxyz(3)
integer(kind=irg),INTENT(INOUT):: ierr
real(kind=sgl)			:: res(3)

real(kind=sgl)			:: XYZ(3), sXYZ(3), T1, T2, c, s, q, LamXYZ(3), edge
integer(kind=irg)		:: p

ierr = 0
edge = 0.50 * (sngl(cPi))**(2.0/3.0)

if (maxval(dabs(lxyz)).gt.edge) then
  res = (/ 0.0, 0.0, 0.0 /)
  ierr = 1
  return
end if

! determine which pyramid pair the point lies in and copy coordinates in correct order (see paper)
p = GetPyramidSingle(lxyz)
select case (p)
 case (1,2)
  sXYZ = lxyz
 case (3,4)
  sXYZ = (/ lxyz(2), lxyz(3), lxyz(1) /)
 case (5,6)
  sXYZ = (/ lxyz(3), lxyz(1), lxyz(2) /)
end select

! scale by grid parameter ratio sc
XYZ = LPs%sc * sXYZ

! transform to the sphere grid via the curved square, and intercept the zero point
if (maxval(abs(XYZ)).eq.0.0) then
  LamXYZ = (/ 0.0, 0.0, 0.0 /)
else
! intercept all the points along the z-axis
  if (maxval(abs(XYZ(1:2))).eq.0.0) then
    LamXYZ = (/ 0.0, 0.0, sngl(LPs%pref) * XYZ(3) /)
  else  ! this is a general grid point
    if (abs(XYZ(2)).le.abs(XYZ(1))) then
      q = LPs%pi12 * XYZ(2)/XYZ(1)
      c = cos(q)
      s = sin(q)
      q = LPs%prek * XYZ(1) / sqrt(LPs%r2-c)
      T1 = (LPs%r2*c - 1.0) * q
      T2 = LPs%r2 * s * q
    else
      q = LPs%pi12 * XYZ(1)/XYZ(2)
      c = cos(q)
      s = sin(q)
      q = LPs%prek * XYZ(2) / sqrt(LPs%r2-c)
      T1 = LPs%r2 * s * q
      T2 = (LPs%r2*c - 1.0) * q
    end if

! transform to sphere grid (inverse Lambert)
! [note that there is no need to worry about dividing by zero, since XYZ(3) can not become zero]
    c = T1**2+T2**2
    s = LPs%Pi * c/(24.0*XYZ(3)**2)
    c = LPs%sPi * c / LPs%r24 / XYZ(3)
    q = sqrt( 1.0 - s )
    LamXYZ = (/ T1 * q, T2 * q, sngl(LPs%pref) * XYZ(3) - c /)
  end if
end if

! reverse the coordinates back to the regular order according to the original pyramid number
select case (p)
 case (1,2)
  res = LamXYZ
 case (3,4)
  res = (/ LamXYZ(3), LamXYZ(1), LamXYZ(2) /)
 case (5,6)
  res = (/ LamXYZ(2), LamXYZ(3), LamXYZ(1) /)
end select

end function Lambert3DCubeForwardSingle

!--------------------------------------------------------------------------
!
! FUNCTION: Lambert3DCubeForwardDouble
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief map from 3D cubic grid to 3D ball 
!
!> @param lxyz 3D coordinates to be considered (double precision)  
!> @param ierr error flag 0 = OK, 1 = outside of unit cube
! 
!> @date 7/12/13    MDG 1.0 original
!--------------------------------------------------------------------------
recursive function Lambert3DCubeForwardDouble(lxyz,ierr) result(res)

use local
use constants

IMPLICIT NONE

real(kind=dbl),INTENT(IN)	:: lxyz(3)
integer(kind=irg),INTENT(INOUT):: ierr
real(kind=dbl)			:: res(3)

real(kind=dbl)			:: XYZ(3), sXYZ(3), T1, T2, c, s, q, LamXYZ(3), edge
integer(kind=irg)		:: p

ierr = 0
edge = 0.5D0 * (cPi)**(2.D0/3.D0)
if (maxval(dabs(lxyz)).gt.edge) then
  res = (/ 0.D0, 0.D0, 0.D0 /)
  ierr = 1
  return
end if

! determine which pyramid pair the point lies in and copy coordinates in correct order (see paper)
p = GetPyramidDouble(lxyz)
select case (p)
 case (1,2)
  sXYZ = lxyz
 case (3,4)
  sXYZ = (/ lxyz(2), lxyz(3), lxyz(1) /)
 case (5,6)
  sXYZ = (/ lxyz(3), lxyz(1), lxyz(2) /)
end select

! scale by grid parameter ratio sc
XYZ = LPs%sc * sXYZ

! transform to the sphere grid via the curved square, and intercept the zero point
if (maxval(dabs(XYZ)).eq.0.D0) then
  LamXYZ = (/ 0.D0, 0.D0, 0.D0 /)
else
! intercept all the points along the z-axis
  if (maxval(dabs(XYZ(1:2))).eq.0.D0) then
    LamXYZ = (/ 0.D0, 0.D0, LPs%pref * XYZ(3) /)
  else  ! this is a general grid point
    if (dabs(XYZ(2)).le.dabs(XYZ(1))) then
      c = dcos(LPs%pi12 * XYZ(2)/XYZ(1))
      s = dsin(LPs%pi12 * XYZ(2)/XYZ(1))
      q = LPs%prek * XYZ(1) / dsqrt(LPs%r2-c)
      T1 = (LPs%r2*c - 1.D0) * q
      T2 = LPs%r2 * s * q
    else
      c = dcos(LPs%pi12 * XYZ(1)/XYZ(2))
      s = dsin(LPs%pi12 * XYZ(1)/XYZ(2))
      q = LPs%prek * XYZ(2) / dsqrt(LPs%r2-c)
      T1 = LPs%r2 * s * q
      T2 = (LPs%r2*c - 1.D0) * q
    end if

! transform to sphere grid (inverse Lambert)
! [note that there is no need to worry about dividing by zero, since XYZ(3) can not become zero]
    c = T1**2+T2**2
    s = LPs%Pi * c/(24.D0*XYZ(3)**2)
    c = LPs%sPi * c / LPs%r24 / XYZ(3)
    q = dsqrt( 1.0 - s )
    LamXYZ = (/ T1 * q, T2 * q, LPs%pref * XYZ(3) - c /)
  end if
end if

! reverse the coordinates back to the regular order according to the original pyramid number
select case (p)
 case (1,2)
  res = LamXYZ
 case (3,4)
  res = (/ LamXYZ(3), LamXYZ(1), LamXYZ(2) /)
 case (5,6)
  res = (/ LamXYZ(2), LamXYZ(3), LamXYZ(1) /)
end select

end function Lambert3DCubeForwardDouble

!--------------------------------------------------------------------------
!
! FUNCTION: Lambert3DCubeInverseSingle
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief map from 3D ball to 3D cubic grid  
!
!> @param xyz 3D coordinates to be considered (single precision)  
!> @param ierr error flag 0 = OK, 1 = 
! 
!> @date 7/12/13    MDG 1.0 original
!--------------------------------------------------------------------------
recursive function Lambert3DCubeInverseSingle(xyz,ierr) result(res)

use local

IMPLICIT NONE

real(kind=sgl),INTENT(IN)	:: xyz(3)
integer(kind=irg),INTENT(INOUT):: ierr
real(kind=sgl)			:: res(3)

real(kind=sgl)			:: rs, xyz3(3), xyz2(3), qxy, q2xy, sq2xy, q, ac, T1inv, T2inv, &
				   xyz1(3), sx, sy, qx2y, sqx2y
	
integer(kind=irg)		:: p

ierr = 0

rs = sqrt(sum(xyz*xyz))
if (rs.gt.LPs%R1) then 
  res = (/ 0.0,0.0,0.0 /) 
  return
end if

if (maxval(abs(xyz)).eq.0.0) then 
  res = (/ 0.0,0.0,0.0 /) 
  return
end if

! determine pyramid
p = GetPyramidSingle(xyz)
select case (p)
 case (1,2)
  xyz3 = xyz
 case (3,4)
  xyz3 = (/ xyz(2), xyz(3), xyz(1) /)
 case (5,6)
  xyz3 = (/ xyz(3), xyz(1), xyz(2) /)
end select

! inverse M_3
q = sqrt( 2.0*rs/(rs+abs(xyz3(3))) )
xyz2 = (/ xyz3(1) * q, xyz3(2) * q, (abs(xyz3(3))/xyz3(3)) * rs / sngl(LPs%pref) /)

! inverse M_2
qxy = xyz2(1)*xyz2(1)+xyz2(2)*xyz2(2)
sx = 1.0
if (xyz2(1).ne.0.0)  sx = abs(xyz2(1))/xyz2(1) 
sy = 1.0
if (xyz2(2).ne.0.0)  sy = abs(xyz2(2))/xyz2(2)

if (abs(xyz2(2)).le.abs(xyz2(1))) then 
  q2xy = qxy + xyz2(1)*xyz2(1)
  sq2xy = sqrt(q2xy)
  q = (LPs%beta/LPs%r2/LPs%R1) * sqrt(q2xy*qxy/(q2xy-abs(xyz2(1))*sq2xy))
  ac = acos( (xyz2(2)*xyz2(2)+abs(xyz2(1))*sq2xy)/LPs%r2/qxy )
  T1inv = q * sx
  T2inv = q * sy * ac/LPs%pi12
else 
  qx2y = qxy + xyz2(2)*xyz2(2)
  sqx2y = sqrt(qx2y)
  q = (LPs%beta/LPs%r2/LPs%R1) * sqrt(qx2y*qxy/(qx2y-abs(xyz2(2))*sqx2y))
  ac = acos( (xyz2(1)*xyz2(1)+abs(xyz2(2))*sqx2y)/LPs%r2/qxy )
  T1inv = q * sx * ac/LPs%pi12
  T2inv = q * sy
end if
xyz1 = (/ T1inv, T2inv, xyz2(3) /)

! inverse M_1
xyz1 = xyz1 / LPs%sc

! reverst the coordinates back to the regular order according to the original pyramid number
select case (p)
 case (1,2)
  res = xyz1
 case (3,4)
  res = (/ xyz1(3), xyz1(1), xyz1(2) /)
 case (5,6)
  res = (/ xyz1(2), xyz1(3), xyz1(1) /)
end select

end function Lambert3DCubeInverseSingle

!--------------------------------------------------------------------------
!
! FUNCTION: Lambert3DCubeInverseDouble
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief map from 3D ball to 3D cubic grid  
!
!> @param xyz 3D coordinates to be considered (double precision)  
!> @param ierr error flag 0 = OK, 1 = 
! 
!> @date 7/12/13    MDG 1.0 original
!--------------------------------------------------------------------------
recursive function Lambert3DCubeInverseDouble(xyz,ierr) result(res)

use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN)	:: xyz(3)
integer(kind=irg),INTENT(INOUT):: ierr
real(kind=dbl)			:: res(3)

real(kind=dbl)			:: rs, xyz3(3), xyz2(3), qxy, q2xy, sq2xy, q, ac, T1inv, T2inv, &
				   xyz1(3), sx, sy, qx2y, sqx2y
integer(kind=irg)		:: p

ierr = 0

rs = dsqrt(sum(xyz*xyz))
if (rs.gt.LPs%R1) then 
  res = (/ 0.D0,0.D0,0.D0 /) 
  return
end if

if (maxval(dabs(xyz)).eq.0.D0) then 
  res = (/ 0.D0,0.D0,0.D0 /) 
  return
end if

! determine pyramid
p = GetPyramidDouble(xyz)
select case (p)
 case (1,2)
  xyz3 = xyz
 case (3,4)
  xyz3 = (/ xyz(2), xyz(3), xyz(1) /)
 case (5,6)
  xyz3 = (/ xyz(3), xyz(1), xyz(2) /)
end select

! inverse M_3
q = dsqrt( 2.D0*rs/(rs+dabs(xyz3(3))) )
xyz2 = (/ xyz3(1) * q, xyz3(2) * q, (dabs(xyz3(3))/xyz3(3)) * rs / LPs%pref /)

! inverse M_2
qxy = xyz2(1)*xyz2(1)+xyz2(2)*xyz2(2)
sx = 1.D0
if (xyz2(1).ne.0.D0)  sx = dabs(xyz2(1))/xyz2(1) 
sy = 1.D0
if (xyz2(2).ne.0.D0)  sy = dabs(xyz2(2))/xyz2(2)

if (dabs(xyz2(2)).le.dabs(xyz2(1))) then 
  q2xy = qxy + xyz2(1)*xyz2(1)
  sq2xy = dsqrt(q2xy)
  q = (LPs%beta/LPs%r2/LPs%R1) * dsqrt(q2xy*qxy/(q2xy-dabs(xyz2(1))*sq2xy))
  ac = dacos( (xyz2(2)*xyz2(2)+dabs(xyz2(1))*sq2xy)/LPs%r2/qxy )
  T1inv = q * sx
  T2inv = q * sy * ac/LPs%pi12
else 
  qx2y = qxy + xyz2(2)*xyz2(2)
  sqx2y = dsqrt(qx2y)
  q = (LPs%beta/LPs%r2/LPs%R1) * dsqrt(qx2y*qxy/(qx2y-dabs(xyz2(2))*sqx2y))
  ac = dacos( (xyz2(1)*xyz2(1)+dabs(xyz2(2))*sqx2y)/LPs%r2/qxy )
  T1inv = q * sx * ac/LPs%pi12
  T2inv = q * sy
end if
xyz1 = (/ T1inv, T2inv, xyz2(3) /)

! inverse M_1
xyz1 = xyz1 / LPs%sc

! reverst the coordinates back to the regular order according to the original pyramid number
select case (p)
 case (1,2)
  res = xyz1
 case (3,4)
  res = (/ xyz1(3), xyz1(1), xyz1(2) /)
 case (5,6)
  res = (/ xyz1(2), xyz1(3), xyz1(1) /)
end select

end function Lambert3DCubeInverseDouble


!--------------------------------------------------------------------------
!
! FUNCTION: GetPyramidSingle
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief determine to which pyramid a point in a cubic grid belongs
!
!> @param xyz 3D coordinates to be considered (single precision)  
! 
!> @date 11/21/12    MDG 1.0 original
!--------------------------------------------------------------------------
recursive function GetPyramidSingle(xyz) result(res)

use local

IMPLICIT NONE

real(kind=sgl),INTENT(IN)	:: xyz(3) 
integer(kind=irg)		:: res, p

logical				:: next

next = .TRUE.
if ((abs(xyz(1)).le.xyz(3)).and.(abs(xyz(2)).le.xyz(3))) then
  p = 1				! pyramid 1
  next = .FALSE.
end if  
if (next) then 
 if ((abs(xyz(1)).le.-xyz(3)).and.(abs(xyz(2)).le.-xyz(3))) then
  p = 2				! pyramid 2
  next = .FALSE.
 end if   
end if

if (next) then
 if ((abs(xyz(3)).le.xyz(1)).and.(abs(xyz(2)).le.xyz(1))) then
  p = 3				! pyramid 3
  next = .FALSE.
 end if  
end if
if (next) then
 if ((abs(xyz(3)).le.-xyz(1)).and.(abs(xyz(2)).le.-xyz(1))) then
  p = 4				! pyramid 4
  next = .FALSE.
 end if  
end if

if (next) then
 if ((abs(xyz(1)).le.xyz(2)).and.(abs(xyz(3)).le.xyz(2))) then
  p = 5				! pyramid 5
  next = .FALSE.
 end if  
end if
if (next) then
 if ((abs(xyz(1)).le.-xyz(2)).and.(abs(xyz(3)).le.-xyz(2))) then
  p = 6				! pyramid 6
  next = .FALSE.
 end if  
end if
res = p

end function GetPyramidSingle

!--------------------------------------------------------------------------
!
! FUNCTION: GetPyramidDouble
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief determine to which pyramid a point in a cubic grid belongs
!
!> @param xyz 3D coordinates to be considered (double precision)  
! 
!> @date 11/21/12    MDG 1.0 original
!--------------------------------------------------------------------------
recursive function GetPyramidDouble(xyz) result(res)

use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN)	:: xyz(3) 
integer(kind=irg)		:: res, p

logical				:: next

next = .TRUE.
if ((dabs(xyz(1)).le.xyz(3)).and.(dabs(xyz(2)).le.xyz(3))) then
  p = 1				! pyramid 1
  next = .FALSE.
end if  
if (next) then
 if ((dabs(xyz(1)).le.-xyz(3)).and.(dabs(xyz(2)).le.-xyz(3))) then
  p = 2				! pyramid 2
  next = .FALSE.
 end if   
end if

if (next) then
 if ((dabs(xyz(3)).le.xyz(1)).and.(dabs(xyz(2)).le.xyz(1))) then
  p = 3				! pyramid 3
  next = .FALSE.
 end if 
end if 
if (next) then
 if ((dabs(xyz(3)).le.-xyz(1)).and.(dabs(xyz(2)).le.-xyz(1))) then
  p = 4				! pyramid 4
  next = .FALSE.
 end if  
end if

if (next) then
 if ((dabs(xyz(1)).le.xyz(2)).and.(dabs(xyz(3)).le.xyz(2))) then
  p = 5				! pyramid 5
  next = .FALSE.
 end if  
end if
if (next) then
 if ((dabs(xyz(1)).le.-xyz(2)).and.(dabs(xyz(3)).le.-xyz(2))) then
  p = 6				! pyramid 6
  next = .FALSE.
 end if  
end if

res = p

end function GetPyramidDouble

!---------------
! and finally, here are the 3D cube and ball to 4D unit quaternion sphere mappings
!---------------

!--------------------------------------------------------------------------
!
! FUNCTION: Lambert3DBallToQuaternionSingle
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief map from 3D ball to unit quaternion sphere 
!
!> @param xyz 3D coordinates to be considered (single precision)  
!> @param ierr error flag 0 = OK, 1 = outside of ball
! 
!> @date 7/12/13    MDG 1.0 original
!--------------------------------------------------------------------------
recursive function Lambert3DBallToQuaternionSingle(xyz, ierr) result(res)

use local
use quaternions

IMPLICIT NONE

real(kind=sgl),INTENT(IN)	:: xyz(3)
integer(kind=irg),INTENT(INOUT):: ierr
real(kind=sgl)			:: res(4)

real(kind=sgl)			:: q, x(7), ft, t
integer(kind=irg)		:: j

ierr = 0
q = sqrt(sum(xyz**2))
! is the input point inside the ball that maps onto the unit quaternion sphere ?
if (q.gt.LPs%R1) then 
  res = (/ 0.0,0.0,0.0,0.0 /)
  ierr = 1
  return
end if

if (maxval(abs(xyz)).eq.0.0) then
  res = (/ 1.0,0.0,0.0,0.0 /)
else
! get the value of t
  x(1) = q**2
  do j=2,7
    x(j) = x(j-1) * x(1)
  end do
  t = 1.0 + sum( x * LPs%tfit )

! and get f(t)
  q = sqrt(1.0-t**2)
  ft = (1.5*(acos(t)-t*q))**(1.0/3.0) / q

! and finally create the quaternion
  res = (/ t, xyz(1)/ft, xyz(2)/ft, xyz(3)/ft /)
end if

end function Lambert3DBallToQuaternionSingle

!--------------------------------------------------------------------------
!
! FUNCTION: Lambert3DBallToQuaternionDouble
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief map from 3D ball to unit quaternion sphere 
!
!> @param xyz 3D coordinates to be considered (double precision)  
!> @param ierr error flag 0 = OK, 1 = outside of ball
! 
!> @date 7/12/13    MDG 1.0 original
!--------------------------------------------------------------------------
recursive function Lambert3DBallToQuaternionDouble(xyz, ierr) result(res)

use local
use quaternions

IMPLICIT NONE

real(kind=dbl),INTENT(IN)	:: xyz(3)
integer(kind=irg),INTENT(INOUT):: ierr
real(kind=dbl)			:: res(4)

real(kind=dbl)			:: q, x(7), ft, t
integer(kind=irg)		:: j

ierr = 0
q = dsqrt(sum(xyz**2))
! is the input point inside the ball that maps onto the unit quaternion sphere ?
if (q.gt.LPs%R1) then 
  res = (/ 0.D0,0.D0,0.D0,0.D0 /)
  ierr = 1
  return
end if

if (maxval(dabs(xyz)).eq.0.D0) then
  res = (/ 1.D0,0.D0,0.D0,0.D0 /)
else
! get the value of t
  x(1) = q**2
  do j=2,7
    x(j) = x(j-1) * x(1)
  end do
  t = 1.D0 + sum( x * LPs%tfit )

! and get f(t)
  q = dsqrt(1.D0-t**2)
  ft = (1.5D0*(dacos(t)-t*q))**(1.D0/3.D0) / q

! and finally create the quaternion
  res = (/ t, xyz(1)/ft, xyz(2)/ft, xyz(3)/ft /)
end if

end function Lambert3DBallToQuaternionDouble


!--------------------------------------------------------------------------
!
! FUNCTION: Lambert3DCubeToQuaternionSingle
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief map from 3D cube to unit quaternion sphere 
!
!> @param xyz 3D coordinates to be considered (single precision)  
!> @param ierr error flag 0 = OK, 1 = outside of cube
! 
!> @date 7/12/13    MDG 1.0 original
!--------------------------------------------------------------------------
recursive function Lambert3DCubeToQuaternionSingle(xyz, ierr) result(res)

use local
use quaternions

IMPLICIT NONE

real(kind=sgl),INTENT(IN)	:: xyz(3)
integer(kind=irg),INTENT(INOUT):: ierr
real(kind=sgl)			:: res(4)

real(kind=sgl)			:: q(3)

q = Lambert3DCubeForwardSingle(xyz,ierr)

res = Lambert3DBallToQuaternionSingle(q,ierr)

end function Lambert3DCubeToQuaternionSingle

!--------------------------------------------------------------------------
!
! FUNCTION: Lambert3DCubeToQuaternionDouble
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief map from 3D cube to unit quaternion sphere 
!
!> @param xyz 3D coordinates to be considered (double precision)  
!> @param ierr error flag 0 = OK, 1 = outside of cube
! 
!> @date 7/12/13    MDG 1.0 original
!--------------------------------------------------------------------------
recursive function Lambert3DCubeToQuaternionDouble(xyz, ierr) result(res)

use local
use quaternions

IMPLICIT NONE

real(kind=dbl),INTENT(IN)	:: xyz(3)
integer(kind=irg),INTENT(INOUT):: ierr
real(kind=dbl)			:: res(4)

real(kind=dbl)			:: q(3)

q = Lambert3DCubeForwardDouble(xyz,ierr)

res = Lambert3DBallToQuaternionDouble(q,ierr)


end function Lambert3DCubeToQuaternionDouble


!---------
! the next couple of routine implement the forward and inverse stereographic projection
!---------


!--------------------------------------------------------------------------
!
! FUNCTION: StereoGraphicForwardSingle
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Forward stereographic projection (from unit sphere to plane) 
!
!> @param xyz 3D coordinates to be considered (single precision)  
!> @param ierr error flag 0 = OK, 1 = coordinates not normalized
! 
!> @date 7/12/13    MDG 1.0 original
!--------------------------------------------------------------------------
recursive function StereoGraphicForwardSingle(xyz, ierr, Radius) result(res)

use local

IMPLICIT NONE

real(kind=sgl),INTENT(IN)		:: xyz(3)
integer(kind=irg),INTENT(INOUT)	:: ierr
real(kind=sgl),INTENT(IN),OPTIONAL	:: Radius
real(kind=sgl)				:: res(2)

real(kind=sgl),parameter		:: eps = 1.E-7

! input point must be on unit sphere
ierr = 0
if ( abs(1.0-sum(xyz**2)).gt.eps) then
  res = (/ 0.0, 0.0 /)
  ierr = 1
  return
end if

! projection
res = (/ xyz(1)/(1.0+xyz(3)), xyz(2)/(1.0+xyz(3)) /)

! scale if necessary
if (present(Radius)) res = Radius * res
 
end function StereoGraphicForwardSingle

!--------------------------------------------------------------------------
!
! FUNCTION: StereoGraphicForwardDouble
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Forward stereographic projection (from unit sphere to plane) 
!
!> @param xyz 3D coordinates to be considered (double precision)  
!> @param ierr error flag 0 = OK, 1 = coordinates not normalized
! 
!> @date 7/12/13    MDG 1.0 original
!--------------------------------------------------------------------------
recursive function StereoGraphicForwardDouble(xyz, ierr, Radius) result(res)

use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN)		:: xyz(3)
integer(kind=irg),INTENT(INOUT)	:: ierr
real(kind=dbl),INTENT(IN),OPTIONAL	:: Radius
real(kind=dbl)				:: res(2)

real(kind=dbl),parameter		:: eps = 1.E-12

! input point must be on unit sphere
ierr = 0
if ( dabs(1.D0-sum(xyz**2)).gt.eps) then
  res = (/ 0.D0, 0.D0 /)
  ierr = 1
  return
end if

! projection
res = (/ xyz(1)/(1.D0+xyz(3)), xyz(2)/(1.D0+xyz(3)) /)

! scale if necessary
if (present(Radius)) res = Radius * res
 
end function StereoGraphicForwardDouble

!--------------------------------------------------------------------------
!
! FUNCTION: StereoGraphicInverseSingle
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Inverse stereographic projection (from plane to upper unit hemisphere) 
!
!> @param xy 2D coordinates to be considered (single precision)  
!> @param ierr error flag 0 = OK, 1 = input cordinates outside projection circle
! 
!> @date 7/12/13    MDG 1.0 original
!--------------------------------------------------------------------------
recursive function StereoGraphicInverseSingle(xy, ierr, Radius) result(res)

use local

IMPLICIT NONE

real(kind=sgl),INTENT(IN)		:: xy(2)
integer(kind=irg),INTENT(INOUT)	:: ierr
real(kind=sgl),INTENT(IN)		:: Radius
real(kind=sgl)				:: res(3)

real(kind=sgl)				:: q, qq

ierr = 0
if (maxval(abs(xy)).eq.0.0) then
  res = (/ 0.0, 0.0, 1.0 /)
else
  qq = sum(xy**2)
  if (qq.gt.Radius**2) then
    res = (/ 0.0, 0.0, 0.0 /)
    ierr = 1
  else
    q = 1.0/(Radius**2 + qq)
    res = (/ 2.0*xy(1), 2.0*xy(2), 1.0-qq /) 
    res = res * q
  end if
end if

end function StereoGraphicInverseSingle

!--------------------------------------------------------------------------
!
! FUNCTION: StereoGraphicInverseDouble
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Inverse stereographic projection (from plane to upper unit hemisphere) 
!
!> @param xy 2D coordinates to be considered (double precision)  
!> @param ierr error flag 0 = OK, 1 = input cordinates outside projection circle
! 
!> @date 7/12/13    MDG 1.0 original
!--------------------------------------------------------------------------
recursive function StereoGraphicInverseDouble(xy, ierr, Radius) result(res)

use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN)		:: xy(2)
integer(kind=irg),INTENT(INOUT)	:: ierr
real(kind=dbl),INTENT(IN)		:: Radius
real(kind=dbl)				:: res(3)

real(kind=dbl)				:: q, qq

ierr = 0
if (maxval(dabs(xy)).eq.0.D0) then
  res = (/ 0.D0, 0.D0, 1.D0 /)
else
  qq = sum(xy**2)
  if (qq.gt.Radius**2) then
    res = (/ 0.D0, 0.D0, 0.D0 /)
    ierr = 1
  else
    q = 1.D0/(Radius**2 + qq)
    res = (/ 2.D0*xy(1), 2.D0*xy(2), 1.D0-qq /) 
    res = res * q
  end if
end if
  
end function StereoGraphicInverseDouble


!---------
! the next couple of routines implement the standard forward and inverse Lambert projection
!---------


!--------------------------------------------------------------------------
!
! FUNCTION: LambertForwardSingle
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Forward Lambert projection (from sphere South pole to plane) 
!
!> @param xyz 3D coordinates to be considered (single precision)  
!> @param ierr error flag 0 = OK, 1 = coordinates not on sphere, 2 = degenerate North pole projection
! 
!> @date 8/01/13    MDG 1.0 original
!--------------------------------------------------------------------------
recursive function LambertForwardSingle(xyz, ierr, Radius) result(res)

use local

IMPLICIT NONE

real(kind=sgl),INTENT(IN)		:: xyz(3)
integer(kind=irg),INTENT(INOUT)	:: ierr
real(kind=sgl),INTENT(IN),OPTIONAL	:: Radius
real(kind=sgl)				:: res(2), q

real(kind=sgl),parameter		:: eps = 1.E-7

! input point must be on sphere with radius Radius
ierr = 0
if ( abs(Radius**2-sum(xyz**2)).gt.eps) then
  res = (/ 0.0, 0.0 /)
  ierr = 1
  return
end if

! if the point is the North pole, then we have a degenerate projection
! onto the entire circle of radius 2*Radius; we'll signal that case
! with ierr=2
if (xyz(3).eq.Radius) then
  res = (/ 0.0, 0.0 /)
  ierr = 2
  return
end if

! otherwise we apply the forward projection
q = sqrt(2.0*Radius/(Radius - xyz(3)))
res = (/ q * xyz(1), q * xyz(2) /)
 
end function LambertForwardSingle

!--------------------------------------------------------------------------
!
! FUNCTION: LambertForwardDouble
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Forward Lambert projection (from sphere South pole to plane) 
!
!> @param xyz 3D coordinates to be considered (single precision)  
!> @param ierr error flag 0 = OK, 1 = coordinates not on sphere, 2 = degenerate North pole projection
! 
!> @date 8/01/13    MDG 1.0 original
!--------------------------------------------------------------------------
recursive function LambertForwardDouble(xyz, ierr, Radius) result(res)

use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN)		:: xyz(3)
integer(kind=irg),INTENT(INOUT)	:: ierr
real(kind=dbl),INTENT(IN),OPTIONAL	:: Radius
real(kind=dbl)				:: res(2), q

real(kind=dbl),parameter		:: eps = 1.E-12

! input point must be on sphere with radius Radius
ierr = 0
if ( dabs(Radius**2-sum(xyz**2)).gt.eps) then
  res = (/ 0.D0, 0.D0 /)
  ierr = 1
  return
end if

! if the point is the North pole, then we have a degenerate projection
! onto the entire circle of radius 2*Radius; we'll signal that case
! with ierr=2
if (xyz(3).eq.Radius) then
  res = (/ 0.D0, 0.D0 /)
  ierr = 2
  return
end if

! otherwise we'll apply the forward projection
q = dsqrt(2.D0*Radius/(Radius - xyz(3)))
res = (/ q * xyz(1), q * xyz(2) /)
 
end function LambertForwardDouble



!--------------------------------------------------------------------------
!
! FUNCTION: LambertInverseSingle
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Inverse Lambert projection (from plane to sphere tangent at South pole) 
!
!> @param xy 2D coordinates to be considered (single precision)  
!> @param ierr error flag 0 = OK, 1 = coordinates outside circle of radius 2R
! 
!> @date 8/01/13    MDG 1.0 original
!--------------------------------------------------------------------------
recursive function LambertInverseSingle(xy, ierr, Radius) result(res)

use local

IMPLICIT NONE

real(kind=sgl),INTENT(IN)		:: xy(2)
integer(kind=irg),INTENT(INOUT)	:: ierr
real(kind=sgl),INTENT(IN),OPTIONAL	:: Radius
real(kind=sgl)				:: res(3), q, tr

! input point must be inside circle with radius 2*Radius
ierr = 0
q = sum(xy**2)
tr = 4.0*Radius**2
if (q.gt.tr) then
  res = (/ 0.0, 0.0, 0.0 /)
  ierr = 1
  return
end if

! projection (automatically takes care of projection from degenerate outer circle)
tr = sqrt(1.0 - q/tr)
res = (/ tr * xy(1), tr * xy(2), -Radius + q/(2.0*Radius) /)
 
end function LambertInverseSingle


!--------------------------------------------------------------------------
!
! FUNCTION: LambertInverseDouble
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Inverse Lambert projection (from plane to sphere tangent at South pole) 
!
!> @param xy 2D coordinates to be considered (single precision)  
!> @param ierr error flag 0 = OK, 1 = coordinates outside circle of radius 2R
! 
!> @date 8/01/13    MDG 1.0 original
!--------------------------------------------------------------------------
recursive function LambertInverseDouble(xy, ierr, Radius) result(res)

use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN)		:: xy(2)
integer(kind=irg),INTENT(INOUT)	:: ierr
real(kind=dbl),INTENT(IN),OPTIONAL	:: Radius
real(kind=dbl)				:: res(3), q, tr

! input point must be inside circle with radius 2*Radius
ierr = 0
q = sum(xy**2)
tr = 4.D0*Radius**2
if (q.gt.tr) then
  res = (/ 0.D0, 0.D0, 0.D0 /)
  ierr = 1
  return
end if

! projection (automatically takes care of projection from degenerate outer circle)
tr = dsqrt(1.D0 - q/tr)
res = (/ tr * xy(1), tr * xy(2), -Radius + q/(2.D0*Radius) /)
 
end function LambertInverseDouble



!--------------------------------------------------------------------------
!
! FUNCTION: Apply2DLaueSymmetry
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Apply the Laue group symmetry to a pair of coordinates 
!
!> @detail This is not a complete symmetry implementation, but rather one
!> that is restricted to the integer coordinate equivalents; for instance,
!> in cubic m-3m symmetry, we should in principle have 24 equivalent points
!> for the Northern hemisphere, but the three-fold axes will be ignored.
!> This means that this is just a placeholder routine until a full symmetry
!> handling can be implemented.
!
!> @param ipx x-coordinate
!> @param ipy y-coordinate
!> @param isym Laue group number
!> @param iequiv array with equivalent coordinates
!> @param nequiv number of equivalent coordinates
!
!> @todo implement full symmetry handling; this may need some more verification tests
! 
!> @date  10/13/98 MDG 1.0 original
!> @date   7/04/01 MDG 2.0 f90
!> @date  12/05/12 MDG 2.1 adapted as a subroutine
!> @date  09/20/13 MDG 3.0 re-adapted to generate 2D equivalent point list
!--------------------------------------------------------------------------
subroutine Apply2DLaueSymmetry(ipx,ipy,isym,iequiv,nequiv)

use local
use io
use diffraction
use multibeams
use dynamical

IMPLICIT NONE

integer(kind=irg),INTENT(IN)	:: ipx
integer(kind=irg),INTENT(IN)	:: ipy
integer(kind=irg),INTENT(IN)	:: isym
integer(kind=irg),INTENT(OUT)	:: iequiv(2,12)
integer(kind=irg),INTENT(OUT)	:: nequiv
 
integer(kind=irg)	:: fly
real(kind=sgl)		:: my

iequiv = 0

! we should always include the input point in the output
iequiv(1:2,1) = (/ ipx, ipy /)
nequiv = 1

if (maxval(abs( (/ ipx, ipy /) )).ne.0.0) then
 select case (isym)
  case (1)
! do nothing (triclinic symmetry)
  case (2)
	iequiv(1:2,2) = (/ ipx, -ipy /)
	nequiv = 2
  case (4)
	iequiv(1:2,2) = (/ -ipy, ipx /)
	iequiv(1:2,3) = (/ ipy, -ipx /)
	iequiv(1:2,4) = (/ -ipx, -ipy /)
	nequiv = 4
  case (3,10)
	iequiv(1:2,2) = (/ -ipx, ipy /)
	iequiv(1:2,3) = (/ ipx, -ipy /)
	iequiv(1:2,4) = (/ -ipx, -ipy /)
	nequiv = 4
  case (5,11)
	iequiv(1:2,2) = (/ -ipx, ipy /)
	iequiv(1:2,3) = (/ ipx, -ipy /)
	iequiv(1:2,4) = (/ -ipx, -ipy /)
	iequiv(1:2,5) = (/ ipy, ipx /)
	iequiv(1:2,6) = (/ -ipy, ipx /)
	iequiv(1:2,7) = (/ ipy, -ipx /)
	iequiv(1:2,8) = (/ -ipy, -ipx /)
	nequiv = 8
! from here on, the coordinate grid is assumed to be hexagonal !!!
  case (6,7)  ! the transformations were computed and verified using Mathematica on 8/28/12
    	fly = floor( float(ipy)/2.0)
    	iequiv(1:2,2) = (/ ipy-2*fly -floor( float(ipx+3*ipy-3*fly)/2.0 ), ipx-ipy+fly /)
	iequiv(1:2,3) = (/ ipy-fly-floor( float(ipx-fly)/2.0 ), -ipx-fly /)
	nequiv = 3
  case (13)  ! the transformations were computed and verified using Mathematica on 9/14/12
  ! mirror is parallel to the a_1 axis
    	my = float(mod(ipy,2))/2.0
   	fly = floor( float(ipy)/2.0)
	iequiv(1:2,2) = (/ ipy-floor(float(ipx+fly)/2.0),int(ipx+ipy/2-my) /)
	iequiv(1:2,3) = (/ ipy-2*fly -floor( float(ipx+3*ipy-3*fly)/2.0), ipx-ipy+fly /)
	iequiv(1:2,4) = (/ ipx, -ipy /)
	iequiv(1:2,5) = (/ ipy-fly-floor( float(ipx-fly)/2.0), -ipx-fly /)
	iequiv(1:2,6) = (/ -fly-floor(float(ipx+ipy-fly)/2.0), int(-ipx+ipy/2+my) /)
	nequiv = 6
   case (8)  ! the transformations were computed and verified using Mathematica on 9/10/12
    	my = float(mod(ipy,2))
    	fly = floor( float(ipy)/2.0)
	iequiv(1:2,2) = (/ ipx-fly-floor( float(ipx+2*ipy-3*fly)/2.0), ipx+fly /)
	iequiv(1:2,3) = (/ ipy-2*fly -floor( float(ipx+3*ipy-3*fly)/2.0), ipx-ipy+fly /)
	iequiv(1:2,4) = (/ int(-ipx+my), -ipy /)
	iequiv(1:2,5) = (/ ipy-fly-floor( float(ipx-fly)/2.0), -ipx-fly /)
	iequiv(1:2,6) = (/ ipx+ipy-floor( float(ipx+ipy-fly)/2.0), int(-ipx+ipy/2+my/2) /)
	nequiv = 6
  case (9)  ! the transformations were computed and verified using Mathematica on 8/28/12
    	my = float(mod(ipy,2))/2.0
    	fly = floor( float(ipy)/2.0)
	iequiv(1:2,2) = (/ fly + floor( float(ipx+ipy -fly)/2.0 ), int(ipx-ipy/2-my) /)
	iequiv(1:2,3) = (/ ipx-fly-floor( float(ipx+2*ipy-3*fly)/2.0), ipx+fly /)
	iequiv(1:2,4) = (/ ipx,-ipy /)
	iequiv(1:2,5) = (/ ipy-2*fly -floor( float(ipx+3*ipy-3*fly)/2.0), ipx-ipy+fly /)
	iequiv(1:2,6) = (/ ipx-floor( float(ipx+2*ipy-fly)/2.0), -ipx-fly /)
	iequiv(1:2,7) = (/ int(-ipx+2*my), -ipy /)
	iequiv(1:2,8) = (/ -fly-floor( float(ipx+ipy-fly)/2.0), int(-ipx+ipy/2+my) /)
	iequiv(1:2,9) = (/ ipy-fly-floor( float(ipx-fly)/2.0), -ipx-fly /)
	iequiv(1:2,10) = (/ int(-ipx+2*my), ipy /)
	iequiv(1:2,11) = (/ ipx+ipy-floor( float(ipx+ipy-fly)/2.0), int((-2*ipx+ipy+2*my)/2) /)
	iequiv(1:2,12) = (/ ipy-2*fly -floor( float(ipx-3*fly)/2.0), ipx+fly /)
	nequiv = 12
! finally, take the special case of the -3m Laue group with mirror at 30 degrees from a_1.... 
 case (12)  ! the transformations were computed and verified using Mathematica on 9/6/12
    	my = float(mod(ipy,2))
    	fly = floor( float(ipy)/2.0)
	iequiv(1:2,2) = (/ ipx+fly-floor( float(ipx-ipy+fly)/2.0), ipx-ipy+fly /)
	iequiv(1:2,3) = (/ ipy-2*fly -floor( float(ipx+3*ipy-3*fly)/2.0), ipx-ipy+fly /)
	iequiv(1:2,4) = (/ ipx-floor( float(ipx+2*ipy-fly)/2.0), -ipx-fly /)
	iequiv(1:2,5) = (/ ipy-fly-floor( float(ipx-fly)/2.0), -ipx-fly /)
	iequiv(1:2,6) = (/ int(-ipx+my), ipy /)
	nequiv = 6
 end select
end if
    
end subroutine Apply2DLaueSymmetry

!--------------------------------------------------------------------------
!
! FUNCTION: Apply2DPGSymmetry
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Apply the 2D point group symmetry to a pair of coordinates 
!
!> @details This routine returns a set of unique equivalent coordinates for the
!> input point.  Note that for special points, the number nequiv must be reduced
!> from its point group order value.
!
!> @param ipx x-coordinate
!> @param ipy y-coordinate
!> @param isym 2D point group number
!> @param iequiv array with equivalent coordinates
!> @param nequiv number of equivalent coordinates
!! 
!> @date  10/02/13 MDG 1.0 first version
!--------------------------------------------------------------------------
subroutine Apply2DPGSymmetry(ipx,ipy,isym,iequiv,nequiv)

use local
use io
use symmetryvars
use symmetry

IMPLICIT NONE

integer(kind=irg),INTENT(IN)	:: ipx
integer(kind=irg),INTENT(IN)	:: ipy
integer(kind=irg),INTENT(IN)	:: isym
integer(kind=irg),INTENT(OUT)	:: iequiv(2,12)
integer(kind=irg),INTENT(OUT)	:: nequiv
 
integer(kind=irg)		:: i, j, pequiv(2,12), mequiv
real(kind=sgl),parameter	:: eps = 1.0E-6 
real(kind=sgl)			:: diff
logical				:: newp

! make sure that the symmetry matrices have been predefined; if not, then
! compute them first
if (TDPG%SYM_pgnum.ne.isym) call Generate2DSymmetry(isym)

! set the order;  note that this may need to reduced for special points
mequiv = TDPG%SYM_MATnum
iequiv = 0
pequiv = 0

! compute the transformed coordinates
do i=1,mequiv
  pequiv(1:2,i) = matmul( (/ ipx, ipy /), TDPG%SYM_direc(i,1:2,1:2) )
end do

! the first point is always unique, so simply copy it
nequiv = 1
iequiv(1:2,1) = pequiv(1:2,1)

! next, identify double entries and remove them from the list
do i=2,mequiv
  newp = .TRUE.
  do j=1,nequiv
   diff= sum(abs(pequiv(1:2,i)-iequiv(1:2,j)))
   if (diff.le.eps) then
     newp  = .FALSE.
   endif
  end do

! yes, it is a new point
  if (newp) then
   nequiv=nequiv+1
   iequiv(1:2,nequiv) = pequiv(1:2,i)
  end if
end do

! that's it.
end subroutine Apply2DPGSymmetry


end module Lambert
