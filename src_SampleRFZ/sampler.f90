! ###################################################################
! Copyright (c) 2014, Marc De Graef/Carnegie Mellon University
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

! ###################################################################
! Copyright (c) 2014, Marc De Graef/Carnegie Mellon University
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
! so3 module
!--------------------------------------------------------------------------
!
! MODULE: so3
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief
!> routines needed for uniform sampling of so3 via the cubochoric parameterization
!
!> @details This standalone module is used as follows:
!>
!> call SampleRFZ(pgnum, nsteps)
!>
!> where pgnum is the point group number (International Tables) and nsteps is the
!> required number of sampling steps along the cubic semi-edge.  The SampleRFZ
!> routine creates a linked list named FZlist which has FZcnt entries in it.  The
!> individual entries are the Rodrigues vector components and a link to the next
!> entry.  So, going through the list afterwards can be done as follows:
!>
!> FZtmp => FZlist                        ! point to the top of the list
!> do i = 1, FZcnt                        ! loop over all entries
!>   do something with FZtmp%rod          ! anything, really ...
!>   FZtmp => FZtmp%next                  ! point to the next entry
!> end do
!>
!
!> @date 5/12/14 MDG 1.0 original, extracted from 5/12/14 version of CTEMsoft2013 main library
!--------------------------------------------------------------------------

module so3

! The module should be processor independent.  This can
! be accomplished by the use of the "kind" parameters.

! Define the "kind" parameters for single and double precision reals,
  integer,parameter        		:: sgl = SELECTED_REAL_KIND(p=6,r=37)
  integer,parameter        		:: dbl = SELECTED_REAL_KIND(p=13,r=200)

! Define the "kind" parameters for short and regular integers,
  integer,parameter        		:: ish = SELECTED_INT_KIND(3)
  integer,parameter        		:: irg = SELECTED_INT_KIND(9)

! these are a bunch of constants used for the mapping; they are all in double precision
type LambertParameters
	real(kind=dbl)		:: Pi    	!  pi
	real(kind=dbl)		:: iPi    	!  1/pi
	real(kind=dbl)		:: sPi 		!  sqrt(pi)
	real(kind=dbl)		:: srt    	!  sqrt(3)/2
	real(kind=dbl)		:: pref		!  sqrt(6/pi)
! the following constants are used for the cube to quaternion hemisphere mapping
	real(kind=dbl)		:: a		! pi^(5/6)/6^(1/6)
	real(kind=dbl)		:: ap		! pi^(2/3)
	real(kind=dbl)		:: sc		! a/ap
	real(kind=dbl)		:: beta		! pi^(5/6)/6^(1/6)/2
	real(kind=dbl)		:: R1		! (3pi/4)^(1/3)
	real(kind=dbl)		:: r2		! sqrt(2)
	real(kind=dbl)		:: r22		! 1/sqrt(2)
	real(kind=dbl)		:: pi12		! pi/12
	real(kind=dbl)		:: pi8		! pi/8
	real(kind=dbl)		:: prek		! R1 2^(1/4)/beta
	real(kind=dbl)		:: r24		! sqrt(24)
        real(kind=dbl)	        :: rvmax2      ! square of max rodrigues vector length
	real(kind=dbl)         :: BP(6)       ! used for Fundamental Zone determination
end type LambertParameters

type(LambertParameters)	:: LPs         ! this variable needs to be initialized using InitLambertParameters

! type definition for linked list of Rodrigues points
type FZpointd
        real(kind=dbl)        :: rod(3)        ! Rodrigues point
        type(FZpointd),pointer	:: next	         ! link to next point
end type FZpointd


! The following two arrays are used to determine the FZtype (FZtarray) and primary rotation axis order (FZoarray)
! for each of the 32 crystallographic point group symmetries (in the order of the International Tables)
!
! 1 (C1), -1 (Ci), [triclinic]
! 2 (C2), m (Cs), 2/m (C2h), [monoclinic]
! 222 (D2), mm2 (C2v), mmm (D2h), [orthorhombic]
! 4 (C4), -4 (S4), 4/m (C4h), 422 (D4), 4mm (C4v), -42m (D2d), 4/mmm (D4h), [tetragonal]
! 3 (C3), -3 (C3i), 32 (D3), 3m (C3v), -3m (D3d), [trigonal]
! 6 (C6), -6 (C3h), 6/m (C6h), 622 (D6), 6mm (C6v), -6m2 (D3h), 6/mmm (D6h), [hexagonal]
! 23 (T), m3 (Th), 432 (O), -43m (Td), m-3m (Oh) [cubic]
!
! FZtype
! 0        no symmetry at all
! 1        cyclic symmetry
! 2        dihedral symmetry
! 3        tetrahedral symmetry
! 4        octahedral symmetry
!
integer(kind=irg),parameter	:: FZtarray(32) = (/ 0,0,1,1,1,2,2,2,1,1,1,2,2,2,2,1,1,2,2,2,1,1,1,2,2,2,2,3,3,4,3,4 /)
integer(kind=irg),parameter	:: FZoarray(32) = (/ 0,0,2,2,2,2,2,2,4,4,4,4,4,4,4,3,3,3,3,3,6,6,6,6,6,6,6,0,0,0,0,0 /)

!-------------------------------
! function/subroutine interfaces
!-------------------------------
! note that only one subroutines is public in this version of the module

! sampler routine
public :: SampleRFZ

! initialize the Lambert projection constants
private :: InitLambertParameters

! convert cubochoric to Rodrigues
private :: cu2ro

! convert homochoric to Rodrigues
private :: ho2ro

! mappings from the 3D cubic grid to the 3D spherical grid
private :: LambertCubeToBall

! auxiliary private functions for the cube to sphere mappings
private :: GetPyramid

! logical functions to determine if point is inside specific FZ
private :: insideCyclicFZ, insideDihedralFZ, insideCubicFZ


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

IMPLICIT NONE

real(kind=dbl)	:: dpi

dpi = 4.D0*datan(1.D0)

! do not change any of the following constants !
LPs%Pi = dpi
LPs%iPi = 1.D0/dpi
LPs%sPi = dsqrt(dpi)
LPs%srt = dsqrt(3.D0)/2.D0
LPs%pref = dsqrt(6.D0*LPs%iPi)

LPs%a = dpi**(5.D0/6.D0)/6.D0**(1.D0/6.D0)
LPs%ap = dpi**(2.D0/3.D0)
LPs%sc = LPs%a/LPs%ap
LPs%beta = 0.5D0*LPs%a
LPs%R1 = (3.D0*dpi/4.D0)**(1.D0/3.D0)
LPs%r2 = dsqrt(2.D0)
LPs%r22 = dsqrt(2.D0)/2.D0
LPs%pi12 = dpi/12.D0
LPs%pi8 = dtan(dpi/8.D0)
LPs%prek =  LPs%R1*2.D0**(0.25D0)/LPs%beta
LPs%r24 = dsqrt(24.D0)

! we need to make sure that we do not consider Rodrigues vectors of infinite length,
! which would correspond to 180 degree rotations; we'll put the max length so that
! we can deal with 179.999 degree rotations, which should be good enough for most
! practical cases (this may be changed if deemed necessary)
LPs%rvmax2 = (dtan(179.999D0*dpi/2.D0/180.D0))**2  ! square of max rodrigues vector length

! truncation values needed for the Cyclic fundamental zones of order 2, 3, 4, and 6
LPs%BP(1:6) = (/ 0.D0, dtan(dpi/4.D0), dtan(dpi/6.D0), LPs%pi8, 0.D0, dtan(LPs%pi12) /)

end subroutine InitLambertParameters


!-------------------------------------------------------------
!-------------------------------------------------------------
!
! Check the paper by Rosca, Morawiec, and De Graef
! for details on this cube-to-ball mapping, in particular
! the labeling of the pyramids that make up the cube.
!
! Citation:  "A new method of constructing a grid in 
! the space of 3D rotations and its applications to
! texture analysis," D. Rosca, A. Morawiec, and M. De Graef,
! submitted to Modeling and Simulations in Materials Science
! and Engineering (April 2014).
!
!-------------------------------------------------------------
!-------------------------------------------------------------


!--------------------------------------------------------------------------
!
! FUNCTION: GetPyramid
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief determine to which pyramid a point in a cubic grid belongs
!
!> @details check the first figure in the cube-to-ball paper for the pyramid labels
!
!> @param xyz 3D coordinates to be considered (double precision)
!
!> @date 11/21/12    MDG 1.0 original
!--------------------------------------------------------------------------
recursive function GetPyramid(xyz) result(res)

IMPLICIT NONE

real(kind=dbl),INTENT(IN)	:: xyz(3)
integer(kind=irg)		:: res

if ((dabs(xyz(1)).le.xyz(3)).and.(dabs(xyz(2)).le.xyz(3))) then
  res = 1				! pyramid 1
  return
end if

if ((dabs(xyz(1)).le.-xyz(3)).and.(dabs(xyz(2)).le.-xyz(3))) then
  res = 2				! pyramid 2
  return
end if

if ((dabs(xyz(3)).le.xyz(1)).and.(dabs(xyz(2)).le.xyz(1))) then
  res = 3				! pyramid 3
  return
end if

if ((dabs(xyz(3)).le.-xyz(1)).and.(dabs(xyz(2)).le.-xyz(1))) then
  res = 4				! pyramid 4
  return
end if

if ((dabs(xyz(1)).le.xyz(2)).and.(dabs(xyz(3)).le.xyz(2))) then
  res = 5				! pyramid 5
  return
end if

! if ((dabs(xyz(1)).le.-xyz(2)).and.(dabs(xyz(3)).le.-xyz(2))) then
  res = 6				! pyramid 6

end function GetPyramid


!--------------------------------------------------------------------------
!
! FUNCTION: LambertCubeToBall
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
recursive function LambertCubeToBall(lxyz,ierr) result(res)

IMPLICIT NONE

real(kind=dbl),INTENT(IN)	:: lxyz(3)
integer(kind=irg),INTENT(INOUT):: ierr
real(kind=dbl)			:: res(3)

real(kind=dbl)			:: XYZ(3), sXYZ(3), T1, T2, c, s, q, LamXYZ(3), edge
integer(kind=irg)		:: p

ierr = 0
edge = 0.5D0 * LPs%ap
if (maxval(dabs(lxyz)).gt.edge) then
  res = (/ 0.D0, 0.D0, 0.D0 /)
  ierr = 1
  return
end if

! determine which pyramid pair the point lies in and copy coordinates in correct order (see paper)
p = GetPyramid(lxyz)
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

end function LambertCubeToBall


!--------------------------------------------------------------------------
!
! FUNCTION: ho2ro
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief convert homochoric to Rodrigues
!
!> @param h homochoric coordinates (double precision)
!
!> @date 8/12/13   MDG 1.0 original
!--------------------------------------------------------------------------
function ho2ro(h) result (res)

IMPLICIT NONE

real(kind=dbl),INTENT(IN)	:: h(3)		!< homochoric coordinates
real(kind=dbl)			:: res(3)

integer(kind=irg)		:: i
real(kind=dbl)			:: hn(3), hmag, s, hm

! fit parameters determined with Mathematica
real(kind=dbl),parameter	:: c(7) = (/ -0.5000096149170321D0, -0.02486606148871731D0, &
              				-0.004549381779362819D0, 0.0005118668366387526D0, &
              				-0.0016500827333575548D0, 0.0007593352203388718D0, &
              				-0.0002040422502566876D0 /)

! normalize h and store the magnitude
hmag = sum(h * h)
if (hmag.eq.0.D0) then
  res = (/ 0.D0, 0.D0, 0.D0 /)
else
  hm = hmag
  hn = h / dsqrt(hmag)

! convert the magnitude to the rotation angle
  s = c(1) * hmag
  do i=2,7
    hm = hm * hmag
    s = s + c(i) * hm
  end do

  s = dtan( dacos(1.D0+s) )
  res = (/ hn(1), hn(2), hn(3) /) * s
end if

end function ho2ro

!--------------------------------------------------------------------------
!
! FUNCTION: cu2ro
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
function cu2ro(c) result (res)

real(kind=dbl), intent(in) 		:: c(3)		!< input coordinates
real(kind=dbl) 				:: res(3), tmp(3), cc(3)

integer(kind=irg)			:: ierr

cc = c
! calling program must have initialized the Lambert parameters!!!!
tmp = LambertCubeToBall(cc,ierr)

! if ierr=1, then the input point does not lie inside the sampling cube.
! the calling program should make sure that this never happens, but if
! it does, we need to alert the user and abort the program right here...
if (ierr.eq.1) then
  write (*,*) 'Fatal Error: the sampling point coordinates are outside sampling cube...'
  write (*,*) c
  write (*,*) '   Sampling cube has semi edge length ',0.5D0*LPs%ap
  STOP
end if

res = ho2ro(tmp)

end function cu2ro



!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! next, we define a number of logical routines, that decide whether or not
! a point in Rodrigues representation lies inside the fundamental zone (FZ)
! for a given crystal symmetry. This follows the Morawiec@Field paper:
!
! A. Morawiec & D. P. Field (1996) Rodrigues parameterization for orientation
! and misorientation distributions, Philosophical Magazine A, 73:4, 1113-1130,
! DOI: 10.1080/01418619608243708
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------


!--------------------------------------------------------------------------
!
! FUNCTION: insideCyclicFZ
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief does rodrigues point lie inside cyclic FZ (for 2, 3, 4, and 6-fold)?
!
!> @param rod Rodrigues coordinates  (double precision)
!> @param order depending on main symmetry axis
!
!> @note calling program MUST have initialized the Lambert parameters first!!!
!
!> @date 5/12/14   MDG 1.0 original
!--------------------------------------------------------------------------
function insideCyclicFZ(rod,order) result(res)

IMPLICIT NONE

real(kind=dbl), INTENT(IN)                :: rod(3)
integer(kind=irg), INTENT(IN)             :: order

logical	                                  :: res

! check the z-component vs. tan(pi/2n)
res = dabs(rod(3)).le.LPs%BP(order)

end function insideCyclicFZ

!--------------------------------------------------------------------------
!
! FUNCTION: insideDihedralFZ
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief does Rodrigues point lie inside dihedral FZ (for 2, 3, 4, and 6-fold)?
!
!> @param rod Rodrigues coordinates (double precision)
!> @param order depending on main symmetry axis
!
!> @note calling program MUST have initialized the Lambert parameters first!!!
!
!> @todo we ignore here the fact that, among others, the 3m point group can be oriented in two ways;
!> @todo this should be fixable in the future with an additional optional argument
!
!> @date 5/12/14   MDG 1.0 original
!--------------------------------------------------------------------------
function insideDihedralFZ(rod,order) result(res)

IMPLICIT NONE

real(kind=dbl), INTENT(IN)                :: rod(3)
integer(kind=irg), INTENT(IN)             :: order

logical	                                  :: res, c1, c2
real(kind=dbl),parameter                  :: r1 = 1.0D0

! first, check the z-component vs. tan(pi/2n)  (same as insideCyclicFZ)
c1 = dabs(rod(3)).le.LPs%BP(order)
res = .FALSE.

! check the square boundary planes if c1=.TRUE.
if (c1) then
  select case (order)
    case (2)
      c2 = (dabs(rod(1)).le.r1).and.(dabs(rod(2)).le.r1)
    case (3)
      c2 =          dabs( LPs%srt*rod(1)+0.5*rod(2)).le.r1
      c2 = c2.and.( dabs( LPs%srt*rod(1)-0.5*rod(2)).le.r1 )
      c2 = c2.and.( dabs(rod(2)).le.r1 )
    case (4)
      c2 = (dabs(rod(1)).le.r1).and.(dabs(rod(2)).le.r1)
      c2 = c2.and.((LPs%r22*dabs(rod(1)+rod(2)).le.r1).and.(LPs%r22*dabs(rod(1)-rod(2)).le.r1))
    case (6)
      c2 =          dabs( 0.5*rod(1)+LPs%srt*rod(2)).le.r1
      c2 = c2.and.( dabs( LPs%srt*rod(1)+0.5*rod(2)).le.r1 )
      c2 = c2.and.( dabs( LPs%srt*rod(1)-0.5*rod(2)).le.r1 )
      c2 = c2.and.( dabs( 0.5*rod(1)-LPs%srt*rod(2)).le.r1 )
      c2 = c2.and.( dabs(rod(2)).le.r1 )
      c2 = c2.and.( dabs(rod(1)).le.r1 )
  end select
  res = c2
end if

end function insideDihedralFZ

!--------------------------------------------------------------------------
!
! FUNCTION: insideCubicFZ
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief does rodrigues point lie inside cubic FZ (octahedral or tetrahedral)?
!
!> @param rod Rodrigues coordinates  (double precision)
!> @param ot 'oct' or 'tet', depending on symmetry
!
!> @note calling program MUST have initialized the Lambert parameters first!!!
!
!> @date 5/12/14   MDG 1.0 original
!--------------------------------------------------------------------------
function insideCubicFZ(rod,ot) result(res)

IMPLICIT NONE

real(kind=dbl), INTENT(IN)                :: rod(3)
character(3), INTENT(IN)                  :: ot

logical	                                  :: res, c1, c2
real(kind=dbl)                            :: rx, ry, rz
real(kind=dbl),parameter	           :: r1 = 1.0D0

rx = rod(1)
ry = rod(2)
rz = rod(3)

res = .FALSE.

! primary cube planes (only needed for octahedral case)
if (ot.eq.'oct') then
  c1 = (maxval(dabs(rod)).le.LPs%pi8)
else
  c1 = .TRUE.
end if

! octahedral truncation planes, both for tetrahedral and octahedral point groups
c2 = (dabs(rx+ry+rz).le.r1)
c2 = c2.and.(dabs(-rx+ry+rz).le.r1)
c2 = c2.and.(dabs(rx-ry+rz).le.r1)
c2 = c2.and.(dabs(rx+ry-rz).le.r1)

! if both c1 and c2, then the point is inside
if (c1.and.c2) res = .TRUE.

end function insideCubicFZ

!--------------------------------------------------------------------------
!
! FUNCTION: SampleRFZ
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Generate a uniform sampling of a Rodriguess FZ
!
!> @note This routine fills in a linked list FZlist of Rodrigues points that
!> are inside a specific fundamental zone determined by the sample point group;
!> this list can then be further dealt with in the calling program.
!
!> @param nsteps number of steps along semi-edge in cubochoric grid
!> @param pgnum point group number to determine the appropriate Rodrigues fundamental zone
!> @param FZlist pointer to start of list
!> @param FZcnt number of entries in list
!
!> @date 5/12/14   MDG 1.0 original
!--------------------------------------------------------------------------
subroutine SampleRFZ(nsteps, pgnum, FZlist, FZcnt)

IMPLICIT NONE

integer(kind=irg), INTENT(IN)        :: nsteps
integer(kind=irg), INTENT(IN)        :: pgnum
type(FZpointd),pointer, INTENT(INOUT):: FZlist
integer(kind=irg), INTENT(OUT)       :: FZcnt

real(kind=dbl)                       :: x, y, z, s, rod(3), delta, rval
type(FZpointd), pointer	      :: FZtmp2, FZtmp
integer(kind=irg)                    :: FZtype, FZorder
logical	                             :: insideFZ

! initialize all special constants
call InitLambertParameters

! cube semi-edge length
s = 0.5D0 * LPS%ap

! step size for sampling of grid; total number of samples = (2*nsteps+1)**3
delta = s/dble(nsteps)

! set the counter to zero
FZcnt = 0

! make sure the linked lists are empty
if (associated(FZlist)) then
  FZtmp => FZlist%next
  FZtmp2 => FZlist
  do
    deallocate(FZtmp2)
    if (.not. associated(FZtmp) ) EXIT
    FZtmp2 => FZtmp
    FZtmp => FZtmp%next
  end do
  nullify(FZlist)
end if

! determine which function we should call for this point group symmetry
FZtype = FZtarray(pgnum)
FZorder = FZoarray(pgnum)

! loop over the cube of volume pi^2; note that we do not want to include
! the opposite edges/facets of the cube, to avoid double counting rotations
! with a rotation angle of 180 degrees.
x = -s
do while (x.lt.s)
  y = -s
  do while (y.lt.s)
    z = -s
    do while (z.lt.s)

! convert to Rodrigues representation
      rod = cu2ro( (/ x, y, z /) )

! is this point inside the selected Rodrigues FZ ?
! limit the divergent Rodrigues space to 179.999 degrees, which
! corresponds to a length of the rod vector of 113924.
! (actually, LPs%rvmax2 is the square of this value, to avoid having to take square roots)
! this only applies to FZtypes 0 and 1; the other FZs are always finite.
       select case (FZtype)
        case (0)
          rval = sum(rod**2)
          if (rval.lt.LPs%rvmax2) insideFZ = .TRUE.
        case (1)
          rval = sum(rod**2)
          if (rval.lt.LPs%rvmax2) insideFZ = insideCyclicFZ(rod,FZorder)
        case (2)
          insideFZ = insideDihedralFZ(rod,FZorder)
        case (3)
          insideFZ = insideCubicFZ(rod,'tet')
        case (4)
          insideFZ = insideCubicFZ(rod,'oct')
       end select

! If insideFZ=.TRUE., then add this point to the linked list FZlist and keep
! track of how many points there are on this list
       if (insideFZ) then
        if (.not.associated(FZlist)) then ! do this only for the first point
          allocate(FZlist)
          FZtmp => FZlist
          nullify(FZtmp%next)
          FZtmp%rod = rod
        else
          allocate(FZtmp%next)
          FZtmp => FZtmp%next
          nullify(FZtmp%next)
          FZtmp%rod = rod
        end if
        FZcnt = FZcnt + 1
       end if

    z = z + delta
  end do
  y = y + delta
 end do
 x = x + delta
end do

! that's it.

end subroutine SampleRFZ



end module so3



!--------------------------------------------------------------------------
!
! PROGRAM: sampler
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Example program illustrating how to use the Rodrigues fundamental zone sampler
!
!> @date 5/12/14   MDG 1.0 original
!--------------------------------------------------------------------------
 MODULE m
   IMPLICIT NONE

   ! Define interface of call-back routine.
   ABSTRACT INTERFACE
	 SUBROUTINE callbackInt (x)
	   USE, INTRINSIC :: ISO_C_BINDING
	   INTEGER(KIND=4), INTENT(IN), VALUE :: x
	 END SUBROUTINE callbackInt

	 SUBROUTINE callbackStr	(s)
	   USE, INTRINSIC :: ISO_C_BINDING
 	   CHARACTER*(*) s
!  CHARACTER(LEN=255), INTENT(IN), VALUE :: s
	 END SUBROUTINE callbackStr

	 SUBROUTINE callbackPtrArray(x, p)
	   USE, INTRINSIC :: ISO_C_BINDING
	   INTEGER(KIND=4), INTENT(IN), VALUE :: x
	   TYPE(C_PTR), INTENT(OUT) :: p
	   END SUBROUTINE callbackPtrArray
     END INTERFACE
end module m




subroutine sampler(pgnum, nsteps, cstruct, cproc)

	use m
	use so3

	USE, INTRINSIC :: ISO_C_BINDING
	IMPLICIT NONE
	TYPE(C_FUNPTR), INTENT(IN), VALUE :: cproc

  	integer :: memaddress

	integer(kind=irg)	        :: i, pgnum, nsteps, j

	type(FZpointd),pointer         :: FZlist, FZtmp        ! pointers
	integer(kind=irg)              :: FZcnt                ! counts number of entries in linked list
	!real(kind=dbl), ALLOCATABLE :: FZarray(:,:)
	real(C_DOUBLE), pointer :: out(:)
	type(C_PTR) :: outp
	integer ierr
	
	TYPE, BIND(C) :: MYFTYPE
		INTEGER(C_INT) :: I, J
		REAL(C_FLOAT) :: S
		real(C_DOUBLE) :: D
	END TYPE MYFTYPE

	TYPE (MYFTYPE) :: cstruct


	PROCEDURE(callbackPtrArray), POINTER :: proc
	CALL C_F_PROCPOINTER (cproc, proc)
	! point group number
	write (*,*) 'Point group number : ', pgnum
	!read (*,*) pgnum

	! point group number
	write (*,*) 'Number of intervals along the cube semi-edge length : ', nsteps
	!read (*,*) nsteps

	write(*,*) cstruct%I, cstruct%J, cstruct%S, cstruct%D 
	! get the linked list for the FZ for point group symmetry pgnum for 100 steps along the cubic semi-edge
	call SampleRFZ(nsteps, pgnum, FZlist, FZcnt)
	call proc(FZcnt, outp)
	CALL C_F_POINTER(outp, out, [3*FZcnt])

	memaddress = loc(out)
  print *, memaddress

 ! write (*,*) 'out memory addres ', out
	write (*,*) 'Total number of orientations inside Rodrigues Fundamental Zone ',FZcnt
	!allocate(FZarray(3, FZcnt))
	!allocate(out(3*FZcnt), stat=ierr)
	write(*,*) 'Error Number', ierr



	! now we have the linked list so we can do anything we want with it.
	! in this test program, we create a VTK file so that we can visualize the RFZ with ParaView
	close(UNIT=22,STATUS='keep')
	open (UNIT=22,FILE='test.vtk',FORM='formatted',STATUS='unknown')

	write (22,"(A)") '# vtk DataFile Version 2.0'
	write (22,"(A)") 'Uniform sampling of Rodrigues FZ'
	write (22,"(A)") 'ASCII'
	write (22,"(A)") 'DATASET POLYDATA'
	write (22,"(A,I8,A)") 'POINTS ',FZcnt,' float'

	!write(*,*) FZarray(:,400)

	! scan through the list
	FZtmp => FZlist
	do i = 1, FZcnt
	  write (22,"(3F16.4)") FZtmp%rod
	  !write(*,*) i, FZtmp%rod
	  !write(*,*) i, FZarray(:,i)
	  !FZarray(:,i) = FZtmp%rod
	  j = 3*i-2
		  out(j) = FZtmp%rod(1)
		  out(j+1) = FZtmp%rod(2)
		  out(j+2) = FZtmp%rod(3)

	  FZtmp => FZtmp%next
	end do

	!out = reshape(outTemp, (/3*FZcnt/))

	outp = C_LOC(out)

	write (22,"(A)") ' '
	write (22,"(A,I8)") ' POINT_DATA ',FZcnt
	write (22,"(A)") ' SCALARS radii float 1'
	write (22,"(A)") ' LOOKUP_TABLE default'
	do i=1,FZcnt
	  write (22,"(F13.6)") 0.1
	end do
	close(UNIT=22,STATUS='keep')

end subroutine sampler
