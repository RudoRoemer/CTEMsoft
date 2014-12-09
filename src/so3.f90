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
! CTEMsoft2013:so3.f90
!--------------------------------------------------------------------------
!
! MODULE: so3
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief everything that has to do with sampling of rotation space SO(3)
!
!> @todo verify that this is correct for point groups with multiple settings, eg, 3m
!
!> @date 05/29/14 MDG 1.0 original
!> @date 10/02/14 MDG 2.0 removed globals + rewrite
!--------------------------------------------------------------------------
module so3

use local
use typedefs
use constants
use rotations

IMPLICIT NONE

! sampler routine
public :: SampleRFZ

! logical functions to determine if point is inside specific FZ
private :: insideCyclicFZ, insideDihedralFZ, insideCubicFZ


contains

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
!> @date 05/12/14  MDG 1.0 original
!> @date 10/02/14  MDG 2.0 rewrite
!--------------------------------------------------------------------------
recursive function insideCyclicFZ(rod,order) result(res)

IMPLICIT NONE

real(kind=dbl), INTENT(IN)                :: rod(4)
integer(kind=irg), INTENT(IN)             :: order

logical                                   :: res

res = .FALSE.

if (rod(4).ne.infty) then
! check the z-component vs. tan(pi/2n)
  res = dabs(rod(3)*rod(4)).le.LPs%BP(order)
else
  if (rod(3).eq.0.D0) res = .TRUE.
endif

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
!> @todo for now, we ignore here the fact that, among others, the 3m point group can be oriented in two ways;
!> @todo this should be fixable in the future with an additional optional argument
!
!> @date 05/12/14  MDG 1.0 original
!> @date 10/02/14  MDG 2.0 rewrite
!--------------------------------------------------------------------------
recursive function insideDihedralFZ(rod,order) result(res)

IMPLICIT NONE

real(kind=dbl), INTENT(IN)                :: rod(4)
integer(kind=irg), INTENT(IN)             :: order

logical                                   :: res, c1, c2
real(kind=dbl)                            :: r(3)
real(kind=dbl),parameter                  :: r1 = 1.0D0

r(1:3) = rod(1:3) * rod(4)

! first, check the z-component vs. tan(pi/2n)  (same as insideCyclicFZ)
c1 = dabs(r(3)).le.LPs%BP(order)
res = .FALSE.

! check the square boundary planes if c1=.TRUE.
if (c1) then
  select case (order)
    case (2)
      c2 = (dabs(r(1)).le.r1).and.(dabs(r(2)).le.r1)
    case (3)
      c2 =          dabs( LPs%srt*r(1)+0.5*r(2)).le.r1
      c2 = c2.and.( dabs( LPs%srt*r(1)-0.5*r(2)).le.r1 )
      c2 = c2.and.( dabs(r(2)).le.r1 )
    case (4)
      c2 = (dabs(r(1)).le.r1).and.(dabs(r(2)).le.r1)
      c2 = c2.and.((LPs%r22*dabs(r(1)+r(2)).le.r1).and.(LPs%r22*dabs(r(1)-r(2)).le.r1))
    case (6)
      c2 =          dabs( 0.5*r(1)+LPs%srt*r(2)).le.r1
      c2 = c2.and.( dabs( LPs%srt*r(1)+0.5*r(2)).le.r1 )
      c2 = c2.and.( dabs( LPs%srt*r(1)-0.5*r(2)).le.r1 )
      c2 = c2.and.( dabs( 0.5*r(1)-LPs%srt*r(2)).le.r1 )
      c2 = c2.and.( dabs(r(2)).le.r1 )
      c2 = c2.and.( dabs(r(1)).le.r1 )
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
!> @date 05/12/14  MDG 1.0 original
!> @date 10/02/14  MDG 2.0 rewrite
!--------------------------------------------------------------------------
recursive function insideCubicFZ(rod,ot) result(res)

IMPLICIT NONE

real(kind=dbl), INTENT(IN)                :: rod(4)
character(3), INTENT(IN)                  :: ot

logical                                   :: res, c1, c2
real(kind=dbl)                            :: rx, ry, rz
real(kind=dbl),parameter                  :: r1 = 1.0D0

rx = rod(1)*rod(4)
ry = rod(2)*rod(4)
rz = rod(3)*rod(4)

res = .FALSE.

! primary cube planes (only needed for octahedral case)
if (ot.eq.'oct') then
  c1 = (maxval(dabs( (/ rx, ry, rz /) )).le.LPs%pi8) 
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
!>
!> Here's how you would use this routine in a main program:
!>
!> use so3
!>
!> integer(kind=irg)       :: FZcnt, nsteps, pgnum
!> type(FZpointd),pointer  :: FZlist, FZtmp
!> 
!> nullify(FZlist)
!> FZcnt = 0
!> nsteps = 10
!> pgnum = 32
!> call sampleRFZ(nsteps, pgnum, FZcnt, FZlist)
!> 
!> Then you can access all the entries in the list and, for instance, convert them to Euler angles...
!>
!> FZtmp => FZlist                        ! point to the top of the list
!> do i = 1, FZcnt                        ! loop over all entries
!>   eu = ro2eu(FZtmp%rod)                ! convert to Euler angles
!>   do something with eu                 ! for instance, write eu to a file
!>   FZtmp => FZtmp%next                  ! point to the next entry
!> end do
!>
!> If you just want to look at the first 10 entries on the list and show all other orientation representations:
!>
!> type(orientationtyped):: ot
!> 
!> FZtmp => FZlist
!> do i = 1,10
!>   ot = init_orientation(FZtmp%rod,'ro')
!>   call print_orientation(ot)
!>   FZtmp => FZtmp%next
!> end do
!
!> @param nsteps number of steps along semi-edge in cubochoric grid
!> @param pgnum point group number to determine the appropriate Rodrigues fundamental zone
!> @param FZcnt (output) number of points inside fundamental zone
!> @param FZlist (output) linked list of points inside fundamental zone
!
!> @date 05/12/14  MDG 1.0 original
!> @date 10/02/14  MDG 2.0 rewrite, removed all globals, added function arguments
!--------------------------------------------------------------------------
recursive subroutine SampleRFZ(nsteps,pgnum,FZcnt,FZlist)

IMPLICIT NONE

integer(kind=irg), INTENT(IN)        :: nsteps
integer(kind=irg), INTENT(IN)        :: pgnum
integer(kind=irg),INTENT(OUT)        :: FZcnt                ! counts number of entries in linked list
type(FZpointd),pointer,INTENT(OUT)   :: FZlist               ! pointers

real(kind=dbl)                       :: x, y, z, s, rod(4), delta, rval
type(FZpointd), pointer              :: FZtmp, FZtmp2
integer(kind=irg)                    :: FZtype, FZorder
logical                              :: insideFZ

! cube semi-edge length
s = 0.5D0 * LPs%ap

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

! we always want the identity rotation to be the first one on the list
! it is automatically skipped later on...
allocate(FZlist)
FZtmp => FZlist
nullify(FZtmp%next)
FZtmp%rod = (/ 0.D0, 0.D0, 0.D0, 0.D0 /)
FZcnt = 1

! determine which function we should call for this point group symmetry
FZtype = FZtarray(pgnum)
FZorder = FZoarray(pgnum)

! loop over the cube of volume pi^2; note that we do not want to include
! the opposite edges/facets of the cube, to avoid double counting rotations
! with a rotation angle of 180 degrees.  This only affects the cyclic groups.
x = -s
do while (x.lt.s)
  y = -s
  do while (y.lt.s)
    z = -s
    do while (z.lt.s)

! convert to Rodrigues representation
      rod = cu2ro( (/ x, y, z /) )

! is this point inside the selected Rodrigues FZ ?  
! dealing with 180 rotations is needed only for 
! FZtypes 0 and 1; the other FZs are always finite.
       select case (FZtype)
        case (0)
          insideFZ = .TRUE.   ! all points are inside the FZ
        case (1)
          insideFZ = insideCyclicFZ(rod,FZorder)        ! infinity is checked inside this function
        case (2)
          if (rod(4).ne.infty) insideFZ = insideDihedralFZ(rod,FZorder)
        case (3)
          if (rod(4).ne.infty) insideFZ = insideCubicFZ(rod,'tet')
        case (4)
          if (rod(4).ne.infty) insideFZ = insideCubicFZ(rod,'oct')
       end select

! If insideFZ=.TRUE., then add this point to the linked list FZlist and keep
! track of how many points there are on this list
       if (insideFZ) then 
        allocate(FZtmp%next)
        FZtmp => FZtmp%next
        nullify(FZtmp%next)
        FZtmp%rod = rod
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
