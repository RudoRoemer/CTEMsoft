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
!> @date 5/29/14   MDG 1.0 original
!--------------------------------------------------------------------------
module so3

use local
use constants
use Lambert
use rotations

IMPLICIT NONE

! type definition for linked list of Rodrigues points
type FZpointd
        real(kind=dbl)        :: rod(3)        ! Rodrigues point
        type(FZpointd),pointer	:: next	         ! link to next point
end type FZpointd

type(FZpointd),pointer         :: FZlist, FZtmp        ! pointers
integer(kind=irg)              :: FZcnt                ! counts number of entries in linked list

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
!> @note calling program MUST have initialized the Lambert parameters first!!!
!
!> @date 5/12/14   MDG 1.0 original
!--------------------------------------------------------------------------
function insideCyclicFZ(rod,order) result(res)

use local
use Lambert

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

use local
use Lambert

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

use local
use Lambert

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
!
!> @note calling program MUST have initialized the Lambert parameters first!!!
!
!> @date 5/12/14   MDG 1.0 original
!--------------------------------------------------------------------------
subroutine SampleRFZ(nsteps,pgnum)

use local
use Lambert
use rotations

IMPLICIT NONE

integer(kind=irg), INTENT(IN)        :: nsteps
integer(kind=irg), INTENT(IN)        :: pgnum

real(kind=dbl)                       :: x, y, z, s, rod(3), delta, rval
type(FZpointd), pointer	       :: FZtmp2
integer(kind=irg)                    :: FZtype, FZorder
logical	                             :: insideFZ

! initialize all special constants
call InitLambertParameters

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
