
this is a new version of the basic library, in which everything has been simplified
and duplication of code has been reduced significantly.  The new code is set up so
that it should be much easier to deal with multi-phase materials.

As an example, here is how to deal with two crystal structures:

use local
use initializers
use crystalvars
use crystal
use files
use gvectors
use dynamical
use diffraction

IMPLICIT NONE

real(kind=sgl)		        :: voltage, dmin, k(3), dp
character(fnlen)		:: xtalname
type(unitcell)	                :: cellA, cellB
logical	                        :: verbose
type(reflisttype),pointer	:: rl
integer(kind=irg)               :: ii, skip, gg(3)

! general variable inits
verbose = .FALSE.
dmin = 0.01
voltage = 20000.0

! load and init structure A
xtalname = 'Cu.xtal'
call Initialize_Cell(xtalname, dmin, voltage, verbose)
call CopyFromCell(cellA)

! init and load structure B
xtalname = 'GaP.xtal'
call Initialize_Cell(xtalname, dmin, voltage, verbose)
call CopyFromCell(cellB)


DynFN = (/ 1.0, 1.0, 1.0 /) 
k = DynFN
call NormVec(k,'r')
k = k/sngl(mLambda)
call Initialize_ReflectionList(k,verbose)

! switching to crystal A 
call CopyToCell(cellA)
write (*,*) 'Lattice parameters '
write (*,*) cell%a, cell%b, cell%c
write (*,*) cell%alpha, cell%beta, cell%gamma
rl => cell%reflist%next
do
  if (.not.associated(rl)) EXIT
  dp = dot_product(rl%hkl,DynFN)
  write (*,*) rl%hkl, dp
  rl => rl%next
end do
!
! switching to crystal B 
call CopyToCell(cellB)
write (*,*) 'Lattice parameters '
write (*,*) cell%a, cell%b, cell%c
write (*,*) cell%alpha, cell%beta, cell%gamma
rl => cell%reflist%next
do
  if (.not.associated(rl)) EXIT
  dp = dot_product(rl%hkl,DynFN)
  write (*,*) rl%hkl, dp
  rl => rl%next
end do

etc .... 




