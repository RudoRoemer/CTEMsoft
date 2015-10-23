!######################################################################
! PROGRAM :: main program
!
!
!	contains :: main loops (1D)
!
!
!######################################################################

program main

use local
use constants
use io
use error
use quaternions
use Lambert
use typedefs
use rotations

real*4,allocatable :: eu(:)
real*4,allocatable :: b(:)
real*4,allocatable :: qu(:)
real*4,allocatable :: om(:,:)
real*4,allocatable :: ax(:)
real*4,allocatable :: ro(:)

allocate(eu(3))
allocate(b(3))
allocate(qu(4))
allocate(om(3,3))
allocate(ax(4))
allocate(ro(4))
!0.0000000000000000  1.5707963705062866  3.1415927410125732

eu(1) = 3.1415927410125732D0      
eu(2) = 0.0D0
eu(3) = 0.0D0
print*,"Euler IN:",eu(1),eu(2),eu(3)

ax = eu2ax(eu);
print*,"ax OUT:",ax(1),ax(2),ax(3),ax(4)

ro = eu2ro(eu);
print*,"Ro OUT:",ro(1),ro(2),ro(3),ro(4)



qu = eu2qu(eu);
print*,"qu OUT:",qu(2),qu(3),qu(4),qu(1)


! eu = qu2om(om)
! print*,"Euler Out:",eu(1),eu(2),eu(3)

! qu = eu2qu(eu);
! print*,"qu OUT:",qu(2),qu(3),qu(4),qu(1)


om = eu2om(eu)
print*,"OM"
print*,om(1,1),om(1,2),om(1,3)
print*,om(2,1),om(2,2),om(2,3)
print*,om(3,1),om(3,2),om(3,3)



! om = qu2om_d(qu)
! print*,"OM"
! print*,om(1,1),om(1,2),om(1,3)
! print*,om(2,1),om(2,2),om(2,3)
! print*,om(3,1),om(3,2),om(3,3)




end program main
