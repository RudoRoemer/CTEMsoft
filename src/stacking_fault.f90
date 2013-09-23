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
! CTEMsoft2013:stacking_fault.f90
!--------------------------------------------------------------------------
!
! MODULE: stacking_fault
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief Provides routines to compute the displacement vector for various stacking faults.
! 
!> @date   04/29/11 MDG 1.0 original
!> @date   06/04/13 MDG 2.0 rewrite + quaternions instead of rotations
!--------------------------------------------------------------------------
module stacking_fault

use local  

type stackingfaulttype
  real(kind=sgl)     :: lpu(3),tpu(3),lpb(3),lpbc(3),tpb(3),plane(3),sep,id,jd, &
                                lptop(3),lpbot(3),tptop(3),tpbot(3),thetan,a_if(3,3), &
                                lpr(3),tpr(3)
  real(kind=sgl),allocatable     :: zpos(:,:)
end type stackingfaulttype

type (stackingfaulttype), allocatable  :: SF(:)

contains


!--------------------------------------------------------------------------
!
! FUNCTION: point_inside_triangle
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief  determines whether or not a point lies inside or outside a triangle (in 2D)
!
!> @details based on http://www.blackpawn.com/texts/pointinpoly/default.html
! 
!> @param v0 vertex coordinates
!> @param v1 vertex coordinates
!> @param v2 vertex coordinates
! 
!> @date 1/5/99   MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   03/25/13 MDG 3.0 updated IO
!--------------------------------------------------------------------------
logical function point_inside_triangle(v0,v1,v2)

IMPLICIT NONE

real(kind=sgl),INTENT(IN) 	:: v0(2), v1(2), v2(2) 
real(kind=sgl)				:: dot00, dot01, dot02, dot11, dot12, invdenom, u, v
logical 					:: inside

dot00 = v0(1)*v0(1) + v0(2)*v0(2) ! matmul(v0,v0)
dot01 = v0(1)*v1(1) + v0(2)*v1(2) ! matmul(v0,v1)
dot02 = v0(1)*v2(1) + v0(2)*v2(2) ! matmul(v0,v2)
dot11 = v1(1)*v1(1) + v1(2)*v1(2) ! matmul(v1,v1)
dot12 = v1(1)*v2(1) + v1(2)*v2(2) ! matmul(v1,v2)

! use barycentric coordinates
invdenom = 1.0 / (dot00*dot11 - dot01*dot01)
u = (dot11*dot02 - dot01*dot12)*invdenom
v = (dot00*dot12 - dot01*dot02)*invdenom

inside = .FALSE.
if ((u.ge.0.0).and.(v.ge.0.0).and.((u+v).le.1.0)) inside=.TRUE.

point_inside_triangle = inside

end function point_inside_triangle

!--------------------------------------------------------------------------
!
! SUBROUTINE: rank_points
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief rank points for inside/outside analysis
!
!> @details based on http://www.blackpawn.com/texts/pointinpoly/default.html
! 
!> @param p1 vertex coordinates
!> @param p2 vertex coordinates
!> @param p3 vertex coordinates
!> @param p4 vertex coordinates
!> @param xx sorted x coordinates
!> @param yy sorted y coordinates
! 
!> @date 1/5/99   MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   03/25/13 MDG 3.0 updated IO
!--------------------------------------------------------------------------
subroutine rank_points(p1,p2,p3,p4,xx,yy)

IMPLICIT NONE

real(kind=sgl),INTENT(IN) 	:: p1(2), p2(2), p3(2), p4(2) 
real(kind=sgl),INTENT(OUT)	:: xx(4), yy(4) 
real(kind=sgl)				::a(6), ma 

! set the first point
xx(1) = p1(1)
yy(1) = p1(2)

! compute the areas of six permutations
a(1) = p1(1)*p2(2)-p1(2)*p2(1) + p2(1)*p3(2)-p2(2)*p3(1) + p3(1)*p4(2)-p3(2)*p4(1) + p4(1)*p1(2)-p4(2)*p1(1)
a(2) = p1(1)*p2(2)-p1(2)*p2(1) + p2(1)*p4(2)-p2(2)*p4(1) + p4(1)*p3(2)-p4(2)*p3(1) + p3(1)*p1(2)-p3(2)*p1(1)
a(3) = p1(1)*p3(2)-p1(2)*p3(1) + p3(1)*p2(2)-p3(2)*p2(1) + p2(1)*p4(2)-p2(2)*p4(1) + p4(1)*p1(2)-p4(2)*p1(1)
a(4) = p1(1)*p3(2)-p1(2)*p3(1) + p3(1)*p4(2)-p3(2)*p4(1) + p4(1)*p2(2)-p4(2)*p2(1) + p2(1)*p1(2)-p2(2)*p1(1)
a(5) = p1(1)*p4(2)-p1(2)*p4(1) + p4(1)*p2(2)-p4(2)*p2(1) + p2(1)*p3(2)-p2(2)*p3(1) + p3(1)*p1(2)-p3(2)*p1(1)
a(6) = p1(1)*p4(2)-p1(2)*p4(1) + p4(1)*p3(2)-p4(2)*p3(1) + p3(1)*p2(2)-p3(2)*p2(1) + p2(1)*p1(2)-p2(2)*p1(1)
a = abs(a)
ma = maxval(a)

if (a(1).eq.ma) then 
  xx(2:4) = (/ p2(1),p3(1),p4(1) /)
  yy(2:4) = (/ p2(2),p3(2),p4(2) /)
  return
end if
if (a(2).eq.ma) then 
  xx(2:4) = (/ p2(1),p4(1),p3(1) /)
  yy(2:4) = (/ p2(2),p4(2),p3(2) /)
  return
end if
if (a(3).eq.ma) then 
  xx(2:4) = (/ p3(1),p2(1),p4(1) /)
  yy(2:4) = (/ p3(2),p2(2),p4(2) /)
  return
end if
if (a(4).eq.ma) then 
  xx(2:4) = (/ p3(1),p4(1),p2(1) /)
  yy(2:4) = (/ p3(2),p4(2),p2(2) /)
  return
end if
if (a(5).eq.ma) then 
  xx(2:4) = (/ p4(1),p2(1),p3(1) /)
  yy(2:4) = (/ p4(2),p2(2),p3(2) /)
  return
end if
if (a(6).eq.ma) then 
  xx(2:4) = (/ p4(1),p3(1),p2(1) /)
  yy(2:4) = (/ p4(2),p3(2),p2(2) /)
end if

end subroutine rank_points

!--------------------------------------------------------------------------
!
! FUNCTION: point_inside_polygon
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief  determines whether or not a point lies inside or outside a polygon (in 2D)
!
!> @details based on http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html (Fortran version)
! 
!> @param px point coordinate x
!> @param py point coordinate y
!> @param xx vertex coordinates
!> @param yy vertex coordinates
! 
!> @date 1/5/99   MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   03/25/13 MDG 3.0 updated IO
!--------------------------------------------------------------------------
logical function point_inside_polygon(px,py,xx,yy)

IMPLICIT NONE

real(kind=sgl),INTENT(IN) 	:: px, py, xx(4), yy(4)
real(kind=sgl) 				::  x(4), y(4), z
integer(kind=irg) 			:: i, j, inorout
logical 					:: mx,my,nx,ny

x = xx-px
y = yy-py
inorout = -1

do i=1,4
  j = 1+mod(i,4)
  mx = x(i).ge.0.0
  nx = x(j).ge.0.0
  my = y(i).ge.0.0
  ny = y(j).ge.0.0
  if (.not.((my.or.ny).and.(mx.or.nx)).or.(mx.and.nx)) cycle
  if (my.and.ny.AND.(mx.or.nx).and..not.(mx.and.nx)) then
    inorout = -inorout
    cycle
  else
    z  = (y(i)*x(j)-x(i)*y(j))/(x(j)-x(i))
    if (z.lt.0.0) cycle
    if (z.eq.0.0) then 
      inorout = 0
      exit
    end if
    if (z.gt.0.0) inorout=-inorout
  endif
end do

point_inside_polygon = inorout.ge.0

end function point_inside_polygon


!--------------------------------------------------------------------------
!
! SUBROUTINE: makestackingfault
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief compute parameters for a stacking fault
!
!> @details  This subroutine computes the geometrical parameters for a 
!> stacking fault.  It computes, among others, the coordinates of the 
!> centers of the partial dislocations, the intersections of each dislocation
!> line with the top and bottom foil surface, and an array that indicates, for
!> each image pixel, whether or not the corresponding integration column 
!> contains this fault; if it does not, the  value in the array is equal to 
!> -10000; if it does, then the value is equal to the point where the fault
!> plane intersects with the column, measured from the top surface.
!> In short, anything that is needed in the CalcR routine and can be computed 
!> ahead of time, is computed here.  The routine also calls the makedislocation 
!> routine to create the partials.
! 
!> @param inum
!> @param DF_L column edge length
!> @param nx
!> @param ny
!> @param DF_g
!> @param ndl
!> @param dinfo trigger for verbose output
! 
!> @todo This routine seems to have errors in it, so it will need tobe debugged entirely !
!
!> @date 1/5/99   MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   03/25/13 MDG 3.0 updated IO
!--------------------------------------------------------------------------
subroutine makestackingfault(inum,DF_L,nx,ny,DF_g,ndl,dinfo)
 
use local
use math
use constants
use files
use foilmodule
use dislocation
use crystal
use crystalvars
use symmetry
use symmetryvars

IMPLICIT NONE

real(kind=sgl),INTENT(IN)	:: DF_L
integer(kind=irg),INTENT(IN)	:: inum, nx, ny, DF_g(3), ndl, dinfo

real(kind=sgl)  			:: fpn(3),am(4,4),midpoint(3), &
                  				 lptopi(3),lpboti(3),tptopi(3),tpboti(3),gg,det,A(4), xx(4), yy(4), tmp(3), &
						 tmp2(3), cc
!real(kind=sgl)  			:: a_ir(3,3),a_id(3,3),fpn(3),SFlpuf(3),SFtpuf(3),alpha,am(4,4),midpoint(3),Fi(3), &
!                  				 lptopi(3),lpboti(3),tptopi(3),tpboti(3),gg,det,A(4),dp(3), xx(4), yy(4), tmp(3), &
!						 tmp2(3), cc
integer(kind=irg) 			:: i,j,info,ipiv,minx,maxx,miny,maxy
!logical 					:: inside1,inside2, inside

! begin by computing the angle between the projected fault normal vector 
! in the image plane

! first transform the plane normal into the direct space reference frame
call TransSpace(SF(inum)%plane,tmp,'r','d')
call NormVec(tmp,'d')
tmp2 = foil%F
call NormVec(tmp2,'d')
cc = CalcAngle(tmp2,tmp,'d')
! then take the double cross product Fx(pxF) with p the fault normal
call CalcCross(tmp,tmp2,fpn,'d','d',0)
call NormVec(fpn,'d')
call CalcCross(tmp2,fpn,tmp,'d','d',0)
gg = sqrt(tmp(1)**2+tmp(2)**2)
write (*,*) 'projected fault plane normal : ',tmp
write (*,*) 'angle to foil normal : ',cc*180.0/cPi
fpn(3) = cos(cc)
fpn(1:2) = tmp(1:2)*sqrt(1.0-fpn(3)**2)/gg
write (*,*) 'image space fault plane normal : ',fpn
! the resulting vector lies in the plane normal to the foil normal, and we
! need its components to compute the angle thetan (normalize the first two components)

! and get the angle
SF(inum)%thetan = atan2(fpn(2),fpn(1))

if (dinfo.eq.1) write (*,*) ' thetan = ',SF(inum)%thetan

write (*,*) 'scale factor  = ',DF_L

! next compute the location of the partial centers in the image space
SF(inum)%lpr(1:3) = (/ SF(inum)%id-0.5*sin(SF(inum)%thetan)*SF(inum)%sep, &
                                   SF(inum)%jd+0.5*cos(SF(inum)%thetan)*SF(inum)%sep,0.0 /) / DF_L
SF(inum)%tpr(1:3) = (/ SF(inum)%id+0.5*sin(SF(inum)%thetan)*SF(inum)%sep, &
                                   SF(inum)%jd-0.5*cos(SF(inum)%thetan)*SF(inum)%sep,0.0 /) / DF_L

if (dinfo.eq.0) write (*,*) 'lpr_i = ',SF(inum)%lpr(1:3)
if (dinfo.eq.0) write (*,*) 'tpr_i = ',SF(inum)%tpr(1:3)

! call makedislocation for each of the partials
  DL(ndl+1)%id = SF(inum)%lpr(1) 
  DL(ndl+1)%jd = SF(inum)%lpr(2) 
  DL(ndl+1)%u = SF(inum)%lpu
  DL(ndl+1)%burg = SF(inum)%lpb
  DL(ndl+1)%g = float(DF_g)
  call makedislocation(ndl+1,dinfo,DF_L)
  if (dinfo.eq.0) write (*,*) 'Leading Partial Position ',DL(ndl+1)%id,DL(ndl+1)%jd
    
  DL(ndl+2)%id = SF(inum)%tpr(1) 
  DL(ndl+2)%jd = SF(inum)%tpr(2) 
  DL(ndl+2)%u = SF(inum)%tpu
  DL(ndl+2)%burg = SF(inum)%tpb
  DL(ndl+2)%g = float(DF_g)
  call makedislocation(ndl+2,dinfo,DF_L)
  if (dinfo.eq.0)  write (*,*) 'Trailing Partial Position ',DL(ndl+2)%id,DL(ndl+2)%jd

! copy the top and bottom dislocation intersections (computed in make_dislocation) 
! into the corresponding variables of the SF record
SF(inum)%lpbot = DL(ndl+1)%top
SF(inum)%lptop = DL(ndl+1)%bottom
SF(inum)%tpbot = DL(ndl+2)%top
SF(inum)%tptop = DL(ndl+2)%bottom

! obviously, these four points need to lie in a single plane; at this point, we check that this is indeed the case
! by computing the volume of the tetrahedron formed by these four points; if the volume is zero, then the 
! points are co-planar.  (Use LAPACK's LU-decomposition and compute the product of the diagonal elements of U)
am(1:4,1) = (/ SF(inum)%lptop(1:3),1.0 /)
am(1:4,2) = (/ SF(inum)%lpbot(1:3),1.0 /) 
am(1:4,3) = (/ SF(inum)%tptop(1:3),1.0 /) 
am(1:4,4) = (/ SF(inum)%tpbot(1:3),1.0 /) 
call sgetrf(4,4,am,4,ipiv,info)
det = abs(am(1,1)*am(2,2)*am(3,3)*am(4,4))
if (dinfo.eq.0) write (*,*) 'determinant (should be zero) = ',det

! ok, next we need to figure out which image pixels lie on the projection of the stacking fault plane.
! these lines need to be replaced with lines that transform the coordinates to the tilted foil case !!!
lptopi = SF(inum)%lptop
lpboti = SF(inum)%lpbot
tptopi = SF(inum)%tptop
tpboti = SF(inum)%tpbot
write (*,*) 'SF parameters :'
write (*,*) lptopi,' <> ',lpboti
write (*,*) tptopi,' <> ',tpboti

! define the array that will contain the zpos values
allocate(SF(inum)%zpos(nx,ny))
SF(inum)%zpos = -10000.0    ! set all points to be outside the projected SF

! first determine the smaller box
minx = nint(min( lptopi(1),lpboti(1),tptopi(1),tpboti(1) )) -2
maxx = nint(max( lptopi(1),lpboti(1),tptopi(1),tpboti(1) )) +2
miny = nint(min( lptopi(2),lpboti(2),tptopi(2),tpboti(2) )) -2
maxy = nint(max( lptopi(2),lpboti(2),tptopi(2),tpboti(2) )) +2

! the fault edges may fall outside of the viewing frame (origin at the center !!!)
if (minx.lt.(-nx/2+1)) minx=-nx/2+1
if (maxx.gt.nx/2) maxx=nx/2
if (miny.lt.(-ny/2+1)) miny=-ny/2+1
if (maxy.gt.ny/2) maxy=ny/2

write (*,*) 'Integer fault box = ',minx,maxx,miny,maxy

! get the equation of the stacking fault plane in the image reference frame
! first the unit plane normal in image space
A(1:3) = fpn(1:3)
midpoint = 0.25*(lptopi+lpboti+tptopi+tpboti)
A(4) = sum(fpn*midpoint)
if (dinfo.eq.0) write (*,*) 'fault plane parameters : ',A, midpoint

! rank the corner points so that the polygon is convex
! call rank_points(tpboti(1:2),lpboti(1:2),lptopi(1:2),tptopi(1:2),xx,yy)
xx = (/ lptopi(1), tptopi(1), tpboti(1),lpboti(1) /)
yy = (/ lptopi(2), tptopi(2), tpboti(2),lpboti(2) /)

! for all of the points inside this box:
do i=minx,maxx
  do j=miny,maxy
    if (point_inside_polygon( float(i), float(j), xx, yy )) then ! the point lies inside the projected region, so we need the depth of the SF plane at this position
      SF(inum)%zpos(i+nx/2,j+ny/2) =   ( DF_L * ( A(4) - A(1)*float(i) - A(2)*float(j) )/A(3) )
    end if
!    if (i.eq.0) write (*,*) j, DF_L * ( A(4) - A(1)*float(i) - A(2)*float(j) )/A(3), &
!         point_inside_polygon( float(i), float(j), xx, yy )
  end do
end do
if (dinfo.eq.1) write (*,*) 'fault plane pixels determined'

! let's also make sure that the leading partial Burgers vector is translated to the 
! cartesian reference frame, so that it can be used directly by the CalcR routine
SF(inum)%lpbc = matmul(cell%dsm,SF(inum)%lpb)

! that should do it for the stacking fault...  The rest 
! takes place in the CalcRLocal routine.
end subroutine makestackingfault


!--------------------------------------------------------------------------
!
! SUBROUTINE: read_stacking_fault_data
!
!> @author Marc De Graef, Carnegie Melon University
!
!> @brief  read stacking fault namelist files
!
!> @param numsf number of stacking faults
!> @param numdisl number of dislocations
!> @param sfname name of staking fault namelist file (string array)
!> @param DF_L 
!> @param DF_npix number of x-pixels
!> @param DF_npiy number of y-pixels
!> @param DF_gf 
!> @param dinfo logical to trigger verbose output
! 
!> @date    1/5/99   MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   06/04/13 MDG 3.0 rewrite
!--------------------------------------------------------------------------
subroutine read_stacking_fault_data(numsf,numdisl,sfname,DF_L,DF_npix,DF_npiy,DF_g,dinfo)

use local
use io
use files
use dislocation

IMPLICIT NONE

integer(kind=irg),INTENT(IN)		:: numsf, DF_npix, DF_npiy, DF_g(3), dinfo
integer(kind=irg),INTENT(INOUT)	:: numdisl
character(fnlen),INTENT(IN)		:: sfname(maxdefects)
real(kind=sgl),INTENT(IN)		:: DF_L

integer(kind=irg) 			:: i,SFplane(3)
real(kind=sgl)     			:: SFi,SFj,SFsep,SFlpu(3),SFlpb(3),SFtpu(3),SFtpb(3)

namelist / SFdata / SFi, SFj, SFsep, SFplane, SFlpu, SFlpb, SFtpu, SFtpb

! allocate the necessary memory
allocate(SF(numsf))

! if the dislocation memory has not yet been allocated, do it here...
if ( .not.allocated(DL) ) then
  allocate(DL(2*numsf))
endif

! read the namelist files for all of the stacking faults
 do i=1,numsf
    mess = 'opening '//sfname(i); call Message("(/A)")
    OPEN(UNIT=dataunit,FILE=sfname(i),DELIM='APOSTROPHE')
    READ(UNIT=dataunit,NML=SFdata)
    CLOSE(UNIT=dataunit)
! transform the fault fractional coordinates to nm in the image reference frame
    SF(i)%id = SFi * 0.5 * float(DF_npix) ! * DF_L  (zooming is done later in the image reference frame)
    SF(i)%jd = SFj * 0.5 * float(DF_npiy) ! * DF_L
    SF(i)%sep = SFsep
    SF(i)%plane = SFplane
    SF(i)%lpu = SFlpu
    SF(i)%lpb = SFlpb
    SF(i)%tpu = SFtpu
    SF(i)%tpb = SFtpb
! initialize the stacking fault variables and both partial dislocations
    call makestackingfault(i,DF_L,DF_npix,DF_npiy,DF_g,numdisl,dinfo)
    numdisl = numdisl + 2
 end do
 
end subroutine read_stacking_fault_data

end module stacking_fault

