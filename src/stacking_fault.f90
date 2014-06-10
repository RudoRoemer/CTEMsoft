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
! CTEMsoft2013:stacking_fault.f90
!--------------------------------------------------------------------------
!
! MODULE: stacking_fault
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Provides routines to compute the displacement vector for various stacking faults.
!
!> @note Remember that all defects are defined with respect to the foil in the zero-tilt
!> orientation !!!!  This means that the foil reference frame then coincides with the 
!> image reference frame, apart from a potential scaling factor.  
!>
!> Remember also that the quat_mult routine for quaternion multiplication expects 
!> the subscripts for the a_xx rotations to be read from right to left, whereas the 
!> quat_rotate_vector routine expects them the other way around...  So, to combine
!> the quaternions a_fc (from crystal to foil) with a_if (from foil to image) we
!> need to do a_ic=quat_mult(a_if,a_fc); if we want to rotate the vector v from crystal
!> to image frame, then we need quat_rotate_vector(a_ci,v).  This is really important
!> and has caused some confusion/errors in the past; it is a consequence of my choice for
!> the quat_rotate_vector order of operations... 
! 
!> @date   04/29/11 MDG 1.0 original
!> @date   06/04/13 MDG 2.0 rewrite + quaternions instead of rotations
!--------------------------------------------------------------------------
module stacking_fault

use local  
use typedefs

contains


!--------------------------------------------------------------------------
!
! FUNCTION: point_inside_triangle
!
!> @author Marc De Graef, Carnegie Mellon University
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
real(kind=sgl)			:: dot00, dot01, dot02, dot11, dot12, invdenom, u, v
logical 			:: inside

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
!> @author Marc De Graef, Carnegie Mellon University
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
real(kind=sgl)			::a(6), ma 

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
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  determines whether or not a point lies inside or outside a polygon (in 2D)
!
!> @details based on http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html (Fortran version)
! 
!> @note Copyright (c) 1970-2003, Wm. Randolph Franklin
!>
!> Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
!> and associated documentation files (the "Software"), to deal in the Software without restriction, 
!> including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
!> and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, 
!> subject to the following conditions:

!> Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimers.
!> Redistributions in binary form must reproduce the above copyright notice in the documentation and/or 
!> other materials provided with the distribution.
!> The name of W. Randolph Franklin may not be used to endorse or promote products derived from this Software 
!> without specific prior written permission.
!> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED 
!> TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
!> THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
!> CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
!> DEALINGS IN THE SOFTWARE.
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
function point_inside_polygon(px,py,xx,yy) result(inorout)

IMPLICIT NONE

real(kind=sgl),INTENT(IN) 	:: px, py, xx(4), yy(4)
real(kind=sgl) 			:: x(4), y(4), z
integer(kind=irg) 		:: i, j, inorout
logical 			:: mx,my,nx,ny

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
  if (.not.(my.and.ny.AND.(mx.or.nx).and..not.(mx.and.nx))) then
    z  = (y(i)*x(j)-x(i)*y(j))/(x(j)-x(i))
    if (z.lt.0.0) cycle
    if (z.eq.0.0) then 
      inorout = 0
      exit
    end if
    if (z.gt.0.0) inorout=-inorout
  else
    inorout = -inorout
    cycle
  endif
end do

end function point_inside_polygon

!--------------------------------------------------------------------------
!
! SUBROUTINE: makestackingfault
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute parameters for a stacking fault
!
!> @details  This subroutine computes the geometrical parameters for a 
!> stacking fault.  It computes, among others, the coordinates of the 
!> centers of the partial dislocations, the intersections of each dislocation
!> line with the top and bottom foil surface, and an array that indicates, for
!> each image pixel, whether or not the corresponding integration column 
!> contains this fault; if it does not, the  value in the array is set to 
!> -10000; if it does, then the value is equal to the point where the fault
!> plane intersects with the column, measured from the top surface.
!> In short, anything that is needed in the CalcR routine and can be computed 
!> ahead of time, is computed here.  The routine also calls the makedislocation 
!> routine to create the partials.
! 
!> @param cell unit cell pointer
!> @param defects defect structure
!> @param foil foil structure
!> @param inum
!> @param DF_L column edge length
!> @param nx
!> @param ny
!> @param DF_g
!> @param ndl
!> @param dinfo trigger for verbose output
!
!> @date   11/05/13 MDG 1.0 new attempt to replace faulty original routine
!> @date   11/13/13 MDG 1.1 traced error to problem with transformations in defectmodule
!> @date   11/13/13 MDG 1.2 changed SF normal transformation for zpos array computation (to be tested)
!> @date   06/09/14 MDG 2.0 added defects and cell as arguments
!> @date   06/10/14 MDG 2.1 added foil as argument 
!--------------------------------------------------------------------------
subroutine makestackingfault(defects,cell,foil,inum,DF_L,nx,ny,DF_g,dinfo)
 
use math
use constants
use files
use foilmodule
use dislocation
use quaternions
use rotations
use crystal
use symmetry

IMPLICIT NONE

type(defecttype),INTENT(INOUT)              :: defects
type(unitcell),pointer	                     :: cell
type(foiltype),INTENT(INOUT)                :: foil
real(kind=sgl),INTENT(IN)	             :: DF_L
integer(kind=irg),INTENT(IN)	             :: inum, nx, ny, DF_g(3), dinfo

real(kind=sgl)  		:: fpn(3),am(4,4),midpoint(3), ex(3), ey(3),&
                  		 lptopi(3),lpboti(3),tptopi(3),tpboti(3),det,A(4), xx(4), yy(4), tmp(3), &
				 planenormal(3), rzero(3), unita(3)

integer(kind=irg) 		:: i,j,info,ipiv,minx,maxx,miny,maxy, ndl

ndl = defects%numdisl

! we begin by computing the geometry in the foil reference frame, which is the cartesian frame 
! for zero sample tilt;  sample tilts are applied once we known the partial dislocation geometry
call TransSpace(cell,defects%SF(inum)%plane,tmp,'r','c')
call NormVec(cell,tmp,'c')
planenormal =  tmp

call CalcCross(cell, planenormal, (/ 0.0,0.0,1.0 /),  unita, 'c', 'c', 0)
call NormVec(cell,unita,'c')
fpn = planenormal
if (dinfo.eq.1) write (*,*) ' unita should have zero third component ',unita, fpn

! the fault plane goes through the point rzero in the foil center plane
rzero = (/ defects%SF(inum)%id, defects%SF(inum)%jd, 0.0 /)

! this leads to the location of the partial dislocation centers
defects%SF(inum)%lpr = rzero + 0.5*defects%SF(inum)%sep*unita / DF_L
defects%SF(inum)%tpr = rzero - 0.5*defects%SF(inum)%sep*unita / DF_L

if (dinfo.eq.1) write (*,*) 'lpr_i = ',defects%SF(inum)%lpr(1:3)
if (dinfo.eq.1) write (*,*) 'tpr_i = ',defects%SF(inum)%tpr(1:3)

! call makedislocation for each of the partials
  defects%DL(ndl+1)%id = defects%SF(inum)%lpr(1) 
  defects%DL(ndl+1)%jd = defects%SF(inum)%lpr(2) 
  defects%DL(ndl+1)%u =  defects%SF(inum)%lpu
  defects%DL(ndl+1)%burg = defects%SF(inum)%lpb
  defects%DL(ndl+1)%g = float(DF_g)
  call makedislocation(defects,cell,foil,ndl+1,dinfo,DF_L)
  if (dinfo.eq.1) write (*,*) 'Leading Partial Position ',defects%DL(ndl+1)%id,defects%DL(ndl+1)%jd
    
  defects%DL(ndl+2)%id = defects%SF(inum)%tpr(1) 
  defects%DL(ndl+2)%jd = defects%SF(inum)%tpr(2) 
  defects%DL(ndl+2)%u = defects%SF(inum)%tpu
  defects%DL(ndl+2)%burg = defects%SF(inum)%tpb
  defects%DL(ndl+2)%g = float(DF_g)
  call makedislocation(defects,cell,foil,ndl+2,dinfo,DF_L)
  if (dinfo.eq.1)  write (*,*) 'Trailing Partial Position ',defects%DL(ndl+2)%id,defects%DL(ndl+2)%jd

! copy the top and bottom dislocation intersections (computed in make_dislocation) 
! into the corresponding variables of the SF record

defects%SF(inum)%lpbot = defects%DL(ndl+1)%bottom
defects%SF(inum)%lptop = defects%DL(ndl+1)%top
defects%SF(inum)%tpbot = defects%DL(ndl+2)%bottom
defects%SF(inum)%tptop = defects%DL(ndl+2)%top

! obviously, these four points need to lie in a single plane; at this point, we check that this is indeed the case
! by computing the volume of the tetrahedron formed by these four points; if the volume is zero, then the 
! points are co-planar.  (Use LAPACK's LU-decomposition and compute the product of the diagonal elements of U)
am(1:4,1) = (/ defects%SF(inum)%lptop(1:3),1.0 /)
am(1:4,2) = (/ defects%SF(inum)%lpbot(1:3),1.0 /) 
am(1:4,3) = (/ defects%SF(inum)%tptop(1:3),1.0 /) 
am(1:4,4) = (/ defects%SF(inum)%tpbot(1:3),1.0 /) 
call sgetrf(4,4,am,4,ipiv,info)
det = abs(am(1,1)*am(2,2)*am(3,3)*am(4,4))
if (dinfo.eq.1) write (*,*) 'determinant (should be zero) = ',det

! ok, next we need to figure out which image pixels lie on the projection of the stacking fault plane.
! We need to transform the corner points into the image reference frame !!!
lptopi = quat_rotate_vector( conjg(foil%a_fi), dble(defects%SF(inum)%lptop))
lpboti = quat_rotate_vector( conjg(foil%a_fi), dble(defects%SF(inum)%lpbot))
tptopi = quat_rotate_vector( conjg(foil%a_fi), dble(defects%SF(inum)%tptop))
tpboti = quat_rotate_vector( conjg(foil%a_fi), dble(defects%SF(inum)%tpbot))
if (dinfo.eq.1) then
  write (*,*) 'SF parameters :'
  write (*,*) lptopi,' <> ',lpboti
  write (*,*) tptopi,' <> ',tpboti
end if

! define the array that will contain the zpos values
allocate(defects%SF(inum)%zpos(nx,ny))
defects%SF(inum)%zpos = -10000.0    ! set all points to be outside the projected SF

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

if (dinfo.eq.1) write (*,*) 'Integer fault box = ',minx,maxx,miny,maxy

! get the equation of the stacking fault plane in the image reference frame
! first the unit plane normal in image space
! we'll take two vectors: ex = from ltop to ttop; ey = from ltop to lbot
ex = tptopi - lptopi
ey = lpboti - lptopi
call NormVec(cell,ex,'c')
call NormVec(cell,ey,'c')
call CalcCross(cell,ex,ey,fpn,'c','c',0)

A(1:3) = fpn ! quat_rotate_vector( conjg(foil%a_fi), dble(fpn(1:3)) )
midpoint = 0.25*(lptopi+lpboti+tptopi+tpboti) ! quat_rotate_vector( conjg(foil%a_fi), dble(0.25*(lptopi+lpboti+tptopi+tpboti) ))
A(4) = sum(A(1:3)*midpoint)
if (dinfo.eq.1) write (*,*) 'fault plane parameters : ',A, midpoint

! rank the corner points so that the polygon is convex
! call rank_points(tpboti(1:2),lpboti(1:2),lptopi(1:2),tptopi(1:2),xx,yy)
xx = (/ lptopi(1), tptopi(1), tpboti(1),lpboti(1) /)
yy = (/ lptopi(2), tptopi(2), tpboti(2),lpboti(2) /)

! for all of the points inside this box:
do i=minx,maxx
  do j=miny,maxy
    if (point_inside_polygon( float(i), float(j), xx, yy ).gt.0) then 
! the point lies inside the projected region, 
! so we need the depth of the SF plane at this position, taking into account the 
! proper coordinate transformation (depth must be expressed in image reference frame)
        defects%SF(inum)%zpos(i+nx/2,j+ny/2) =  DF_L * ( A(4) - A(1)*float(i) - A(2)*float(j) )/A(3)
    end if
  end do
end do
if (dinfo.eq.1) write (*,*) 'fault plane pixels determined'

! let's also make sure that the SF displacement vector is translated to the 
! cartesian reference frame, so that it can be used directly by the CalcR routine
defects%SF(inum)%lpbc = matmul(cell%dsm,defects%SF(inum)%Rdisp)

! that should do it for the stacking fault...  The rest 
! takes place in the CalcR routine.
end subroutine makestackingfault


!--------------------------------------------------------------------------
!
! SUBROUTINE: makestackingfaultECCI
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute parameters for a stacking fault in ECCI mode
!
!> @details  This subroutine computes the geometrical parameters for a 
!> stacking fault.  It computes, among others, the coordinates of the surface
!> intersections of the partial dislocations, and an array that indicates, for
!> each image pixel, whether or not the corresponding integration column 
!> contains this fault; if it does not, the  value in the array is set to 
!> -10000; if it does, then the value is equal to the point where the fault
!> plane intersects with the column, measured from the top surface.
!> In short, anything that is needed in the CalcR routine and can be computed 
!> ahead of time, is computed here.  The routine also calls the makedislocation 
!> routine to create the partials.
! 
!> @param defects defect structure
!> @param cell unit cell pointer
!> @param foil foil structure
!> @param inum
!> @param DF_L column edge length
!> @param nx
!> @param ny
!> @param DF_g
!> @param ndl
!> @param dinfo trigger for verbose output
!
!> @date   11/05/13 MDG 1.0 new attempt to replace faulty original routine
!> @date   11/13/13 MDG 1.1 traced error to problem with transformations in defectmodule
!> @date   11/13/13 MDG 1.2 changed SF normal transformation for zpos array computation (to be tested)
!> @date   12/17/13 MDG 1.3 branch from original routine to deal with different ECCI geometry
!> @date   12/18/13 MDG 1.4 debug of stacking fault location array
!> @date   06/09/14 MDG 2.0 added cell, SF, YD arguments
!> @date   06/10/14 MDG 2.1 added foil argument
!--------------------------------------------------------------------------
subroutine makestackingfaultECCI(defects,cell,foil,inum,DF_L,nx,ny,DF_g,dinfo)
 
use math
use constants
use files
use foilmodule
use YSHmodule
use dislocation
use quaternions
use rotations
use crystal
use symmetry

IMPLICIT NONE

type(defecttype),INTENT(INOUT) :: defects
type(unitcell),pointer	        :: cell
type(foiltype),INTENT(INOUT)   :: foil
real(kind=sgl),INTENT(IN)	:: DF_L
integer(kind=irg),INTENT(IN)	:: inum, nx, ny, DF_g(3), dinfo

real(kind=sgl)  		:: fpn(3),am(4,4),midpoint(3), ex(3), ey(3),&
                  		 lptopi(3),lpboti(3),tptopi(3),tpboti(3),det,A(4), xx(4), yy(4), tmp(3), &
				 planenormal(3), rzero(3), unita(3), lun(3), tun(3), zang, zu, zz

integer(kind=irg) 		:: i,j,info,ipiv,minx,maxx,miny,maxy, ndl

ndl = defects%numYdisl

! we begin by computing the geometry in the foil reference frame, which is the cartesian frame 
! for zero sample tilt;  sample tilts are applied once we known the partial dislocation geometry
call TransSpace(cell,defects%SF(inum)%plane,tmp,'r','c')
call NormVec(cell,tmp,'c')
planenormal =  tmp

call CalcCross(cell, planenormal, (/ 0.0,0.0,1.0 /),  unita, 'c', 'c', 0)
call NormVec(cell,unita,'c')
fpn = planenormal
if (dinfo.eq.1) write (*,*) ' unita should have zero third component ',unita, fpn


! the fault plane goes through the point rzero in the foil top surface
rzero = (/ defects%SF(inum)%id, defects%SF(inum)%jd, sngl(foil%z0)*0.5 /)

! this leads to the location of the partial dislocation intersections, and these must
! be Yoffe dislocations !!!
defects%SF(inum)%lpr = rzero + 0.5*defects%SF(inum)%sep*unita / DF_L
defects%SF(inum)%tpr = rzero - 0.5*defects%SF(inum)%sep*unita / DF_L

if (dinfo.eq.1) write (*,*) 'lpr_i = ',defects%SF(inum)%lpr(1:3)
if (dinfo.eq.1) write (*,*) 'tpr_i = ',defects%SF(inum)%tpr(1:3)


! convert line directions to the Cartesian crystal reference frame
call TransSpace(cell,defects%SF(inum)%lpu,lun,'d','c')
call TransSpace(cell,defects%SF(inum)%tpu,tun,'d','c')
! normalize both vectors
call NormVec(cell,lun,'c')
call NormVec(cell,tun,'c')


! call makeYSHdislocation for each of the partials [must be done in main program]
!  if (.not. allocated(YD)) allocate(defects%YD(3*maxdefects))

  defects%YD(ndl+1)%id = defects%SF(inum)%lpr(1) 
  defects%YD(ndl+1)%jd = defects%SF(inum)%lpr(2)
  defects%YD(ndl+1)%u(1:3) = dble(defects%SF(inum)%lpu(1:3))
  defects%YD(ndl+1)%burg(1:3) = dble(defects%SF(inum)%lpb(1:3))
  defects%YD(ndl+1)%g = float(DF_g)
  defects%YD(ndl+1)%sig = defects%SF(inum)%poisson
  
  defects%YD(ndl+2)%id = defects%SF(inum)%tpr(1) 
  defects%YD(ndl+2)%jd = defects%SF(inum)%tpr(2)
  defects%YD(ndl+2)%u(1:3) = dble(defects%SF(inum)%tpu(1:3))
  defects%YD(ndl+2)%burg(1:3) = dble(defects%SF(inum)%tpb(1:3))
  defects%YD(ndl+2)%g = float(DF_g)
  defects%YD(ndl+2)%sig = defects%SF(inum)%poisson

  call makeYSHdislocation(defects,cell,foil,ndl+1,dinfo,DF_L)
  if (dinfo.eq.1) write (*,*) 'Leading Partial Position ',defects%YD(ndl+1)%id,defects%YD(ndl+1)%jd
  call makeYSHdislocation(defects,cell,foil,ndl+2,dinfo,DF_L)
  if (dinfo.eq.1) write (*,*) 'Trailing Partial Position ',defects%YD(ndl+2)%id,defects%YD(ndl+2)%jd



! first find the length of the dislocation line inside the foil along the
! dislocation z-axis, which is the line direction; also, compute the intersection
! points of the line with the top and bottom surfaces  (all components in [nm])
zang = CalcAngle(cell, defects%YD(ndl+1)%u,foil%F,'d')
zz = cos(zang)
if (abs(zz).gt.0.00001) then 
  zu = abs(foil%z0/zz)
else
  zu = 100000.0         ! this is when the dislocation is nearly parallel to the foil
end if

! transform the line direction to the foil reference frame
  tmp = quat_rotate_vector( foil%a_fc, dble(lun) ) / DF_L

! determine the top and bottom intersection coordinates 
  defects%YD(ndl+1)%top = (/ defects%YD(ndl+1)%id, defects%YD(ndl+1)%jd, foil%z0*0.5 /)
  if (zz.gt.0.0) then  ! u points to the top of the foil
    defects%YD(ndl+1)%bottom = (/ defects%YD(ndl+1)%id - tmp(1)*zu, defects%YD(ndl+1)%jd - tmp(2)*zu, -0.5D0*foil%z0 /)
  else                 ! u points to the bottom of the foil
    defects%YD(ndl+1)%bottom = (/ defects%YD(ndl+1)%id + tmp(1)*zu, defects%YD(ndl+1)%jd + tmp(2)*zu, -0.5D0*foil%z0 /)
  end if  


! first find the length of the dislocation line inside the foil along the
! dislocation z-axis, which is the line direction; also, compute the intersection
! points of the line with the top and bottom surfaces  (all components in [nm])
zang = CalcAngle(cell, defects%YD(ndl+2)%u,foil%F,'d')
zz = cos(zang)
if (abs(zz).gt.0.00001) then 
  zu = abs(foil%z0/zz)
else
  zu = 100000.0         ! this is when the dislocation is nearly parallel to the foil
end if
! transform the line direction to the foil reference frame
  tmp = quat_rotate_vector( foil%a_fc, dble(tun) ) / DF_L

! determine the top and bottom intersection coordinates 
  defects%YD(ndl+2)%top = (/ defects%YD(ndl+2)%id, defects%YD(ndl+2)%jd, foil%z0*0.5 /)
  if (zz.gt.0.0) then  ! u points to the top of the foil
    defects%YD(ndl+2)%bottom = (/ defects%YD(ndl+2)%id - tmp(1)*zu, defects%YD(ndl+2)%jd - tmp(2)*zu, -0.5D0*foil%z0 /)
  else                 ! u points to the bottom of the foil
    defects%YD(ndl+2)%bottom = (/ defects%YD(ndl+2)%id + tmp(1)*zu, defects%YD(ndl+2)%jd + tmp(2)*zu, -0.5D0*foil%z0 /)
  end if  


! copy the top and bottom dislocation intersections
! into the corresponding variables of the SF record

defects%SF(inum)%lpbot = defects%YD(ndl+1)%bottom
defects%SF(inum)%lptop = defects%YD(ndl+1)%top
defects%SF(inum)%tpbot = defects%YD(ndl+2)%bottom
defects%SF(inum)%tptop = defects%YD(ndl+2)%top

! obviously, these four points need to lie in a single plane; at this point, we check that this is indeed the case
! by computing the volume of the tetrahedron formed by these four points; if the volume is zero, then the 
! points are co-planar.  (Use LAPACK's LU-decomposition and compute the product of the diagonal elements of U)
am(1:4,1) = (/ defects%SF(inum)%lptop(1:3),1.0 /)
am(1:4,2) = (/ defects%SF(inum)%lpbot(1:3),1.0 /) 
am(1:4,3) = (/ defects%SF(inum)%tptop(1:3),1.0 /) 
am(1:4,4) = (/ defects%SF(inum)%tpbot(1:3),1.0 /) 
call sgetrf(4,4,am,4,ipiv,info)
det = abs(am(1,1)*am(2,2)*am(3,3)*am(4,4))
if (dinfo.eq.1) write (*,*) 'determinant (should be zero) = ',det

! ok, next we need to figure out which image pixels lie on the projection of the stacking fault plane.
! We need to transform the corner points into the image reference frame !!!
lptopi = quat_rotate_vector( conjg(foil%a_fi), dble(defects%SF(inum)%lptop))
lpboti = quat_rotate_vector( conjg(foil%a_fi), dble(defects%SF(inum)%lpbot))
tptopi = quat_rotate_vector( conjg(foil%a_fi), dble(defects%SF(inum)%tptop))
tpboti = quat_rotate_vector( conjg(foil%a_fi), dble(defects%SF(inum)%tpbot))
if (dinfo.eq.1) then
  write (*,*) 'SF parameters :'
  write (*,*) lptopi,' <> ',lpboti
  write (*,*) tptopi,' <> ',tpboti
end if

! define the array that will contain the zpos values
allocate(defects%SF(inum)%zpos(nx,ny))
defects%SF(inum)%zpos = -10000.0    ! set all points to be outside the projected SF

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

if (dinfo.eq.1) write (*,*) 'Integer fault box = ',minx,maxx,miny,maxy

! get the equation of the stacking fault plane in the image reference frame
! first the unit plane normal in image space
! we'll take two vectors: ex = from ltop to ttop; ey = from ltop to lbot
ex = tptopi - lptopi
ey = lpboti - lptopi
call NormVec(cell,ex,'c')
call NormVec(cell,ey,'c')
call CalcCross(cell,ex,ey,fpn,'c','c',0)

A(1:3) = fpn ! quat_rotate_vector( conjg(foil%a_fi), dble(fpn(1:3)) )
midpoint = 0.25*(lptopi+lpboti+tptopi+tpboti) ! quat_rotate_vector( conjg(foil%a_fi), dble(0.25*(lptopi+lpboti+tptopi+tpboti) ))
A(4) = sum(A(1:3)*midpoint)
if (dinfo.eq.1) write (*,*) 'fault plane parameters : ',A, midpoint

! rank the corner points so that the polygon is convex
! call rank_points(tpboti(1:2),lpboti(1:2),lptopi(1:2),tptopi(1:2),xx,yy)
xx = (/ lptopi(1), tptopi(1), tpboti(1),lpboti(1) /)
yy = (/ lptopi(2), tptopi(2), tpboti(2),lpboti(2) /)

! for all of the points inside this box:
do i=minx,maxx
  do j=miny,maxy
    if (point_inside_polygon( float(i), float(j), xx, yy ).gt.0) then 
! the point lies inside the projected region, 
! so we need the depth of the SF plane at this position, taking into account the 
! proper coordinate transformation (depth must be expressed in image reference frame)
        defects%SF(inum)%zpos(i+nx/2,j+ny/2) = ( A(4) - A(1)*float(i) - A(2)*float(j) )/A(3)
    end if
  end do
end do
if (dinfo.eq.1) write (*,*) 'fault plane pixels determined'

!open(unit=20,file='sf.data',status='unknown',form='unformatted')
!write(20) nx, ny
!write (20) defects%SF(inum)%zpos
!close(unit=20,status='keep')

! let's also make sure that the SF displacement vector is translated to the 
! cartesian reference frame, so that it can be used directly by the CalcR routine
defects%SF(inum)%lpbc = matmul(cell%dsm,defects%SF(inum)%Rdisp)

! that should do it for the stacking fault...  The rest 
! takes place in the CalcR routine.
end subroutine makestackingfaultECCI


!--------------------------------------------------------------------------
!
! SUBROUTINE: read_stacking_fault_data
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  read stacking fault namelist files
!
!> @param defects defect structure
!> @param cell unit cell pointer
!> @param foil foil structure
!> @param DF_L 
!> @param DF_npix number of x-pixels
!> @param DF_npiy number of y-pixels
!> @param DF_gf 
!> @param dinfo logical to trigger verbose output
!> @param ECCI logical optional to indicate ECCI formatting rather than regular TEM
! 
!> @date    1/5/99  MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   06/04/13 MDG 3.0 rewrite
!> @date   12/17/13 MDG 3.1 added ECCI mode
!> @date   06/09/14 MDG 4.0 added cell, defects arguments
!> @date   06/10/14 MDG 4.1 added foil argument
!--------------------------------------------------------------------------
subroutine read_stacking_fault_data(defects,cell,foil,DF_L,DF_npix,DF_npiy,DF_g,dinfo,ECCI)

use YSHmodule
use io
use files
use dislocation

IMPLICIT NONE

type(unitcell),pointer	                :: cell
type(defecttype),INTENT(INOUT)         :: defects
type(foiltype),INTENT(INOUT)           :: foil
integer(kind=irg),INTENT(IN)		:: DF_npix, DF_npiy, DF_g(3), dinfo
real(kind=sgl),INTENT(IN)		:: DF_L
logical,INTENT(IN),OPTIONAL	        :: ECCI

integer(kind=irg) 			:: i,SFplane(3)
real(kind=sgl)     			:: SFi,SFj,SFsep,SFlpu(3),SFlpb(3),SFtpu(3),SFtpb(3),SFR(3), poisson

namelist / SFdata / SFi, SFj, SFsep, SFplane, SFlpu, SFlpb, SFtpu, SFtpb, SFR, poisson

! allocate the necessary memory 
 allocate(defects%SF(defects%numsf))

! read the namelist files for all of the stacking faults
 do i=1,defects%numsf
    call Message('opening '//trim(defects%sfname(i)), frm = "(/A)")
    SFR = (/ 0.0, 0.0, 0.0 /)
    poisson = 0.0
    OPEN(UNIT=dataunit,FILE=trim(defects%sfname(i)),DELIM='APOSTROPHE')
    READ(UNIT=dataunit,NML=SFdata)
    CLOSE(UNIT=dataunit)
! transform the fault fractional coordinates to nm in the image reference frame
    defects%SF(i)%id = SFi * 0.5 * float(DF_npix) ! * DF_L  (zooming is done later in the image reference frame)
    defects%SF(i)%jd = SFj * 0.5 * float(DF_npiy) ! * DF_L
    defects%SF(i)%sep = SFsep
    defects%SF(i)%plane = SFplane
    defects%SF(i)%poisson = poisson
    defects%SF(i)%lpu = SFlpu
    defects%SF(i)%lpb = SFlpb
    defects%SF(i)%tpu = SFtpu
    defects%SF(i)%tpb = SFtpb
    if (sum(abs(SFR)).eq.0.0) then  
      defects%SF(i)%Rdisp = SFlpb
    else
      defects%SF(i)%Rdisp = SFR
    end if
! initialize the stacking fault variables and both partial dislocations; this might depend
! on the imaging mode (TEM vs. ECCI)
    if (present(ECCI)) then
      call makestackingfaultECCI(defects,cell,foil,i,DF_L,DF_npix,DF_npiy,DF_g,dinfo)
      defects%numYdisl = defects%numYdisl + 2
    else
      call makestackingfault(defects,cell,foil,i,DF_L,DF_npix,DF_npiy,DF_g,dinfo)
      defects%numdisl = defects%numdisl + 2
    end if
 end do
 
end subroutine read_stacking_fault_data

end module stacking_fault

