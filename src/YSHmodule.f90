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
! CTEMsoft2013:YSHModule.f90
!--------------------------------------------------------------------------
!
! MODULE: YSHModule
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Provides routines to compute the displacement vector for surface intersecting dislocations.
! 
!> @date   04/29/11 MDG 1.0 original
!> @date   06/04/13 MDG 2.0 rewrite + quaternions instead of rotations
!--------------------------------------------------------------------------
module YSHModule

use local
use quaternions

IMPLICIT NONE

type YDtype
  real(kind=dbl)     	:: burg(3), burgd(3), u(3), un(3), g(3), gn(3), id, jd, zu, bs, be, bx, beta
  real(kind=dbl)     	:: alpha, ca, sa, ta, cota,  top(3), bottom(3), sig
  real(kind=dbl)	:: a_dc(4), a_id(4), a_di(4)
end type YDtype

type (YDtype), allocatable  :: YD(:)    


contains



!--------------------------------------------------------------------------
!
! FUNCTION: YSHDisp
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  compute the displacement field of an inclined dislocation intersecting the foil surface
!
!> @details compute the displacement field of an inclined dislocation intersection the top surface of 
!> the foil, taking into account surface relaxations for the isotropic elastic case (cubic only) ... 
!>
!> equations are based on the Shaibani&Hazzledine 1981 paper, along with special limits for 
!> the alpha->0 case, which were derived by MDG using Mathematica.
!
!> @param x dislocation x-coordinate
!> @param y dislocation y-coordinate
!> @param z dislocation z-coordinate
!> @param ii dislocation number
!
!> @todo There may be a problem with dislocations normal to the foil surface, likely a typographical error
!> in the SH paper; this needs to be resolved further, which may require explicit repetition of all 
!> analytical computations!
! 
!> @date    1/5/99   MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   06/04/13 MDG 3.0 rewrite
!--------------------------------------------------------------------------
function YSHDisp(x,y,z,ii) result(res)

use local
use constants

IMPLICIT NONE

real(kind=dbl),INTENT(IN) 	:: x,y,z
integer(kind=irg),INTENT(IN)  	::  ii

real(kind=dbl)    			:: eta, zeta, etap, zetap, r, oms, omts, xx, sgn, om, omp, AA, BB, BBp, th, &
                                 		k, lam, alA, alB, u, v, w, ms, S, du, dv, dw, qe, me, De, qx, mx, Dx, rr, eps
real(kind=dbl) 				:: res(3) 

! initialize the geometrical parameters
eta  = y*YD(ii)%ca - z*YD(ii)%sa
zeta = y*YD(ii)%sa + z*YD(ii)%ca
etap  = -y*YD(ii)%ca - z*YD(ii)%sa
zetap = y*YD(ii)%sa - z*YD(ii)%ca
r = sqrt(x**2+y**2+z**2) 
oms = 1.D0-YD(ii)%sig
omts = 1.D0-2.D0*YD(ii)%sig

! cover the special case of negative x values (based on IDL tests)
xx = x
if (xx.lt.0.D0) then
  xx = dabs(x)
  sgn = -1.D0
else 
  sgn = 1.D0
end if

! more parameters
om =  (datan2(y,x)-datan2(eta,x)+datan2(x*r*YD(ii)%sa,eta*y+x**2*YD(ii)%ca))
omp= (datan2(y,x)-datan2(etap,x)+datan2(x*r*YD(ii)%sa,etap*y-x**2*YD(ii)%ca))

AA = r-z
BB  = r-zeta
BBp = r-zetap 
th = 2.D0*oms*(omp-om)
lam = omts*dlog(BBp/BB)
alA = dlog(AA)
alB = dlog(BB)


u = 0.D0
v = 0.D0
w = 0.D0
eps = 1.0D-6

! screw component first
if (abs(YD(ii)%bs).gt.eps) then 
  ms = x*sin(2.D0*YD(ii)%alpha)/r/BB
  S = YD(ii)%bs/(4.D0*cPi)
  if (YD(ii)%alpha.gt.0.01) then 
    du = x*ms+2.D0*eta*YD(ii)%ca**2/BB+2.D0*omts*YD(ii)%cota*(-1.D0+YD(ii)%ca+&
            YD(ii)%ca*alA-y*YD(ii)%sa/AA-alB)-sin(2.D0*YD(ii)%alpha)
    dv = y*ms-2.D0*x*YD(ii)%ca/BB-YD(ii)%sa*(omp-om)+2.D0*omts*YD(ii)%cota*(x*YD(ii)%sa/AA-om*YD(ii)%ca)
    dw = z*ms+YD(ii)%ca*(omp-om)-2.D0*omts*om*YD(ii)%ca
  else 
    du = 2.D0*y/(r-z)
    dv = -2.D0*x*(r+z)/(x**2+y**2)
    dw = cPi + datan2(y,x) - datan2(-y,x)
  end if
  u = u+du*S
  v = v-sgn*dv*S
  w = w+sgn*dw*S
end if

! then the edge component in the y-z plane
if (abs(YD(ii)%be).gt.eps) then 
  qe = x*(1.D0/BBp-1.D0/BB+2.D0*z*YD(ii)%ca/BB**2)
  me = -qe/r-4.D0*oms*x*YD(ii)%ca**2/r/BB
  De = YD(ii)%be/(8.D0*cPi*oms)
  if (YD(ii)%alpha.gt.0.01) then 
    k = 4.D0*oms*omts*YD(ii)%cota**2
    du = x*me+lam+2.D0*YD(ii)%ca*(z+2.D0*oms*eta*YD(ii)%sa)/BB-4.D0*oms*YD(ii)%sa**2+k*(1.D0-YD(ii)%ca-YD(ii)%ca*alA+&
             y*YD(ii)%sa/AA+alB)
    dv = y*me+qe*YD(ii)%sa+th*YD(ii)%ca+k*(-x*YD(ii)%sa/AA+om*YD(ii)%ca)
    dw = z*me+qe*YD(ii)%ca+th*YD(ii)%sa-2.D0*x*YD(ii)%ca*(1.D0/BBp+omts/BB)+k*om*YD(ii)%sa
!    write (*,*) du,dv,dw
  else 
    rr = x**2+y**2
    du = 2.D0*z/(r-z)+4.D0*x**2*(YD(ii)%sig*rr-r**2)/r/AA**2/(r+z)+2.D0*omts*oms*((x**2+z*(z-r))/AA**2+alA)+omts*dlog((r+z)/AA)
    dv = 4.D0*x*y*(rr*YD(ii)%sig-r**2)/r/AA**2/(r+z)+2.D0*x*y*(rr+2.D0*z*(r+z))*oms*omts/rr**2+&
            2.D0*oms*(cPi + datan2(y,x) - datan2(-y,x))
    dw = 4.D0*x*rr*YD(ii)%sig*(z-2.D0*r*oms)/r/AA**2/(r+z)
  end if
  u = u+du*De
  v = v+sgn*dv*De
  w = w+sgn*dw*De
end if

! and finally the bx edge component
if (abs(YD(ii)%bx).gt.eps) then 
  qx = etap/BBp-eta/BB-2.D0*z*eta*YD(ii)%ca/BB**2
  mx = -qx/r+2.D0*omts*y*YD(ii)%ca/r/BB
  Dx = YD(ii)%bx/(8.D0*cPi*oms)
  if (YD(ii)%alpha.gt.0.01) then 
    k = 4.D0*oms*omts*YD(ii)%cota**2
    du = x*mx+th+k*(x*YD(ii)%ta/AA-om)
    dv = y*mx+qx*YD(ii)%sa-lam*YD(ii)%ca-2.D0*YD(ii)%ca*(z*YD(ii)%ca+omts*y*YD(ii)%sa)/BB+k*(-1.D0+YD(ii)%ca-alA+y*YD(ii)%ta/AA+ &
          YD(ii)%ca*alB)
    dw = z*mx+qx*YD(ii)%ca-lam*YD(ii)%sa-2.D0*etap*YD(ii)%ca/BBp+4.D0*YD(ii)%ca*(oms*y*YD(ii)%ca-omts*z*YD(ii)%sa)/BB+ &
           k*YD(ii)%ta*(YD(ii)%ca-alA+YD(ii)%ca*alB)+4.D0*oms*YD(ii)%ca*YD(ii)%cota
 else 
    rr = x**2+y**2
    du = -4.D0*x*y*(rr*YD(ii)%sig-r**2)/r/AA**2/(r+z)-2.D0*x*y*(rr+2.D0*z*(r+z))*oms*omts/rr**2+&
            2.D0*oms*(cPi + datan2(y,x) - datan2(-y,x))
    dv = 2.D0*z/(r-z)-4.D0*y**2*(YD(ii)%sig*rr-r**2)/r/AA**2/(r+z)+2.D0*omts*oms*(-1.D0+(z*(r+z)-y**2)/AA**2+alA)- &
            omts*dlog((r+z)/AA)
    dw = 0.D0     ! not sure if this limit is correct ... Mathematica gives a directedinfinity value for the limit, which mght mean that the 
    ! original YSH expression in the paper is incorrect for the w component ... 
  end if

  u = u+sgn*du*Dx
  v = v+dv*Dx
  w = w+dw*Dx
end if

! and return the displacement components
!write (*,*) u,v,w
res = (/ u,v,w /)

end function YSHDisp


!--------------------------------------------------------------------------
!
! FUNCTION: makeYSHdislocation
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief pre-compute geometrical parametersfor the Yoffe&Shaibani&Hazzledine (YSH) surface-relaxed dislocation in an elastically
!> isotropic matrix. 
!
!> @detailsThese parameters are then used in the CalcR routine.
!>
!> We implemented the YSH expressions instead of Yoffe's 
!> since the former are more easily handled for numerical computations.
!>
!> SH have redefined the x-y-z reference frame used by Yoffe to fall along
!> the dislocation line itself.  As a result, the Burgers vector must be decomposed
!> into a screw component and two edge components, one in the plane of the 
!> discontinuity, the other normal to that plane (which is by definition the x-axis).
!> Check the SH paper for more details.
!
!> @param i dislocation number
!> @param dinfo triggers verbose output
!> @param L column edge length
!
!> @todo Convert IO to Write_Value calls
! 
!> @date    1/5/99   MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   06/04/13 MDG 3.0 rewrite+added quaternions
!--------------------------------------------------------------------------
subroutine makeYSHdislocation(i,dinfo, L)    

use local
use foilmodule
use constants
use crystal
use rotations

IMPLICIT NONE

real(kind=sgl),INTENT(IN)	:: L
real(kind=dbl)  	:: alpha, beta, tu(3), tx(3), ty(3), te(3), tb(3), bl, fx(3), fy(3), fz(3), dx, dy, a_di(3,3)
integer(kind=irg) 	:: i,dinfo

! first, determine the alpha angle between the 
! negative z-axis, which is really the negative foil normal, and the line direction
! (make sure to reduce the angle to [0,90] interval).
! Each YSH dislocation must have a line direction that points INTO the foil !!!
alpha = CalcAngle(foil%F,dble(YD(i)%u),'d')*180.0/cPi
if (alpha.ge.90.0) then 
  alpha = 180.0-alpha
  write (*,*) 'Foil normal = ',foil%F
  write (*,*) 'line direction = ',YD(i)%u
  write (*,*) '  --> alpha angle = ',alpha
  alpha = alpha*cPi/180.0
else 
  write (*,*) 'YSH dislocations must have line directions pointing into the foil ! '
  stop
end if

! normalize the line direction
call TransSpace(YD(i)%u,tu,'d','c')
call NormVec(tu,'c')

! consider the case of alpha=0 separately
if (alpha.gt.0.0) then
  call TransSpace(foil%F,ty,'d','c')
  call NormVec(ty,'c')   !  F
  call CalcCross(tu,ty,tx,'c','c',0)       ! x = u x F
  call NormVec(tx,'c')
  call CalcCross(tx,tu,te,'c','c',0)             ! e = x x u
  call NormVec(te,'c')
  call CalcCross(ty,tx,ty,'c','c',0)
  call NormVec(ty,'c')
else
  tx = foil%qn
  call CalcCross(tx,tu,te,'c','c',0)             ! e = x x u
  call NormVec(te,'c')
end if  
bl = CalcLength(YD(i)%burg,'d')
write (*,*) ' tx = ',tx
write (*,*) ' te = ',te
write (*,*) ' tu = ',tu
write (*,*) ' ty = ',ty
write (*,*) ' bl = ',bl
call TransSpace(YD(i)%burg,tb,'d','c')
call NormVec(tb,'c')
YD(i)%bx = bl * CalcDot(tb,tx,'c')   ! edge component normal to cut plane
YD(i)%be = bl * CalcDot(tb,te,'c')   ! edge component in cut plane
YD(i)%bs = bl * CalcDot(tb,tu,'c')  ! screw component

write (*,*) 'Burgers vector components (bx,be,bs) ',YD(i)%bx,YD(i)%be,YD(i)%bs
! verified MDG 7/31/11


! we will also need to know the quaternion rotation between the dislocation reference frame 
! and the foil reference frame, so that we can transform the foil coordinates to defect 
! coordinates...  We need the angle beta between the defect x axis (tx) and the foil x axis,
! which is the first column of the foil%a_fc matrix ...   We must make sure that this angle
! is measured in a CCW sense.
! projection of defect x axis onto foil x and y axes
call TransSpace(foil%q,fx,'d','c')
call TransSpace(foil%F,fz,'d','c')
call NormVec(fx,'c')
call NormVec(fz,'c')
call CalcCross(fz,fx,fy,'c','c',0)
dx = CalcDot(tx,fx,'c')
dy = CalcDot(tx,fy,'c')
! use the arctan function to get the angle with correct quadrant computation
YD(i)%beta = atan2(dy,dx)
write (*,*) ' dx = ',dx
write (*,*) ' dy = ',dy
write (*,*) ' beta = ',YD(i)%beta
beta = YD(i)%beta
a_di(1,1:3) = (/ cos(beta), sin(beta), 0.D0 /)
a_di(2,1:3) = (/ -sin(beta), cos(beta), 0.D0 /)
a_di(3,1:3) = (/ 0.D0, 0.D0, 1.D0 /)
YD(i)%a_di = om2qu(a_di)
YD(i)%a_id = conjg(YD(i)%a_di)

write (*,*) YD(i)%a_di

! finally some geometrical parameters needed for the displacement field computation...
YD(i)%alpha =  alpha
YD(i)%ca = cos(alpha)
YD(i)%sa = sin(alpha)
YD(i)%ta = tan(alpha)
YD(i)%cota = 1.0/YD(i)%ta

write (*,*) 'angulars = ',cos(alpha),sin(alpha),tan(alpha)

! that's it! the rest is handled in the CalcR routine.

end subroutine makeYSHdislocation


!--------------------------------------------------------------------------
!
! SUBROUTINE: read_YSH_dislocation_data
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  prints quaternion and equivalent rotation matrix for a given input quaternion
! 
!> @param dislYname name of dislocation namelist file (string array)
!> @param numYdisl number of dislocations
!> @param DF_npix number of x-pixels
!> @param DF_npiy number of y-pixels
!> @param DF_gf 
!> @param L 
!> @param dinfo logical to trigger verbose output
! 
!> @date 1/5/99   MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   03/25/13 MDG 3.0 updated IO
!--------------------------------------------------------------------------
subroutine read_YSH_dislocation_data(dislYname,numYdisl,DF_npix,DF_npiy,DF_gf,L,dinfo)

use local
use io
use files

IMPLICIT NONE

character(50),INTENT(IN)	 	:: dislYname(3*maxdefects)
integer(kind=irg),INTENT(IN)	:: numYdisl, DF_npix, DF_npiy, dinfo
real(kind=sgl),INTENT(IN)	:: DF_gf(3), L

integer(kind=irg) 			:: i
real(kind=sgl) 				:: id,jd,u(3),bv(3),poisson

namelist / dislocationdata / id, jd, u, bv, poisson

! allocate the memory for the dislocation parameters
  allocate(YD(numYdisl))

! these are just the individual dislocations; the ones that belong to 
! stacking faults are handled separately
   do i=1,numYdisl
    mess = 'opening '//dislYname(i); call Message("(/A)")
    OPEN(UNIT=dataunit,FILE=dislYname(i),DELIM='APOSTROPHE')
    READ(UNIT=dataunit,NML=dislocationdata)
    CLOSE(UNIT=dataunit)
    
! top-of-the-foil intersection of dislocation line is transformed to foil coordinates [nm] with DL(i)%kd=0 (center of foil) [verified 4/23/11]
! the point (0,0) is at the center of the image ... hence the factor of 0.5
    YD(i)%id = id * 0.5 * float(DF_npix) ! * L   scaling (zooming) is done later in the image reference frame...
    YD(i)%jd = jd * 0.5 * float(DF_npiy) ! * L
    YD(i)%u = u
    YD(i)%burg = bv
    YD(i)%g = DF_gf
    YD(i)%sig = poisson
    
! and pre-compute the dislocation displacement field parameters
       call makeYSHdislocation(i,dinfo, L)    
  end do
  
end subroutine read_YSH_dislocation_data


end module YSHModule


