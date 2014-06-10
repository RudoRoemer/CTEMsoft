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
!> @date   11/21/13 MDG 3.0 verification with new libraries and new ECCI program
!> @date   06/09/14 MDG 4.0 removed global variable YD
!--------------------------------------------------------------------------
module YSHModule

use local
use quaternions
use typedefs

IMPLICIT NONE 

contains



!--------------------------------------------------------------------------
!
! FUNCTION: YSHDisp
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  compute the displacement field of an inclined dislocation intersecting the foil surface
!
!> @details compute the displacement field of an inclined dislocation intersecting the top surface of 
!> the foil, taking into account surface relaxations for the isotropic elastic case (cubic only) ... 
!>
!> equations are based on the Shaibani&Hazzledine 1981 paper, along with special limits for 
!> the alpha->0 case, which were derived by MDG using Mathematica. 
!
!> @paraqm defects defects structure
!> @param x dislocation x-coordinate
!> @param y dislocation y-coordinate
!> @param z dislocation z-coordinate
!> @param ii dislocation number
!
!> @todo There is a problem with dislocations normal to the foil surface, likely a typographical error
!> in the SH paper; this needs to be resolved further, which may require explicit repetition of all 
!> analytical computations! Mathematica gives an infinite limit for the bx edge case when normal
!> to the foil surface.
! 
!> @date    1/5/99  MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   06/04/13 MDG 3.0 rewrite
!> @date   11/21/13 MDG 3.1 verification
!> @date   06/09/14 MDG 4.0 added defects as argument
!--------------------------------------------------------------------------
function YSHDisp(defects,x,y,z,ii) result(res)

use constants

IMPLICIT NONE

type(defecttype),INTENT(INOUT) :: defects
real(kind=dbl),INTENT(IN) 	:: x,y,z
integer(kind=irg),INTENT(IN)  	:: ii

real(kind=dbl)    		:: eta, zeta, etap, zetap, r, oms, omts, xx, sgn, om, omp, AA, BB, BBp, th, &
                                 k, lam, alA, alB, u, v, w, ms, S, du, dv, dw, qe, me, De, qx, mx, Dx, rr, eps
real(kind=dbl) 			:: res(3) 

! initialize the geometrical parameters
eta  = y*defects%YD(ii)%ca - z*defects%YD(ii)%sa
zeta = y*defects%YD(ii)%sa + z*defects%YD(ii)%ca
etap  = -y*defects%YD(ii)%ca - z*defects%YD(ii)%sa
zetap = y*defects%YD(ii)%sa - z*defects%YD(ii)%ca
r = sqrt(x**2+y**2+z**2) 
oms = 1.D0-defects%YD(ii)%sig
omts = 1.D0-2.D0*defects%YD(ii)%sig

! cover the special case of negative x values (based on IDL tests)
xx = x
sgn = 1.D0
if (xx.lt.0.D0) then
  xx = dabs(x)
  sgn = -1.D0
else 
  sgn = 1.D0
end if

! more parameters
om =  (datan2(y,xx)-datan2(eta,xx)+datan2(xx*r*defects%YD(ii)%sa,eta*y+xx**2*defects%YD(ii)%ca))
omp= (datan2(y,xx)-datan2(etap,xx)+datan2(xx*r*defects%YD(ii)%sa,etap*y-xx**2*defects%YD(ii)%ca))

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
if (abs(defects%YD(ii)%bs).gt.eps) then 
  ms = xx*sin(2.D0*defects%YD(ii)%alpha)/r/BB
  S = defects%YD(ii)%bs/(4.D0*cPi)
  if (defects%YD(ii)%alpha.gt.0.01) then 
    du = xx*ms+2.D0*eta*defects%YD(ii)%ca**2/BB+2.D0*omts*defects%YD(ii)%cota*(-1.D0+defects%YD(ii)%ca+&
            defects%YD(ii)%ca*alA-y*defects%YD(ii)%sa/AA-alB)-sin(2.D0*defects%YD(ii)%alpha)
    dv = y*ms-2.D0*xx*defects%YD(ii)%ca/BB-defects%YD(ii)%sa*(omp-om)+2.D0*omts*defects%YD(ii)%cota* &
         (xx*defects%YD(ii)%sa/AA-om*defects%YD(ii)%ca)
    dw = z*ms+defects%YD(ii)%ca*(omp-om)-2.D0*omts*om*defects%YD(ii)%ca
  else 
    du = 2.D0*y/(r-z)
    dv = -2.D0*xx*(r+z)/(xx**2+y**2)
    dw = cPi + datan2(y,xx) - datan2(-y,xx)
  end if
  u = u+du*S
  v = v-sgn*dv*S
  w = w+sgn*dw*S
end if

! then the edge component in the y-z plane
if (abs(defects%YD(ii)%be).gt.eps) then 
  qe = xx*(1.D0/BBp-1.D0/BB+2.D0*z*defects%YD(ii)%ca/BB**2)
  me = -qe/r-4.D0*oms*xx*defects%YD(ii)%ca**2/r/BB
  De = defects%YD(ii)%be/(8.D0*cPi*oms)
  if (defects%YD(ii)%alpha.gt.0.01) then 
    k = 4.D0*oms*omts*defects%YD(ii)%cota**2
    du = xx*me+lam+2.D0*defects%YD(ii)%ca*(z+2.D0*oms*eta*defects%YD(ii)%sa)/BB-4.D0*oms*defects%YD(ii)%sa**2+&
           k*(1.D0-defects%YD(ii)%ca-defects%YD(ii)%ca*alA+y*defects%YD(ii)%sa/AA+alB)
    dv = y*me+qe*defects%YD(ii)%sa+th*defects%YD(ii)%ca+k*(-xx*defects%YD(ii)%sa/AA+om*defects%YD(ii)%ca)
    dw = z*me+qe*defects%YD(ii)%ca+th*defects%YD(ii)%sa-2.D0*xx*defects%YD(ii)%ca*(1.D0/BBp+omts/BB)+k*om*defects%YD(ii)%sa
!    write (*,*) du,dv,dw
  else 
    rr = xx**2+y**2
    du = 2.D0*z/(r-z)+4.D0*xx**2*(defects%YD(ii)%sig*rr-r**2)/r/AA**2/(r+z)+2.D0*omts*oms*((xx**2+z*(z-r))/AA**2+alA)+&
         omts*dlog((r+z)/AA)
    dv = 4.D0*xx*y*(rr*defects%YD(ii)%sig-r**2)/r/AA**2/(r+z)+2.D0*xx*y*(rr+2.D0*z*(r+z))*oms*omts/rr**2+&
            2.D0*oms*(cPi + datan2(y,xx) - datan2(-y,xx))
    dw = 4.D0*xx*rr*defects%YD(ii)%sig*(z-2.D0*r*oms)/r/AA**2/(r+z)
  end if
  u = u+du*De
  v = v+sgn*dv*De
  w = w+sgn*dw*De
end if

! and finally the bx edge component
if (abs(defects%YD(ii)%bx).gt.eps) then 
  qx = etap/BBp-eta/BB-2.D0*z*eta*defects%YD(ii)%ca/BB**2
  mx = -qx/r+2.D0*omts*y*defects%YD(ii)%ca/r/BB
  Dx = defects%YD(ii)%bx/(8.D0*cPi*oms)
  if (defects%YD(ii)%alpha.gt.0.01) then 
    k = 4.D0*oms*omts*defects%YD(ii)%cota**2
    du = xx*mx+th+k*(xx*defects%YD(ii)%ta/AA-om)
    dv = y*mx+qx*defects%YD(ii)%sa-lam*defects%YD(ii)%ca-2.D0*defects%YD(ii)%ca*&
        (z*defects%YD(ii)%ca+omts*y*defects%YD(ii)%sa)/BB+ &
        k*(-1.D0+defects%YD(ii)%ca-alA+y*defects%YD(ii)%ta/AA+defects%YD(ii)%ca*alB)
    dw = z*mx+qx*defects%YD(ii)%ca-lam*defects%YD(ii)%sa-2.D0*etap*defects%YD(ii)%ca/BBp+&
        4.D0*defects%YD(ii)%ca*(oms*y*defects%YD(ii)%ca-omts*z*defects%YD(ii)%sa)/BB+ &
        k*defects%YD(ii)%ta*(defects%YD(ii)%ca-alA+defects%YD(ii)%ca*alB)+4.D0*oms*defects%YD(ii)%ca*defects%YD(ii)%cota
 else 
    rr = xx**2+y**2
    du = -4.D0*xx*y*(rr*defects%YD(ii)%sig-r**2)/r/AA**2/(r+z)-2.D0*xx*y*(rr+2.D0*z*(r+z))*oms*omts/rr**2+&
            2.D0*oms*(cPi + datan2(y,xx) - datan2(-y,xx))
    dv = 2.D0*z/(r-z)-4.D0*y**2*(defects%YD(ii)%sig*rr-r**2)/r/AA**2/(r+z)+2.D0*omts*oms*(-1.D0+(z*(r-z)-y**2)/AA**2-alA)- &
            omts*dlog((r+z)/AA)
    dw = 0.D0     ! not sure if this limit is correct ... Mathematica gives a directedinfinity value for the limit, which might mean that the 
    ! original YSH expression in the paper is incorrect for the w component ... this needs to be rederived and verified !!!
  end if

  u = u+sgn*du*Dx
  v = v+dv*Dx
  w = w+dw*Dx
end if

! and return the displacement components
res = (/ u,v,w /)

end function YSHDisp


!--------------------------------------------------------------------------
!
! FUNCTION: makeYSHdislocation
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief pre-compute geometrical parametersf or the Yoffe&Shaibani&Hazzledine (YSH) 
!> surface-relaxed dislocation in an elastically isotropic matrix. 
!
!> @details These parameters are then used in the CalcR routine.
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
!> @param defects defects structure
!> @param cell unit cell pointer
!> @param foil foil structure
!> @param i dislocation number
!> @param dinfo triggers verbose output
!> @param L column edge length
!
!> @todo Convert IO to Write_Value calls
! 
!> @date  1/5/99  MDG 1.0 original
!> @date  5/19/01 MDG 2.0 f90 version
!> @date 11/27/01 MDG 2.1 added kind support
!> @date 06/04/13 MDG 3.0 rewrite+added quaternions
!> @date 11/21/13 MDG 3.1 verification + rewrite of output handling
!> @date 06/09/14 MDG 4.0 added cell and defects arguments
!> @date 06/10/14 MDG 4.1 added foil argument
!--------------------------------------------------------------------------
subroutine makeYSHdislocation(defects,cell,foil,i,dinfo, L)    

use foilmodule
use constants
use crystal
use io
use error
use rotations

IMPLICIT NONE

type(unitcell),pointer	              :: cell
type(defecttype),INTENT(INOUT)       :: defects
type(foiltype),INTENT(INOUT)         :: foil
integer(kind=irg),INTENT(IN)        :: i
integer(kind=irg),INTENT(IN)        :: dinfo
real(kind=sgl),INTENT(IN)           :: L

real(kind=dbl)  	             :: alpha, beta, tu(3), tx(3), ty(3), te(3), tb(3), bl, fx(3), fy(3), fz(3), &
                                       dx, dy, a_di(3,3), io_real(3)

! first, determine the alpha angle between the 
! negative z-axis, which is really the negative foil normal, and the line direction
! (make sure to reduce the angle to [0,90] interval).
! Each YSH dislocation must have a line direction that points INTO the foil !!!
alpha = CalcAngle(cell,foil%F,dble(defects%YD(i)%u),'d')*180.0/cPi
if (alpha.ge.90.0) then 
  alpha = 180.0-alpha
  if (dinfo.eq.1) then 
    io_real(1:3) = foil%F(1:3)
    call WriteValue('Foil normal = ', io_real, 3, "('[',3F5.1,']')")
    io_real(1:3) = defects%YD(i)%u(1:3)
    call WriteValue('line direction = ', io_real, 3, "('[',3F5.1,']')")
    io_real(1) = alpha
    call WriteValue(' --> alpha angle = ', io_real, 1, "(F5.1)")
  end if
  alpha = alpha*cPi/180.0
else 
  call FatalError('makeYSHdislocation','YSH dislocations must have line directions pointing into the foil ! ')
end if

! normalize the line direction
call TransSpace(cell,defects%YD(i)%u,tu,'d','c')
call NormVec(cell,tu,'c')

! consider the case of alpha=0 separately
if (alpha.gt.0.0) then
  call TransSpace(cell,foil%F,ty,'d','c')
  call NormVec(cell,ty,'c')                     !  F
  call CalcCross(cell,tu,ty,tx,'c','c',0)       ! x = u x F
  call NormVec(cell,tx,'c')
  call CalcCross(cell,tx,tu,te,'c','c',0)       ! e = x x u
  call NormVec(cell,te,'c')
  call CalcCross(cell,ty,tx,ty,'c','c',0)
  call NormVec(cell,ty,'c')
else
  tx = foil%qn
  call CalcCross(cell,tx,tu,te,'c','c',0)       ! e = x x u
  call NormVec(cell,te,'c')
end if  
bl = CalcLength(cell,defects%YD(i)%burg,'d')

if (dinfo.eq.1) then 
  io_real(1:3) = tx(1:3)
  call WriteValue(' tx = ',io_real, 3, "(3F5.1)")
  io_real(1:3) = te(1:3)
  call WriteValue(' te = ',io_real, 3, "(3F5.1)")
  io_real(1:3) = tu(1:3)
  call WriteValue(' tu = ',io_real, 3, "(3F5.1)")
  io_real(1:3) = ty(1:3)
  call WriteValue(' ty = ',io_real, 3, "(3F5.1)")
  io_real(1) = bl
  call WriteValue(' bl = ',io_real, 1, "(F8.3)")
end if

call TransSpace(cell,defects%YD(i)%burg,tb,'d','c')
call NormVec(cell,tb,'c')
defects%YD(i)%bx = bl * CalcDot(cell,tb,tx,'c')   ! edge component normal to cut plane
defects%YD(i)%be = bl * CalcDot(cell,tb,te,'c')   ! edge component in cut plane
defects%YD(i)%bs = bl * CalcDot(cell,tb,tu,'c')   ! screw component

if (dinfo.eq.1) then 
  io_real(1:3) = (/ defects%YD(i)%bx,defects%YD(i)%be,defects%YD(i)%bs /)
  call WriteValue('Burgers vector components (bx,be,bs) ', io_real, 3, "(3F12.6)") 
end if
! verified MDG 7/31/11


! we will also need to know the quaternion rotation between the dislocation reference frame 
! and the foil reference frame, so that we can transform the foil coordinates to defect 
! coordinates...  We need the angle beta between the defect x axis (tx) and the foil x axis,
! which is the first column of the foil%a_fc matrix ...   We must make sure that this angle
! is measured in a CCW sense.

! projection of defect x axis onto foil x and y axes
call TransSpace(cell,foil%q,fx,'d','c')
call TransSpace(cell,foil%F,fz,'d','c')
call NormVec(cell,fx,'c')
call NormVec(cell,fz,'c')
call CalcCross(cell,fz,fx,fy,'c','c',0)
dx = CalcDot(cell,tx,fx,'c')
dy = CalcDot(cell,tx,fy,'c')

! use the arctan function to get the angle with correct quadrant computation
defects%YD(i)%beta = atan2(dy,dx) !+ cPi*0.5

if (dinfo.eq.1) then 
  io_real(1) = dx
  call WriteValue(' dx = ', io_real, 1, "(F8.3)")
  io_real(1) = dy
  call WriteValue(' dy = ', io_real, 1, "(F8.3)")
  io_real(1) = defects%YD(i)%beta
  call WriteValue(' beta = ', io_real, 1, "(F8.3)")
end if

! convert to a quaternion
beta = defects%YD(i)%beta
a_di(1,1:3) = (/ cos(beta), sin(beta), 0.D0 /)
a_di(2,1:3) = (/ -sin(beta), cos(beta), 0.D0 /)
a_di(3,1:3) = (/ 0.D0, 0.D0, 1.D0 /)
defects%YD(i)%a_di = om2qu(a_di)
defects%YD(i)%a_id = conjg(defects%YD(i)%a_di)

if (dinfo.eq.1) then 
  write (*,*) 'beta = ',beta
  write (*,*) defects%YD(i)%a_di
  write (*,*) defects%YD(i)%a_id
end if


! finally some geometrical parameters needed for the displacement field computation...
defects%YD(i)%alpha =  alpha
defects%YD(i)%ca = cos(alpha)
defects%YD(i)%sa = sin(alpha)
defects%YD(i)%ta = tan(alpha)
defects%YD(i)%cota = 1.0/defects%YD(i)%ta

! that's it! the rest is handled in the CalcR routine.

end subroutine makeYSHdislocation


!--------------------------------------------------------------------------
!
! SUBROUTINE: read_YSH_dislocation_data
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read Yoffe dislocation input files
! 
!> @param defects defects structure
!> @param cell unit cell pointer
!> @param foil foil structure
!> @param DF_npix number of x-pixels
!> @param DF_npiy number of y-pixels
!> @param DF_gf 
!> @param L 
!> @param dinfo logical to trigger verbose output
! 
!> @date  1/5/99  MDG 1.0 original
!> @date  5/19/01 MDG 2.0 f90 version
!> @date 11/27/01 MDG 2.1 added kind support
!> @date 03/25/13 MDG 3.0 updated IO
!> @date 11/21/13 MDG 3.1 verification
!> @date 06/10/14 MDG 4.0 added defects, cell and foil arguments
!--------------------------------------------------------------------------
subroutine read_YSH_dislocation_data(defects,cell,foil,DF_npix,DF_npiy,DF_gf,L,dinfo)

use io
use files

IMPLICIT NONE

type(defecttype),INTENT(INOUT) :: defects
type(unitcell),pointer	        :: cell
type(foiltype),INTENT(INOUT)   :: foil
integer(kind=irg),INTENT(IN)	:: DF_npix, DF_npiy, dinfo
real(kind=sgl),INTENT(IN)	:: DF_gf(3), L

integer(kind=irg) 	        :: i
real(kind=sgl) 			:: id,jd,u(3),bv(3),poisson

namelist / dislocationdata / id, jd, u, bv, poisson

! allocate the memory for the dislocation parameters
allocate(defects%YD(3*maxdefects))

! these are just the individual dislocations; the ones that belong to 
! stacking faults are handled separately
   do i=1,defects%numYdisl
    call Message('opening '//trim(defects%dislYname(i)), frm = "(/A)")
    OPEN(UNIT=dataunit,FILE=trim(defects%dislYname(i)),DELIM='APOSTROPHE')
    READ(UNIT=dataunit,NML=dislocationdata)
    CLOSE(UNIT=dataunit)
    
! top-of-the-foil intersection of dislocation line is transformed to foil coordinates [nm] with DL(i)%kd=0 (center of foil) [verified 4/23/11]
! the point (0,0) is at the center of the image ... hence the factor of 0.5
    defects%YD(i)%id = id * 0.5 * float(DF_npix) ! * L   scaling (zooming) is done later in the image reference frame...
    defects%YD(i)%jd = jd * 0.5 * float(DF_npiy) ! * L
    defects%YD(i)%u = u
    defects%YD(i)%burg = bv
    defects%YD(i)%g = DF_gf
    defects%YD(i)%sig = poisson
    
! and pre-compute the dislocation displacement field parameters
    call makeYSHdislocation(defects,cell,foil,i,dinfo, L)    
  end do
  
end subroutine read_YSH_dislocation_data


end module YSHModule


