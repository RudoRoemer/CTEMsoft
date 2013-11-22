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
! CTEMsoft2013:dislocation.f90
!--------------------------------------------------------------------------
!
! MODULE: dislocation
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Provides routines to compute the displacement vector for various dislocations.
! 
!> @date   04/29/11 MDG 1.0 original
!> @date   06/04/13 MDG 2.0 rewrite + quaternions instead of rotations
!--------------------------------------------------------------------------
module dislocation

use local
use quaternions

! define a dislocation, along with all matrices and such, but using quaternions instead of rotations !
type dislocationtype
  real(kind=dbl)     		:: burg(3),burgd(3),u(3),un(3),g(3),gn(3),id,jd, zfrac, zu
  real(kind=dbl)     		:: top(3), bottom(3)
  real(kind=dbl)		:: a_dc(4), a_id(4), a_di(4), a_df(4)
  complex(kind=dbl)  		:: dismat(3,3),pa(3)
end type dislocationtype

type (dislocationtype), allocatable  :: DL(:)    


contains

!--------------------------------------------------------------------------
!
! SUBROUTINE: TransFourthRankTensor 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  transform a fourth rank tensor using a given transformation matrix
! 
!> @note This is one of the very few places in this package where we still use the 
!> good old-fashioned rotation matrix instead of quaternions... Note also that we
!> use the 6x6 notation for the tensors, so we need to convert them to real tensor
!> notation before carrying out the rotations.
!
!> @param al rotation matrix
!> @param cin unrotated tensor
!> @param cout rotated tensor
! 
!> @date 1/5/99   MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   06/04/13 MDG 3.0 rewrite
!--------------------------------------------------------------------------
subroutine TransFourthRankTensor(al,cin,cout)

use local

IMPLICIT NONE

real(kind=dbl),INTENT(IN)	:: al(3,3)
real(kind=sgl),INTENT(IN)	:: cin(6,6)
real(kind=sgl),INTENT(OUT)	:: cout(6,6) 

real(kind=sgl)				:: cold(3,3,3,3), cnew(3,3,3,3)
integer(kind=irg)			:: i,j,k,l,p,q,r,s,delta(3,3),gamma(6,2)

! initalize a bunch of variables
cold = 0.0
cnew = 0.0
cout = 0.0
delta(1,1) = 1; delta(1,2) = 6; delta(1,3) = 5
delta(2,1) = 6; delta(2,2) = 2; delta(2,3) = 4
delta(3,1) = 5; delta(3,2) = 4; delta(3,3) = 3
gamma(1,1) = 1; gamma(1,2) = 1
gamma(2,1) = 2; gamma(2,2) = 2
gamma(3,1) = 3; gamma(3,2) = 3
gamma(4,1) = 2; gamma(4,2) = 3
gamma(5,1) = 1; gamma(5,2) = 3
gamma(6,1) = 1; gamma(6,2) = 2

! convert to real tensor indices
do i=1,3
 do j=1,3
  do k=1,3
   do l=1,3
    cold(i,j,k,l) = cin(delta(i,j),delta(k,l))
   end do
  end do
 end do
end do

! and transform
do i=1,3
 do j=1,3
  do k=1,3
   do l=1,3
    do p=1,3
     do q=1,3
      do r=1,3
       do s=1,3
        cnew(i,j,k,l) = cnew(i,j,k,l) + al(i,p)*al(j,q)*al(k,r)*al(l,s)*cold(p,q,r,s)
       end do
      end do
     end do
    end do
   end do
  end do
 end do
end do

! and convert to 6-index notation again
do i=1,6
 do j=1,6
  cout(i,j) = cnew(gamma(i,1),gamma(i,2),gamma(j,1),gamma(j,2))
 end do
end do
! That's it.

end subroutine TransFourthRankTensor
 

!--------------------------------------------------------------------------
!
! SUBROUTINE: laguer
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  Laguerre polynomial, based on Numerical Recipes routine of the same name
!
!> @details but adapted to the sextic equation case considered here.
!
!> @param a
!> @param m
!> @param x
!> @param eps
!> @param polish
! 
!> @date 1/5/99   MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   06/04/13 MDG 3.0 rewrite
!--------------------------------------------------------------------------
subroutine laguer(a,m,x,eps,polish)

use local

IMPLICIT NONE

complex(kind=dbl),INTENT(IN)	:: a(*)
integer(kind=irg),INTENT(IN)		:: m
complex(kind=dbl),INTENT(INOUT)	:: x
real(kind=sgl),INTENT(IN)		:: eps
logical,INTENT(IN)				:: polish

complex(kind=dbl)				:: dx,x1,b,d,f,g,h,sq,gp,gm,g2,zero
real(kind=sgl),parameter    		:: epss=6.E-8
integer(kind=irg),parameter 		:: maxit=1000
integer(kind=irg)           			:: iter,j
real(kind=dbl)              			:: dxold,errr,abx,cdx

      zero = cmplx(0.0,0.0,dbl)
      dxold = abs(x)
      
      do iter=1,maxit
        b=a(m+1)
        errr=abs(b)
        d=zero
        f=zero
        abx=abs(x)
        do j=m,1,-1
          f=x*f+d
          d=x*d+b
          b=x*b+a(j)
          errr=abs(b)+abx*errr
	end do
        errr=epss*errr
        if (abs(b).le.errr) then
          dx=zero
          return
        else
          g=d/b
          g2=g*g
          h=g2-2.D0*f/b
          sq=sqrt((m-1)*(m*h-g2))
          gp=g+sq
          gm=g-sq
          if (abs(gp).lt.abs(gm)) gp=gm
          dx=m/gp
        end if
        x1=x-dx
        if (x.eq.x1) return
        x=x1
        cdx=abs(dx)
        if ((iter.gt.6).and.(cdx.ge.dxold)) return
	dxold=cdx
        if (.NOT.polish) then
          if (abs(dx).le.eps*abs(x)) return
	end if
     end do
     pause 'too many iterations'
     return
end subroutine laguer

!--------------------------------------------------------------------------
!
! SUBROUTINE: zroots
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  Polynomial roots, based on Numerical Recipes routine of the same name
!
!> @details but adapted to the sextic equation case considered here.
!
!> @param a
!> @param roots
! 
!> @date 1/5/99   MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   06/04/13 MDG 3.0 rewrite
!--------------------------------------------------------------------------
subroutine zroots(a,roots)

use local

IMPLICIT NONE

complex(kind=dbl),INTENT(IN)           	:: a(6)
complex(kind=dbl),INTENT(OUT)		::roots(6)
integer(kind=irg)           				:: i,j,jj,m
real,parameter              				:: eps = 1.E-6
complex(kind=dbl)           			:: ad(7), x, b, c, czero
 
 m=6
 czero = cmplx(0.0,0.0,dbl)
 do j=1,7
  ad(j) = a(j)
 end do

! get the roots
 do j=6,1,-1
  x = czero
  call laguer(ad,j,x,eps,.FALSE.)
  if (abs(aimag(x)).le.2.*eps**2*abs(real(x))) x=cmplx(real(x),0.D0,dbl)
  roots(j)=x
  b=ad(j+1)
  do jj=j,1,-1
    c=ad(jj)
    ad(jj)=b
    b=x*b+c
  end do
 end do

! polish the roots
 do j=1,6
  call laguer(a,m,roots(j),eps,.TRUE.)
 end do

! and prepare to return to the calling routine
 outerloop: do j=2,6
  x=roots(j)
  innerloop: do i=j-1,1,-1
   if (real(roots(i)).le.real(x)) then
     roots(i+1) = x
     cycle outerloop
   end if
   roots(i+1)=roots(i)
  end do innerloop
  i=0
  roots(i+1)=x
 end do outerloop
 return
end subroutine

!------------------------
! these following three routines should be replaced by something else in the io.f90 file
! using a proper interface for all possible formats
!------------------------
subroutine PrintMatrix(s,a)

use local
use io

IMPLICIT NONE

real(kind=sgl)   :: a(3,3)
integer(kind=irg):: i,j
character(4)     :: s

write (*,"(A/)") s
do i=1,3
  write (*,"(3(F12.5,2x))") (a(i,j),j=1,3)
end do
write (*,"(/)")

end subroutine

subroutine PrintMatrixd(s,a)

use local
use io

IMPLICIT NONE

real(kind=dbl)   :: a(3,3)
integer(kind=irg):: i,j
character(4)     :: s

write (*,"(A/)") s
do i=1,3
  write (*,"(3(F12.5,2x))") (a(i,j),j=1,3)
end do
write (*,"(/)")

end subroutine

subroutine PrintMatrixcd(s,a)

use local
use io

IMPLICIT NONE

complex(kind=dbl)   :: a(3,3)
integer(kind=irg):: i,j
character(4)     :: s

write (*,"(A/)") s
do i=1,3
  write (*,*) (a(i,j),j=1,3)
end do
write (*,"(/)")

end subroutine
! 

!--------------------------------------------------------------------------
!
! SUBROUTINE: makedislocation
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  Compute the dismat displacement matrix for a given dislocation
!
!> @details This subroutine computes the matrix dismat that describes the displacement field
!> of a dislocation.  The routine needs the elastic moduli tensor, the transformation
!> matrix between the crystal and dislocation reference frames, and the dislocation
!> Burgers vector.  The routine computes the arrays dismat and pa, which should be used as follows:
!> 
!> R_k = 2.0*real([ sum_a=1^3 (dismat(k,a)*log(Z_a)) ]),
!> 
!> with Z_a = x_1 + pa(a)*x_2
!>
!>  [see CalcR subroutine for more information]
!>
!> We must also make sure that the x=0 plane of the defect reference frame contains the 
!> incident beam direction, to avoid getting stacking-fault fringes in the wrong plane...
!> Actual stacking faults are added in using a different module (stacking_fault.f90).
!>
!> @param inum
!> @param dinfo
!> @param DF_L column width
! 
!> @date   1/ 5/99 MDG 1.0 original
!> @date   5/19/01 MDG 2.0 f90 version
!> @date  11/27/01 MDG 2.1 added kind support
!> @date  06/04/13 MDG 3.0 rewrite
!> @date  10/30/13 MDG 3.1 debug of all rotation parts
!--------------------------------------------------------------------------
subroutine makedislocation(inum,dinfo,DF_L)

use local
use math
use constants
use foilmodule
use crystal
use crystalvars
use symmetry
use symmetryvars
use quaternions
use rotations

IMPLICIT NONE

integer(kind=irg),INTENT(IN)		:: inum
integer(kind=irg),INTENT(IN)		:: dinfo
real(kind=sgl),INTENT(IN)		:: DF_L

real(kind=dbl)          		:: zz,zang,zmin
real(kind=sgl)          		:: ec(6,6),lec(6,6)
real(kind=dbl)				:: a_dc(3,3),tmp(3),ex(3),ey(3)
real(kind=dbl)          		:: Bij(3,3),Hij(3,3)
complex(kind=dbl)       		:: a(0:6),b(0:6),c(0:6),d(0:6),e(0:6),ff(0:6),tt(5,0:6),s(0:6),roots(6), &
                                         zero,pasq(3),mat(3,3),aka(3,3),Lia(3,3),Mai(3,3),v(3),pas
integer(kind=irg)       		:: i,j,k,l,imin,ind(3),jnd(3)


! convert line direction and g-vector to the Cartesian crystal reference frame
call TransSpace(DL(inum)%u,DL(inum)%un,'d','c')
call TransSpace(DL(inum)%g,DL(inum)%gn,'r','c')
! normalize both vectors
call NormVec(DL(inum)%un,'c')
call NormVec(DL(inum)%gn,'c')

! first find the length of the dislocation line inside the foil along the
! dislocation z-axis, which is the line direction; also, compute the intersection
! points of the line with the top and bottom surfaces  (all components in [nm])
zang = CalcAngle(DL(inum)%u,dble(foil%F),'d')
zz = cos(zang)
if (abs(zz).gt.0.00001) then 
  DL(inum)%zu = 0.5*foil%z0/zz
else
  DL(inum)%zu = 100000.0         ! this is when the dislocation is nearly parallel to the foil
end if

! transform the line direction to the foil reference frame
tmp = quat_rotate_vector( foil%a_fc, dble(DL(inum)%un) ) / DF_L

if (dinfo.eq.1) then
  write (*,*) 'transformed line direction ', tmp, zang, zz
end if

! determine the top and bottom intersection coordinates 
if (zz.gt.0.0) then  ! u points to the top of the foil
    DL(inum)%top = (/ DL(inum)%id + tmp(1)*DL(inum)%zu, DL(inum)%jd + tmp(2)*DL(inum)%zu, 0.5D0*foil%z0 /)
    DL(inum)%bottom = (/ DL(inum)%id - tmp(1)*DL(inum)%zu, DL(inum)%jd - tmp(2)*DL(inum)%zu, -0.5D0*foil%z0 /)
else                       ! u points to the bottom of the foil
    DL(inum)%top = (/ DL(inum)%id - tmp(1)*DL(inum)%zu, DL(inum)%jd - tmp(2)*DL(inum)%zu, -0.5D0*foil%z0 /)
    DL(inum)%bottom = (/ DL(inum)%id + tmp(1)*DL(inum)%zu, DL(inum)%jd + tmp(2)*DL(inum)%zu, 0.5D0*foil%z0 /)
end if  

if (dinfo.eq.1) then
  write (*,*) DL(inum)%id,DL(inum)%jd
  write (*,*) 'dislocation top intersection at ',DL(inum)%top
  write (*,*) 'dislocation bottom intersection at ',DL(inum)%bottom
end if 

! a_dc (crystal to defect)  matrix corrected on 11/29/10 to put defect x-axis in the plane of u and B
if (dinfo.eq.1) then
  write (*,*) 'cartesian quantities'
  write (*,*) 'unit line direction = ',DL(inum)%un
  write (*,*) 'unit beam direction = ',foil%Bn
end if

! transform beam direction (currently in foil frame) to cartesian 
tmp = quat_rotate_vector(conjg(foil%a_fc), dble(foil%Bn))
!tmp = quat_rotate_vector(conjg(foil%a_fc), (/ 0.0D0, 0.0D0, -1.0D0/) )
call NormVec(tmp,'c')

! the defect z axis is the line direction and x is in the plane of u and B to avoid the intrinsic discontinuity (cut plane)
a_dc(3,1:3) = DL(inum)%un(1:3)
call CalcCross(dble(DL(inum)%un),tmp,ex,'c','c',0)
call NormVec(ex,'c')
a_dc(1,1:3) = ex(1:3)
call CalcCross(dble(DL(inum)%un),ex,ey,'c','c',0)
call NormVec(ey,'c')
a_dc(2,1:3) = ey(1:3)
DL(inum)%a_dc = om2qu(a_dc)

if (dinfo.eq.1) then
  call PrintMatrixd('a_dc',a_dc)
end if

! a_di (image to defect)
DL(inum)%a_di = quat_mult( DL(inum)%a_dc, conjg(foil%a_ic) )
DL(inum)%a_id = conjg(DL(inum)%a_di)

if (dinfo.eq.1) then
  call print_orientation(init_orientation(DL(inum)%a_di,'qu'),'om','a_di: ')
  call print_orientation(init_orientation(DL(inum)%a_id,'qu'),'om','a_id: ')
end if

! finally, get the foil to defect transformation (used in defect module)
DL(inum)%a_df = quat_mult( DL(inum)%a_di, conjg(foil%a_fi) )

! Burgers vector (in the defect reference frame !!!)
! first transform Burgers vector to crystal cartesian reference frame
call TransSpace(dble(DL(inum)%burg),tmp,'d','c')
! then convert this to the defect reference frame
DL(inum)%burgd(1:3) = quat_rotate_vector(DL(inum)%a_dc,dble(tmp))

if (dinfo.eq.1) then
  write (*,*) 'rotated burgers vector  = ', DL(inum)%burgd(1:3) 
end if

! transform the elastic moduli
lec = foil%elmo

! transform lec to defect reference frame
a_dc = qu2om(DL(inum)%a_dc)
call TransFourthRankTensor(a_dc,lec,ec)
if (dinfo.eq.1)  then 
  write (*,*) 'Elasticity tensor in defect reference frame'
  do i=1,6 
    write (*,"(6(F8.4,2x))") (ec(i,j),j=1,6)
  end do
  write (*,*) '----'
end if

! next, create the sextic polynomial
zero = cmplx(0.0,0.0,dbl)
a=zero; b=zero; c=zero; d=zero; e=zero; ff=zero
a(0:2) = (/ cmplx(ec(1,1),0.0,dbl), cmplx(ec(1,6)*2.0,0.0,dbl),     cmplx(ec(6,6),0.0,dbl) /)
b(0:2) = (/ cmplx(ec(6,6),0.0,dbl), cmplx(ec(2,6)*2.0,0.0,dbl),     cmplx(ec(2,2),0.0,dbl) /)
c(0:2) = (/ cmplx(ec(5,5),0.0,dbl), cmplx(ec(4,5)*2.0,0.0,dbl),     cmplx(ec(4,4),0.0,dbl) /)
d(0:2) = (/ cmplx(ec(5,6),0.0,dbl), cmplx(ec(4,6)+ec(2,5),0.0,dbl), cmplx(ec(2,4),0.0,dbl) /)
e(0:2) = (/ cmplx(ec(1,5),0.0,dbl), cmplx(ec(1,4)+ec(5,6),0.0,dbl), cmplx(ec(4,6),0.0,dbl) /)
ff(0:2) = (/ cmplx(ec(1,6),0.0,dbl), cmplx(ec(1,2)+ec(6,6),0.0,dbl), cmplx(ec(2,6),0.0,dbl) /)
tt = zero
s = zero

! matrix elements
do j=0,6 
 do i=0,j 
  tt(1,j) = tt(1,j) + a(j-i)*b(i)
  tt(2,j) = tt(2,j) + d(j-i)*e(i)
  tt(3,j) = tt(3,j) + a(j-i)*d(i)
  tt(4,j) = tt(4,j) + b(j-i)*e(i)
  tt(5,j) = tt(5,j) + c(j-i)*ff(i)
 end do
end do

! determinant leading to the sextic equation
do j=0,6 
 do i=0,j 
  s(j) = s(j) + tt(1,j-i)*c(i) + 2.0*tt(2,j-i)*ff(i) - tt(3,j-i)*d(i) - tt(4,j-i)*e(i) - tt(5,j-i)*ff(i)
 end do
end do

! get the complex root pairs
call zroots(s,roots)

! then, solve the equation for the vector A_k using the roots with positive imaginary part.
k=1
do j=1,6
  if (aimag(roots(j)).gt.0.0) then
    DL(inum)%pa(k) = roots(j)
    k=k+1
  end if
end do

! renumber them to avoid the symmetry degeneracy (see page 328 Head et al.)
v(1:3) = ec(5,5) + 2.0*DL(inum)%pa(1:3)*ec(4,5) + ec(4,4)*DL(inum)%pa(1:3)**2
zmin = 100.0
imin = 0

! where is the smallest value ?
do i=1,3
  if (abs(v(i)).lt.zmin) then
    imin=i
    zmin = abs(v(i))
  end if
end do

! is the 3rd one the smallest ? if not, then swap with the current 3rd one.
if (imin.ne.3) then
  pas = DL(inum)%pa(imin)
  DL(inum)%pa(imin)=DL(inum)%pa(3)
  DL(inum)%pa(3)=pas
end if
  pas = DL(inum)%pa(1)
  DL(inum)%pa(1)=DL(inum)%pa(2)
  DL(inum)%pa(2)=pas

! eliminate really small numbers
do i=1,3
  if (abs(aimag(DL(inum)%pa(i))).lt.1.0e-8)  DL(inum)%pa(i)=cmplx(real(DL(inum)%pa(i)),0.0,dbl)
  if (abs(real(DL(inum)%pa(i))).lt.1.0e-8)   DL(inum)%pa(i)=cmplx(0.0,aimag(DL(inum)%pa(i)),dbl)
end do
if (dinfo.eq.1) then
  write (*,*) ' sextic roots'
  do i=1,3
    write (*,*) DL(inum)%pa(i)
  end do
  write (*,*) '---'
end if

!  compute the A_ka vectors (see description on page 328 of Head et al.)
pasq = DL(inum)%pa**2
if (dinfo.eq.1) write (*,*) 'Aka vectors'
do k=1,3
  mat = zero
  mat(1,1) = ec(1,1)+2.D0*DL(inum)%pa(k)*ec(1,6)+ec(6,6)*pasq(k)
  mat(2,2) = ec(6,6)+2.D0*DL(inum)%pa(k)*ec(2,6)+ec(2,2)*pasq(k)
  mat(3,3) = ec(5,5)+2.D0*DL(inum)%pa(k)*ec(4,5)+ec(4,4)*pasq(k)
  mat(2,3) = ec(5,6)+DL(inum)%pa(k)*(ec(4,6)+ec(2,5))+ec(2,4)*pasq(k)
  mat(1,3) = ec(1,5)+DL(inum)%pa(k)*(ec(1,4)+ec(5,6))+ec(4,6)*pasq(k)
  mat(1,2) = ec(1,6)+DL(inum)%pa(k)*(ec(1,2)+ec(6,6))+ec(2,6)*pasq(k)
  if (k.eq.1) then
    aka(1,1) = mat(2,2)*mat(3,3)-mat(2,3)*mat(2,3)
    aka(1,2) = mat(1,3)*mat(2,3)-mat(1,2)*mat(3,3)
    aka(1,3) = mat(1,2)*mat(2,3)-mat(1,3)*mat(2,2)
  end if
  if (k.eq.2) then
    aka(2,1) = mat(1,3)*mat(2,3)-mat(1,2)*mat(3,3)
    aka(2,2) = mat(1,1)*mat(3,3)-mat(1,3)*mat(1,3)
    aka(2,3) = mat(1,3)*mat(1,2)-mat(1,1)*mat(2,3)
  end if
  if (k.eq.3) then
    aka(3,1) = mat(1,2)*mat(2,3)-mat(1,3)*mat(2,2)
    aka(3,2) = mat(1,3)*mat(1,2)-mat(1,1)*mat(2,3)
    aka(3,3) = mat(1,1)*mat(2,2)-mat(1,2)*mat(1,2)
  end if  
  if (dinfo.eq.1) write (*,*) k,(aka(k,j),j=1,3)
end do
aka = transpose(aka)

! next, create the L_ialpha matrix
ind = (/ 6, 2, 4 /)
jnd = (/ 1, 6, 5 /)
Lia = zero
do i=1,3 
 do j=1,3
  do k=1,3
   Lia(i,j) = Lia(i,j)+(ec(ind(i),jnd(k))+DL(inum)%pa(j)*ec(ind(i),ind(k)))*aka(k,j)
  end do
 end do
end do
if (dinfo.eq.1)  call PrintMatrixcd('Lia',Lia)

! and invert it
call cInvert(Lia,Mai)
if (dinfo.eq.1)  call PrintMatrixcd('Mai',Mai)

! compute Bij ( real matrix )
Bij = 0.D0
do i=1,3
 do j=1,3
  do k=1,3
   Bij(i,j) = Bij(i,j) - aimag(aka(i,k))*real(Mai(k,j)) - real(aka(i,k))*aimag(Mai(k,j))
  end do
 end do
end do
if (dinfo.eq.1)  call PrintMatrixd('Bij',Bij)

! and invert to get Hij
call mInvert(Bij,Hij,.FALSE.)
if (dinfo.eq.1) call PrintMatrixd('Hij',Hij)

! compute matrix (this is what actually gets to be used for the 
! displacement field); needs to know the Burgers vector.
DL(inum)%dismat = zero
do k=1,3
 do l=1,3
   do i=1,3
    do j=1,3
     DL(inum)%dismat(k,l) = DL(inum)%dismat(k,l) + DL(inum)%burgd(i)*Hij(j,i)*Mai(l,j)*aka(k,l)
    end do
   end do
 end do
end do

! scale by 1/4pi
DL(inum)%dismat = DL(inum)%dismat*0.25D0/cPi
if (dinfo.eq.1)  call PrintMatrixcd('dismat',DL(inum)%dismat)

! and return to calling routine
end subroutine makedislocation



!--------------------------------------------------------------------------
!
! SUBROUTINE: read_dislocation_data
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  read dislocation namelist files
!
!> @param dislname name of dislocation namelist file (string array)
!> @param numdisl number of dislocations
!> @param numsf number of stacking faults
!> @param DF_npix number of x-pixels
!> @param DF_npiy number of y-pixels
!> @param DF_gf 
!> @param L 
!> @param dinfo logical to trigger verbose output
! 
!> @date    1/5/99   MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   06/04/13 MDG 3.0 rewrite
!--------------------------------------------------------------------------
subroutine read_dislocation_data(dislname,numdisl,numsf,DF_npix,DF_npiy,DF_gf,L,dinfo)

use local
use io
use files

IMPLICIT NONE

character(fnlen),INTENT(IN)	 	:: dislname(3*maxdefects)
integer(kind=irg),INTENT(IN)		:: numdisl, numsf, DF_npix, DF_npiy, dinfo
real(kind=sgl),INTENT(IN)		:: DF_gf(3), L

integer(kind=irg) 			:: i
real(kind=dbl) 				:: id,jd,u(3),bv(3),zfrac, poisson

namelist / dislocationdata / id, jd, u, bv, zfrac, poisson


! zfrac goes between -0.5 and +0.5, with -0.5 being the top surface and +0.5 the bottom
! this only really matters for dislocations that are parallel to the foil surfaces

! allocate the memory for the dislocation parameters
  allocate(DL(numdisl+2*numsf))

if (dinfo.eq.1) then
  do i=1,numdisl
    write (*,*) i,'->',trim(dislname(i)),'<-'
  end do
end if

! these are just the individual dislocations; the ones that belong to 
! stacking faults are handled separately
   do i=1,numdisl
    zfrac = 0.0   ! unless the dislocation line is parallel to the foil surface, zfrac is always 0.0 
    
    mess = 'opening '//trim(dislname(i)); call Message("(/A)")
    OPEN(UNIT=dataunit,FILE=trim(dislname(i)),DELIM='APOSTROPHE')
    READ(UNIT=dataunit,NML=dislocationdata)
    CLOSE(UNIT=dataunit)
    
! center of dislocation inside the foil is transformed to foil coordinates [nm] with DL(i)%kd=0 (center of foil) [verified 4/23/11]
! the point (0,0) is at the center of the image ... hence the factor of 0.5
    DL(i)%id = id * 0.5 * float(DF_npix) ! * L   scaling (zooming) is done later in the image reference frame...
    DL(i)%jd = jd * 0.5 * float(DF_npiy) ! * L
    DL(i)%u = u
    DL(i)%burg = bv
    DL(i)%g = DF_gf
    DL(i)%zfrac = zfrac - 0.5
     
! and pre-compute the dislocation displacement field parameters
       call makedislocation(i,dinfo, L)
  end do
  
end subroutine read_dislocation_data


end module dislocation
