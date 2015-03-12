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
! CTEMsoft2013:defectmodule.f90
!--------------------------------------------------------------------------
!
! MODULE: defectmodule
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Provides a routine to compute the displacement vector for an array of defects.
! 
!> @date   04/29/11 MDG 1.0 original
!> @date   06/04/13 MDG 2.0 rewrite + quaternions instead of rotations
!> @date   11/13/13 MDG 2.1 fixed error with coordinate transformations (after very long bug search!)
!--------------------------------------------------------------------------
module defectmodule

use local
use typedefs

! we have to get rid of ALL these global variables !!!
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

! complex(kind=dbl),allocatable    	:: DF_Sarray(:,:,:),theta(:),DF_Svoid(:,:)
!real(kind=sgl),allocatable       	:: images(:,:,:),DF_inclusion(:,:)
!
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------


contains


!--------------------------------------------------------------------------
!
! SUBROUTINE: CalcR
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief returns the total displacement vector for each slice in a column
!
!> @details Note that the end result MUST be expressed in the cartesian reference frame !
!
!> @param defects defect structure
!> @param cell unit cell pointer
!> @param foil foil structure
!> @param i integer x coordinate 
!> @param j integer y coordinate
!
!> @note This entire routine was thoroughly verified after the quaternion conversion !
!>
!> General comment for those who wish to add other defects...
!>
!> the general procedure to implement a defect displacement field is as follows:
!> - if the defect has its own reference frame, then transform (xpos,ypos,zpos) to
!>   that frame (see the dislocation section below for an example), and then do
!>   the computation of the displacement vector and express it in the cartesian frame.
!>
!> - if the defect uses the foil reference frame (e.g., voids, inclusions), then use tmpf
!>  as the current position vector.

!> @date  10/20/98 MDG 1.0 original
!> @date   5/22/01 MDG 2.0 f90
!> @date  11/27/01 MDG 2.1 added kind support
!> @date  03/26/13 MDG 3.0 updated IO
!> @date  10/30/13 MDG 3.1 debug of coordinate rotations
!> @date  11/13/13 MDG 3.2 finally, the bug has been found!  
!> @date  06/09/14 MDG 4.0 introduced defects argument and simplified routine
!> @date  06/10/14 MDG 4.1 added foil argument
!--------------------------------------------------------------------------
subroutine CalcR(defects,cell,foil,i,j)

use local
use constants
use crystal
use dislocation
use YSHModule
use foilmodule
use void
use stacking_fault
use inclusion
use quaternions
use rotations

IMPLICIT NONE

type(defecttype),INTENT(INOUT)        :: defects
type(unitcell),pointer	               :: cell
type(foiltype),INTENT(INOUT)          :: foil
integer(kind=irg),INTENT(IN)    	:: i,j

integer(kind=irg)			:: k, islice, ii
real(kind=dbl)        			:: dis,xpos,ypos,zpos,sumR(3),thick,tmp(3),tmp2(3), &
					   tmpf(3),u(3),zaamp,zaphase,zar,zai,zr(3),zi(3), &
                                 	   zt,fx,fy,fz,a_fm(3,3)    !,&
!                                nu,x,y,z,zn,t,pre,r1,r2,r3,th,rn 
                         			 
complex(kind=dbl)     			:: za(3)
complex(kind=sgl)     			:: zero
logical               			:: isvoid

! scale the image coordinates with respect to the origin at the center of the image
 xpos = float(i-defects%DF_npix/2)*defects%DF_L
 ypos = float(j-defects%DF_npiy/2)*defects%DF_L

! determine the starting point of the z-integration for the tilted foil
! this depends on the foil normal components which give the equation
! of the top foil plane as F . r = z0/2, from which we get zt...
 a_fm = qu2om(foil%a_fm)
 fx = a_fm(3,1)
 fy = a_fm(3,2)
 fz = a_fm(3,3) 
 zt = foil%zb*0.5 - (fx*xpos + fy*ypos)/fz
 
! initialize some other variables
 thick = foil%zb
 zero = cmplx(0.0,0.0)

! loop over all slices (this is the main loop)
 sliceloop: do islice = 1,defects%DF_nums 

! zpos is the position down the column, starting at zt (in image coordinates)
    zpos = zt - float(islice)*defects%DF_slice
    
! set the displacements to zero
    sumR = 0.0
        
! set the position in the foil reference frame
    tmpf = quat_LPstar( foil%a_fi, dble( (/ xpos, ypos, zpos /)) )

!------------
!----VOIDS---
!------------
! voids are easy to deal with; we simply return -10000 for each point tmpf that lies inside
! one of the voids; the calling routine then knows to use the void scattering matrix.
if (defects%numvoids.ne.0) then 
! are we inside a void ?
    isvoid = .FALSE.
    voidloop: do ii=1,defects%numvoids
! subtract the void position from the current slice position to get the relative position vector
     tmp = tmpf -  (/ defects%voids(ii)%xpos, defects%voids(ii)%ypos, defects%voids(ii)%zpos /)
     dis = CalcLength(cell,tmp,'c')
     if (dis.lt.defects%voids(ii)%radius) then ! inside void
       isvoid = .TRUE.
       exit voidloop
     end if
    end do voidloop
! skip the rest of the computation for this slice if we are inside a void
    if (isvoid.eqv..TRUE.) then 
      defects%DF_R(islice,1) = -10000.0
      cycle sliceloop
    end if
 end if 

! ok, if we get here, then we're not inside a void...

!------------------
!----CURVED FOIL---
!------------------
! first we take the foil shape into account using equations (8.28) and (8.29)
 sumR = sumR + float(islice)*defects%DF_slice*foil%sg(i,j)*defects%DF_gstar

!-----------------
!--DISLOCATIONS--
!-----------------
! let's put a few dislocations in ... (see section 8.4.2)
do ii=1,defects%numdisl
! compute the difference vector between the current (xpos,ypos,zpos) in the foil reference frame
! and the defect center coordinate
  tmp2 =  tmpf - dble((/ defects%DF_L*defects%DL(ii)%id, defects%DF_L*defects%DL(ii)%jd, defects%DL(ii)%zfrac*foil%z0 /))

! then convert the difference vector to the defect reference frame for this dislocation (we will only need the x and y coordinates)
  tmp = quat_LPstar( defects%DL(ii)%a_df, tmp2 ) 


! compute x1 + p_alpha x2  (eq. 8.38)
  za(1:3) = tmp(1) + defects%DL(ii)%pa(1:3)*tmp(2)
! compute the displacement vector u (eq. 8.38) [this expands the log of a complex number and takes the real part only,
! taking proper care of the branch cut] 
   if (tmp(1).gt.0.0) then
   do k=1,3
    zar =  real(za(k))
    zai = aimag(za(k))
    zaamp = abs(za(k))
    zaphase = abs(zai/zar)
    zr(k) = log(zaamp)
    zi(k) = atan(zaphase)
    if (zar.le.0.0) then
      if (zai.lt.0.0) zi(k) = -cPi+zi(k)
      if (zai.eq.0.0) zi(k) = cPi
      if (zai.gt.0.0) zi(k) = cPi-zi(k)
    else
      if (zai.lt.0.0) zi(k) = -zi(k)
    end if
   end do
  else
   do k=1,3
    zar =  real(za(k))
    zai = aimag(za(k))
    zaamp = abs(za(k))
    zaphase = abs(zai/zar)
    zr(k) = log(zaamp)
    zi(k) = atan(zaphase)
    if (zar.le.0.0) then
      if (zai.gt.0.0) zi(k) = cPi-zi(k)
      if (zai.eq.0.0) zi(k) = cPi
      if (zai.lt.0.0) zi(k) = cPi+zi(k)
    else
      if (zai.lt.0.0) zi(k) = 2.0*cPi-zi(k)
      if (zai.eq.0.0) zi(k) = 0.0
    end if
   end do  
  end if
  u = 2.0*real(matmul(defects%DL(ii)%dismat,cmplx(zr,zi)))
! transform displacement vector u to the Cartesian crystal reference frame
  u = quat_LPstar( conjg(defects%DL(ii)%a_dc), dble(u) )  
  sumR = sumR + u 
end do

!-------------------------------------
!--SURFACE INTERSECTING DISLOCATIONS--
!-------------------------------------
! this part is mostly used for ECCI-type image simulations, not for CTEM or STEM,
! although it could probably be used there as well; we would need to extend it 
! to incorporate both top and bottom foil surfaces

! do we have any dislocations with surface relaxations ?  YSH model
if (defects%numYdisl.gt.0) then 
   do ii=1,defects%numYdisl
! first, figure out what the coordinates are in the YSH reference frame for this dislocation ... 
! translate to the defect origin
     tmp =  tmpf -  (/ defects%DF_L*defects%YD(ii)%id, defects%DF_L*defects%YD(ii)%jd, foil%z0*0.5 /)

! rotate into the defect reference frame
     tmp = quat_LPstar( defects%YD(ii)%a_di, tmp )   

! compute the displacement vector
!     u = sngl(YSHDisp(dble(tmp(2)),-dble(tmp(1)),dble(tmp(3)),ii))
     u = sngl(YSHDisp(defects,dble(tmp(1)),dble(tmp(2)),dble(tmp(3)),ii))

! and rotate back to the image reference frame
     u = quat_LPstar( defects%YD(ii)%a_id, u )
     u = quat_LPstar( conjg(foil%a_ic), u ) 

! that should do it !
     sumR = sumR + u
   end do
end if

!--------------------
!--STACKING FAULTS--
!--------------------
! stacking faults (this is easy because we've already done all the work in the stacking_fault module)
! all we need is the z-value at which the stacking fault plane is crossed in this particular image
! column; from that point on, we simply add the leading partial Burgers vector to the total displacement.
do ii=1,defects%numsf
  if ((zpos.lt.defects%SF(ii)%zpos(i,j)).and.(defects%SF(ii)%zpos(i,j).ne.-10000.0)) then 
    sumR = sumR + defects%SF(ii)%lpbc
  end if
end do

!--------------------
!--LARGE INCLUSIONS--  currently commented out 
!--------------------
! Mader's expression for the displacement field of a large inclusion
!   if (0.eq.1.) then 
!    nu = 0.25
!    ce = 0.005
!    rn = 25.0*DF_L
!    x = (float(i-DF_npix/2)-0.5)*DF_L
!    y = (float(j-DF_npiy/2)-0.5)*DF_L
!    z = float(k)*DF_slice
!    zn = 100.5*DF_slice
!    t = DF_slice * DF_nums
!    pre = (1.0+nu)/(3.0*(1.0-nu))*ce*rn**3
! 
!    r1 = sqrt(x**2+y**2+(z-zn)**2)
!    r2 = sqrt(x**2+y**2+(z+zn)**2)
!    r3 = sqrt(x**2+y**2+(2.0*t-z-zn)**2)
! 
!    if (((r1.eq.0.0).or.(r2.eq.0.0)).or.(r3.eq.0.0)) then
!      return
!    else
!     dis = (1.0/r1**3+(3.0-4.0*nu)/r2**3-6.0*z*(z+zn)/r2**5+(3.0-4.0*nu)/r3**3-6.0*(t-z)*(2.0*t-z-zn)/r3**5)
!     rx = x*dis
!     ry = y*dis
!     rz = (z-zn)/r1**3-(3.0-4.0*nu)*((z+zn)/r2**3+(2.0*t-z-zn)/r3**3)-6.0*z*(z+zn)**2/r2**5 + &
!          2.0*z/r2**3+6.0*(t-z)*(2.0*t-z-zn)**2/r3**5-2.0*(t-z)/r3**3
! 
!     sumR = pre*(/ rx, ry, rz /)
!     return
!    end if
!   end if

!--------------------
!--SMALL INCLUSIONS--
!--------------------
! then the coherent precipitates, using the model in section 8.4.1
  if (defects%numinc.gt.0) then
   do ii=1,defects%numinc
! subtract the inclusion position from the current slice position to get the relative position vector
     tmp = tmpf - (/ defects%inclusions(ii)%xpos, defects%inclusions(ii)%ypos, defects%inclusions(ii)%zpos /)
     dis = CalcLength(cell,tmp,'c')
     if (dis.ge.defects%inclusions(ii)%radius) then ! outside particle
       tmp = tmp*(defects%inclusions(ii)%radius/dis)**3
     end if
     sumR = sumR + defects%inclusions(ii)%C*tmp
   end do
  end if

! TO BE IMPLEMENTED FOR RICHARD LESAR'S Discrete Dislocation Dynamics ! 
! finally any displacement fields defined by the user routine UserDisp
! sumR = sumR + UserDisp()


! TO BE IMPLEMENTED FOR YUNZHI WANG's Dislocation Simulations ! 
! finally any displacement fields defined by the user routine UserDisp
! sumR = sumR + UserDisp()        

   defects%DF_R(islice,1:3) = sumR(1:3)
  end do sliceloop ! main loop over the slices

end subroutine CalcR


end module defectmodule
