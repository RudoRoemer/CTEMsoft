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
! EMsoft:foilmodule.f90
!--------------------------------------------------------------------------
!
! MODULE: foilmodule
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief everything that has to do with the sample foil
! 
!> @date 1/5/99   MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   06/03/13 MDG 3.0 updated IO, introduced quaternions instead of rotation matrices
!> @data   11/14/13 MDG 3.1 changed F and q to integer indices from float
!> @date   06/10/14 MDG 4.0 removed global variable foil
!--------------------------------------------------------------------------

module foilmodule

use local
use typedefs
use quaternions

! type (foiltype)        :: foil

contains

!--------------------------------------------------------------------------
!
! SUBROUTINE: initialize_foil_geometry
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  Initializes the foil geometry
!
!> @details This new implementation uses quaternions for all rotations. 
! 
!> @param cell unit cell pointer
!> @param dinfo
! 
!> @date  1/ 5/99 MDG 1.0 original
!> @date  1/11/10 MDG 2.0 rewrite of beam direction part
!> @date  3/28/11 MDG 2.1 code verified
!> @date  4/23/11 MDG 2.2 redefined origin to be at center of image
!> @date  6/03/13 MDG 3.0 replaced rotation matrices by quaternions throughout
!> @date 10/30/13 MDG 3.1 complete debug of quaternion and rotation implementation 
!> @date 06/09/14 MDG 4.0 added cell and foil as argument
!--------------------------------------------------------------------------
subroutine initialize_foil_geometry(cell,foil,dinfo)

use local
use typedefs
use math
use constants
use crystal
use symmetry
use io
use rotations

IMPLICIT NONE

type(unitcell),pointer          :: cell
type(foiltype),INTENT(INOUT)   :: foil
integer(kind=sgl),INTENT(IN)    :: dinfo

real(kind=dbl)                  :: ey(3), ex(3), tt, dx, dy 
real(kind=sgl)                  :: io_real(3)
real(kind=dbl)                  :: cp,sp,cs,ss,cr,sr, ca, sa, a_fc(3,3)
integer(kind=irg)               :: i,j
type(orientationtyped)          :: ot
character(10)                   :: pret

! determine the foil-to-microscope transformations [verified on 4/23/11, converted to quaternions on 6/4/13, 
! verified 10/30/13, and again on 11/11/13 after some changes elsewhere]
  if (foil%alR.eq.0.D0) then 
! the double tilt holder transformation a_fm; note quaternions, hence we need the half-angles !
! a_fm transforms a vector v FROM the microscope reference frame To the foil reference frame
! using the quat_rotate_vector routine.
    cp = dcos(foil%alP*0.5D0)
    sp = dsin(foil%alP*0.5D0)
    ca = dcos(foil%alP)
    sa = dsin(foil%alP)
    cs = dcos(foil%alS*0.5D0)
    ss = dsin(foil%alS*0.5D0)
    foil%a_fm = conjg( quat_mult( (/ cs, 0.D0, ss*ca, ss*sa /), (/ cp, sp, 0.D0, 0.D0 /) ) )
  else
! the rotation tilt holder transformation a_fm [verified on 4/23/11, converted to quaternions on 6/4/13,
! and again on 11/11/13 after changes elsewhere]
    cp = dcos(foil%alP*0.5D0)
    sp = dsin(foil%alP*0.5D0)
    cr = dcos(foil%alR*0.5D0)
    sr = dsin(foil%alR*0.5D0)
    ca = dcos(foil%alP)
    sa = dsin(foil%alP)
    foil%a_fm = conjg( quat_mult( (/ cr,0.D0, -sr*sa, sr*ca /), (/ cp, sp,0.D0,0.D0  /) ) )
  end if  
  if (dinfo.eq.1) then
    pret = 'a_fm: '
    ot = init_orientation(foil%a_fm,'qu')
    call print_orientation(ot,'om',pret) 
  end if

! a_mi (image to microscope) apart from a scale factor, these two are identical 
! The EM book uses a beta rotation angle between the image and the microscope,
! but that is really not necessary because we already fix the image with respect to
! the microscope by defining q (the horizontal image direction) to point to the 
! airlock. [verified 4/23/11, converted to quaternions on 6/4/13]
! So we'll keep this transformation equal to the identity at all times.
  foil%a_mi = (/ 1.D0,0.D0,0.D0,0.D0 /)   ! identity quaternion 
  if (dinfo.eq.1) then
    pret = 'a_mi: '
    ot = init_orientation(foil%a_mi,'qu')
    call print_orientation(ot,'om',pret) 
  end if
  
! This allows us to get the beam direction, since we know the foil normal and tilt angles
! The beam direction is the inverse transform of the microscope e_z-axis to the foil reference frame [verified 11/12/13]
  foil%B = quat_Lp( conjg(foil%a_fm), (/ 0.D0,0.D0,-1.D0 /) )
  foil%Bn = foil%B
  call NormVec(cell,foil%Bn,'c')
  if (dinfo.eq.1) then
    io_real(1:3) = foil%B(1:3)
    call WriteValue('  Beam direction (foil reference frame) = ',io_real,3,"('[',F12.5,',',F12.5,',',F12.5,']')")
  end if

! transform both the foil normal and the q-vector to the crystal cartesian reference frame (eq. 8.8) [verified 4/23/11,
! and again on 11/12/13 afterchanges elsewhere]
  call TransSpace(cell,foil%F,foil%Fn,'d','c')
  call TransSpace(cell,foil%q,foil%qn,'r','c')
  call NormVec(cell,foil%Fn,'c')
  call NormVec(cell,foil%qn,'c')
! a_fc (crystal to foil)  
  a_fc(3,1:3) = foil%Fn(1:3)
  a_fc(1,1:3) = foil%qn(1:3)
  call CalcCross(cell,foil%Fn,foil%qn,ey,'c','c',0)
  call NormVec(cell,ey,'c')
  a_fc(2,1:3) = ey(1:3)
  foil%a_fc = om2qu(a_fc)
  if (dinfo.eq.1) then
    pret = 'a_fc: '
    ot = init_orientation(foil%a_fc,'qu')
    call print_orientation(ot,'om',pret) 
  end if
  
! a_mc (crystal to microscope)
  foil%a_mc = quat_mult( conjg(foil%a_fm), foil%a_fc )
  if (dinfo.eq.1) then
    pret = 'a_mc: '
    ot = init_orientation(foil%a_mc,'qu')
    call print_orientation(ot,'om',pret) 
  end if
  
! a_ic (crystal to image)
  foil%a_ic = quat_mult( conjg(foil%a_mi) , foil%a_mc)
  if (dinfo.eq.1) then
    pret = 'a_ic: '
    ot = init_orientation(foil%a_ic,'qu')
    call print_orientation(ot,'om',pret) 
  end if
  
! a_fi (image to foil)
  foil%a_fi = quat_mult( foil%a_fc , conjg(foil%a_ic) )
  if (dinfo.eq.1) then
    pret = 'a_fi: '
    ot = init_orientation(foil%a_fi,'qu')
    call print_orientation(ot,'om',pret) 
  end if
  
! express the beam direction in the Bravais reference frame [verified 4/23/11, and again on 11/12/13
! after changes elsewhere]
  ex = quat_Lp( conjg(foil%a_fc), dble(foil%Bn) ) 
  call TransSpace(cell,ex,ey,'c','d')
  call NormVec(cell,ey,'c')
  if (dinfo.eq.1) then
    io_real(1:3) = ey(1:3)
    call WriteValue('  Beam direction (crystal reference frame) = ', io_real, 3, "('[',F12.5,',',F12.5,',',F12.5,']'/)")
  end if
  
! define the foil shape (for now as an elliptic paraboloid z = brx * (x-xc)^2 + bry * (y-yc)^2) 
if (.not.allocated(foil%sg)) allocate(foil%sg(foil%npix,foil%npiy))
! if the foil is not bent, then we set this array to zero, otherwise we compute the elliptical paraboloid
if ((foil%brx.eq.0.0).and.(foil%bry.eq.0.0)) then
  if (dinfo.eq.1) then
    call Message(' Initializing a flat foil ', frm = "(A)")
  end if
  foil%sg = 0.0
else
  dx = foil%npix*0.5
  dy = foil%npiy*0.5
  do i=1,foil%npix
   tt = foil%brx * (float(i)-dx-foil%cpx)**2
    do j=1,foil%npiy
! initialize the foil shape function; we assume that the center of the elliptic paraboloid is at location (cpx,cpy)
! presumably, this surface could also be a saddle point if the brx and bry values have opposite sign ...
      foil%sg(i,j) = tt + foil%bry * (float(j)-dy-foil%cpy)**2+ 2.0*foil%brxy * (float(j)-dy-foil%cpy)*(float(i)-dx-foil%cpx)
    end do
  end do
  if (dinfo.eq.1) then
    call Message(' Initializing a bent foil ', frm = "(A)")
    io_real(1)=minval(foil%sg); io_real(2)=maxval(foil%sg); 
    call WriteValue('Range of local excitation error deviations : ', io_real, 2, "(F10.6,',',F10.6/)")
  end if
end if


end subroutine initialize_foil_geometry

!--------------------------------------------------------------------------
!
! SUBROUTINE: read_foil_data
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief  Reads the foil geometry data
!
!> @param cell unit cell pointer
!> @param foil foil structure
!> @param foilnmlfile name of foil namelist file
!> @param npix number of x image pixels
!> @param npiy number of y image pixels
!> @param L pixel size for column approximation
!> @param dinfo flag to print information
! 
!> @date 010/5/11 MDG 1.0 original
!> @date 06/04/13 MDG 2.0 general rewrite
!> @date 06/10/14 MDG 3.0 added cell and foil as arguments
!--------------------------------------------------------------------------
subroutine read_foil_data(cell,foil,foilnmlfile,npix,npiy,L,dinfo)

use crystal
use io
use files
use rotations

IMPLICIT NONE

type(unitcell),pointer          :: cell
type(foiltype),INTENT(INOUT)   :: foil
character(fnlen),INTENT(IN)     :: foilnmlfile
integer(kind=irg),INTENT(IN)    :: npix, npiy, dinfo
real(kind=sgl),INTENT(IN)       :: L

integer(kind=irg)               :: i,j
real(kind=sgl)                  :: foilF(3), foilq(3), foilalP, foilalS, foilalR, foilbeP, &
                                foilz0, foilelmo(6,6),brx,bry,brxy,cpx,cpy,x, io_real(1)
real(kind=dbl)                  :: amat(3,3)
real(kind=sgl), parameter       :: cPi=3.1415926536

namelist / foildata / foilF, foilq, foilalP, foilalS, foilalR, foilbeP, foilz0, foilelmo, brx, bry, brxy, cpx, cpy

! set the default values
 foilelmo = 0.0                         ! elastic moduli
 foilF = (/ 0.0,0.0,1.0 /)              ! foil normal in direct space Bravais reference frame 
 foilq = (/ 1.0,0.0,0.0 /)              ! reciprocal space vector along primary tilt axis towards airlock
 foilalP = 0.0                          ! primary tilt angle in degrees
 foilalS = 0.0                          ! secondary tilt angle (for double tilt holder)
 foilalR = 0.0                          ! secondary tilt angle (for rotation tilt holder)
 foilbeP = 0.0                          ! angle of primary tilt axis w.r.t. image bottom edge
 foilz0 = 100.0                         ! foil thickness in nm
 brx = 0.0                              ! parameters to describe the foil shape as a quadratic surface 
 bry = 0.0
 brxy = 0.0
 cpx = 0.0                              ! center of the foil quadratic surface within [-1,1] range in pixel coordinates
 cpy = 0.0
 
 if (dinfo.eq.1) then 
  call Message('opening '//foilnmlfile, frm = "(/A/)")
 end if

 OPEN(UNIT=dataunit,FILE=trim(foilnmlfile),DELIM='APOSTROPHE')
 READ(UNIT=dataunit,NML=foildata)
 CLOSE(UNIT=dataunit)

! verify that the foil normal (in real space) and q (in reciprocal space) are orthogonal
x = CalcDot(cell,foilF,foilq,'c')
if (abs(x).gt.0.005) then
  call Message('Foil normal must be orthogonal to q', frm = "(A)")
  stop
end if

! and assign these values to the appropriate slots in foil%   [verified 4/23/11]
foil%F = foilF                          ! foil normal (w.r.t. Bravais lattice)
foil%q = foilq                          ! reciprocal vector pointing to airlock and normal to foilB
foil%alP = foilalP*cPi/180.0            ! convert the tilt and rotation angles to radians
foil%alS = foilalS*cPi/180.0
foil%alR = foilalR*cPi/180.0
foil%z0 = foilz0                        ! foil thickness in nanometers
foil%npix = npix                        ! image size (this duplicates some values, but it's easier this way)
foil%npiy = npiy
foil%brx = brx                          ! foil shape parameter for the x-direction
foil%bry = bry                          ! same for the y-direction
foil%brxy = brxy                        ! and a cross term
foil%cpx = cpx * float(npix) * 0.5 * L  ! we'll define the foil shape center w.r.t. to the center of the image in [nm] coordinates
foil%cpy = cpy * float(npiy) * 0.5 * L  ! 

! copy the elastic moduli to the proper array, taking into account the array symmetry
! (the units of the moduli do not matter, as long as they are all the same; only the upper
! triangle needs to be defined in the namelist file)
do i=1,6
  do j=i,6
    foil%elmo(i,j) = foilelmo(i,j)
    if (i.ne.j) foil%elmo(j,i) = foilelmo(i,j)
  end do
end do

! initialize a bunch of foil related quantities, using quaternions for all rotations
call initialize_foil_geometry(cell,foil,dinfo)

! compute the projected thickness
amat = qu2om(foil%a_fm)
foil%zb = foil%z0/amat(3,3)
if (dinfo.eq.1) then
  io_real(1)=foil%z0
  call WriteValue('Nominal foil thickness = ', io_real, 1, "(F8.3)")
  io_real(1)=foil%zb
  call WriteValue('Effective foil thickness = ', io_real, 1, "(F8.3/)")
end if

end subroutine read_foil_data



end module foilmodule

