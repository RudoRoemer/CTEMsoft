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
! EMsoft:EMdymod.f90
!--------------------------------------------------------------------------
!
! MODULE: EMdymod
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief routines that can be called by external code
!
!> @date  10/16/15 MDG 1.0 original
!--------------------------------------------------------------------------
module EMdymod


contains

!--------------------------------------------------------------------------
!
! SUBROUTINE:SingleEBSDPattern
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief This function can be called as a standalone function to compute an EBSD pattern
!
!> @etails The main purpose of this routine and its accompanying wrapper routine is to
!> provide a way for an external program to compute a channeling pattern.  The idea is that 
!> all the necessary arrays and variables are passed in by reference as arguments, without
!> the need for the routine to fetch any other data from files etc...  The initial goal is
!> to have a function that can be called with the CALL_EXTERNAL mechanism in IDL, but 
!> in the long run this will also be the approach for calling the routine from C/C++, which
!> is an essential part of integration with DREAM.3D.  This routine is a simplified version
!> of the core of the EMECP program. 
!>
!> This routine will first compute the detector arrays rgx etc. if necessary, and then perform
!> the usual interpolation from the square Lambert projection. The pattern will be a basic pattern,
!> without any intensity scaling or binning etc; the calling program should take care of those 
!> operations.
!
!> @param nipar number of entries in ipar
!> @param nfpar number of entries in fpar
!> @param n1, n2 dimension parameters for accum_e
!> @param m1 dimension parameter for mLPNH and mLPSH
!> @param o1,o2 dimension parameters for EBSDpattern
!> @param ipar array with integer input parameters
!> @param fpar array with float input parameters
!> @param accum_e array with Monte Carlo histogram
!> @param mLPNH Northern hemisphere master pattern
!> @param mLPSH Southern hemisphere master pattern
!> @param EBSDpattern output array
!
!> @date 10/16/15 SS 1.0 original
!--------------------------------------------------------------------------
subroutine SingleEBSDPattern(nipar, nfpar, n1, n2, m1, o1, o2, ipar, fpar, accum_e, mLPNH, mLPSH, EBSDpattern)

! the input parameters are all part of a ipar and fpar input arrays instead of the usual namelist structures.
! The following is the mapping:
!
! ipar(1) = 1 if rgx, rgy, rgz detector arrays need to be computed, 0 if not (arrays will have save status)
! ipar(2) = enl%numsx
! ipar(3) = enl%numsy
! ipar(4) = enl%numEbins
! ipar(5) = enl%nsx
! ipar(6) = enl%nsy
! ipar(7) = nx = (enl%nsx-1)/2
! ipar(8) = enl%npx

! fpar(1) = enl%xpc
! fpar(2) = enl%ypc
! fpar(3) = enl%delta
! fpar(4) = enl%MCsig
! fpar(5) = enl%omega
! fpar(6) = enl%thetac
! fpar(7) = enl%L
! fpar(8) = enl%beamcurrent
! fpar(9) = enl%dwelltime
! fpar(10-13) = quaternion for requested Euler angles

use local
use constants
use Lambert
use quaternions

IMPLICIT NONE

integer(8),INTENT(IN)                   :: nipar
integer(8),INTENT(IN)                   :: nfpar
integer(8),INTENT(IN)                   :: n1
integer(8),INTENT(IN)                   :: n2
integer(8),INTENT(IN)                   :: ipar(nipar)
integer(8),INTENT(IN)                   :: m1
integer(8),INTENT(IN)                   :: o1
integer(8),INTENT(IN)                   :: o2
real(kind=sgl),INTENT(IN)               :: fpar(nfpar)
real(kind=sgl),INTENT(IN)               :: accum_e(n1,-n2:n2,-n2:n2)
real(kind=sgl),INTENT(IN)               :: mLPNH(-m1:m1, -m1:m1, n1)
real(kind=sgl),INTENT(IN)               :: mLPSH(-m1:m1, -m1:m1, n1)
real(kind=sgl),INTENT(OUT)              :: EBSDpattern(o1, o2)

! variables that must be saved for the next time this function is called
real(kind=irg),allocatable,save         :: accum_e_detector(:,:,:)
real(kind=sgl),allocatable,save         :: rgx(:,:), rgy(:,:), rgz(:,:)
real(kind=sgl),save                     :: prefactor

! other variables
real(kind=sgl),allocatable              :: scin_x(:), scin_y(:)                 ! scintillator coordinate ararays [microns]
real(kind=sgl),parameter                :: dtor = 0.0174533  ! convert from degrees to radians
real(kind=sgl)                          :: alp, ca, sa, cw, sw, quat(4)
real(kind=sgl)                          :: L2, Ls, Lc     ! distances
integer(kind=irg)                       :: nix, niy, binx, biny,  nixp, niyp, i, j, Emin, Emax, istat, k      ! various parameters
real(kind=sgl)                          :: dc(3), scl           ! direction cosine array
real(kind=sgl)                          :: sx, dx, dxm, dy, dym, rhos, x, bindx         ! various parameters
real(kind=sgl)                          :: ixy(2)
real(kind=dbl),parameter                :: nAmpere = 6.241D+18 



!====================================
! ------ generate the detector rgx, rgy, rgz arrays if needed (calling program must decide this via ipar(1))
!====================================
if (ipar(1).eq.1) then

! This needs to be done only once for a given detector geometry
  allocate(scin_x(ipar(2)),scin_y(ipar(3)),stat=istat)
  
  scin_x = - ( fpar(1) - ( 1.0 - float(ipar(2)) ) * 0.5 - (/ (i-1, i=1,ipar(2)) /) ) * fpar(3)
  scin_y = ( fpar(2) - ( 1.0 - float(ipar(3)) ) * 0.5 - (/ (i-1, i=1,ipar(3)) /) ) * fpar(3)

! auxiliary angle to rotate between reference frames
  alp = 0.5 * cPi - (fpar(4) - fpar(6)) * dtor
  ca = cos(alp)
  sa = sin(alp)

  cw = cos(fpar(5) * dtor)
  sw = sin(fpar(5) * dtor)

! compute auxilliary interpolation arrays
  if (allocated(rgx)) deallocate(rgx, rgy, rgz)

  allocate(rgx(ipar(2),ipar(3)), rgy(ipar(2),ipar(3)), rgz(ipar(2),ipar(3)))

  L2 = fpar(7) * fpar(7)
  do j=1,ipar(2)
    sx = L2 + scin_x(j) * scin_x(j)
    Ls = -sw * scin_x(j) + fpar(7) * cw
    Lc = cw * scin_x(j) + fpar(7) * sw
    do i=1,ipar(3)
     rhos = 1.0/sqrt(sx + scin_y(i)**2)
     rgx(j,i) = (scin_y(i) * ca + sa * Ls) * rhos
     rgy(j,i) = Lc * rhos
     rgz(j,i) = (-sa * scin_y(i) + ca * Ls) * rhos
    end do
  end do

  deallocate(scin_x, scin_y)

!====================================
! ------ create the equivalent detector energy array
!====================================
! from the Monte Carlo energy data, we need to extract the relevant
! entries for the detector geometry defined above.  Once that is 
! done, we can get rid of the larger energy array
!

! determine the scale factor for the Lambert interpolation; the square has
! an edge length of 2 x sqrt(pi/2)
  scl = float(ipar(5)) 

! energy summation will go over all energy bins
  Emin = 1
  Emax = ipar(4)

  if (allocated(accum_e_detector)) deallocate(accum_e_detector)

  allocate(accum_e_detector(ipar(4),ipar(2),ipar(3)))

  do i=1,ipar(2)
    do j=1,ipar(3)
! do the coordinate transformation for this detector pixel
       dc = (/ rgx(i,j), rgy(i,j), rgz(i,j) /)
! make sure the third one is positive; if not, switch all 
       if (dc(3).lt.0.0) dc = -dc
! convert these direction cosines to coordinates in the Rosca-Lambert projection
        ixy = scl * LambertSphereToSquare( dc, istat )
        x = ixy(1)
        ixy(1) = ixy(2)
        ixy(2) = -x
! four-point interpolation (bi-quadratic)
        nix = int(ipar(5)+ixy(1))-ipar(5)
        niy = int(ipar(6)+ixy(2))-ipar(6)
        dx = ixy(1)-nix
        dy = ixy(2)-niy
        dxm = 1.0-dx
        dym = 1.0-dy
! interpolate the intensity 
        do k=1,ipar(4)
          accum_e_detector(k,i,j) = accum_e(k,nix,niy) * dxm * dym + &
                                    accum_e(k,nix+1,niy) * dx * dym + &
                                    accum_e(k,nix,niy+1) * dxm * dy + &
                                    accum_e(k,nix+1,niy+1) * dx * dy
        end do
    end do
  end do 
  accum_e_detector = accum_e_detector * 0.25
  prefactor = 0.25D0 * nAmpere * fpar(8) * fpar(9)  * 1.0D-15 / sum(accum_e_detector)
end if   ! end of ipar(1)=1 test

! from here on, we simply compute the EBSD pattern by interpolation, using the saved arrays from above
! no intensity scaling or anything else...
scl = dble(ipar(8)) 

EBSDpattern = 0.0

quat(1:4) = fpar(10:13)

do i=1,ipar(2)
  do j=1,ipar(3)
! do the active coordinate transformation for this euler angle
    dc = quat_Lp(quat,  (/ rgx(i,j), rgy(i,j), rgz(i,j) /) )
! normalize dc
    dc = dc/sqrt(sum(dc*dc))
! convert these direction cosines to coordinates in the Rosca-Lambert projection (always square projection !!!)
    ixy = scl * LambertSphereToSquare( dc, istat )

    if (istat.eq.0) then 
! four-point interpolation (bi-quadratic)
      nix = int(ipar(8)+ixy(1))-ipar(8)
      niy = int(ipar(8)+ixy(2))-ipar(8)
      nixp = nix+1
      niyp = niy+1
      if (nixp.gt.ipar(8)) nixp = nix
      if (niyp.gt.ipar(8)) niyp = niy
      if (nix.lt.-ipar(8)) nix = nixp
      if (niy.lt.-ipar(8)) niy = niyp
      dx = ixy(1)-nix
      dy = ixy(2)-niy
      dxm = 1.0-dx
      dym = 1.0-dy
      if (dc(3).gt.0.0) then 
        do k=1,ipar(4) 
          EBSDpattern(i,j) = EBSDpattern(i,j) + accum_e_detector(k,i,j) * ( mLPNH(nix,niy,k) * dxm * dym +&
                                        mLPNH(nixp,niy,k) * dx * dym + mLPNH(nix,niyp,k) * dxm * dy + &
                                        mLPNH(nixp,niyp,k) * dx * dy )
        end do
      else
        do k=1,ipar(4) 
          EBSDpattern(i,j) = EBSDpattern(i,j) + accum_e_detector(k,i,j) * ( mLPSH(nix,niy,k) * dxm * dym +&
                                        mLPSH(nixp,niy,k) * dx * dym + mLPSH(nix,niyp,k) * dxm * dy + &
                                        mLPSH(nixp,niyp,k) * dx * dy )
        end do
      end if
    end if
  end do
end do
EBSDpattern = prefactor * EBSDpattern

return
end subroutine SingleEBSDPattern

!--------------------------------------------------------------------------
!
! SUBROUTINE:SingleECPattern
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief This function can be called as a standalone function to compute an electron channeling pattern
!> based on Marc's code above
!
!> @etails The main purpose of this routine and its accompanying wrapper routine is to
!> provide a way for an external program to compute an EBSD pattern.  The idea is that 
!> all the necessary arrays and variables are passed in by reference as arguments, without
!> the need for the routine to fetch any other data from files etc...  The initial goal is
!> to have a function that can be called with the CALL_EXTERNAL mechanism in IDL, but 
!> in the long run this will also be the approach for calling the routine from C/C++, which
!> is an essential part of integration with DREAM.3D.  This routine is a simplified version
!> of the core of the EMEBSD program. 
!>
!> This routine will first compute the incident cone vectors etc. if necessary, and then perform
!> the usual interpolation from the square Lambert projection. The pattern will be a basic pattern,
!> without any intensity scaling or binning etc; the calling program should take care of those 
!> operations.
!
!> @param nipar number of entries in ipar
!> @param nfpar number of entries in fpar
!> @param n1, n2 dimension parameters for accum_e
!> @param m1 dimension parameter for mLPNH and mLPSH 
!> @param m2 dimension parameter for mLPNH and mLPSH 
!> @param o1,o2 dimension parameters for ECpattern
!> @param ipar array with integer input parameters
!> @param fpar array with float input parameters
!> @param accum_e array with Monte Carlo histogram
!> @param mLPNH Northern hemisphere master pattern
!> @param mLPSH Southern hemisphere master pattern
!> @param EBSDpattern output array
!
!> @date 10/16/15  SS 1.0 original
!> @date 11/02/14 MDG 1.1 put all integer parameters inside ipar and fixed size of ipar/fpar
!--------------------------------------------------------------------------
!subroutine SingleECPattern(nipar, nfpar, n1, n2, m1, m2, o1, o2, ipar, fpar, accum_e, mLPNH, mLPSH, ECpattern)
subroutine SingleECPattern(ipar, fpar, ECpattern, accum_e, mLPNH, mLPSH)

! the input parameters are all part of a ipar and fpar input arrays instead of the usual namelist structures.
! The following is the mapping:
!
! ipar(1) = 1 if GetVectorsCone detector arrays need to be computed, 0 if not (arrays will have save status)
! ipar(2) = detnumsx
! ipar(3) = detnumsy
! ipar(4) = numEbins
! ipar(5) = mcnsx
! ipar(6) = numset
! ipar(7) = mpnpx

! fpar(1) = enl%thetac
! fpar(2) = enl%sampletilt
! fpar(3) = enl%workingdistance
! fpar(4) = enl%Rin
! fpar(5) = enl%Rout
! fpar(6) = enl%sigstart
! fpar(7) = enl%sigend
! fpar(8) = enl%sigstep
! fpar(9-12) =  quaternion for requested Euler angles

use local
use constants
use Lambert
use quaternions

IMPLICIT NONE

integer(8),PARAMETER                    :: nipar=7
integer(8),PARAMETER                    :: nfpar=12
integer(8),INTENT(IN)                   :: ipar(nipar)
real(kind=sgl),INTENT(IN)               :: fpar(nfpar)
real(kind=sgl),INTENT(OUT)              :: ECpattern(ipar(2),ipar(3))
real(kind=sgl),INTENT(IN)               :: accum_e(ipar(4),-ipar(5):ipar(5),-ipar(5):ipar(5))
real(kind=sgl),INTENT(IN)               :: mLPNH(-ipar(7):ipar(7), -ipar(7):ipar(7), ipar(6))
real(kind=sgl),INTENT(IN)               :: mLPSH(-ipar(7):ipar(7), -ipar(7):ipar(7), ipar(6))

real(kind=sgl),allocatable,save         :: klist(:,:,:), rgx(:,:), rgy(:,:), rgz(:,:), weightfact(:)

real(kind=dbl),parameter                :: Rtod = 57.2957795131D0
real(kind=dbl),parameter                :: dtoR = 0.01745329251D0

real(kind=sgl)                          :: kk(3), thetacr, ktmax, delta, wf, quat(4)
real(kind=sgl)                          :: mLPNH_sum(-ipar(7):ipar(7), -ipar(7):ipar(7)) 
real(kind=sgl)                          :: mLPSH_sum(-ipar(7):ipar(7), -ipar(7):ipar(7))
integer(kind=irg)                       :: istat, imin, imax, jmin, jmax, ii ,jj, nazimuth, npolar, nsig
integer(kind=irg)                       :: ipolar, iazimuth, isig, isampletilt, nix, niy, nixp, niyp, isigp
real(kind=sgl)                          :: thetain, thetaout, polar, azimuthal, delpolar, delazimuth, om(3,3)
real(kind=sgl)                          :: dc(3), scl, deltheta, acc_sum, MCangle, ixy(2), dx, dy, dxm, dym, dp


write (*,*) 'f90 accum_e        : ',shape(accum_e)
write (*,*) 'f90 mLPNH          : ',shape(mLPNH)
write (*,*) 'f90 mLPSH          : ',shape(mLPSH)
write (*,*) 'f90 ECpattern      : ',shape(ECpattern)

!==================================================================================
! ------ generate the detector klist, rgx, rgy, rgz, weightfactors arrays if needed 
!------- (calling program must decide this via ipar(1))
!==================================================================================

if (ipar(1) .eq. 1) then

    if (allocated(klist)) deallocate(klist)
    allocate(klist(1:3,-ipar(2):ipar(2),-ipar(2):ipar(2)), stat=istat)
    kk = (/0.0,0.0,1.0/)
    thetacr = DtoR*fpar(1)
    ktmax = tan(thetacr)
    delta = 2.0*ktmax/dble(ipar(2)-1)

    imin = 1
    imax = ipar(2)
    jmin = 1
    jmax = ipar(3)
     
    do ii = imin, imax
        do jj = jmin, jmax
            klist(1:3,ii,jj) = (/-ktmax+delta*(ii-1),-ktmax+delta*(jj-1),0.0/) + kk(1:3)
            klist(1:3,ii,jj) =  klist(1:3,ii,jj)/sqrt(sum( klist(1:3,ii,jj)**2))
        end do
    end do

    thetain = atan2(fpar(4),fpar(3))
    thetaout = atan2(fpar(5),fpar(3))

    om(1,:) = (/cos(fpar(2)),0.0,sin(fpar(2))/)
    om(2,:) = (/0.0,1.0,0.0/)
    om(3,:) = (/-sin(fpar(2)),0.0,cos(fpar(2))/)

    npolar = nint((thetaout - thetain)*180.0/cPi) + 1
    delpolar = (thetaout - thetain)/float(npolar-1)

    nazimuth = 361
    delazimuth = 2.0*cPi/float(nazimuth-1)

    if (allocated(rgx)) deallocate(rgx, rgy, rgz)
    allocate(rgx(npolar, nazimuth), rgy(npolar, nazimuth), rgz(npolar, nazimuth), stat=istat)

    do ipolar = 1,npolar
         polar = thetain + float(ipolar-1)*delpolar

         do iazimuth = 1,nazimuth
             azimuthal = float(iazimuth-1)*delazimuth

             dc(1) = cos(azimuthal)*sin(polar)
             dc(2) = sin(azimuthal)*sin(polar)
             dc(3) = cos(polar)

             dc = matmul(om,dc)

             rgx(ipolar,iazimuth) = dc(1)
             rgy(ipolar,iazimuth) = dc(2)
             rgz(ipolar,iazimuth) = dc(3)
        end do
    end do

!===================================================================
! ------ generate the weight factors from the monte carlo histogram
!===================================================================

scl = float(ipar(5))
nsig = nint(fpar(1) + abs(fpar(2))) + 1

deltheta = (fpar(1)+abs(fpar(2)))/float(nsig-1)

if (allocated(weightfact)) deallocate(weightfact)

allocate(weightfact(1:nsig), stat=istat)
weightfact = 0.0

do isig = 1,nsig
    acc_sum = 0.0
    MCangle = (isig - 1)*deltheta
    isampletilt = nint((MCangle - fpar(6))/fpar(8))
    
    if (isampletilt .lt. 1) then
        isampletilt = abs(isampletilt) + 1
    else
        isampletilt = isampletilt + 1
    end if

    do ipolar = 1,npolar
        do iazimuth = 1,nazimuth
            dc(1:3) = (/rgx(ipolar,iazimuth), rgy(ipolar,iazimuth), rgz(ipolar,iazimuth)/)
! convert to Rosca-lambert projection
            ixy = scl *  LambertSphereToSquare( dc, istat )
            nix = int(ipar(5)+ixy(1))-ipar(5)
            niy = int(ipar(5)+ixy(2))-ipar(5)
            nixp = nix+1
            niyp = niy+1
            if (nixp.gt.ipar(5)) nixp = nix
            if (niyp.gt.ipar(5)) niyp = niy
            if (nix.lt.-ipar(5)) nix = nixp
            if (niy.lt.-ipar(5)) niy = niyp
            dx = ixy(1)-nix
            dy = ixy(2)-niy
            dxm = 1.0-dx
            dym = 1.0-dy
            
            acc_sum = 0.25*(accum_e(isampletilt,nix,niy) * dxm * dym + &
                            accum_e(isampletilt,nixp,niy) * dx * dym + &
                            accum_e(isampletilt,nix,niyp) * dxm * dy + &
                            accum_e(isampletilt,nixp,niyp) * dx * dy)
             
            weightfact(isig) = weightfact(isig) + acc_sum

        end do
    end do
end do

weightfact(1:nsig) = weightfact(1:nsig)/weightfact(1)

end if

!===================================================================
! ------ perform interpolation from square lambert map
!===================================================================

quat(1:4) = fpar(9:12)
scl = float(ipar(7))

mLPNH_sum = sum(mLPNH,3)
mLPSH_sum = sum(mLPSH,3)

write (*,*) 'imin/imax : ',imin,imax

do ii = imin, imax
    do jj = jmin, jmax

        dc(1:3) = klist(1:3,ii,jj)
        dp = DOT_PRODUCT(dc(1:3),(/sin(fpar(2)*dtoR),0.D0,cos(fpar(2)*dtoR)/))      
      
        MCangle = acos(dp)*Rtod

        isig = int(MCangle) + 1
        if (isig .gt. nsig) isig = nsig

        isigp = isig + 1
        if (isigp .gt. nsig) isigp = nsig

        dx = MCangle - int(MCangle)
        dxm =  1.0 - dx
 
        wf = weightfact(isig) * dxm + weightfact(isigp) * dx
 
        dc = quat_LP(quat, dc)
        dc = dc/sqrt(sum(dc*dc))

        ixy = scl * LambertSphereToSquare( dc, istat )
        nix = int(ipar(7)+ixy(1))-ipar(7)
        niy = int(ipar(7)+ixy(2))-ipar(7)
        nixp = nix+1
        niyp = niy+1
        if (nixp.gt.ipar(7)) nixp = nix
        if (niyp.gt.ipar(7)) niyp = niy
        if (nix.lt.-ipar(7)) nix = nixp
        if (niy.lt.-ipar(7)) niy = niyp
        dx = ixy(1)-nix
        dy = ixy(2)-niy
        dxm = 1.0-dx
        dym = 1.0-dy

        if (dc(3).gt.0.0) then 
            ECpattern(ii,jj) = wf * ( mLPNH_sum(nix,niy) * dxm * dym + &
                         mLPNH_sum(nixp,niy) * dx * dym + &
                         mLPNH_sum(nix,niyp) * dxm * dy + &
                         mLPNH_sum(nixp,niyp) * dx * dy )

        else
            ECpattern(ii,jj) =  wf * ( mLPSH_sum(nix,niy) * dxm * dym + &
                         mLPSH_sum(nixp,niy) * dx * dym + &
                         mLPSH_sum(nix,niyp) * dxm * dy + &
                         mLPSH_sum(nixp,niyp) * dx * dy )
        end if

    end do
end do

write (*,*) 'max = ',maxval(ECpattern), shape(ECpattern)
write (*,*) 'max(NH) ',maxval(mLPNH_sum), maxval(mLPSH_sum)

open(unit=10,file='test.out',form='unformatted',status='unknown')
write(10) ipar(2), ipar(3)
write(10) ECpattern
close(unit=10,status='keep')


end subroutine SingleECPattern

!--------------------------------------------------------------------------
!
! SUBROUTINE:SingleEBSDPatternWrapper
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief wrapper routine for SingleEBSDPattern
!>
!> see https://groups.google.com/forum/#!topic/comp.lang.idl-pvwave/Gk0xxVFbW8E
!>
!
!> @param argc number of argument
!> @param argv pointers to subroutine parameters
!
!> @date 10/16/15 MDG 1.0 original
!--------------------------------------------------------------------------
function SingleEBSDPatternWrapper(argc, argv) bind(c, name='SingleEBSDPatternWrapper') 

use,INTRINSIC :: ISO_C_BINDING

IMPLICIT NONE

INTEGER(c_size_t), VALUE, INTENT(IN)            :: argc 
type(c_ptr), dimension(argc), INTENT(INOUT)     :: argv
REAL(c_float)                                   :: SingleEBSDPatternWrapper

! wrapper function dependent declarations; they are all pointers 
! since we pass everything by reference from IDL 
integer(c_size_t), pointer                      :: nipar, nfpar, n1, n2, m1, o1, o2 
integer(c_size_t),dimension(:), pointer         :: ipar
real(c_float), dimension(:), pointer            :: fpar
real(c_float), dimension(:,:), pointer          :: EBSDpattern
real(c_float), dimension(:,:,:), pointer        :: accum_e, mLPNH, mLPSH

! the following line just helps in identifying the correct order of the subroutine arguments...
!                             1      2      3   4   5   6   7   8     9     10       11     12     13
!subroutine SingleEBSDPattern(nipar, nfpar, n1, n2, m1, o1, o2, ipar, fpar, accum_e, mLPNH, mLPSH, EBSDpattern)
!
! transform the C pointers above to fortran pointers, and use them in the regular function call
call c_f_pointer(argv(1),nipar) 
call c_f_pointer(argv(2),nfpar) 
call c_f_pointer(argv(3),n1) 
call c_f_pointer(argv(4),n2) 
call c_f_pointer(argv(5),m1) 
call c_f_pointer(argv(6),o1) 
call c_f_pointer(argv(7),o2) 
call c_f_pointer(argv(8),ipar,(/nipar/)) 
call c_f_pointer(argv(9),fpar,(/nfpar/)) 
call c_f_pointer(argv(10),accum_e,(/n1,2*n2+1,2*n2+1/)) 
call c_f_pointer(argv(11),mLPNH,(/2*m1+1,2*m1+1,n1/)) 
call c_f_pointer(argv(12),mLPSH,(/2*m1+1,2*m1+1,n1/)) 
call c_f_pointer(argv(13),EBSDpattern,(/o1,o2/)) 

call SingleEBSDPattern(nipar, nfpar, n1, n2, m1, o1, o2, ipar, fpar, accum_e, mLPNH, mLPSH, EBSDpattern)

SingleEBSDPatternWrapper = 1._c_float
end function SingleEBSDPatternWrapper


!--------------------------------------------------------------------------
!
! SUBROUTINE:SingleECPatternWrapper
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief wrapper routine for SingleECPPattern; based on Marc's routine above
!>
!> see https://groups.google.com/forum/#!topic/comp.lang.idl-pvwave/Gk0xxVFbW8E
!>
!
!> @param argc number of argument
!> @param argv pointers to subroutine parameters
!
!> @date 10/28/15  SS 1.0 original
!> @date 11/02/15 MDG 1.1 simplified parameters
!--------------------------------------------------------------------------
function SingleECPatternWrapper(argc, argv) bind(c, name='SingleECPatternWrapper') 

use,INTRINSIC :: ISO_C_BINDING

IMPLICIT NONE

INTEGER(c_size_t), VALUE, INTENT(IN)            :: argc 
type(c_ptr), dimension(argc), INTENT(INOUT)     :: argv
REAL(c_float)                                   :: SingleECPatternWrapper

! wrapper function dependent declarations; they are all pointers 
! since we pass everything by reference from IDL 
integer(c_size_t)                               :: nipar, nfpar
integer(c_size_t),dimension(:), pointer         :: ipar
real(c_float), dimension(:), pointer            :: fpar
real(c_float), dimension(:,:,:), pointer        :: accum_e, mLPNH, mLPSH
real(c_float), dimension(:,:), pointer          :: ECpattern

! the following line just helps in identifying the correct order of the subroutine arguments...
!                             1      2     3       4       5      6
!subroutine SingleECPattern(ipar, fpar, accum_e, mLPNH, mLPSH, ECPattern)
!
! transform the C pointers above to fortran pointers, and use them in the regular function call
nipar = 7
nfpar = 12
call c_f_pointer(argv(1),ipar,(/nipar/)) 
call c_f_pointer(argv(2),fpar,(/nfpar/)) 
call c_f_pointer(argv(3),ECpattern,(/ipar(2),ipar(3)/))
call c_f_pointer(argv(4),accum_e,(/ipar(4),2*ipar(5)+1,2*ipar(5)+1/))
call c_f_pointer(argv(5),mLPNH,(/2*ipar(7)+1, 2*ipar(7)+1, ipar(6)/))
call c_f_pointer(argv(6),mLPSH,(/2*ipar(7)+1, 2*ipar(7)+1, ipar(6)/))

call SingleECPattern(ipar, fpar, ECpattern, accum_e, mLPNH, mLPSH)

write (*,*) 'maxval(ECpattern) ',maxval(ECpattern),shape(ECpattern)

SingleECPatternWrapper = 1._c_float
end function SingleECPatternWrapper


end module EMdymod
