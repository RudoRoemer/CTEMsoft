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
! CTEMsoft2013:STEMmodule.f90
!--------------------------------------------------------------------------
!
! MODULE: STEMmodule
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Provides routines to handle the STEM detector geometry and weightfactors.
! 
!> @date   04/29/11 MDG 1.0 original
!> @date   06/12/13 MDG 2.0 rewrite 
!--------------------------------------------------------------------------
module STEMmodule

use local

type STEMtype
	character(fnlen)		:: weightoutput
	character(2)  			:: geometry
	integer(kind=irg) 		:: numberofsvalues, numk, numCL
	real(kind=sgl) 			:: BFradius,ADFinnerradius,ADFouterradius,kt,beamconvergence,cameralength, &
	                          	   BFmrad,ADFimrad,ADFomrad, diffapmrad, diffapmcenter, CLarray(20)
	logical,allocatable  		:: ZABFweightsarray(:,:,:),ZAADFweightsarray(:,:,:)       ! only used for the zone axis case
	real(kind=sgl),allocatable  	:: sgarray(:,:),BFweightsarray(:,:,:),ADFweightsarray(:,:,:)   ! only used for the systematic row case
end type STEMtype

type(STEMtype) :: STEM


contains

!--------------------------------------------------------------------------
!
! SUBROUTINE: init_STEM
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief initialize the weight factors for the systematic row case.
!
!> @param nn number of beams
!> @param g systematic row basic vector
! 
!> @date   04/29/11 MDG 1.0 original
!> @date   06/12/13 MDG 2.0 rewrite 
!--------------------------------------------------------------------------
subroutine init_STEM(nn,g)

use local
use io
use crystal
use diffraction

IMPLICIT NONE

integer(kind=irg),INTENT(IN)		:: nn
integer(kind=irg),INTENT(IN)		:: g(3)

integer(kind=irg) 			:: i,j,n,ira,jj,k,kk, iCL
real(kind=sgl) 				:: glen, thb, alp, omega_c, omega_min, omega_max,omega,a,b,c,th,dom,p,q,dr,dx
real(kind=sgl),parameter 		:: cPi=3.141592654

! these are only used to debug this routine
real(kind=sgl),allocatable 		:: thetar(:),outar(:,:,:)
logical    				:: debug = .FALSE., diffappresent = .FALSE., apinBF=.FALSE. , apinADF = .FALSE.

! this routine initializes the excitation error arrays and the weight-factor arrays for systematic row STEM signals
! we'll assume that the central beam is centered on the BF detector; then we can 
! compute the complete geometry by working in mrad units throughout.

! allocate the excitation error array areal(1..nn,1..STEM%numberofsvalues)
  allocate(STEM%sgarray(nn,STEM%numberofsvalues))

! determine the lower and upper bounds of the excitation error for the fundamental reflection G
  thb = CalcDiffAngle(g(1),g(2),g(3))*0.5  ! Bragg angle in radians

! convert k_t to the alp and omega angles (in radians)
  glen = CalcLength(float(g),'r')
  alp = -2.0*STEM%kt*thb
  omega_c = cPi*0.5+alp
  omega_min = omega_c - STEM%beamconvergence/1000.0
  omega_max = omega_c + STEM%beamconvergence/1000.0

! step size
  dom = (omega_max - omega_min)/float(STEM%numberofsvalues-1)

! and for each value in between, compute each reflection's excitation error
ira = (nn+1)/2

do j=1,STEM%numberofsvalues
! set omega angle
   omega = omega_min+float(j-1)*dom
   do i=1,nn
    n = -ira+i
! excitation error
    STEM%sgarray(nn+1-i,j) = -n*glen*cos(omega)-(1.0-sqrt(1.0-(n*mLambda*glen*sin(omega))**2))/mLambda
   end do
end do

if (debug) then
  allocate(thetar(25))
  thetar = 0.0
  thetar(1:7) = (/ a,b,c,th,dom,dr,dx /)
end if 

! next, we compute the weightfactors, i.e., how much does each excitation error value contribute
! to the BF or ADF signal?  The weight factor is basically the length of the chord across the overlap
! area of the diffraction disk and the detector, which requires a little bit of math to figure out;
! the math employs the concept of the radical line (see mathworld.com section on circle-circle intersections)

! this computation is carried out in mrad units !
allocate(STEM%BFweightsarray(nn,STEM%numberofsvalues,STEM%numCL),STEM%ADFweightsarray(nn,STEM%numberofsvalues,STEM%numCL))

STEM%BFweightsarray = 0.0
STEM%ADFweightsarray = 0.0


outerCLloop: do iCL=1,STEM%numCL    ! this is the outer loop over the microscope camera lengths (very long loop !!!)

! fist, convert the detector parameters to mrad units
STEM%BFmrad = atan(STEM%BFradius/STEM%CLarray(iCL))*1000.0
STEM%ADFimrad = atan(STEM%ADFinnerradius/STEM%CLarray(iCL))*1000.0
STEM%ADFomrad = atan(STEM%ADFouterradius/STEM%CLarray(iCL))*1000.0
if (STEM%diffapmrad.ne.0.0) diffappresent = .TRUE.

! then, for each point inside each diffraction disk, determine where it falls with respect to
! the BD and ADF detectors ... Also, look for disk overlaps as they might require amplitudes
! to be added in instead of intensities (for starters, we could just not allow that to happen...)

! rename some variables to short symbols
a = STEM%ADFimrad
b = STEM%ADFomrad
c = STEM%BFmrad
th = STEM%beamconvergence
n = STEM%numberofsvalues
if (diffappresent) then
  dr = STEM%diffapmrad
  dx = STEM%diffapmcenter
end if
omega_min =  -th
omega_max = th
dom = 2.0*th/float(n-1)

if (.not.diffappresent) then   ! there is no diffraction aperture, so compute the regular weight factors
! first, do the math for the g=0 disk  (we're dropping common factors of 2 for the weights) 
i = ira
do j=(n+1)/2,n
  omega = omega_min+float(j-1)*dom
  if (th.gt.c) then       ! the zero disk is larger than the BF detector, so it (potentially) gives two signals
    if (omega.le.c) STEM%BFweightsarray(i,j,iCL) = sqrt(c**2-omega**2)
    if (th.ge.a) then    ! there's overlap with the ADF detector
      if (omega.le.a) then  ! the first part needs to have a bit subtracted
        STEM%ADFweightsarray(i,j,iCL) = sqrt((th**2-omega**2)) - sqrt((a**2-omega**2))
      else   ! the second part does not
        STEM%ADFweightsarray(i,j,iCL) = sqrt((th**2-omega**2)) 
      end if
    end if
  else                         ! the zero disk is smaller than the BF detector, so only a BF signal
   STEM%BFweightsarray(i,j,iCL) = sqrt((th**2-omega**2))
  end if
! then apply symmetry for the other half of the g=0 disk
  if (j.ne.(n+1)/2) then
    jj = n+1 - j
    STEM%BFweightsarray(i,jj,iCL) = STEM%BFweightsarray(i,j,iCL)
    STEM%ADFweightsarray(i,jj,iCL) = STEM%ADFweightsarray(i,j,iCL)
  end if  
end do  ! that completes the central disk weight factors

! the other disks are quite a bit more difficult to deal with ... there are a lot of possible cases to consider ...
do i=ira+1,nn      ! loop over the positive reflections of the systematic row (the rest follows by symmetry)
! redefine a couple of parameters
  j = i - ira
  thb = CalcDiffAngle(j*g(1),j*g(2),j*g(3))*1000.0  ! diffraction angle in mrad
  omega_min = thb - th
  omega_max = thb + th
! only used for debugging
 if (debug)  thetar(7+j) = thb
! first check if a part of this disk lies inside the BF detector
  if (omega_min.lt.c) then     ! yes, it does, so determine the BF weight factors for this disk
    if (omega_max.le.c) then  ! does it lie completely inside the BF detector?
      do j=1,n   ! yes it does, so compute all weight factors
        omega = omega_min + float(j-1)*dom
        STEM%BFweightsarray(i,j,iCL) = sqrt(th**2 - (thb-omega)**2)
        STEM%BFweightsarray(2*ira-i,n+1-j,iCL) = STEM%BFweightsarray(i,j,iCL)
      end do
    else  ! no, there's some potential overlap with the ADF detector 
      do j=1,n   ! once again, there are a few cases
        omega = omega_min + float(j-1)*dom    ! this is the position
        p = (thb**2-th**2+a**2)*0.5/thb             ! this is the location of the radical line for the ADF detector
        q = (thb**2-th**2+c**2)*0.5/thb             ! this is the location of the radical line for the BF detector
        if (omega.le.q) then   ! this point contributes to the BF detector
          STEM%BFweightsarray(i,j,iCL) = sqrt(th**2 - (thb-omega)**2)
        end if
        if ((omega.gt.q).and.(omega.le.c)) then   ! this point contributes to the BF detector
          STEM%BFweightsarray(i,j,iCL) = sqrt(c**2 - omega**2)
        end if
        if ((omega_max.ge.a).and.(omega.ge.p).and.(omega.le.a)) then ! this point contributes to the ADF detector (using radical line position)
          STEM%ADFweightsarray(i,j,iCL) = sqrt(th**2 -  (thb-omega)**2)  - sqrt(a**2 - omega**2)
        end if
         if ((omega_max.ge.a).and.(omega.gt.a)) then ! this point lies on the ADF detector 
          STEM%ADFweightsarray(i,j,iCL) = sqrt(th**2 -  (thb-omega)**2)
        end if
       STEM%BFweightsarray(2*ira-i,n+1-j,iCL) = STEM%BFweightsarray(i,j,iCL)       
       STEM%ADFweightsarray(2*ira-i,n+1-j,iCL) = STEM%ADFweightsarray(i,j,iCL)
      end do
    end if 
  else    ! no, it does not intersect the BF detector, so this disk can only contribute to the ADF weight factors
! once more there are several cases which we'll treat in increasing value of the position...
      do j=1,n 
        omega = omega_min + float(j-1)*dom    ! this is the position
        p = (thb**2-th**2+a**2)*0.5/thb             ! this is the location of the radical line for the inner ADF detector edge
        q = (thb**2-th**2+b**2)*0.5/thb             ! this is the location of the radical line for the outer ADF detector edge
        if ((omega.lt.a).and.(omega.ge.p))  then    ! inside the inner ADF edge, but close enough to contribute
          STEM%ADFweightsarray(i,j,iCL) = sqrt(th**2 -  (thb-omega)**2)  - sqrt(a**2 - omega**2)
        end if
         if ((omega.ge.a).and.(omega_max.le.b)) then ! this point lies on the ADF detector 
          STEM%ADFweightsarray(i,j,iCL) = sqrt(th**2 -  (thb-omega)**2)
        end if
        if ((omega_max.gt.b).and.(omega.le.q)) then   ! this point lies on the ADF detector
          STEM%ADFweightsarray(i,j,iCL) = sqrt(th**2 - (thb-omega)**2)
        end if
        if ((omega_max.gt.b).and.(omega.gt.q).and.(omega.le.b))  then   ! this point contributes to the ADF detector
          STEM%ADFweightsarray(i,j,iCL) = sqrt(b**2 - omega**2)
        end if
        STEM%ADFweightsarray(2*ira-i,n+1-j,iCL) = STEM%ADFweightsarray(i,j,iCL)
      end do
  end if

end do
end if ! end of regular weight factors without a diffraction aperture



if (diffappresent) then   ! there is a diffraction aperture, so revisit the weight factors.
! once again, there are many different cases that need to be addressed...
  
! we do not allow for a diffraction aperture that overlaps the boundary between BF and ADF detectors,
! nor an aperture that lies entirely beyond the ADF detector
! first the BF test
  if ( ((dx-dr).gt.-c).and.((dx+dr).lt.c))  apinBF = .TRUE.

! then the ADF detector
  if ( (((dx-dr).gt.-b).and.((dx+dr).lt.-c)) .or.(((dx-dr).gt.a).and.((dx+dr).lt.b)) )  apinADF = .TRUE. 

! if the aperture is outside the ADF detector, or it overlaps the space between the detectors, then abort
  if ( .not.apinBF .and. .not.apinADF ) then
    mess = 'Please fix input: Diffraction aperture outside BF detector disk or ADF ring !'  
     call Message("(A)")
     write (*,*) apinBF, apinADF, a,b,c, dx-dr, dx+dr, dx,dr
     stop
  end if



if (apinBF) then
! figure out which diffraction disk(s) contribute to the BF detector
 do i=1,nn      ! loop over all reflections of the systematic row
! redefine a couple of parameters
  j = -(nn-1)/2-1+i
  if (j.ne.0) then 
    thb = (j/abs(j)) * CalcDiffAngle(j*g(1),j*g(2),j*g(3))*1000.0  ! diffraction angle in mrad
  else
    thb = 0.0
  end if  
 ! only used for debugging
 if (debug)  thetar(7+i) = thb
  omega_min = thb - th
  omega_max = thb + th
! check whether or not there is any overlap between this disk and the diffraction aperture opening
  if ((omega_max.lt.(dx-dr)).or.(omega_min.gt.(dx+dr))) then  ! this disks does not fall inside the diffraction aperture 
    STEM%BFweightsarray(i,1:n,iCL) = 0.0
  else
! case 1: dx-dr < omega_min < omega_max < dx+dr
    if ((omega_max.lt.(dx+dr)).and.(omega_min.gt.(dx-dr))) then
      do k=1,n
        omega = omega_min+float(k-1)*dom
        kk = k
        if (j.lt.0) kk = n+1-k 
        STEM%BFweightsarray(i,kk,iCL) =  sqrt((th**2-(thb-omega)**2))
     end do 
    end if
! case 2: omega_min < dx-dr  < dx+dr < omega_max 
    if ((omega_max.gt.(dx+dr)).and.(omega_min.lt.(dx-dr))) then
      do k=1,n
        omega = omega_min+float(k-1)*dom
        kk = k
        if (j.lt.0) kk = n+1-k 
        if ((omega.gt.dx-dr).and.(omega.lt.dx+dr)) then
           STEM%BFweightsarray(i,kk,iCL) =  sqrt((dr**2-(omega-dx)**2))
        end if
      end do 
    end if
! case 3: omega_min < dx-dr   < omega_max < dx+dr
    if ((omega_min.lt.(dx-dr)).and.(omega_max.lt.(dx+dr))) then
       p = ((dx-thb)**2-dr**2+th**2)/2.0/(dx-thb)+thb
       do k=1,n
         omega = omega_min+float(k-1)*dom
         kk = k
!         if (j.lt.0) kk = n+1-k 
        if ((omega.gt.dx-dr).and.(omega.le.p)) then
          STEM%BFweightsarray(i,kk,iCL) =  sqrt(dr**2-(omega-dx)**2)
        end if
        if ((omega.gt.p).and.(omega.le.omega_max)) then
          STEM%BFweightsarray(i,kk,iCL) =  sqrt((th**2-(thb-omega)**2))       
        end if
      end do 
    end if
 ! case 4:  dx-dr   < omega_min < dx+dr < omega_max
    if ((omega_min.gt.(dx-dr)).and.(omega_max.gt.(dx+dr))) then
       p = ((dx-thb)**2-th**2+dr**2)/2.0/(thb-dx) + dx
        do k=1,n
        omega = omega_min+float(k-1)*dom
        kk = k
!        if (j.lt.0) kk = n+1-k 
        if ((omega.gt.p).and.(omega.le.dx+dr)) then
          STEM%BFweightsarray(i,kk,iCL) = sqrt((dr**2-(omega-dx)**2)) 
        end if
        if ((omega.gt.omega_min).and.(omega.le.p)) then
           STEM%BFweightsarray(i,kk,iCL) =  sqrt((th**2-(thb-omega)**2))
        end if
      end do 
    end if 
  end if
 end do  ! this completes the BF weight factors when a diffraction aperture is present and apinBF=.TRUE.
end if 


! next determine the ADF weight factors in the presence of an aperture
if (apinADF) then
! figure out which diffraction disk(s) contribute to the ADF detector
 do i=1,nn      ! loop over all reflections of the systematic row
! redefine a couple of parameters
    j = -(nn-1)/2-1+i
    if (j.ne.0) then 
      thb = (j/abs(j)) * CalcDiffAngle(j*g(1),j*g(2),j*g(3))*1000.0  ! diffraction angle in mrad
    else
      thb = 0.0
   end if  
 ! only used for debugging
 if (debug)  thetar(7+i) = thb
  omega_min = thb - th
  omega_max = thb + th
! check whether or not there is any overlap between this disk and the diffraction aperture opening
  if ((omega_max.lt.(dx-dr)).or.(omega_min.gt.(dx+dr))) then  ! this disks does not fall inside the diffraction aperture 
    STEM%ADFweightsarray(i,1:n,iCL) = 0.0
  else
! case 1: dx-dr < omega_min < omega_max < dx+dr
    if ((omega_max.lt.(dx+dr)).and.(omega_min.gt.(dx-dr))) then
      do k=1,n
        omega = omega_min+float(k-1)*dom
        kk = k
        if (j.lt.0) kk = n+1-k 
        STEM%ADFweightsarray(i,kk,iCL) =  sqrt((th**2-(thb-omega)**2))
     end do 
    end if
! case 2: omega_min < dx-dr  < dx+dr < omega_max 
    if ((omega_max.gt.(dx+dr)).and.(omega_min.lt.(dx-dr))) then
      do k=1,n
        omega = omega_min+float(k-1)*dom
        kk = k
        if (j.lt.0) kk = n+1-k 
        if ((omega.gt.dx-dr).and.(omega.lt.dx+dr)) then
           STEM%ADFweightsarray(i,kk,iCL) =  sqrt((dr**2-(omega-dx)**2))
        end if
      end do 
    end if
! case 3: omega_min < dx-dr   < omega_max < dx+dr
    if ((omega_min.lt.(dx-dr)).and.(omega_max.lt.(dx+dr))) then
       p = ((dx-thb)**2-dr**2+th**2)/2.0/(dx-thb)+thb
       do k=1,n
         omega = omega_min+float(k-1)*dom
         kk = k
!         if (j.lt.0) kk = n+1-k 
        if ((omega.gt.dx-dr).and.(omega.le.p)) then
          STEM%ADFweightsarray(i,kk,iCL) =  sqrt(dr**2-(omega-dx)**2)
        end if
        if ((omega.gt.p).and.(omega.le.omega_max)) then
          STEM%ADFweightsarray(i,kk,iCL) =  sqrt((th**2-(thb-omega)**2))       
        end if
      end do 
    end if
 ! case 4:  dx-dr   < omega_min < dx+dr < omega_max
    if ((omega_min.gt.(dx-dr)).and.(omega_max.gt.(dx+dr))) then
       p = ((dx-thb)**2-th**2+dr**2)/2.0/(thb-dx) + dx
        do k=1,n
        omega = omega_min+float(k-1)*dom
        kk = k
!        if (j.lt.0) kk = n+1-k 
        if ((omega.gt.p).and.(omega.le.dx+dr)) then
          STEM%ADFweightsarray(i,kk,iCL) = sqrt((dr**2-(omega-dx)**2)) 
        end if
        if ((omega.gt.omega_min).and.(omega.le.p)) then
           STEM%ADFweightsarray(i,kk,iCL) =  sqrt((th**2-(thb-omega)**2))
        end if
      end do 
    end if 
  end if
 end do  ! this completes the ADF weight factors when a diffraction aperture is present and apinADF=.TRUE.
 
end if

end if ! if aperture is present


end do outerCLloop   ! see line 146


! and the rest is also only used for debugging purposes
if (debug) then 
  allocate(outar(2*nn,STEM%numberofsvalues,STEM%numCL))
  outar(1:nn,1:STEM%numberofsvalues,1:STEM%numCL) = STEM%BFweightsarray(1:nn,1:STEM%numberofsvalues,1:STEM%numCL)
  outar(nn+1:2*nn,1:STEM%numberofsvalues,1:STEM%numCL) = STEM%ADFweightsarray(1:nn,1:STEM%numberofsvalues,1:STEM%numCL)
! to make sure that everything is correct, let's export this array so that we can display it in IDL
  open(unit=dataunit,file='STEMprofiles.data',status='unknown',form='unformatted')
  write(unit=dataunit) nn,STEM%numberofsvalues,STEM%numCL
  write(unit=dataunit) thetar
  write(unit=dataunit) outar
  close(unit=dataunit,status='keep')
end if

end subroutine init_STEM



!--------------------------------------------------------------------------
!
! SUBROUTINE: init_STEM_ZA
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief initialize weight factors for zone-axis STEM case
! 
!> @note This will need to be reconsidered when we implement sectored detectors ... 
!
!> @param nn number of reflections
! 
!> @date   04/29/11 MDG 1.0 original
!> @date   06/12/13 MDG 2.0 rewrite 
!--------------------------------------------------------------------------
subroutine init_STEM_ZA(nn)

use local
use crystal
use diffraction
use dynamical
use foilmodule
use kvectors
use gvectors

IMPLICIT NONE

integer(kind=irg),INTENT(IN) :: nn

integer(kind=irg)                   :: ik,ig, iCL
real(kind=sgl)                      :: ll(3), lpg(3), gg(3), glen, gplen, kpg

! this routine initializes the excitation error arrays and the weight-factor arrays for zone axis STEM signals
! the weightfactors are quite a bit different from the ones for the systematic row case;
! they are simpler in principle, since each point in the diffracted disk can only lie in one
! place, and hence only contributes to one detector.  However, not all points in a disk
! contribute to the same detector...  The length of the vector k_t+g, expressed in mrad,
! is what needs to be compared to the radii of the BF and ADF detectors.  For each incident 
! beam direction, we take the tangential component of the wave vector and loop over all
! reflections to compute the relevant angle; this then allows us to assign the weight factors
! which are now either 1 or 0 (so they can be stored as logicals).

! allocate the excitation error array areal(1..nn,1..STEM%numk)
  allocate(STEM%sgarray(nn,STEM%numk))
  
! transform the foil normal to real space and normalize
  call TransSpace(sngl(foil%F),DynFN,'d','r')
  call NormVec(DynFN,'r')

! allocate the weight factor arrays, one entry for each beam direction, reflection, and camera length
  allocate(STEM%ZABFweightsarray(nn,STEM%numk,STEM%numCL),STEM%ZAADFweightsarray(nn,STEM%numk,STEM%numCL))
  STEM%ZABFweightsarray = .FALSE.
  STEM%ZAADFweightsarray = .FALSE.

! loop over the wave vector linked list
  ktmp => khead
  beamloopCL: do ik=1,STEM%numk
    ll = ktmp%kt        ! this is the tangential component of the wave vector
! and loop over all reflections
    rltmpa => reflist%next
    reflectionloopCL: do ig=1,nn
      gg = float(rltmpa%hkl)
      glen = CalcLength(gg,'r')
      lpg = ll + gg                ! Laue + g
      gplen = CalcLength(lpg,'r')
      kpg = 2000.0*asin(0.50*sngl(mLambda)*gplen)    ! 2theta in mrad
      do iCL=1,STEM%numCL
        STEM%BFmrad = atan(STEM%BFradius/STEM%CLarray(iCL))*1000.0
        STEM%ADFimrad = atan(STEM%ADFinnerradius/STEM%CLarray(iCL))*1000.0
        STEM%ADFomrad = atan(STEM%ADFouterradius/STEM%CLarray(iCL))*1000.0
        if (kpg.le.STEM%BFmrad) STEM%ZABFweightsarray(ig,ik,iCL) = .TRUE.
        if ((kpg.ge.STEM%ADFimrad).AND.(kpg.le.STEM%ADFomrad)) STEM%ZAADFweightsarray(ig,ik,iCL) = .TRUE.
      end do  ! loop over camera lengths
      STEM%sgarray(ig,ik) = Calcsg(gg,sngl(ktmp%k),DynFN)
 ! and we move to the next reflection in the list
      rltmpa => rltmpa%next
    end do reflectionloopCL  
    ktmp => ktmp%next
  end do beamloopCL

!  open(unit=dataunit,file='ZAbfprofiles.data',status='unknown',form='unformatted')
!  write(unit=dataunit) nn,STEM%numk
!  write(unit=dataunit) int(STEM%ZABFweightsarray)
!  close(unit=dataunit,status='keep')
!  open(unit=dataunit,file='ZAadfprofiles.data',status='unknown',form='unformatted')
!  write(unit=dataunit) nn,STEM%numk
!  write(unit=dataunit)  int(STEM%ZAADFweightsarray)
!  close(unit=dataunit,status='keep')

! that's it folks!
end subroutine init_STEM_ZA



!--------------------------------------------------------------------------
!
! SUBROUTINE: read_STEM_data
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read detector and other parameters for the STEM case
! 
!> @param STEMnmlfile filename of the namelist file
!> @param geometry 'SR' for systematic row or 'ZA' for zone axis
!> @param nn number of reflections
!> @param g fundamental g-vector for systematic row
!> @param kt tangential wave vector component
!> @param numk number of distinct wave vectors (optional)
!> @param beamdiv beam divergence parameter (optional)
! 
!> @date   04/29/11 MDG 1.0 original
!> @date   06/12/13 MDG 2.0 rewrite 
!> @date   11/26/13 MDG 2.1 made geometry an input parameter instead of part of the STEMdata namelist
!--------------------------------------------------------------------------
subroutine read_STEM_data(STEMnmlfile,geometry,nn,g,kt,numk,beamdiv)

use local
use io
use files

IMPLICIT NONE

character(fnlen),INTENT(IN)			:: STEMnmlfile
character(2),INTENT(IN) 			:: geometry  ! 'SR' or 'ZA'
integer(kind=irg),INTENT(IN)			:: nn
integer(kind=irg),INTENT(IN)			:: g(3)
real(kind=sgl),INTENT(IN)			:: kt
integer(kind=irg),INTENT(IN),OPTIONAL		:: numk
real(kind=sgl),INTENT(OUT),OPTIONAL		:: beamdiv

integer(kind=irg) 				:: numberofsvalues,numCL
real(kind=sgl) 					:: BFradius,ADFinnerradius,ADFouterradius,beamconvergence,cameralength, &
                          				diffaprad, diffapcenter,CLarray(20)
character(fnlen) 				:: weightoutput

namelist / STEMdata / BFradius, ADFinnerradius, ADFouterradius, cameralength, numCL, &
                      beamconvergence, numberofsvalues, diffaprad, diffapcenter, weightoutput, CLarray

! set default values (based on OSU Tecnai F20)
BFradius = 3.5			! mm
ADFinnerradius = 3.5		! mm
ADFouterradius = 10.0		! mm
cameralength = 100.0 		! mm
numCL = 0	 		! number of camera lengths to be used  (if zero, then use cameralength instead)
CLarray = 0.0 			! values for the camera lengths
beamconvergence = 2.0		! mrad
numberofsvalues = 33		! integer
diffaprad = 0.0			! diffraction aperture radius in mrad, 0.0 if no aperture is present
diffapcenter = 0.0		! position of center of diffraction aperture in mrad along systematic row
weightoutput = '' 		! string with filename root for graphical output of weight profiles, empty if not needed

! read the namelist file
mess = 'opening '//trim(STEMnmlfile); call Message("(/A)")
OPEN(UNIT=dataunit,FILE=trim(STEMnmlfile),DELIM='APOSTROPHE')
READ(UNIT=dataunit,NML=STEMdata)
CLOSE(UNIT=dataunit)

if (PRESENT(beamdiv)) beamdiv=beamconvergence

! set the parameters
  STEM%BFradius = BFradius
  STEM%ADFinnerradius = ADFinnerradius
  STEM%ADFouterradius = ADFouterradius
  STEM%kt = kt
  STEM%cameralength = cameralength
  STEM%beamconvergence = beamconvergence
  if (numCL.ne.0) then 
    STEM%numCL = numCL
    STEM%CLarray = CLarray
  else
    STEM%numCL = 1
    STEM%CLarray = cameralength  
  end if
! make sure the number of s values is an odd number
  if (mod(numberofsvalues,2).eq.0) numberofsvalues=numberofsvalues+1
  STEM%numberofsvalues = numberofsvalues
  STEM%diffapmrad = diffaprad
  STEM%diffapmcenter = diffapcenter
  STEM%weightoutput = weightoutput
  STEM%geometry = geometry
  if (PRESENT(numk)) then
    STEM%numk = numk
  else
    STEM%numk = STEM%numberofsvalues
  end if

! and initialize all other STEM related arrays 
if (.not.PRESENT(beamdiv)) then
  if (geometry.eq.'SR') then
    call init_STEM(nn,g)
  else
    call init_STEM_ZA(nn)
  end if
end if 

end subroutine read_STEM_data

end module STEMmodule

