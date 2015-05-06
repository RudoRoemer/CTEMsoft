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
! EMsoft:EMoverlap.f90
!--------------------------------------------------------------------------
!
! PROGRAM: EMoverlap
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief dynamical diffraction patterns for overlapping crystals
!
!> @date  12/20/13 MDG 1.0 new program, partially based on EMorient
!--------------------------------------------------------------------------
program EMoverlap

use local
use files
use io

IMPLICIT NONE

character(fnlen)	:: nmldeffile

! deal with the command line arguments, if any
nmldeffile = 'EMoverlap.nml'
progname = 'EMoverlap.f90'
progdesc = 'Dynamical diffraction for overlapping crystals'

call EMsoft(progname, progdesc)

call Interpret_Program_Arguments(nmldeffile,2,(/ 0, 50 /) )

! perform the zone axis computations
call ComputeOverlapPattern(nmldeffile)

end program EMoverlap

!--------------------------------------------------------------------------
!
! SUBROUTINE:ComputeOverlapPattern
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute a dynamical diffraction pattern for two overlapping crystals
!
!> @detail The first part of the program is very much like that of EMorient, namely
!> the computation of the transformation matrix for the specified orientation relation.
!> This is then followed by the computattion of scattered amplitudes for the top
!> crystal for a number of different foil thickness values.  Then, for each of these
!> scattered beams and each sample thickness, we compute the scattered amplitudes for 
!> crystal B which results inthe intensities of all the double diffracted beams.
!
!> @param nmlfile namelist file name
!
!< @date 12/20/13  MDG 1.0 initial version
!--------------------------------------------------------------------------
subroutine ComputeOverlapPattern(nmlfile)

use symmetryvars
use symmetry
use crystalvars
use crystal
use constants
use gvectors
use kvectors
use error
use io
use math
use local
use files
use diffraction
use multibeams
use dynamical
use timing


IMPLICIT NONE

character(fnlen),INTENT(IN)    :: nmlfile

type(unitcell)    		:: cellA, cellB
type(orientation) 		:: orel
real(kind=sgl)                 :: tA(3), tB(3), gA(3), gB(3), c1, c2, io_real(6), TT(3,3), convergence, thetac, &
                                  galen, ktmax, bragg, thickinc, maxthick, voltage, dmin, maxg, gl, kB(3), fnB(3), &
                                  kc(3), hklA(3), klaue(2), newkA(3)
real(kind=dbl)			:: kstar(3), TTinv(3,3)
integer(kind=irg)              :: i, j, nkt, skip, dgn, isym, pgnum, io_int(6), gga(3), ggb(3),  npx, npy, numk, ijmax, &
                                  sLUT(3), imh, imk, iml, ix, iy, iz, numthick, istat, dp, numg, sdm(2), maxHOLZ, kA(3), fnA(3),&
                                  NBeamsA, icnt, NBeamsB, jj, nstart
character(3)                   :: method
character(fnlen)               :: xtalnameA, xtalnameB, outname
complex(kind=dbl)              :: pre, czero, cone, DUgpA

complex(kind=dbl),allocatable  :: LUTA(:,:,:), amplitudeA(:,:), amplitudeB(:,:),ScatMat(:,:), amp(:)
real(kind=sgl),allocatable	 :: thickA(:), thickB(:), inten(:,:), k0pg(:,:)
type(reflisttype),pointer 	 :: reflistA

namelist /overlaplist/ stdout, xtalnameA, xtalnameB, voltage, kA, fnA, dmin, tA, tB, gA, gB, outname, &
                       numthick, thickinc, convergence, nkt, maxg, klaue

! since we are passing a pointer (reflistA) to a subroutine, we must declare that routine here
! with an appropriate interface.
interface
  recursive subroutine Compute_A_DynMat(LUTA, imh, imk, iml, reflistA, kk)

  use local
  use dynamical
  use error
  use constants
  use crystal
  use diffraction
  use io
  use gvectors

  complex(kind=dbl)	        :: LUTA(-imh:imh,-imk:imk,-iml:iml)
  integer(kind=irg)	        :: imh, imk, iml
  type(reflisttype),pointer    :: reflistA
  real(kind=dbl)	        :: kk(3)
  end subroutine Compute_A_DynMat
end interface


! set the input parameters to default values (except for xtalnameA and B, which must be present)
xtalnameA = 'undefined'	         ! initial value to check that the keyword is present in the nml file
xtalnameB = 'undefined'	         ! initial value to check that the keyword is present in the nml file
stdout = 6			         ! standard output
voltage = 200000.0		         ! acceleration voltage [V]
kA = (/ 0, 0, 1 /)		         ! beam direction [direction indices]
fnA = (/ 0, 0, 1 /)		         ! foil normal [direction indices]
klaue = (/ 0.0, 0.0 /)                  ! Laue center position for crystal A
dmin = 0.025			         ! smallest d-spacing to include in dynamical matrix [nm]
tA = (/ 0.0, 0.0, 1.0 /)                ! direction in crystal A
tB = (/ 0.0, 0.0, 1.0 /)                ! parallel direction in crystal B
gA = (/ 1.0, 0.0, 0.0 /)                ! reciprocal vector in crystal A
gB = (/ 1.0, 0.0, 0.0 /)                ! parallel reciprocal vector in crystal B
maxg = 10.0                             ! max g-distance to consider for output
numthick = 10                           ! number of thickness increments, determines total foil thickness
thickinc = 10.0                         ! increment size for crystal A; B will be maxthick - # thickinc
convergence = 5.0                       ! beam convergence angle [mrad] 
nkt = 0                                 ! number of steps along incident cone radius (0 means only one incident beam direction)

outname = 'overlap.data'                ! output filename

pre = dcmplx(0.D0,1.D0)*cPi	! i times pi
czero = dcmplx(0.D0,0.D0)
cone = dcmplx(1.D0,0.D0)

! read the namelist file
open(UNIT=dataunit,FILE=nmlfile,DELIM='apostrophe',STATUS='old')
read(UNIT=dataunit,NML=overlaplist)
close(UNIT=dataunit,STATUS='keep')

if (trim(xtalnameA).eq.'undefined') then
  call FatalError('ComputeOverlapPattern:',' structure file name for crystal A is undefined in '//nmlfile)
end if

if (trim(xtalnameB).eq.'undefined') then
  call FatalError('ComputeOverlapPattern:',' structure file name for crystal B is undefined in '//nmlfile)
end if

! print some information
 progname = 'EMoverlap.f90'
 progdesc = 'Dynamical diffraction pattern for overlapping crystals'
 call EMsoft

! first get the crystal A data 
 SG%SYM_reduce=.TRUE.
 call CrystalData(xtalnameA)
! store Crystal A matrices
 cellA = cell

! first get the crystal B data 
 call CrystalData(xtalnameB)
! store Crystal B matrices
 cellB = cell

! get orientation relation
 orel%tA = tA
 orel%tB = tB
 orel%gA = gA
 orel%gB = gB

! check for orthonormality using zone equation
 c1=sum(orel%tA*orel%gA)
 if (c1.ne.0.0_sgl) then
  mess = 'Plane does not contain direction (crystal A)'; call Message("(A)")
  STOP
 end if
 c2=sum(orel%tB*orel%gB)
 if (c2.ne.0.0_sgl) then
  mess = 'Plane does not contain direction (crystal B)'; call Message("(A)")
  STOP
 end if

! print the OR
 mess = 'Orientation relation between A and B'; call Message("(A)")
 io_real(1:3) = tA(1:3)
 io_real(4:6) = tB(1:3)
 call WriteValue('Parallel directions: ', io_real, 6, "('[',2(F6.3,','),F6.3,'] // [',2(F6.3,','),F6.3,']')")
 io_real(1:3) = gA(1:3)
 io_real(4:6) = gB(1:3)
 call WriteValue('Parallel planes:     ', io_real, 6, "('(',2(F6.3,','),F6.3,') // (',2(F6.3,','),F6.3,')')")

! compute the OR transformation matrix for B the new lattice and A the old
 TT = ComputeOR(orel, cellA, cellB, 'BA')
 call mInvert(dble(TT),TTinv,.FALSE.)
 
!==========================================
!==========================================
! initialize the first crystal and perform 
! a CBED computation for the entire incident
! beam cone.  We'll store the beam amplitudes
! for each beam orientation and thickness; we
! will only do this for ZOLZ reflections for now,
! and we keep the number of beams fixed to the
! number for the exact zone axis orientation.
! To test all of this we will have only one
! incident beam orientation for crystal A.
!==========================================
!==========================================
 cell = cellA

! set the voltage parameter
 skip = 3
 call CalcWaveLength(dble(voltage),skip)
 
! generate all atom positions
 call CalcPositions('v')
 
! get k and f
 DynFN = float(fnA)
 
! determine the point group number
 j=0
 do i=1,32
  if (SGPG(i).le.cell % SYM_SGnum) j=i
 end do

! use the new routine to get the whole pattern 2D symmetry group, since that
! is the one that determines the independent beam directions.
 dgn = GetPatternSymmetry(kA,j,.TRUE.)
 pgnum = j
 isym = WPPG(dgn) ! WPPG lists the whole pattern point group numbers vs. diffraction group numbers

! determine the shortest reciprocal lattice points for this zone
 call ShortestG(kA,gga,ggb,isym)
 io_int(1:3)=gga(1:3)
 io_int(4:6)=ggb(1:3)
 call WriteValue('Reciprocal lattice vectors in A : ', io_int, 6,"('(',3I3,') and (',3I3,')',/)")

! initialize the HOLZ geometry type
 call GetHOLZGeometry(float(gga),float(ggb),kA,fnA) 

! construct the list of all possible reflections
 method = 'ALL'
 thetac = convergence/1000.0
 call Compute_ReflectionList(dmin,kA,gga,ggb,method,.FALSE.,maxHOLZ,thetac)
 galen = CalcLength(float(gga),'r')

! determine range of incident beam directions
 bragg = CalcDiffAngle(gga(1),gga(2),gga(3))*0.5
  
! convert to ktmax along ga
 ktmax = 0.5*thetac/bragg
  
! set parameters for wave vector computation
 ijmax = 1   ! truncation value for beam directions

! there's no point in using symmetry because we use only a single incident beam direction.
 isym = 1
 npx = 0
 npy = 0
! incorporate the Laue center into the incident wave vector...
 newkA = klaue(1)*gga + klaue(2)*ggb + float(kA) 
 call Calckvectors(dble(newkA),dble(gga),dble(ktmax),npx,npy,numk,isym,ijmax,'Standard',.FALSE.)
 io_int(1) = numk
 call WriteValue('Number of incident beam directions : ', io_int, 1, "(I5)")

! we'll have to use a special routine to create the dynamical matrix; no need to use Bethe potentials
! at this point in time, since we're only computing the dynamical matrix a handful of times.  

 sLUT = shape(LUT)
 imh = (sLUT(1)-1)/2
 imk = (sLUT(2)-1)/2
 iml = (sLUT(3)-1)/2
! deallocate(LUT)
 allocate(LUTA(-imh:imh,-imk:imk,-iml:iml),stat=istat)
 do ix=-imh,imh
  do iy=-imk,imk
   do iz=-iml,iml
     call CalcUcg( (/ ix, iy, iz /) ) ! compute potential Fourier coefficient
     LUTA(ix, iy, iz) = pre*rlp%qg
   end do 
  end do 
 end do  
 NBeamsA = DynNbeamsLinked

! in addition we need to keep the reflection linked lists for both phases, so we'll assign the 
! head of the reflection list for the A phase to the pointer reflistA 
 reflistA => reflist%next 
  
! get the absorption coefficient
 call CalcUcg( (/0,0,0/) )
 DUgpA = pre*rlp%qg
 
! determine the thickness list for A and B
 maxthick = numthick * thickinc
 allocate( thickA(numthick), thickB(numthick) )
 do i=1,numthick
   thickA(i) = float(i-1)*thickinc
   thickB(i) = maxthick - thickA(i)
 end do
 
!==========================================
!========================================== 
! we're done with the initializations for A; let's start the image computation.
! we need to first determine how many reflections we're going to keep.
 cell = cellA
! go through the reflistA and count/tag the ones that belong to the ZOLZ   
 rltmpa  => reflistA
 rltmpa%famnum = 1	! we'll use this field to identify the reflections that will be in the output file
 numg = 1
! mess = 'The following reflections are considered in the ZOLZ of crystal A:' 
! call Message("(A)")
! io_int(1:4) = (/numg, rltmpa%hkl(1),rltmpa%hkl(2),rltmpa%hkl(3) /)
! call WriteValue('', io_int, 4, "(4I5)")
 rltmpa => rltmpa%next
 do
   if (.not.associated(rltmpa)) EXIT
   dp = rltmpa%hkl(1) * kA(1) + rltmpa%hkl(2) * kA(2) + rltmpa%hkl(3) * kA(3)
   gl = CalcLength(float(rltmpa%hkl),'r')
   if ((dp.eq.0).and.(gl.le.maxg)) then
     numg = numg+1
     rltmpa%famnum = numg
!     io_int(1:4) = (/numg, rltmpa%hkl(1),rltmpa%hkl(2),rltmpa%hkl(3) /)
!     call WriteValue('', io_int, 4, "(4I5)")
   else
     rltmpa%famnum = 0
   end if
!   if (associated(rltmpa%next)) rltmpa => rltmpa%next 
   rltmpa => rltmpa%next
 end do


! now we can allocate the output amplitude array for crystal A
 allocate(amplitudeA(numthick,numg),stat=istat)
 amplitudeA = czero

 nstart = numg

! so the next step is to fill the amplitudeA array with the correct 
! amplitudes, computed by means of the scattering matrix formalism 
! first, we get the dynamical matrix for the main incident beam direction 
 ktmp => khead
 call Compute_A_DynMat(LUTA, imh, imk, iml, reflistA, ktmp%k)

! next, convert this to the scattering matrix
 sdm = shape(DynMat)
 allocate(ScatMat(sdm(1),sdm(2)))
 allocate(amp(sdm(1)))

! and iterate over the thickness values
 do i=1,numthick
   ScatMat = czero
   call MatrixExponential(DynMat, ScatMat, dble(thickA(i)), 'Pade', sdm(1))
   amp = czero
   amp(1) = cone
   amp = matmul( ScatMat, amp )
   
! extract the amplitudes that we want to use as incident beams for the crystal B
   icnt = 1
   do j=1,NBeamsA
    if (BetheParameter%reflistindex(j).ne.-1) then
       amplitudeA(i,BetheParameter%reflistindex(j)) = amp(icnt)
       icnt = icnt+1
    end if  
   end do 

 end do
 
! before moving on to crystal B, we need to store the geometrical information on the 
! relevant reflections, and the intensities, and then we'll reset the reflection list
  allocate(inten(numthick,numg))
  inten = cdabs(amplitudeA)**2 ! MNS from cabs to cdabs for gfortran
  allocate(k0pg(3,numg))
  rltmpa  => reflistA
  numg = 1
! open the output file
  open(unit=dataunit,file='crystalA.data', status='unknown', form='unformatted')
  write (dataunit) numk, numthick
  write (dataunit) nstart

  do
   if (.not.associated(rltmpa)) EXIT
   if (rltmpa%famnum.ne.0) then
     hklA = float(rltmpa%hkl)
!     call NormVec(hklA,'r')
     k0pg(1:3,numg) = matmul( TT,ktmp%k + hklA - klaue(1)*gga - klaue(2)*ggb )
     hklA = matmul(TTinv, k0pg(1:3,numg) )
     write (dataunit) hklA
     numg = numg+1
   end if
   rltmpa => rltmpa%next 
  end do
  numg = numg-1
  write (dataunit) inten

! and clean up the k-vector list
  call Delete_kvectorlist()

  close(dataunit,status='keep')

!==========================================
!==========================================
! switch to cell B
!==========================================
!==========================================
 mess = 'Starting to work on crystal B'; call Message("(//A/)")
 cell = cellB

! here we need to loop over all the non-zero amplitude beams from crystal A;
! for each such beam gA, there is now an incident beam with wave vector k_0 + gA
! first we need to get the list of such vectors, each of them converted into 
! the reference frame of crystal B by the OR transformation matrix TT; then we
! step through the list and perform a scattering matrix computation as a function
! of the thickness of B.  For each gA, we thus end up with a list of intensities
! for the locations k_0+gA+gB; we'll need to store those projected locations, as well
! as their intensities for each thickness value.

! first we need to redefine the electron wavelength, since the refraction correction
! in the second crystal might be different from that in the first.
! set the voltage parameter
 skip = 3
 call CalcWaveLength(dble(voltage),skip)

! generate all atom positions for crystal B
 call CalcPositions('v')

! we need to generate the list of B reflections, based on the incident vector k0 
! expressed in the B reference frame; we also need the foil normal expressed in B
 kB(1:3) = k0pg(1:3,1)
 io_real(1:3) = kB(1:3)
 call WriteValue('Main incident beam direction in crystal B (reciprocal) : ', io_real, 3, "(3F10.4)")
 fnB = matmul(TT,float(fnA))
 DynFN = fnB
 io_real(1:3) = fnB(1:3)
 call WriteValue('Foil normal direction in crystal B : ', io_real, 3, "(3F10.4)")

! delete the existing reflection list
 call Delete_gvectorlist()
 deallocate(LUTA)
 if (associated(reflist)) nullify(reflist)

! construct the list of all possible reflections
 call Compute_ReflectionListB(dmin,kB)
 NBeamsB = DynNbeamsLinked 

! and create the wave vector list (we'll need to do this in a different 
! way from the usual Calkvectors routine, since now we need the transformed k0+g as 
! incident beam directions. We'll do the pointer allocations explicitly 
! at this point; this may become a new library subroutine later on...
! compute geometrical factors 
 kstar = kB
 kc = (/ 0.0, 0.0, 0.0 /)
 call NormVec(kstar,'r')                       	! normalize reciprocal beam vector

! allocate the head and tail of the linked list
 allocate(khead,stat=istat)   				! allocate new value
 if (istat.ne.0) call FatalError('ComputeOverlapPattern','unable to allocate khead pointer')
 ktail => khead                      			! tail points to new value
 nullify(ktail%next)                			! nullify next in new value
 numk = 1                          			! keep track of number of k-vectors so far
 ktail%k = kstar/cell%mLambda				! divide by wavelength


! ok, we have the double transmitted beam; next we need all the other
! diffracted beams from A that have a non-zero amplitude.
 do i=2,numg
!  if (sum(inten(1:numthick,i)).ne.0.0) then
    allocate(ktail%next,stat=istat)  		! allocate new value
    if (istat.ne.0) call FatalError('ComputeOverlapPattern:',' unable to allocate pointer')
    ktail => ktail%next               		! tail points to new value
    nullify(ktail%next)              		! nullify next in new value
    numk = numk + 1                 		! keep track of number of k-vectors so far
    kc(1:3) = k0pg(1:3,i)
    call NormVec(kc,'r')
    ktail%k(1:3) = kc/cell%mLambda
!  end if
 end do

 io_int(1) = numk
 call WriteValue('Number of incident beam directions for B : ', io_int, 1, "(I5)") 
 
! next we need to loop over all the reflections with non-zero amplitude from crystal A;
! for each of these, we compute the dynamical matrix and integrate the equations to get
! intensities; then we store the tangential components of the diffracted beams with 
! respect to the original ga and gb from crystal A, along with the intensities vs, thickness. 

! open the output file
 open(unit=dataunit,file=trim(outname), status='unknown', form='unformatted')
 write (dataunit) numk, numthick
 write (dataunit) thickinc

 if (allocated(ScatMat)) deallocate(ScatMat)
 if (allocated(inten)) deallocate(inten)
 if (allocated(amp)) deallocate(amp)
 if (allocated(amplitudeB)) deallocate(amplitudeB)
 
! and perform the loop over all wave vectors
 ktmp => khead
 do i=1,numk
   if (mod(i,10).eq.0) then
     io_int(1) = i
     call WriteValue(' Working on pattern ', io_int, 1, "(I5)")
   end if
   
   call Compute_B_DynMat(ktmp%k)
 
! next, convert this to the scattering matrix
   sdm = shape(DynMat)

   allocate(ScatMat(sdm(1),sdm(2)), amplitudeB(numthick,DynNbeams))
   allocate(amp(sdm(1)))
   allocate(inten(numthick,DynNbeams))

! and iterate over the thickness values
   do j=1,numthick
     ScatMat = czero
     call MatrixExponential(DynMat, ScatMat, dble(thickB(j)), 'Pade', sdm(1))
     amp = czero
     amp(1) = amplitudeA(j,i)
     amp = matmul( ScatMat, amp )
     
! extract the amplitudes that we want to use as incident beams for the crystal B
     icnt = 1
     do jj=1,NBeamsB
       if (BetheParameter%stronglist(jj).ne.0) then
        amplitudeB(j,BetheParameter%stronglist(jj)) = amp(icnt)
        icnt = icnt+1
       end if  
     end do 
   end do 

! convert to intensities
   inten = cdabs(amplitudeB)**2 ! MNS from cabs to cdabs for gfortran

! and store the relevant stuff in the output file   
! we need to store all the g-vectors, the coordinates of k0+g+h, and the intensities
! we'll do all this in the reference frame of crystal A, so we need the inverse of TT
   write (dataunit) DynNbeams
   icnt = 1
   do j=1,NBeamsB
     if (BetheParameter%stronglist(j).ne.0) then
       hklA = BetheParameter%stronghkl(1:3,icnt)
!       call NormVec(hklA,'r')
       hklA = matmul(TTinv, ktmp%k + hklA - klaue(1)*gga - klaue(2)*ggb )
       write (dataunit) hklA
       icnt = icnt+1
     end if
   end do
!   write (dataunit) inten
   write (dataunit) amplitudeB

   deallocate(ScatMat, amplitudeB, amp, inten)
   ktmp => ktmp%next
 end do


close (dataunit, status='keep')
 

end subroutine ComputeOverlapPattern


!--------------------------------------------------------------------------
!
! SUBROUTINE:Compute_A_DynMat
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Computes the dynamical matrix for crystal A
!
!> @details This is a reduced and adapted version of the Compute_DynMat routine. 
!
!> @param LUTA reflection lookup table for phase A
!> @param imh h dimension parameter for LUTA
!> @param imk k dimension parameter for LUTA
!> @param iml l dimension parameter for LUTA
!> @param reflistA reflist pointer
!> @param kk incident wave vector
!
!> @date 11/02/13  MDG 1.0 original
!> @date 11/05/13  MDG 1.1 added variant handling
!> @date 12/20/13  MDG 1.2 adapted from EMgamma program
!--------------------------------------------------------------------------
recursive subroutine Compute_A_DynMat(LUTA, imh, imk, iml, reflistA, kk)

use local
use dynamical
use error
use constants
use crystal
use diffraction
use io
use gvectors

IMPLICIT NONE

complex(kind=dbl),INTENT(IN)		:: LUTA(-imh:imh,-imk:imk,-iml:iml)
integer(kind=irg),INTENT(IN)		:: imh, imk, iml
type(reflisttype),pointer 	        :: reflistA
real(kind=dbl),INTENT(IN)		:: kk(3)		!< incident wave vector 

complex(kind=dbl)  			:: czero,pref
integer(kind=irg) 		 	:: istat,ir,ic,nn, istrong, ig, ll(3)
real(kind=sgl)     			:: gg(3), sgp, lUg, cut1

! has the list of reflections been allocated ?
if (.not.associated(reflistA)) call FatalError('Compute_A_DynMat',' reflection list has not been allocated')

! if the dynamical matrices have already been allocated, deallocate them first
! this is partially so that no program will allocate DynMats itself; it must be done
! via this routine only.
if (allocated(DynMat)) deallocate(DynMat)

! initialize some parameters
pref = dcmplx(0.0,1.0)*cPi	! i times pi
czero = dcmplx(0.0,0.0)	! complex zero

! we don't know yet how many strong reflections there are so we'll need to determine this first
! this number depends on some externally supplied parameters, which we will get from a namelist
! file (which should be read only once by the Set_Bethe_Parameters routine), or from default values
! if there is no namelist file in the folder.
if (BetheParameter%cutoff.eq.0.0) call Set_Bethe_Parameters(.TRUE.)
BetheParameter%weakcutoff = BetheParameter%cutoff ! no Bethe potentials in the current approach

! reset the value of DynNbeams in case it was modified in a previous call 
DynNbeams = DynNbeamsLinked
  	
if (.not.allocated(BetheParameter%stronglist)) allocate(BetheParameter%stronglist(DynNbeams))
if (.not.allocated(BetheParameter%reflistindex)) allocate(BetheParameter%reflistindex(DynNbeams))
BetheParameter%stronglist = 0
BetheParameter%reflistindex = -1

rltmpa => reflistA


! deal with the transmitted beam first
    nn = 1		! nn counts all the scattered beams that satisfy the cutoff condition
    rltmpa%sg = 0.D0    
    BetheParameter%stronglist(1) = 1
    BetheParameter%reflistindex(1) = 1

! loop over all reflections in the linked list    
    rltmpa => rltmpa%next
    reflectionloop: do ig=2,DynNbeamsLinked
      gg = float(rltmpa%hkl)        		! this is the reciprocal lattice vector 
      rltmpa%sg = Calcsg(gg,sngl(kk),DynFN)

! use the reflection num entry to indicate whether or not this
! reflection should be used for the dynamical matrix
! We compare |sg| with a multiple of lambda |Ug|
!
!  |sg|>cutoff lambda |Ug|   ->  don't count reflection
!  cutoff lambda |Ug| > |sg|  -> strong reflection
!
      sgp = abs(rltmpa%sg) 
      lUg = cdabs(rltmpa%Ucg) * cell%mLambda !MNS from cabs to cdabs for gfortran
      cut1 = BetheParameter%cutoff * lUg

      if (sgp.le.cut1) then  ! count this beam
	nn = nn+1
	BetheParameter%stronglist(ig) = nn
	if (rltmpa%famnum.ne.0) then
	  BetheParameter%reflistindex(ig) = rltmpa%famnum
	else
	  BetheParameter%reflistindex(ig) = -1
	end if
      end if

! go to the next beam in the list
      rltmpa => rltmpa%next
    end do reflectionloop

! if we don't have any beams in this list (unlikely, but possible if the cutoff
! parameter has an unreasonable value) then we abort the run
! and we report some numbers to the user 
if (nn.eq.0) then
   mess = ' no beams found for the following parameters:'; call Message("(A)")
   write (stdout,*) ' wave vector = ',kk,'  -> number of beams = ',nn
   mess =  '   -> check cutoff and weakcutoff parameters for reasonableness'; call Message("(A)")
   call FatalError('Compute_A_DynMat','No beams in list')
end if

! next, we define nns to be the number of strong beams.
BetheParameter%nns = nn
	 	
! allocate arrays for strong beam information
if (allocated(BetheParameter%stronghkl)) deallocate(BetheParameter%stronghkl)
if (allocated(BetheParameter%strongsg)) deallocate(BetheParameter%strongsg)
if (allocated(BetheParameter%strongID)) deallocate(BetheParameter%strongID)
allocate(BetheParameter%stronghkl(3,BetheParameter%nns),BetheParameter%strongsg(BetheParameter%nns))
allocate(BetheParameter%strongID(BetheParameter%nns))

BetheParameter%stronghkl = 0
BetheParameter%strongsg = 0.0
BetheParameter%strongID = 0


! here's where we extract the relevant information from the linked list (much faster
! than traversing the list each time...)
rltmpa => reflistA    ! reset the a list
istrong = 0
do ir=1,DynNbeamsLinked
     if (BetheParameter%stronglist(ir).gt.0) then
        istrong = istrong+1
        BetheParameter%stronghkl(1:3,istrong) = rltmpa%hkl(1:3)
        BetheParameter%strongsg(istrong) = rltmpa%sg
! make an inverse index list
	BetheParameter%strongID(istrong) = ir		
     end if
   rltmpa => rltmpa%next
end do


! now we are ready to create the dynamical matrix
DynNbeams = BetheParameter%nns

! allocate DynMat and set to complex zero
allocate(DynMat(DynNbeams,DynNbeams),stat=istat)
DynMat = czero

! ir is the row index
do ir=1,BetheParameter%nns
! ic is the column index
  do ic=1,BetheParameter%nns
! compute the Fourier coefficient of the electrostatic lattice potential 
    if (ic.ne.ir) then  ! not a diagonal entry
      ll = BetheParameter%stronghkl(1:3,ir) - BetheParameter%stronghkl(1:3,ic)
      DynMat(ir,ic) = LUTA(ll(1),ll(2),ll(3)) 
    else  ! it is a diagonal entry, so we need the excitation error and the absorption length	
      DynMat(ir,ir) = cmplx(0.0,2.D0*cPi*BetheParameter%strongsg(ir),dbl)+LUTA(0,0,0)
    end if
  end do
end do
! that should do it for the initialization of the dynamical matrix

end subroutine Compute_A_DynMat


!--------------------------------------------------------------------------
!
! SUBROUTINE:Compute_B_DynMat
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Computes the dynamical matrix for crystal B
!
!> @details This is a reduced and adapted version of the Compute_DynMat routine. 
!
!> @param kk incident wave vector
!
!> @date 11/02/13  MDG 1.0 original
!> @date 11/05/13  MDG 1.1 added variant handling
!> @date 12/20/13  MDG 1.2 adapted from EMgamma program
!--------------------------------------------------------------------------
recursive subroutine Compute_B_DynMat(kk)

use local
use dynamical
use error
use constants
use crystal
use diffraction
use io
use gvectors

IMPLICIT NONE

real(kind=dbl),INTENT(IN)		:: kk(3)		!< incident wave vector 

complex(kind=dbl)  			:: czero,pref
integer(kind=irg) 		 	:: istat,ir,ic,nn, istrong, ig, ll(3)
real(kind=sgl)     			:: gg(3), sgp, lUg, cut1

! has the list of reflections been allocated ?
if (.not.associated(reflist)) call FatalError('Compute_B_DynMat',' reflection list has not been allocated')

! if the dynamical matrices have already been allocated, deallocate them first
! this is partially so that no program will allocate DynMats itself; it must be done
! via this routine only.
if (allocated(DynMat)) deallocate(DynMat)

! initialize some parameters
pref = dcmplx(0.0,1.0)*cPi	! i times pi
czero = dcmplx(0.0,0.0)	! complex zero

! we don't know yet how many strong reflections there are so we'll need to determine this first
! this number depends on some externally supplied parameters, which we will get from a namelist
! file (which should be read only once by the Set_Bethe_Parameters routine), or from default values
! if there is no namelist file in the folder.
if (BetheParameter%cutoff.eq.0.0) call Set_Bethe_Parameters(.TRUE.)
BetheParameter%weakcutoff = BetheParameter%cutoff ! no Bethe potentials in the current approach

! reset the value of DynNbeams in case it was modified in a previous call 
DynNbeams = DynNbeamsLinked

if (allocated(BetheParameter%stronglist)) deallocate(BetheParameter%stronglist)
if (allocated(BetheParameter%reflistindex)) deallocate(BetheParameter%reflistindex)
if (.not.allocated(BetheParameter%stronglist)) allocate(BetheParameter%stronglist(DynNbeams))
if (.not.allocated(BetheParameter%reflistindex)) allocate(BetheParameter%reflistindex(DynNbeams))
BetheParameter%stronglist = 0
BetheParameter%reflistindex = -1

rltmpa => reflist%next

! deal with the transmitted beam first
    nn = 1		! nn counts all the scattered beams that satisfy the cutoff condition
    rltmpa%sg = 0.D0    
    BetheParameter%stronglist(1) = 1
!    BetheParameter%reflistindex(1) = 1

! loop over all reflections in the linked list    
    rltmpa => rltmpa%next
    reflectionloop: do ig=2,DynNbeamsLinked
      gg = float(rltmpa%hkl)        		! this is the reciprocal lattice vector 
      rltmpa%sg = Calcsg(gg,sngl(kk),DynFN)

! use the reflection num entry to indicate whether or not this
! reflection should be used for the dynamical matrix
! We compare |sg| with a multiple of lambda |Ug|
!
!  |sg|>cutoff lambda |Ug|   ->  don't count reflection
!  cutoff lambda |Ug| > |sg|  -> strong reflection
!
      sgp = abs(rltmpa%sg) 
      lUg = cdabs(rltmpa%Ucg) * cell%mLambda ! MNS from cabs to cdabs for gfortran
      cut1 = BetheParameter%cutoff * lUg

      if (sgp.le.cut1) then  ! count this beam
	nn = nn+1
	BetheParameter%stronglist(ig) = nn
      end if

! go to the next beam in the list
      if (associated(rltmpa%next)) rltmpa => rltmpa%next
    end do reflectionloop

! if we don't have any beams in this list (unlikely, but possible if the cutoff
! parameter has an unreasonable value) then we abort the run
! and we report some numbers to the user 
if (nn.eq.0) then
   mess = ' no beams found for the following parameters:'; call Message("(A)")
   write (stdout,*) ' wave vector = ',kk,'  -> number of beams = ',nn
   mess =  '   -> check cutoff and weakcutoff parameters for reasonableness'; call Message("(A)")
   call FatalError('Compute_B_DynMat','No beams in list')
end if

! next, we define nns to be the number of strong beams.
BetheParameter%nns = nn

	 	
! allocate arrays for strong beam information
if (allocated(BetheParameter%stronghkl)) deallocate(BetheParameter%stronghkl)
if (allocated(BetheParameter%strongsg)) deallocate(BetheParameter%strongsg)
if (allocated(BetheParameter%strongID)) deallocate(BetheParameter%strongID)
allocate(BetheParameter%stronghkl(3,BetheParameter%nns),BetheParameter%strongsg(BetheParameter%nns))
allocate(BetheParameter%strongID(BetheParameter%nns))

BetheParameter%stronghkl = 0
BetheParameter%strongsg = 0.0
BetheParameter%strongID = 0


! here's where we extract the relevant information from the linked list (much faster
! than traversing the list each time...)
rltmpa => reflist%next    ! reset the a list

istrong = 0
do ir=1,DynNbeamsLinked
     if (BetheParameter%stronglist(ir).gt.0) then
        istrong = istrong+1
        BetheParameter%stronghkl(1:3,istrong) = rltmpa%hkl(1:3)
        BetheParameter%strongsg(istrong) = rltmpa%sg
! make an inverse index list
	BetheParameter%strongID(istrong) = ir		
     end if
     if (associated(rltmpa%next)) rltmpa => rltmpa%next
end do


! now we are ready to create the dynamical matrix
DynNbeams = BetheParameter%nns

! allocate DynMat and set to complex zero
allocate(DynMat(DynNbeams,DynNbeams),stat=istat)
DynMat = czero

! ir is the row index
do ir=1,BetheParameter%nns
! ic is the column index
  do ic=1,BetheParameter%nns
! compute the Fourier coefficient of the electrostatic lattice potential 
    if (ic.ne.ir) then  ! not a diagonal entry
      ll = BetheParameter%stronghkl(1:3,ir) - BetheParameter%stronghkl(1:3,ic)
      DynMat(ir,ic) = LUT(ll(1),ll(2),ll(3)) 
    else  ! it is a diagonal entry, so we need the excitation error and the absorption length	
      DynMat(ir,ir) = cmplx(0.0,2.D0*cPi*BetheParameter%strongsg(ir),dbl)+LUT(0,0,0)
    end if
  end do
end do
! that should do it for the initialization of the dynamical matrix


end subroutine Compute_B_DynMat




!--------------------------------------------------------------------------
!
! SUBROUTINE: Compute_ReflectionListB
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute the entire reflection list for crystal B
!
!> @details also computes the LUT (LookUpTable) that stores all the scattering
!> potential coefficients that are needed to fill the dynamical matrix (only for the 
!> general case).  In this case, we include all reflections and let the DynMat 
!> routine decide which ones to consider explicitly.
!
!> @param dmin minimum d-spacing to allow in the list
!> @param k incident wave vector
!
!> @date 12/23/13 MDG 1.0 original (based on Compute_ReflectionList in gvectors module)
!--------------------------------------------------------------------------
subroutine Compute_ReflectionListB(dmin,k)

use local
use io
use crystal
use constants
use dynamical
use diffraction
use gvectors
use symmetry

IMPLICIT NONE

real(kind=sgl),INTENT(IN)			:: dmin
real(kind=sgl),INTENT(IN)			:: k(3)

integer(kind=irg)				:: imh, imk, iml, gg(3), ix, iy, iz,  istat, numr, ir
real(kind=sgl)					:: dhkl, io_real(9), ddt
integer(kind=irg)				:: io_int(3)
complex(kind=dbl)				:: pre, czero



! set threshold for double diffraction detection
  ddt = 1.0e-10
  pre = dcmplx(0.D0,1.D0)*cPi	! i times pi
  czero = dcmplx(0.D0,0.D0)

! The master list is easily created by brute force
 imh = 1
 do 
   imh = imh + 1
   dhkl = 1.0/CalcLength(  (/float(imh) ,0.0_sgl,0.0_sgl/), 'r')
   if (dhkl.lt.dmin) EXIT
 end do
 imk = 1
 do 
   imk = imk + 1
   dhkl = 1.0/CalcLength( (/0.0_sgl,float(imk),0.0_sgl/), 'r')
  if (dhkl.lt.dmin) EXIT
 end do
 iml = 1
 do 
   iml = iml + 1
   dhkl = 1.0/CalcLength( (/0.0_sgl,0.0_sgl,float(iml)/), 'r')
   if (dhkl.lt.dmin) EXIT
 end do
 io_int = (/ imh, imk, iml /)
 call WriteValue(' Range of reflections along a*, b* and c* = ',io_int,3)

! the LUT array stores all the Fourier coefficients, so that we only need to compute them once...
  if (allocated(LUT)) deallocate(LUT)
  allocate(LUT(-2*imh:2*imh,-2*imk:2*imk,-2*iml:2*iml),stat=istat)
write (*,*) 'zeroing LUT'
  LUT = czero
write (*,*) 'allocated LUT'

! allocate an array to keep track of all families of reflections
  if (allocated(refdone)) then
    write (*,*) 'refdone is already allocated; deallocating now'
    deallocate(refdone)
  end if
  allocate(refdone(-imh:imh,-imk:imk,-iml:iml),stat=istat)
  refdone = 0
write (*,*) 'allocated refdone'
 
! allocate an array that keeps track of potential double diffraction reflections
!if (allocated(dbdiffB)) deallocate(dbdiffB)
  allocate(dbdiffB(-2*imh:2*imh,-2*imk:2*imk,-2*iml:2*iml),stat=istat)
  dbdiffB = .FALSE.
!write (*,*) 'allocated dbdiffb'


  DynNbeams = 1
  gg = (/ 0,0,0 /)
  call AddReflection( gg )   ! this guarantees that 000 is always the first reflection
  rltail%famhkl = gg
  call CalcUcg(gg)   
  DynUpz = rlp%Vpmod
  io_real(1) = rlp%xgp
  call WriteValue(' Normal absorption length = ', io_real, 1)

! and add this reflection to the look-up table
  LUT(gg(1),gg(2),gg(3)) = pre*rlp%qg
  refdone(gg(1),gg(2),gg(3)) = 2 !.TRUE. 

! now do the same for the other allowed reflections
! note that the lookup table must be twice as large as the list of participating reflections,
! since the dynamical matrix uses g-h as its index !!!  However, the linked list of reflections
! should only contain the g, h reflections separately, not their differences. (i.e., smaller box)
ixl: do ix=-2*imh,2*imh
iyl:  do iy=-2*imk,2*imk
izl:   do iz=-2*iml,2*iml
        gg = (/ ix, iy, iz /)
        if (IsGAllowed(gg)) then

! if this g is inside the original box, then add it to the linked list
         if ((abs(ix).le.imh).and.(abs(iy).le.imk).and.(abs(iz).le.iml)) then 

! if we haven't visited this reflection yet
!	   if ( .not.refdone(ix,iy,iz) ) then
	   if ( refdone(ix,iy,iz) .eq. 0) then

! compute the family of reflections
	     call CalcFamily( gg, numr, 'r' )
	     do ir=1,numr
		 refdone(itmp(ir,1),itmp(ir,2),itmp(ir,3)) = 1	! and flag the reflection
 	     end do
! next add those reflections that are inrange
	     do ir=1,numr
		  call AddReflection( (/ itmp(ir,1),itmp(ir,2),itmp(ir,3) /) )
         	  LUT(itmp(ir,1),itmp(ir,2),itmp(ir,3)) = pre*rlp%qg
 		  refdone(itmp(ir,1),itmp(ir,2),itmp(ir,3)) = 2	! and flag the reflection
! flag this reflection as a double diffraction candidate if cabs(Ucg)<ddt threshold
	         if (cabs(rlp%Ucg).le.ddt) then 
         	  dbdiff(itmp(ir,1),itmp(ir,2),itmp(ir,3)) = .TRUE.
         	 end if      	 
	     end do
          end if 
         else ! not inside smaller box

! add the reflection to the look up table
 	   call CalcUcg( gg )
           LUT(ix, iy, iz) = pre*rlp%qg
! flag this reflection as a double diffraction candidate if cabs(Ucg)<ddt threshold
           if (cabs(rlp%Ucg).le.ddt) then 
             dbdiff(ix,iy,iz) = .TRUE.
           end if

         end if ! not inside smaller box

        end if ! IsGAllowed
       end do izl
      end do iyl
    end do ixl
  io_int(1) = DynNbeams
  DynNbeamsLinked = DynNbeams
  call WriteValue(' Length of the master list of reflections : ', io_int, 1, "(I8)")
  
end subroutine Compute_ReflectionListB

