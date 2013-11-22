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
! CTEMsoft2013:CTEMECP.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMECP 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Zone axis electron channeling patterns
!
!> @date 03/18/10 MDG 1.0 f90
!> @date 08/09/10 MDG 2.0 corrected weight factors and g-vector ordering problem
!> @date 11/18/13 MDG 3.0 major rewrite with new libraries 
!--------------------------------------------------------------------------
program CTEMECP

use local
use files
use io

IMPLICIT NONE

character(fnlen)			:: nmldeffile

! deal with the command line arguments, if any
nmldeffile = 'CTEMECP.nml'
progname = 'CTEMECP.f90'
call Interpret_Program_Arguments(nmldeffile,1,(/ 40 /) )

! perform the zone axis computations
call ECpattern(nmldeffile)

end program CTEMECP

!--------------------------------------------------------------------------
!
! SUBROUTINE:ECpattern
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute a large angle zone axis electron channeling pattern
!
!> @note This is really very similar to a LACBED computation, except that 
!> the final intensity computation is somewhat different.  We could in 
!> principle also include the Kossel pattern computation in this program.
!> This program now also includes the Bethe potential approximation, to 
!> hopefully speed things up a little bit...  For now, the Kossel pattern
!> calculation part has been commented out; needs some more work, but no 
!> time available at the moment (end of 2013)...
!
!> @param nmlfile namelist file name
!
!> @date 11/18/13  MDG 1.0 major rewrite from older ECP program
!--------------------------------------------------------------------------
subroutine ECpattern(nmlfile)


use symmetryvars
use symmetry
use crystalvars
use crystal
use constants
use gvectors
use kvectors
use error
use io
use local
use files
use diffraction
use multibeams
use dynamical
use timing

IMPLICIT NONE

character(fnlen),INTENT(IN)	       :: nmlfile

character(fnlen)         :: xtalname, outname
character(3)             :: method
real(kind=sgl)           :: voltage, dmin, convergence, startthick, thickinc, thetac, galen, bragg, klaue(2), io_real(6)
integer(kind=irg)        :: numthick, nt, npix, skip, dgn, pgnum, io_int(6), maxHOLZ, ik, k(3), numk, ga(3), gb(3), &
                            fn(3), nn, npx, npy, isym, numset, it, ijmax, jp

real(kind=dbl)           :: ctmp(192,3),arg, frac, abcdist(3), albegadist(3)
integer                  :: i,j,ir,nat(100),kk(3),&
                            n,ipx,ipy,minbeams, &
                            maxbeams,gzero,ic,ip,ikk     ! counters
real(kind=sgl)           :: ktmax,& ! maximum tangential component of wave vector
                            pre,&   ! prefactors 
                            tpi,Znsq, kkl, &
                            DBWF 
real,allocatable         :: thick(:), sr(:,:,:), EKI(:,:,:) ! thickness array, results
!real(kind=dbl), allocatable           :: Kossel(:)
complex(kind=dbl),allocatable         :: Lgh(:,:,:),Sgh(:,:)
complex(kind=dbl)        :: czero
logical	                 :: distort

namelist /ECPlist/ stdout, xtalname, voltage, k, fn, dmin, distort, abcdist, albegadist, ktmax, &
                   startthick, thickinc, numthick, npix, outname

! set the input parameters to default values (except for xtalname, which must be present)
xtalname = 'undefined'		! initial value to check that the keyword is present in the nml file
stdout = 6			! standard output
voltage = 30000.0		! acceleration voltage [V]
k = (/ 0, 0, 1 /)		! beam direction [direction indices]
fn = (/ 0, 0, 1 /)		! foil normal [direction indices]
dmin = 0.025			! smallest d-spacing to include in dynamical matrix [nm]
distort = .FALSE.                      ! distort the input unit cell ?  
abcdist = (/ 0.4, 0.4, 0.4/)           ! distorted a, b, c [nm]
albegadist = (/ 90.0, 90.0, 90.0 /)    ! distorted angles [degrees]
ktmax = 5.0                   ! beam convergence in units of |g_a|
startthick = 2.0		! starting thickness [nm]
thickinc = 2.0			! thickness increment
numthick = 10			! number of increments
npix = 256			! output arrays will have size npix x npix
outname = 'ecp.data'	! output filename
gzero = 1
frac = 0.05

! read the namelist file
open(UNIT=dataunit,FILE=nmlfile,DELIM='apostrophe',STATUS='old')
read(UNIT=dataunit,NML=ECPlist)
close(UNIT=dataunit,STATUS='keep')

if (trim(xtalname).eq.'undefined') then
  call FatalError('CTEMECP:',' structure file name is undefined in '//nmlfile)
end if

! print some information
 progname = 'CTEMECP.f90'
 progdesc = 'Large angle electron channeling pattern simulation'
 call CTEMsoft

! first get the crystal data and microscope voltage
 SG%SYM_reduce=.TRUE.
 call CrystalData(xtalname)

! do we need to apply a unit cell distortion ?  If so, then recompute the 
! unit cell parameters with the new lattice parameters and then continue 
if (distort) then
!! apply a deformation to this unit cell
    io_real(1:3) = abcdist(1:3)
    call WriteValue(' New lattice parameters a, b, and c [nm] ', io_real, 3, "(2(F8.5,','),F8.5)")
    io_real(1:3) = albegadist(1:3)
    call WriteValue(' New lattice angles alpha, beta, gamma [degrees] ', io_real, 3, "(2(F6.2,','),F6.2)")
    cell%a = abcdist(1); cell%b = abcdist(2); cell%c = abcdist(3)
    cell%alpha = albegadist(1); cell%beta = albegadist(2); cell%gamma = albegadist(3)
    call CalcMatrices
end if

! initialize the wave length and lattice potential computations
 skip = 3
 call CalcWaveLength(dble(voltage),skip)

! generate all atom positions
 call CalcPositions('v')

! transform the foil normal
 call TransSpace(float(fn),DynFN,'d','r')
 call NormVec(DynFN,'r')

! determine the point group number
 j=0
 do i=1,32
  if (SGPG(i).le.cell % SYM_SGnum) j=i
 end do

! use the new routine to get the whole pattern 2D symmetry group, since that
! is the one that determines the independent beam directions.
 dgn = GetPatternSymmetry(k,j,.TRUE.)
 pgnum = j
 isym = WPPG(dgn) ! WPPG lists the whole pattern point group numbers vs. diffraction group numbers

! determine the shortest reciprocal lattice points for this zone
 call ShortestG(k,ga,gb,isym)
 io_int(1:3)=ga(1:3)
 io_int(4:6)=gb(1:3)
 call WriteValue(' Reciprocal lattice vectors : ', io_int, 6,"('(',3I3,') and (',3I3,')',/)")

! initialize the HOLZ geometry type
 call GetHOLZGeometry(float(ga),float(gb),k,fn) 

! construct the list of all possible reflections
 method = 'ALL'
 thetac = (ktmax * 2.0 * bragg)*1000.0
 call Compute_ReflectionList(dmin,k,ga,gb,method,.FALSE.,maxHOLZ,thetac)
 galen = CalcLength(float(ga),'r')

! determine range of incident beam directions
!  bragg = CalcDiffAngle(ga(1),ga(2),ga(3))*0.5
  
! convert to ktmax along ga
!  ktmax = 0.5*thetac/bragg

! the number of pixels across the disk is equal to 2*npix + 1
  npx = npix
  npy = npx
  io_int(1) = 2.0*npx + 1
  call WriteValue('Number of image pixels along diameter of central disk = ', io_int, 1, "(I4)")
  mess=' '; call Message("(A/)")
  
! set parameters for wave vector computation
  klaue = (/ 0.0, 0.0 /)
  ijmax = float(npx)**2   ! truncation value for beam directions

! for now, the solution to the symmetry problem is to do the computation for the entire 
! illumination cone without application of symmetry.  Instead, we'll get the speed up by 
! going to multiple cores later on.
  isym = 1
  call CalckvectorsSymmetry(dble(k),dble(ga),dble(ktmax),npx,npy,numk,isym,ijmax,klaue,.FALSE.)
  io_int(1)=numk
  call WriteValue('Starting computation for # beam directions = ', io_int, 1, "(I8)")

! force dynamical matrix routine to read new Bethe parameters from file
  call Set_Bethe_Parameters(.TRUE.)

! set the thickness array
  nt = numthick
  allocate(thick(nt))
  thick = startthick + thickinc * (/ (i-1,i=1,nt) /)

  nat = 0

  mess = 'Starting main computation'
  call Message("(A/)")
  
!----------------------------MAIN COMPUTATIONAL LOOP-----------------------
! point to the first beam direction
  ktmp => khead
  minbeams = 10000
  maxbeams = 0
  czero = cmplx(0.0,0.0,dbl)
  pre = cmplx(0.0,cPi,dbl)
  call CalcUcg((/0,0,0/))   ! get the normal absorption parameter
  DynUpz = rlp%Vpmod
  numset = cell % ATOM_ntype  ! number of special positions in the unit cell
  tpi = 2.0*cPi
! allocate space for the results
  allocate(sr(2*npx+1,2*npy+1,nt)) ! ,EKI(2*npx+1,2*npy+1,nt))

!  work through the beam direction list
  beamloop: do ik=1,numk
! time the computation
   if (ik.eq.1) call Time_start

	ip = -ktmp%i
 	jp =  ktmp%j

! compute the dynamical matrix using Bloch waves with Bethe potentials; note that the IgnoreFoilNormal flag
! has been set to .FALSE.; if it is set to .TRUE., the computation of the ZOLZ will still be mostly correct,
! but the excitation errors of the HOLZ reflections will be increasingly incorrect with HOLZ order.  This was
! useful during program testing but should probably be removed as an option altogether...
	call Compute_DynMat('BLOCHBETHE', ktmp%k, ktmp%kt, .FALSE.)
        nn = DynNbeams

! then we need to initialize the Sgh array for the strong beams;
! this may need to be modified if we want to include real detector
! geometries and such; check Rossouw's paper for more information.
! this is now modified with respect to the older version to include
! Bethe potential strong reflections only  
        allocate(Sgh(nn,nn),Lgh(nn,nn,nt)) ! ,Kossel(nt))        
        Sgh = cmplx(0.0,0.0)

! for each special position we need to compute this array
        do ip=1,numset
            call CalcOrbit(ip,n,ctmp)
            nat(ip) = n
! get Zn-squared for this special position
            Znsq = float(cell%ATOM_type(ip))**2
! loop over all contributing reflections
! ir is the row index
            do ir=1,nn
! ic is the column index
             do ic=1,nn
              kk = BetheParameter%stronghkl(1:3,ir) - BetheParameter%stronghkl(1:3,ic)
! We'll assume isotropic Debye-Waller factors for now ...
! That means we need the square of the length of s=  kk^2/4
              kkl = 0.25 * CalcLength(float(kk),'r')**2
! here is where we should insert the proper weight factor, Z^2 exp[-M_{h-g}]
! and also the detector geometry...   For now, we do nothing with the detector
! geometry; the Rossouw et al 1994 paper lists a factor A that does not depend
! on anything in particular, so we assume it is 1. 
              do ikk=1,n
! get the argument of the complex exponential
                arg = tpi*sum(kk(1:3)*ctmp(ikk,1:3))
! Debye-Waller exponential
                DBWF = exp(-cell%ATOM_pos(ip,5)*kkl)
! multiply with the prefactor and add
                Sgh(ir,ic) = Sgh(ir,ic) + cmplx(Znsq * DBWF,0.0) * cmplx(cos(arg),sin(arg))
              end do
             end do
            end do  
          end do

! solve the dynamical eigenvalue equation
          allocate(W(nn),CG(nn,nn),alpha(nn))
          call CalcLgh(nn,nt,thick,ktmp%kn,gzero,Lgh) !,Kossel)
          deallocate(W,CG,alpha)

! and store the resulting values
          ipx = ktmp%i + npx + 1
          ipy = ktmp%j + npy + 1
          do it=1,nt
            sr(ipx,ipy,it) = real(sum(Lgh(1:nn,1:nn,it)*Sgh(1:nn,1:nn)))/float(sum(nat))
          end do
!          EKI(ipx,ipy,1:nt) = Kossel(1:nt)

          deallocate(Lgh, Sgh) !, Kossel)
! select next beam direction
          ktmp => ktmp%next
          
! update computation progress
   if (float(ik)/float(numk) .gt. frac) then
    io_int(1) = nint(100.0*frac) 
    call WriteValue('       ', io_int, 1, "(1x,I3,' percent completed')") 
    frac = frac + 0.05
   end if  

  end do beamloop


! stop the clock and report the total time     
  call Time_stop(numk)

! here we probably need to store some additional information for the IDL interface  
  open(unit=dataunit,file=trim(outname),status='unknown',form='unformatted')
  write (dataunit) 2*npx+1,2*npy+1,nt
  write (dataunit) sr
!  write (dataunit) EKI
  close(unit=dataunit,status='keep')

end subroutine ECpattern

!--------------------------------------------------------------------------
!
! SUBROUTINE:CalcLgh
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief integrate the Bloch wave function over the foil thickness
!
!> @param nn number of strong beams
!> @param nt number of thickness values
!> @param thick array of thickness values
!> @param kn normal component of incident wave vector
!> @param gzero index of zero beam (should always be the first one; legacy parameter)
!> @param Lgh output array
!> @param Kossel Kossel intensity array
!
!> @date 11/18/13  MDG 1.0 major rewrite from older ECP program; merged with ECPz
!--------------------------------------------------------------------------
subroutine CalcLgh(nn,nt,thick,kn,gzero,Lgh) !,Kossel)

use local
use io
use files
use diffraction
use dynamical
use constants

IMPLICIT NONE

integer(kind=sgl),INTENT(IN)        :: nn
integer(kind=sgl),INTENT(IN)        :: nt
real(kind=sgl),INTENT(IN)           :: thick(nt)
real(kind=dbl),INTENT(IN)           :: kn
integer(kind=sgl),INTENT(IN)        :: gzero
complex(kind=dbl),INTENT(OUT)       :: Lgh(nn,nn,nt)
!real(kind=dbl),INTENT(OUT)          :: Kossel(nt)

integer                             :: i,j,it,ig,ih,IPIV(nn)
complex(kind=dbl),allocatable       :: CGinv(:,:), Minp(:,:),tmp3(:,:)
complex(kind=dbl)                   :: Ijk(nn,nn),q
real(kind=dbl)                      :: qq, s, t

allocate(CGinv(nn,nn),Minp(nn,nn),tmp3(nn,nn))

! compute the eigenvalues and eigenvectors
! using the LAPACK CGEEV, CGETRF, and CGETRI routines
! 
! then get the eigenvalues and eigenvectors
 Minp = DynMat
 IPIV = 0

 call BWsolve(Minp,W,CG,CGinv,nn,IPIV)

! then compute the integrated intensity matrix
 W = W/cmplx(2.0*kn,0.0)

! first do the Lgh matrices, looping over the thickness
do it=1,nt
! recall that alpha(1:nn) = CGinv(1:nn,gzero)
! first the Ijk matrix
 do i=1,nn
  do j=1,nn
    q = 2.0*cPi*thick(it)*cmplx(aimag(W(i))+aimag(W(j)),real(W(i))-real(W(j)))
    Ijk(i,j) = conjg(CGinv(i,gzero)) * (1.0-exp(-q))/q * CGinv(j,gzero)
  end do
 end do

! then the summations for Lgh
 do ih=1,nn
   do i=1,nn
      tmp3(ih,i) = sum(Ijk(i,1:nn)*CG(ih,1:nn))
   end do
 end do
 do ig=1,nn
  do ih=1,nn
     Lgh(ih,ig,it) = sum(conjg(CG(ig,1:nn))*tmp3(ih,1:nn))
  end do
 end do
end do ! thickness loop

! and finally the Kossel intensities (commented out for now)
! Kossel = 0.D0
! do j=1,nn
!    qq = -4.D0*cPi*abs(aimag(W(j)))
!    s = cabs(CGinv(j,gzero))**2
!    do it=1,nt
!      t = qq*thick(it)
!      if (abs(t).lt.1.D0) Kossel(it) = Kossel(it) +  s * exp(t)
!    end do
! end do

deallocate(CGinv,Minp,tmp3)
end subroutine

