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
!> @date 03/03/14 MDG 3.1 rewrite with scattering matrix formalism
!--------------------------------------------------------------------------
program CTEMECPSM

use local
use files
use io

IMPLICIT NONE

character(fnlen)	:: nmldeffile

! deal with the command line arguments, if any
nmldeffile = 'CTEMECP.nml'
progname = 'CTEMECP.f90'
call Interpret_Program_Arguments(nmldeffile,2,(/ 0, 40 /) )

! perform the zone axis computations
call ECpattern(nmldeffile)

end program CTEMECPSM

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
!< @date 11/22/13  MDG 1.1 output modified for IDL interface
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
real(kind=sgl)           :: voltage, dmin, startthick, thickinc, thetac, galen, bragg, klaue(2), io_real(6), &
                            kstar(3), gperp(3), delta
integer(kind=irg)        :: numthick, nt, npix, skip, dgn, pgnum, io_int(6), maxHOLZ, ik, k(3), numk, ga(3), gb(3), &
                            fn(3), nn, npx, npy, isym, numset, it, ijmax, jp

real(kind=dbl)           :: ctmp(192,3),arg, abcdist(3), albegadist(3)
integer                  :: i,j,ir,nat(100),kk(3),&
                            n,ipx,ipy,minbeams, &
                            maxbeams,gzero,ic,ip,ikk     ! counters
real(kind=sgl)           :: ktmax,& ! maximum tangential component of wave vector
                            pre,&   ! prefactors 
                            tpi,Znsq, kkl, &
                            DBWF, frac
real,allocatable         :: thick(:), sr(:,:,:) ! thickness array, results
!real(kind=dbl), allocatable           :: Kossel(:)
complex(kind=dbl),allocatable         :: Lgh(:,:,:),Sgh(:,:)
complex(kind=dbl)        :: czero
logical	                 :: distort

namelist /ECPlist/ stdout, xtalname, voltage, k, fn, dmin, distort, abcdist, albegadist, ktmax, &
                   startthick, thickinc, numthick, npix, outname, thetac

! set the input parameters to default values (except for xtalname, which must be present)
xtalname = 'undefined'		        ! initial value to check that the keyword is present in the nml file
stdout = 6			        ! standard output
voltage = 30000.0		        ! acceleration voltage [V]
k = (/ 0, 0, 1 /)		        ! beam direction [direction indices]
fn = (/ 0, 0, 1 /)		        ! foil normal [direction indices]
dmin = 0.025			        ! smallest d-spacing to include in dynamical matrix [nm]
distort = .FALSE.                      ! distort the input unit cell ?  
abcdist = (/ 0.4, 0.4, 0.4/)           ! distorted a, b, c [nm]
albegadist = (/ 90.0, 90.0, 90.0 /)    ! distorted angles [degrees]
ktmax = 0.0                            ! beam convergence in units of |g_a|
thetac = 0.0                           ! beam convergence in mrad (either ktmax or thetac must be given)
startthick = 2.0		        ! starting thickness [nm]
thickinc = 2.0			        ! thickness increment
numthick = 10			        ! number of increments
npix = 256			        ! output arrays will have size npix x npix
outname = 'ecp.data'        	        ! output filename

! init some parameters
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
else
    abcdist = (/ cell%a, cell%b, cell%c /)
    albegadist = (/ cell%alpha, cell%beta, cell%gamma /)
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
 bragg = CalcDiffAngle(ga(1),ga(2),ga(3))*0.5
 if (ktmax.ne.0.0) then 
   thetac = (ktmax * 2.0 * bragg)*1000.0
 else
   ktmax = thetac / (2000.0 * bragg)
 end if
 io_real(1) = thetac
 call WriteValue(' Pattern convergence angle [mrad] = ',io_real,1,"(F8.3)")
 io_real(1) = bragg*1000.0
 call WriteValue(' Bragg angle of g_a [mrad] = ',io_real,1,"(F6.3)")
 call Compute_ReflectionList(dmin,k,ga,gb,method,.FALSE.,maxHOLZ,thetac/1000.0)
 galen = CalcLength(float(ga),'r')

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
  call CalckvectorsSymmetry(dble(k),dble(ga),dble(ktmax),npx,npy,numk,isym,ijmax,klaue)
  io_int(1)=numk
  call WriteValue('Starting computation for # beam directions = ', io_int, 1, "(I8)")

! force dynamical matrix routine to read new Bethe parameters from file
  call Set_Bethe_Parameters(.TRUE.)

! set the thickness array
  nt = numthick
  allocate(thick(nt))
  thick = startthick + thickinc * (/ (i-1,i=1,nt) /)

  nat = 0
  
  call Time_report(frac)
  call Time_start

  
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

	ip = -ktmp%i
 	jp =  ktmp%j

! compute the dynamical matrix using Bloch waves with Bethe potentials; note that the IgnoreFoilNormal flag
! has been set to .FALSE.; if it is set to .TRUE., the computation of the ZOLZ will still be mostly correct,
! but the excitation errors of the HOLZ reflections will be increasingly incorrect with HOLZ order.  This was
! useful during program testing but should probably be removed as an option altogether...
	call Compute_DynMat('BLOCHBETHE', ktmp%k, ktmp%kt, .FALSE.)
        nn = DynNbeams
        
! next we need to convert the Bloch dynamical matrix to hte A matrix for the scattering matrix formalism...



!write (*,*) ik, nn

!if (ik.eq.1) then
!  do ir=1,nn
!    write (*,*)  ir,' : ',BetheParameter%stronghkl(1:3,ir) 
!  end do
!end if

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
            Znsq = float(cell%ATOM_type(ip))**2 * cell%ATOM_pos(ip,4)
! loop over all contributing reflections
! ir is the row index
            do ir=1,nn
! ic is the column index
             do ic=1,nn
              kk = BetheParameter%stronghkl(1:3,ir) - BetheParameter%stronghkl(1:3,ic)
! We'll assume isotropic Debye-Waller factors for now ...
! That means we need the square of the length of s=  kk^2/4
              kkl = 0.25 * CalcLength(float(kk),'r')**2
! Debye-Waller exponential
                DBWF = Znsq * exp(-cell%ATOM_pos(ip,5)*kkl)
! here is where we should insert the proper weight factor, Z^2 exp[-M_{h-g}]
! and also the detector geometry...   For now, we do nothing with the detector
! geometry; the Rossouw et al 1994 paper lists a factor A that does not depend
! on anything in particular, so we assume it is 1. 
              do ikk=1,n
! get the argument of the complex exponential
                arg = tpi*sum(kk(1:3)*ctmp(ikk,1:3))
! multiply with the prefactor and add
                Sgh(ir,ic) = Sgh(ir,ic) + cmplx(DBWF,0.0) * cmplx(cos(arg),sin(arg))
              end do
             end do
            end do  
          end do
          
!write (*,*)  'calling CalcLghSM'
! solve the dynamical eigenvalue equation
          call CalcLghSM(nn,nt,thick,ktmp%kn,gzero,Lgh) !,Kossel)

! and store the resulting values
          ipx = ktmp%i + npx + 1
          ipy = ktmp%j + npy + 1
          do it=1,nt
            sr(ipx,ipy,it) = real(sum(Lgh(1:nn,1:nn,it)*Sgh(1:nn,1:nn)))/float(sum(nat))
          end do

          deallocate(Lgh, Sgh) !, Kossel)
! select next beam direction
          ktmp => ktmp%next
          
! update computation progress
   if (float(ik)/float(numk) .gt. frac) then
     call Time_remaining(ik,numk)
     frac = frac + 0.05
   end if  

  end do beamloop

! stop the clock and report the total time     
  call Time_stop(numk)

! store additional information for the IDL interface  
  open(unit=dataunit,file=trim(outname),status='unknown',action='write',form='unformatted')
! write the program identifier
  write (dataunit) trim(progname)
! write the version number
  write (dataunit) scversion
! first write the array dimensions
  write (dataunit) 2*npix+1,2*npix+1,nt
! then the name of the crystal data file
  write (dataunit) xtalname
! altered lattice parameters
  if (distort) then 
    write (dataunit) 1
  else
    write (dataunit) 0
  end if
! new lattice parameters and angles
  write (dataunit) abcdist
  write (dataunit) albegadist
! the accelerating voltage [V]
  write (dataunit) voltage
! convergence angle [mrad]
  write (dataunit) thetac
! max kt value in units of ga
  write (dataunit) ktmax
! the zone axis indices
  write (dataunit) k
! the foil normal indices
  write (dataunit) fn
! number of k-values in disk
  write (dataunit) numk
! dmin value
  write (dataunit) dmin
! horizontal reciprocal lattice vector
  write (dataunit) ga  
! length horizontal reciprocal lattice vector (need for proper Laue center coordinate scaling)
  write (dataunit) galen
! we need to store the gperp vectors
  delta = 2.0*ktmax*galen/float(2*npix+1)        ! grid step size in nm-1 
  call TransSpace(float(k),kstar,'d','r')        ! transform incident direction to reciprocal space
  call CalcCross(float(ga),kstar,gperp,'r','r',0)! compute g_perp = ga x k
  call NormVec(gperp,'r')                        ! normalize g_perp
  write (dataunit) delta
  write (dataunit) gperp
! eight integers with the labels of various symmetry groups
  write (dataunit) (/ pgnum, PGLaue(pgnum), dgn, PDG(dgn), BFPG(dgn), WPPG(dgn), DFGN(dgn), DFSP(dgn) /)
! thickness data
  write (dataunit) startthick, thickinc
! and the actual data array
  write (dataunit) sr
  close(unit=dataunit,status='keep')

  mess = 'Data stored in output file '//trim(outname)
  call Message("(/A/)")

end subroutine ECpattern

!--------------------------------------------------------------------------
!
! SUBROUTINE:CalcLghSM
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
!> @date 03/03/14  MDG 2.0 version that works with the scattering matrix...
!--------------------------------------------------------------------------
recursive subroutine CalcLghSM(nn,nt,thick,kn,gzero,Lgh) !,Kossel)

use local
use io
use files
use diffraction
use dynamical
use constants
use math

IMPLICIT NONE

integer(kind=sgl),INTENT(IN)        :: nn
integer(kind=sgl),INTENT(IN)        :: nt
real(kind=sgl),INTENT(IN)           :: thick(nt)
real(kind=dbl),INTENT(IN)           :: kn
integer(kind=sgl),INTENT(IN)        :: gzero
complex(kind=dbl),INTENT(OUT)       :: Lgh(nn,nn,nt)
!real(kind=dbl),INTENT(OUT)          :: Kossel(nt)

integer                             :: i
complex(kind=dbl),allocatable       :: Minp(:,:),Azz(:,:),ampl(:),ampl2(:)
real(kind=sgl)                      :: dthick
  
  allocate(Minp(nn,nn),Azz(nn,nn),ampl(nn),ampl2(nn))

  dthick = thick(2)-thick(1)
  Minp = DynMat * dcmplx(0.D0,cPi * mLambda)
  call MatrixExponential(Minp, Azz, dble(dthick), 'Pade', nn)  

  ampl = dcmplx(0.D0,0.D0)
  ampl(1) = dcmplx(1.0D0,0.D0)
  Lgh = dcmplx(0.D0,0.D0)

! add some thickness handling here !!!
! thge integration uses a small step size, but there might
! be only a small number of thicknesses for which output is
! requested !

  do i=1,nt
    ampl2 = matmul(Azz,ampl)
    if (i.eq.1) then 
      Lgh(1:nn,1:nn,i) = spread(ampl2(1:nn),dim=2,ncopies=nn)*spread(conjg(ampl2(1:nn)),dim=1,ncopies=nn)
    else
      Lgh(1:nn,1:nn,i) = Lgh(1:nn,1:nn,i-1)+spread(ampl2(1:nn),dim=2,ncopies=nn)*spread(conjg(ampl2(1:nn)),dim=1,ncopies=nn)
    end if
!    do ih=1,nn
!      do ig=1,nn
!        Lgh(ig,ih,i) = Lgh(ig,ih,i) + ampl2(ih) * conjg(ampl2(ig))
!      end do
!    end do
    ampl = ampl2
  end do

  
  deallocate(Minp,Azz,ampl,ampl2)
end subroutine

