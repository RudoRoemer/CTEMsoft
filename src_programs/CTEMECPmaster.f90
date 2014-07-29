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
! CTEMsoft2013:CTEMECPmaster.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMECP
!
!> @author Marc De Graef/ Saransh Singh, Carnegie Mellon University
!
!> @brief Zone axis electron channeling master patterns
!
!> @date 03/18/10 MDG 1.0 f90
!> @date 08/09/10 MDG 2.0 corrected weight factors and g-vector ordering problem
!> @date 11/18/13 MDG 3.0 major rewrite with new libraries
!> @date 06/27/14 MDG 4.0 removal of all globals; separation of namelist handling from computation
!> @date 06/30/14 MDG 4.1 added OpenMP
!> @DATE 07/23/14 SS  4.2 rewrite master pattern
!--------------------------------------------------------------------------

program CTEMECPmaster

use local
use typedefs
use NameListTypedefs
use NameListHandlers
use files
use io

IMPLICIT NONE

character(fnlen)                        :: nmldeffile, progname, progdesc
type(ECPMasterNameListType)                   :: ecpnl

nmldeffile = 'CTEMECPmaster.nml'
progname = 'CTEMECPmaster.f90'
progdesc = 'Master pattern generation for Electron channeling pattern'

! deal with the command line arguments, if any
call Interpret_Program_Arguments(nmldeffile,2,(/ 0, 40 /), progname)

! deal with the namelist stuff
call GetECPMasterNameList(nmldeffile,ecpnl)

! print some information
call CTEMsoft(progname, progdesc)

! perform the zone axis computations
call ECmasterpattern(ecpnl, progname)

end program CTEMECPmaster

!--------------------------------------------------------------------------
!
! SUBROUTINE:ECmasterpattern
!
!> @author Marc De Graef/Saransh Singh, Carnegie Mellon University
!
!> @brief compute a master electron channeling pattern
!
!> @note This is really very similar to a LACBED and EBSDmaster computation, except that
!> the final intensity computation is somewhat different.  We could in
!> principle also include the Kossel pattern computation in this program.
!> This program now also includes the Bethe potential approximation, to
!> hopefully speed things up a little bit...
!
!> @param ecpnl name list structure
!> @param progname program name
!
!> @date 11/18/13  MDG 1.0 major rewrite from older ECP program
!> @date 11/22/13  MDG 1.1 output modified for IDL interface
!> @date 03/04/14  MDG 1.2 added scattering matrix mode
!> @date 06/27/14  MDG 2.0 removal of globals, split of namelist and computation; OpenMP
!> @date 06/30/14  MDG 2.1 debug; found some inconsistent array (de)allocations
!> @date 07/23/14  SS  2.2 conversion to master pattern simulation
!--------------------------------------------------------------------------
subroutine ECmasterpattern(ecpnl, progname)

use local
use typedefs
use NameListTypedefs
use crystal
use symmetry
use Lambert
use initializers
use constants
use gvectors
use kvectors
use error
use io
use files
use diffraction
use omp_lib
use MBModule

IMPLICIT NONE

type(ECPNameListType),INTENT(IN)        :: ecpnl
character(fnlen),INTENT(IN)             :: progname

!character(3)                    :: method
!real(kind=sgl)                  :: galen, bragg, klaue(2), io_real(6), kstar(3), gperp(3), delta, thetac, &
!kk(3), ktmax, FN(3), kn, fnat
!integer(kind=irg)               :: nt, skip, dgn, pgnum, io_int(6), maxHOLZ, ik, numk, ga(3), gb(3), TID, &
!nn, npx, npy, isym, numset, it, ijmax, jp, istat, iequiv(2,12), nequiv, NUMTHREADS, &
!nns, nnw, nref, tots, totw
integer(kind=irg)                :: numzbins, totnum_el, numsx, numsy ! read from energyfile
real(kind=dbl)                   :: EkeV, Ehistmin, Ebinsize, depthmax, depthstep, sig, omega ! read from energyfile
real(kind=dbl), allocatable      :: accum_z(:,:,:,:), accum_e(:,:,:,:)
character(4)                     :: MCmode
real(kind=dbl)                   :: io_real(6) ! auxiliary variables
!real(kind=dbl)                  :: ctmp(192,3),arg
!integer                         :: i,j,ir, n,ipx,ipy,gzero,ic,ip,ikk
!real(kind=sgl)                  :: pre, tpi,Znsq, kkl, DBWF, frac
!real,allocatable                :: thick(:), sr(:,:,:)
!complex(kind=dbl),allocatable   :: Lgh(:,:,:),Sgh(:,:),Sghtmp(:,:,:)
!complex(kind=dbl)               :: czero
!real(kind=sgl),allocatable      :: karray(:,:)
!integer(kind=irg),allocatable   :: kij(:,:), nat(:)
!complex(kind=dbl),allocatable   :: DynMat(:,:)
!logical                         :: verbose

type(unitcell),pointer          :: cell
!type(gnode)                     :: rlp
!type(DynType)                   :: Dyn
type(kvectorlist),pointer       :: khead, ktmp
!type(symdata2D)                 :: TDPG
!type(BetheParameterType)        :: BetheParameters
!type(reflisttype),pointer       :: reflist, firstw,rltmp


! init some parameters
gzero = 1
frac = 0.05

nullify(khead)
nullify(ktmp)

allocate(cell)

!=============================================
!=============================================
! ---------- read Monte Carlo output file and extract necessary parameters
! first, we need to load the data from the MC program.
call Message('opening '//trim(emnl%energyfile), frm = "(A)" )

open(dataunit,file=trim(emnl%energyfile),status='unknown',form='unformatted')

! lines from CTEMMC.f90... these are the things we need to read in...
! write (dataunit) progname
!! write the version number
! write (dataunit) scversion
!! then the name of the crystal data file
! write (dataunit) xtalname
!! energy information etc...
! write (dataunit) numEbins, numzbins, numsx, numsy, num_el*NUMTHREADS, NUMTHREADS
! write (dataunit) EkeV, Ehistmin, Ebinsize, depthmax, depthstep
! write (dataunit) sig, omega
! write (dataunit) MCmode
!! and here are the actual results
! write (dataunit) accum_e
! write (dataunit) accum_z

read (dataunit) oldprogname
read (dataunit) MCscversion
read (dataunit) xtalname

read(dataunit) numzbins, numsx, numsy, totnum_el
numsx = (numsx - 1)/2
numsy = (numsy - 1)/2

read (dataunit) EkeV, Ehistmin, Ebinsize, depthmax, depthstep
io_real(1:5) = (/ EkeV, Ehistmin, Ebinsize, depthmax, depthstep /)
call WriteValue(' EkeV, Ehistmin, Ebinsize, depthmax, depthstep ',io_real,5,"(4F10.5,',',F10.5)")

read (dataunit) sig, omega
read (dataunit) MCmode

!@TODO Modify the MC program to include only accum_z array in energyfile
allocate(accum_e(numEbins,-numsx:numsx,-nsumy:numsy),accum_z(numEbins,numzbins,-numsx/10:numsx/10,-numsy/10:numsy/10),stat=istat)
read(dataunit) accum_e
! actually, we do not yet need the accum_e array for ECP. This will be removed with an updated version of the MC code
! but we need to skip it in this unformatted file so that we can read the accum_z array ...
deallocate(accum_e)

read(dataunit) accum_z    ! we only need this array for the depth integrations

close(dataunit,status='keep')
call Message(' -> completed reading '//trim(emnl%energyfile), frm = "(A)")

!=============================================
!=============================================


!=============================================
!=============================================
! crystallography section

nullify(cell)
allocate(cell)

! load the crystal structure and compute the Fourier coefficient lookup table
verbose = .TRUE.
call Initialize_Cell(cell,Dyn,rlp, xtalname, ecpnl%dmin, sngl(1000.0*EkeV),verbose)

! determine the point group number
j=0
do i=1,32
if (SGPG(i).le.cell % SYM_SGnum) j=i
end do

! use the new routine to get the whole pattern 2D symmetry group, since that
! is the one that determines the independent beam directions.
dgn = GetPatternSymmetry(cell,ecpnl%k,j,.TRUE.)
pgnum = j
isym = WPPG(dgn) ! WPPG lists the whole pattern point group numbers vs. diffraction group numbers

! determine the shortest reciprocal lattice points for this zone
call ShortestG(cell,ecpnl%k,ga,gb,isym)
io_int(1:3)=ga(1:3)
io_int(4:6)=gb(1:3)
call WriteValue(' Reciprocal lattice vectors : ', io_int, 6,"('(',3I3,') and (',3I3,')',/)")

! diffraction geometry
bragg = CalcDiffAngle(cell,ga(1),ga(2),ga(3))*0.5
if (ecpnl%ktmax.ne.0.0) then
thetac = (ecpnl%ktmax * 2.0 * bragg)*1000.0
ktmax = ecpnl%ktmax
io_real(1) = thetac
else
ktmax = ecpnl%thetac / (2000.0 * bragg)
!write (*,*) 'ktmax value = ',ktmax
thetac = ecpnl%thetac
io_real(1) = ecpnl%thetac
end if
call WriteValue(' Pattern convergence angle [mrad] = ',io_real,1,"(F8.3)")
io_real(1) = bragg*1000.0
call WriteValue(' Bragg angle of g_a [mrad] = ',io_real,1,"(F6.3)")

! the number of pixels across the disk is equal to 2*npix + 1
npx = ecpnl%npix
npy = npx
io_int(1) = 2.0*npx + 1
call WriteValue('Number of image pixels along diameter of central disk = ', io_int, 1, "(I4/)")

! for now, the solution to the symmetry problem is to do the computation for the entire
! illumination cone without application of symmetry.  Instead, we'll get the speed up by
! going to multiple cores later on.
isym = 1
! set parameters for wave vector computation
klaue = (/ 0.0, 0.0 /)
ijmax = float(npx)**2   ! truncation value for beam directions

call CalckvectorsSymmetry(khead,cell,TDPG,dble(ecpnl%k),dble(ga),dble(ktmax),npx,npy,numk,isym,ijmax,klaue)
io_int(1)=numk
call WriteValue('Starting computation for # beam directions = ', io_int, 1, "(I8)")

! force dynamical matrix routine to read new Bethe parameters from file
call Set_Bethe_Parameters(BetheParameters,.TRUE.)
!write(*,*) 'BEthe : ',BetheParameters%c1, BetheParameters%c2, BetheParameters%c3

! set the thickness array
nt = ecpnl%numthick
allocate(thick(nt))
thick = ecpnl%startthick + ecpnl%thickinc * (/ (i-1,i=1,nt) /)

!----------------------------MAIN COMPUTATIONAL LOOP-----------------------
czero = cmplx(0.D0,0.D0)
pre = cmplx(0.D0,1.D0) * cPi
tpi = 2.D0*cPi
numset = cell % ATOM_ntype  ! number of special positions in the unit cell
allocate(nat(numset))
nat = 0
fnat = 1.0/float(sum(cell%numat(1:numset)))
!write (*,*) 'asymmetric unit cell # ',numset,fnat

! in preparation for the threaded portion of the program, we need to
! copy the wave vectors into an array rather than a linked list
allocate(karray(4,numk), kij(2,numk),stat=istat)
! point to the first beam direction
ktmp => khead
! and loop through the list, keeping k, kn, and i,j
karray(1:3,1) = sngl(ktmp%k(1:3))
karray(4,1) = sngl(ktmp%kn)
kij(1:2,1) = (/ ktmp%i, ktmp%j /)
do ik=2,numk
ktmp => ktmp%next
karray(1:3,ik) = sngl(ktmp%k(1:3))
karray(4,ik) = sngl(ktmp%kn)
kij(1:2,ik) = (/ ktmp%i, ktmp%j /)
end do
! and remove the linked list
call Delete_kvectorlist(khead)

! allocate space for the results
allocate(sr(2*npx+1,2*npy+1,nt))
sr = 0.0

! set the number of OpenMP threads
call OMP_SET_NUM_THREADS(ecpnl%nthreads)
io_int(1) = ecpnl%nthreads
call WriteValue(' Attempting to set number of threads to ',io_int, 1, frm = "(I4)")

! use OpenMP to run on multiple cores ...
!$OMP PARALLEL default(shared) PRIVATE(DynMat,ik,TID,kk,kn,ipx,ipy,iequiv,nequiv,fnat,ip,jp,reflist,firstw,nns,nnw,nref) &
!$OMP& PRIVATE(Sgh, Lgh, SGHtmp, FN)

NUMTHREADS = OMP_GET_NUM_THREADS()
TID = OMP_GET_THREAD_NUM()

nullify(reflist)
nullify(firstw)

nns = 0
nnw = 0
tots = 0
totw = 0

!$OMP DO SCHEDULE(DYNAMIC,100)

!  work through the beam direction list
beamloop: do ik=1,numk

! generate the reflectionlist
kk(1:3) = karray(1:3,ik)
FN = kk
call Initialize_ReflectionList(cell, reflist, BetheParameters, FN, kk, ecpnl%dmin, nref)

! determine strong and weak reflections
call Apply_BethePotentials(cell, reflist, firstw, BetheParameters, nref, nns, nnw)

! generate the dynamical matrix
allocate(DynMat(nns,nns))
call GetDynMat(cell, reflist, firstw, rlp, DynMat, nns, nnw)

! then we need to initialize the Sgh and Lgh arrays
if (allocated(Sgh)) deallocate(Sgh)
if (allocated(Lgh)) deallocate(Lgh)
if (allocated(Sghtmp)) deallocate(Sghtmp)

allocate(Sghtmp(nns,nns,numset),Lgh(nns,nns,nt),Sgh(nns,nns))
Sgh = czero
Sghtmp = czero
Lgh = czero
nat = 0
call CalcSgh(cell,reflist,nns,numset,Sghtmp,nat)

! sum Sghtmp over the sites
Sgh = sum(Sghtmp,3)

! solve the dynamical eigenvalue equation
kn = karray(4,ik)
call CalcLghECP(DynMat,Lgh,nns,nt,thick,dble(kn),gzero)
deallocate(DynMat,Sghtmp)

! and store the resulting values
ipx = kij(1,ik)
ipy = kij(2,ik)

!$OMP CRITICAL
if (isym.ne.1) then
call Apply2DLaueSymmetry(ipx,ipy,isym,iequiv,nequiv)
iequiv(1,1:nequiv) = iequiv(1,1:nequiv) + npx + 1
iequiv(2,1:nequiv) = iequiv(2,1:nequiv) + npy + 1
do ip=1,nequiv
do jp=1,nt
sr(iequiv(1,ip),iequiv(2,ip),jp) = sr(iequiv(1,ip),iequiv(2,ip),jp) + &
real(sum(Lgh(1:nns,1:nns,jp)*Sgh(1:nns,1:nns)))
end do
end do
else
do jp=1,nt
sr(ipx+npx+1,ipy+npy+1,jp) = real(sum(Lgh(1:nns,1:nns,jp)*Sgh(1:nns,1:nns)))
end do
end if
totw = totw + nnw
tots = tots + nns
!$OMP END CRITICAL

! if (TID.eq.0) write (*,*) ik,sr(ipx,ipy,1:3)

! if (sr(ipx,ipy,1).gt.1000.0) write (*,*) TID, ik, sr(ipx,ipy,1), fnat, maxval(cdabs(Lgh)), maxval(cdabs(Sgh))

deallocate(Lgh, Sgh)

if (mod(ik,2500).eq.0) then
io_int(1) = ik
call WriteValue('  completed beam direction ',io_int, 1, "(I8)")
! write(*,*) minval(sr),maxval(sr)
end if

call Delete_gvectorlist(reflist)

end do beamloop

!$OMP END PARALLEL

sr = sr*fnat

! store additional information for the IDL interface
open(unit=dataunit,file=trim(ecpnl%outname),status='unknown',action='write',form='unformatted')
! write the program identifier
write (dataunit) trim(progname)
! write the version number
write (dataunit) scversion
! first write the array dimensions
write (dataunit) 2*ecpnl%npix+1,2*ecpnl%npix+1,nt
! then the name of the crystal data file
write (dataunit) ecpnl%xtalname
! altered lattice parameters; also combine compmode in this parameter
!  Bloch waves, no distortion: 0
!  Bloch waves, distortion:    1
!  ScatMat, no distortion      2
!  ScatMat, distortion         3
! if (distort) then
!   if (compmode.eq.'Blochwv') then
write (dataunit) 1
!   else
!     write (dataunit) 3
!   end if
! else
!   if (compmode.eq.'Blochwv') then
!     write (dataunit) 0
!   else
!     write (dataunit) 2
!   end if
! end if
! new lattice parameters and angles
write (dataunit) (/ cell%a, cell%b, cell%c /)  ! abcdist
write (dataunit) (/ cell%alpha, cell%beta, cell%gamma /) ! albegadist
! the accelerating voltage [V]
write (dataunit) ecpnl%voltage
! convergence angle [mrad]
write (dataunit) thetac
! max kt value in units of ga
write (dataunit) ktmax
! the zone axis indices
write (dataunit) ecpnl%k
! the foil normal indices
write (dataunit) ecpnl%fn
! number of k-values in disk
write (dataunit) numk
! dmin value
write (dataunit) ecpnl%dmin
! horizontal reciprocal lattice vector
write (dataunit) ga
! length horizontal reciprocal lattice vector (need for proper Laue center coordinate scaling)
write (dataunit) galen
! we need to store the gperp vectors
delta = 2.0*ktmax*galen/float(2*ecpnl%npix+1)        ! grid step size in nm-1
call TransSpace(cell,float(ecpnl%k),kstar,'d','r')        ! transform incident direction to reciprocal space
call CalcCross(cell,float(ga),kstar,gperp,'r','r',0)! compute g_perp = ga x k
call NormVec(cell,gperp,'r')                        ! normalize g_perp
write (dataunit) delta
write (dataunit) gperp
! eight integers with the labels of various symmetry groups
write (dataunit) (/ pgnum, PGLaue(pgnum), dgn, PDG(dgn), BFPG(dgn), WPPG(dgn), DFGN(dgn), DFSP(dgn) /)
! thickness data
write (dataunit) ecpnl%startthick, ecpnl%thickinc
! and the actual data array
write (dataunit) sr
close(unit=dataunit,status='keep')

call Message('Data stored in output file '//trim(ecpnl%outname), frm = "(/A/)")

tots = nint(float(tots)/float(numk))
totw = nint(float(totw)/float(numk))

io_int(1:2) = (/ tots, totw /)
call WriteValue(' Average # strong, weak beams = ',io_int, 2, "(I5,',',I5/)")

end subroutine ECpattern
