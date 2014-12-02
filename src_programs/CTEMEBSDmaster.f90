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
! CTEMsoft2013:CTEMEBSDmaster.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMEBSDmaster
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief CTEMEBSDzE computes the energy-dependent master EBSD pattern for a given structure
!
!> @todo implement full symmetry use; implement multiple reflection output
!> or, easier perhaps, selection of one reflection;
!>
!> implement OpenMP multithreading for the actual computation part; requires modifications
!> in CTEMlib.a routines (mostly THREADPRIVATE commands in several modules)
!
!> @date  03/08/12  MDG 1.0 EBSD program for fundamental zone patterns
!> @date  08/17/12  MDG 1.1 added generalized fundamental zone for all crystal symmetries
!> @date  08/20/12  MDG 1.2 modifid FZ to use the Lambert projection of a square sampling grid
!> @date  09/06/12  MDG 1.3 added support for second setting of Laue group -3m
!> @date  11/21/12  MDG 2.0 added full Lambert projection support for both square and hexagonal grids
!>				   the older code is still available for now, but will be removed after validation
!>				   of the newer code.
!> @date  12/04/12  MDG 2.1 added support for equilateral triangle mapping; needs to be validated.
!>				   also modified the structure of the output file, so that EBSD.f90 will
!>				   know which of the inverse mapping methods it should use.
!> @date  12/10/12  MDG 3.0 expanded EBSDFZ program to include energy-dependencies from Monte Carlo
!> @date  12/12/12  MDG 3.1 test to do an actual numerical integration for the I_jk integrals, using MC profiles
!> @date  08/01/13  MDG 4.0 complete rewrite with Lambert format for MC output and new ctemlib.a routines 
!>                          also, the earlier versions would do only one energy value, whereas this new 
!>                          implementation does the complete energy-dependent master pattern
!> @date  01/27/14  MDG 4.1 continued rewrite, fixed problem with kvector list, replaced gvector routines
!>		   	     with updated routines; changed program name to CTEMEBSDmaster.
!> @date  05/03/14  MDG 4.2 test version to resolve bug in the Sgh matrix part (solved)
!> @date  06/19/14  MDG 4.3 rewrite, removal of all globals, split of namelist handling from computation; add OpenMP
!--------------------------------------------------------------------------
program CTEMEBSDmaster

use local
use NameListTypedefs
use NameListHandlers
use files
use io

IMPLICIT NONE

character(fnlen)                        :: nmldeffile, progname, progdesc
type(EBSDMasterNameListType)            :: emnl

nmldeffile = 'CTEMEBSDmaster.nml'
progname = 'CTEMEBSDmaster.f90'
progdesc = 'EBSD Energy-dependent Master Pattern Simulation'

! deal with the command line arguments, if any
call Interpret_Program_Arguments(nmldeffile,1,(/ 21 /), progname)

! deal with the namelist stuff
call GetEBSDMasterNameList(nmldeffile,emnl)

! print some information
call CTEMsoft(progname, progdesc)

! generate a set of master EBSD patterns
 call ComputeMasterPattern(emnl, progname)

end program CTEMEBSDmaster

!--------------------------------------------------------------------------
!
! SUBROUTINE:ComputeMasterPattern
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute an EBSD master pattern as a function of energy
!
!> @param nmlfile namelist file name
!
!> @date 11/29/01  MDG 1.0 original
!> @date 04/08/13  MDG 2.0 rewrite
!> @date 05/14/13  MDG 2.1 replaced IO by namelist file
!> @date 08/01/13  MDG 3.0 complete rewrite, eliminated old Lambert projection
!> @date 09/25/13  MDG 3.1 replaced k-vector code by kvectors module
!> @date 06/19/14  MDG 4.0 no more globals, nml split, added OpenMP
!--------------------------------------------------------------------------
subroutine ComputeMasterPattern(emnl, progname)

use typedefs
use NameListTypedefs
use initializers
use MBmodule
use symmetry
use crystal
use constants
use error
use gvectors
use kvectors
use io
use local
use files
use diffraction
use multibeams
use timing
use Lambert
use omp_lib

IMPLICIT NONE

type(EBSDMasterNameListType),INTENT(IN) :: emnl
character(fnlen),INTENT(IN)             :: progname


real(kind=dbl)          :: ctmp(192,3), arg, sig, omega
integer(kind=irg)      :: isym,i,j,ik,npy,ipx,ipy,debug,iE,izz, izzmax, iequiv(2,12), nequiv, num_el, MCnthreads, & ! counters
                        numk, & ! number of independent incident beam directions
                        ir,nat(100),kk(3), npyhex, skip, ijmax, one, NUMTHREADS, TID, &
                        numset,n,ix,iy,iz, io_int(6), nns, nnw, nref,  &
                        istat,gzero,ic,ip,ikk, totstrong, totweak     ! counters
real(kind=dbl)         :: tpi,Znsq, kkl, DBWF, kin !
real(kind=sgl)          :: io_real(5), selE, kn, FN(3), kkk(3)
real(kind=sgl),allocatable      :: sr(:,:,:,:), srhex(:,:,:,:), EkeVs(:), svals(:) ! results
complex(kind=dbl)               :: czero
complex(kind=dbl),allocatable   :: Lgh(:,:), Sgh(:,:,:)
logical                 :: usehex, switchmirror, verbose
character(fnlen)        :: xtalname
! the following will need to be moved elsewhere at some point...
integer(kind=irg),parameter     :: LaueTest(11) = (/ 149, 151, 153, 156, 158, 160, 161, 164, 165, 166, 167 /)  ! space groups with 2 or mirror at 30 degrees

! Monte Carlo derived quantities
integer(kind=irg)       :: numEbins, numzbins, nsx, nsy    ! variables used in MC energy file
real(kind=dbl)          :: EkeV, Ehistmin, Ebinsize, depthmax, depthstep, etotal ! enery variables from MC program
integer(kind=irg),allocatable :: accum_e(:,:,:), accum_z(:,:,:,:), thick(:)
real(kind=sgl),allocatable :: lambdaE(:,:)
character(fnlen)        :: oldprogname
character(8)            :: MCscversion
character(4)            :: MCmode

type(unitcell),pointer          :: cell
type(DynType),save              :: Dyn
type(gnode),save                :: rlp
type(reflisttype),pointer       :: reflist,firstw, rltmp
type(BetheParameterType)        :: BetheParameters
type(kvectorlist),pointer       :: khead, ktmp
real(kind=sgl),allocatable      :: karray(:,:)
integer(kind=irg),allocatable   :: kij(:,:)
complex(kind=dbl),allocatable   :: DynMat(:,:)

!$OMP THREADPRIVATE(rlp) 

! note that the CTEMMC.f90 program creates a statistical output file that 
! must be read by the present program, so that things like energy etc are 
! known to the program...  The content of the MC output file is as follows:
! write(dataunit) numEbins, numzbins, numsx, numsy
! write (dataunit) EkeV, Ehistmin, Ebinsize, depthmax, depthstep
! write(dataunit) accum_e
! write (dataunit) accum_z
! So that means that the beam voltage range is already known; we do not need to
! ask for this again...

tpi = 2.D0*cPi
czero = dcmplx(0.D0,0.D0)

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

read(dataunit) numEbins, numzbins, nsx, nsy, num_el ! , MCnthreads
nsx = (nsx - 1)/2
nsy = (nsy - 1)/2

! MCnthreads = 8
io_int(1:5) = (/ numEbins, numzbins, nsx, nsy, num_el /)
call WriteValue(' NumEbins, numzbins, nsx, nsy, num_el',io_int,5,"(4I8,',',I8)")
etotal = num_el 

read (dataunit) EkeV, Ehistmin, Ebinsize, depthmax, depthstep
io_real(1:5) = (/ EkeV, Ehistmin, Ebinsize, depthmax, depthstep /)
call WriteValue(' EkeV, Ehistmin, Ebinsize, depthmax, depthstep ',io_real,5,"(4F10.5,',',F10.5)")

read (dataunit) sig, omega
read (dataunit) MCmode

allocate(accum_e(numEbins,-nsx:nsx,-nsy:nsy),accum_z(numEbins,numzbins,-nsx/10:nsx/10,-nsy/10:nsy/10),stat=istat)
read(dataunit) accum_e
! actually, we do not yet need the accum_e array; that is for the actual EBSD image computation program
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

 verbose = .TRUE.
 call Initialize_Cell(cell,Dyn,rlp,xtalname, emnl%dmin, sngl(1000.0*EkeV), verbose)

! the following line needs to be verified ... 
!  hexset = .TRUE.    ! if hexagonal structure this switch selects between three and four index notation (4 if true)

! determine the point group number
 j=0
 do i=1,32
  if (SGPG(i).le.cell%SYM_SGnum) j=i
 end do
 isym = j

! and convert this to the corresponding  Laue point group number since diffraction 
! patterns are always centrosymmetric (hence, there are only 11 different cases for
! the symmetry of the incident beam).   
  isym = PGLaueinv(isym)  

! If the Laue group is # 7, then we need to determine the orientation of the mirror plane.
! The second orientation of the mirror plane is represented by "Laue group" # 12 in this program.
 switchmirror = .FALSE.
 if (isym.eq.7) then
  do i=1,11
    if (cell%SYM_SGnum.eq.LaueTest(i)) switchmirror = .TRUE.
  end do
 end if
 if (switchmirror) then
  isym = 12
  call Message(' Switching computational wedge to second setting for this space group', frm = "(A)")
 end if
 write (*,*) ' Laue group # ',isym, PGTHD(j)

! if this point group is trigonal or hexagonal, we need to switch usehex to .TRUE. so that
! the program will use the hexagonal sampling method
usehex = .FALSE.
if (((isym.ge.6).and.(isym.le.9)).or.(isym.eq.12)) usehex = .TRUE.
! ---------- end of symmetry and crystallography section
!=============================================
!=============================================

!=============================================
!=============================================
! this is where we determine the value for the thickness integration limit for the CalcLgh3 routine...
allocate(EkeVs(numEbins),thick(numEbins))

do i=1,numEbins
  EkeVs(i) = Ehistmin + float(i-1)*Ebinsize
end do

! then, for each energy determine the 95% histogram thickness
izzmax = 0
do iE = 1,numEbins
 do ix=-nsx/10,nsx/10
  do iy=-nsy/10,nsy/10
   istat = sum(accum_z(iE,:,ix,iy))
   izz = 1
   do while (sum(accum_z(iE,1:izz,ix,iy)).lt.(0.95*istat)) 
    izz = izz+1
   end do
   if (izz.gt.izzmax) izzmax = izz
  end do
 end do
 thick(iE) = dble(izzmax) * depthstep
end do

izz = nint(maxval(thick)/depthstep)
allocate(lambdaE(1:numEbins,1:izz),stat=istat)
do iE=1,numEbins
 do iz=1,izz
  lambdaE(iE,iz) = float(sum(accum_z(iE,iz,-nsx/10:nsx/10,-nsy/10:nsy/10)))/etotal
 end do
end do

! and get rid of the accum_z array
deallocate(accum_z)
! ---------- end of 'read Monte Carlo output file and extract necessary parameters' section
!=============================================
!=============================================


!=============================================
!=============================================
! ---------- a couple of initializations
   numset = cell % ATOM_ntype  
   npy = emnl%npx
   allocate(svals(numset),stat=istat)
   gzero = 1  ! index of incident beam
   debug = 0  ! no longer used
! ----------
!=============================================
!=============================================

!=============================================
!=============================================
! ---------- allocate memory for the master pattern
! we need to sample the stereographic projection Northern hemisphere or a portion
! thereoff, depending on the order of the Laue group.  There are 11 Laue groups, 
! which leads to 9 different shapes for the stereographic asymmetric unit for the 
! independent incident beam directions.  
! allocate space for the results (needs to be altered for general symmetry case)
if (emnl%Esel.eq.-1) then
  allocate(sr(-emnl%npx:emnl%npx,-npy:npy,1:numEbins,1:numset),stat=istat)
else
  allocate(sr(-emnl%npx:emnl%npx,-npy:npy,1,1:numset),stat=istat)
end if 

! in the trigonal/hexagonal case, we need intermediate storage arrays
  if (usehex) then
   npyhex = nint(2.0*float(npy)/sqrt(3.0))
   allocate(srhex(-emnl%npx:emnl%npx,-npyhex:npyhex,1:numEbins,1:numset),stat=istat)
  end if

! set various arrays to zero
   sr = 0.0
   if (usehex) then
     srhex = 0.0
   end if
! ---------- end allocate memory for the master pattern
!=============================================
!=============================================

! force dynamical matrix routine to read new Bethe parameters from file
! this will all be changed with the new version of the Bethe potentials
!  call Set_Bethe_Parameters(BetheParameters)

!=============================================
!=============================================
! ---------- from here on, we need to repeat the entire computation for each energy value
! so this is where we could in principle implement an OpenMP approach; alternatively, 
! we could do the inner loop over the incident beam directions in OpenMP (probably simpler)

energyloop: do iE=numEbins,1,-1
! is this a single-energy run ?
   if (emnl%Esel.ne.-1) then
     if (emnl%Esel.ne.iE) CYCLE energyloop
   end if
   
! print a message to indicate where we are in the computation
   io_int(1)=iE
   call Message('Starting computation for energy bin', frm = "(/A$)")
   call WriteValue(' ',io_int,1,"(I4$)")
   io_real(1) = EkeVs(iE)
   call WriteValue('; energy [keV] = ',io_real,1,"(F6.2/)")
   selE = EkeVs(iE)

! set the accelerating voltage
   skip = 3
   cell%voltage = dble(EkeVs(iE)*1000.0)
   call CalcWaveLength(cell, rlp, skip)

!=============================================
! ---------- create the incident beam directions list
! determine all independent incident beam directions (use a linked list starting at khead)
! numk is the total number of k-vectors to be included in this computation;
! note that this needs to be redone for each energy, since the wave vector changes with energy
   nullify(khead)
   if (usehex) then
    call Calckvectors(khead,cell, (/ 0.D0, 0.D0, 1.D0 /), (/ 0.D0, 0.D0, 0.D0 /),0.D0,emnl%npx,npyhex,numk, &
                isym,ijmax,'RoscaLambert',usehex)
   else 
! Calckvectors(k,ga,ktmax,npx,npy,numk,isym,ijmax,mapmode,usehex)
    call Calckvectors(khead,cell, (/ 0.D0, 0.D0, 1.D0 /), (/ 0.D0, 0.D0, 0.D0 /),0.D0,emnl%npx,npy,numk, &
                isym,ijmax,'RoscaLambert',usehex)
   end if
   io_int(1)=numk
   call WriteValue('# independent beam directions to be considered = ', io_int, 1, "(I8)")

! convert the kvector linked list into arrays for OpenMP
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

  verbose = .FALSE.
  totstrong = 0
  totweak = 0

! ---------- end of "create the incident beam directions list"
!=============================================

! here's where we introduce the OpenMP calls, to spead up the overall calculations...

! set the number of OpenMP threads 
  call OMP_SET_NUM_THREADS(emnl%nthreads)
  io_int(1) = emnl%nthreads
  call WriteValue(' Attempting to set number of threads to ',io_int, 1, frm = "(I4)")

! use OpenMP to run on multiple cores ... 
!!$OMP PARALLEL default(shared) COPYIN(rlp) &
!$OMP PARALLEL COPYIN(rlp) &
!$OMP& PRIVATE(DynMat,Sgh,Lgh,ik,FN,TID,kn,ipx,ipy,ix,iequiv,nequiv,reflist,firstw) &
!$OMP& PRIVATE(kkk,nns,nnw,nref,svals,nat,io_int)

  NUMTHREADS = OMP_GET_NUM_THREADS()
  TID = OMP_GET_THREAD_NUM()

!$OMP DO SCHEDULE(DYNAMIC,100)    
! ---------- and here we start the beam direction loop
   beamloop:do ik = 1,numk

!=============================================
! ---------- create the master reflection list for this beam direction
! Then we must determine the masterlist of reflections (also a linked list);
! This list basically samples a large reciprocal space volume; it does not 
! distinguish between zero and higher order Laue zones, since that 
! distinction becomes meaningless when we consider the complete 
! reciprocal lattice.  
     nullify(reflist)
!    kkk = karray(1:3,ik)
     kkk = karray(1:3,ik)
     FN = kkk
!    call TransSpace(cell,kkk,FN,'r','d')
!    call NormVec(cell,FN,'d')
!    FN = karray(1:3,ik)
     call Initialize_ReflectionList(cell, reflist, BetheParameters, FN, kkk, emnl%dmin, nref, verbose)
! ---------- end of "create the master reflection list"
!=============================================

! lines copied from CTEMKossel, may need to be modified ... 

! determine strong and weak reflections
     nullify(firstw)
     nns = 0
     nnw = 0
     call Apply_BethePotentials(cell, reflist, firstw, BetheParameters, nref, nns, nnw)

! generate the dynamical matrix
     allocate(DynMat(nns,nns))
     call GetDynMat(cell, reflist, firstw, rlp, DynMat, nns, nnw)
     totstrong = totstrong + nns
     totweak = totweak + nnw

! if (TID.eq.0) write (*,*) ik, FN, nref, nns, nnw

! then we need to initialize the Sgh and Lgh arrays
     if (allocated(Sgh)) deallocate(Sgh)
     if (allocated(Lgh)) deallocate(Lgh)
     allocate(Sgh(nns,nns,numset),Lgh(nns,nns))
     Sgh = czero
     Lgh = czero
     nat = 0
     call CalcSgh(cell,reflist,nns,numset,Sgh,nat)

! for now, we're disabling the kinematical part
! solve the dynamical eigenvalue equation for this beam direction  Lgh,thick,kn,nn,gzero,kin,debug
     kn = karray(4,ik)
     call CalcLgh(DynMat,Lgh,dble(thick(iE)),dble(kn),nns,gzero,depthstep,lambdaE(iE,1:izzmax),izzmax)
     deallocate(DynMat)

! dynamical contribution
     svals = 0.0
     do ix=1,numset
       svals(ix) = real(sum(Lgh(1:nns,1:nns)*Sgh(1:nns,1:nns,ix)))
     end do
     svals = svals/float(sum(nat(1:numset)))

! and store the resulting values, applying point group symmetry where needed.
     ipx = kij(1,ik)
     ipy = kij(2,ik)
     call Apply2DLaueSymmetry(ipx,ipy,isym,iequiv,nequiv)
!$OMP CRITICAL
     if (usehex) then
       do ix=1,nequiv
         srhex(iequiv(1,ix),iequiv(2,ix),iE,1:numset) = svals(1:numset)
        end do
     else
        if (emnl%Esel.eq.-1) then
         do ix=1,nequiv
           sr(iequiv(1,ix),iequiv(2,ix),iE,1:numset) = svals(1:numset)
          end do
        else
         do ix=1,nequiv
           sr(iequiv(1,ix),iequiv(2,ix),1,1:numset) = svals(1:numset)
          end do
        endif
     end if
!$OMP END CRITICAL
  
     if (mod(ik,2500).eq.0) then
       io_int(1) = ik
       call WriteValue('  completed beam direction ',io_int, 1, "(I8)")
     end if

     call Delete_gvectorlist(reflist)

    end do beamloop

! end of OpenMP portion
!$OMP END PARALLEL

  deallocate(karray, kij)

! stop the clock and report the total time     
!   call Time_stop(numk)

! write (*,*) 'Some statistics :'
! write (*,*) 'Average number of strong beams : ',float(BetheParameter%totstrong)/float(numk)
! write (*,*) '          (min,max) : ',BetheParameter%minstrong,BetheParameter%maxstrong
! write (*,*) 'Average number of weak beams : ',float(BetheParameter%totweak)/float(numk)
! write (*,*) '          (min,max) : ',BetheParameter%minweak,BetheParameter%maxweak


! Finally, if this was sampled on a hexagonal array, we need to do barycentric interpolation
! to the standard square array for final output and use of the subsequent program.
! [this interpolation scheme must be verified; it is possible that there is an off-by-one error somewhere ...]
!if (usehex) then
!  delta = dsqrt(2.D0)/dble(npx)
!  srt = 2.D0/dsqrt(3.D0)
!! copy the central row without modifications
!  sr(-npx:npx,0) = srhex(-npx:npx,0)
!  srkin(-npx:npx,0) = srkinhex(-npx:npx,0)
!! we'll go through the array with pairs of horizontal rows at a time
!  do j=1,npy-1
!! determine which way the triangle is oriented for this row of the square array
!    jh = floor(j*srt)
!    if (mod(jh,2).eq.0) then ! even numbers mean triangle points down
!      h = delta/srt - (j*delta - float(jh)*delta/srt)
!      lambda = 0.5D0 - h/delta/dsqrt(3.D0)
!      omtl = 1.D0-2.D0*lambda
!      do i=-npx+1,npx-1  ! perform the barycentric interpolation
!! positive row, pay attention to hexagonal coordinate transformation !
!	sr(i,j) = ( srhex(i-1,jh+1) + srhex(i,jh+1) )*lambda + omtl * srhex(i,jh)
!	srkin(i,j) = ( srkinhex(i-1,jh+1) + srkinhex(i,jh+1) )*lambda + omtl * srkinhex(i,jh)
!! negative row
!	sr(i,-j) = ( srhex(i-1,-jh-1) + srhex(i,-jh-1) )*lambda + omtl * srhex(i,-jh)
!	srkin(i,-j) = ( srkinhex(i-1,-jh-1) + srkinhex(i,-jh-1) )*lambda + omtl * srkinhex(i,-jh)
!      end do
!    else
!      h = j*delta - float(jh)*delta/srt
!      lambda = 0.5D0 - h/delta/dsqrt(3.D0)
!      omtl = 1.D0-2.D0*lambda
!      do i=-npx+1,npx-1  ! perform the barycentric interpolation
!! positive row, pay attention to hexagonal coordinate transformation !
!	sr(i,j) = ( srhex(i-1,jh) + srhex(i,jh) )*lambda + omtl * srhex(i,jh+1)
!	srkin(i,j) = ( srkinhex(i-1,jh) + srkinhex(i,jh) )*lambda + omtl * srkinhex(i,jh+1)
!! negative row
!	sr(i,-j) = ( srhex(i-1,-jh) + srhex(i,-jh) )*lambda + omtl * srhex(i,-jh-1)
!	srkin(i,-j) = ( srkinhex(i-1,-jh) + srkinhex(i,-jh) )*lambda + omtl * srkinhex(i,-jh-1)
!      end do
!    end if
!  end do
!end if

! since these computations can take a long time, here we store 
! all the output at the end of each pass through the energyloop.


  io_int(1) = nint(float(totstrong)/float(numk))
  call WriteValue(' -> Average number of strong reflections = ',io_int, 1, "(I5)")
  io_int(1) = nint(float(totweak)/float(numk))
  call WriteValue(' -> Average number of weak reflections   = ',io_int, 1, "(I5)")


  open(unit=dataunit,file=trim(emnl%outname),status='unknown',action='write',form = 'unformatted')
! write the program identifier
  write (dataunit) progname
! write the version number
  write (dataunit) scversion
! then the name of the crystal data file
  write (dataunit) xtalname
! then the name of the corresponding Monte Carlo data file
  write (dataunit) emnl%energyfile
! energy information and array size    
  if (emnl%Esel.eq.-1) then
    write (dataunit) emnl%npx,npy,numEbins,numset
    write (dataunit) EkeVs
  else
    one = 1
    write (dataunit) emnl%npx,npy,one,numset 
    write (dataunit) selE
  end if
! atom type array for asymmetric unit  
  write (dataunit) cell%ATOM_type(1:numset)
! is this a regular (square) or hexagonal projection ?
  if (usehex) then 
    write (dataunit) 'hexago'
  else
    write (dataunit) 'square'
  end if
! and finally the results array
  write (dataunit) sr
  close(unit=dataunit,status='keep')

 if ((emnl%Esel.eq.-1).and.(iE.ne.1)) then 
  call Message('Intermediate data stored in file '//trim(emnl%outname), frm = "(A/)")
 end if

 if ((emnl%Esel.eq.-1).and.(iE.eq.1)) then 
  call Message('Final data stored in file '//trim(emnl%outname), frm = "(A/)")
 end if

end do energyloop

if (emnl%Esel.ne.-1) then
  call Message('Final data stored in file '//trim(emnl%outname), frm = "(A/)")
end if

end subroutine ComputeMasterPattern
