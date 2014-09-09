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
! CTEMsoft2013:CTEMKosselmaster.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMKosselmaster
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief CTEMKosselmaster computes the master Kossel pattern for a given structure
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
!> @date  09/09/14  MDG 5.0 forked from CTEMEBSDmaster to CTEMKosselmaster program
!--------------------------------------------------------------------------
program CTEMKosselmaster

use local
use NameListTypedefs
use NameListHandlers
use files
use io

IMPLICIT NONE

character(fnlen)                        :: nmldeffile, progname, progdesc
type(KosselMasterNameListType)          :: kmnl

nmldeffile = 'CTEMKosselmaster.nml'
progname = 'CTEMKosselmaster.f90'
progdesc = 'Kossel Master Pattern Simulation'

! deal with the command line arguments, if any
call Interpret_Program_Arguments(nmldeffile,1,(/ 13 /), progname)

! deal with the namelist stuff
call GetKosselMasterNameList(nmldeffile,kmnl)

! print some information
call CTEMsoft(progname, progdesc)

! generate a series of master Kossel patterns
 call ComputeKosselMasterPattern(kmnl, progname)

end program CTEMKosselmaster

!--------------------------------------------------------------------------
!
! SUBROUTINE:ComputeKosselMasterPattern
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
!> @date 09/09/14  MDG 5.0 new version for Kossel master pattern
!--------------------------------------------------------------------------
subroutine ComputeKosselMasterPattern(kmnl, progname)

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

type(KosselMasterNameListType),INTENT(IN) :: kmnl
character(fnlen),INTENT(IN)               :: progname


real(kind=dbl)          :: ctmp(192,3), arg
integer(kind=irg)      :: isym,i,j,ik,npy,ipx,ipy,debug,izz, izzmax, iequiv(2,12), nequiv, num_el, MCnthreads, & ! counters
                        numk, & ! number of independent incident beam directions
                        ir,kk(3), npyhex, skip, ijmax, one, NUMTHREADS, TID, &
                        n,ix,iy, io_int(6), nns, nnw, nref,  &
                        istat,gzero,ic,ip,ikk, totstrong, totweak     ! counters
real(kind=dbl)         :: tpi,Znsq, kkl, DBWF, kin !
real(kind=sgl)          :: io_real(5), selE, kn, FN(3), kkk(3)
real(kind=sgl),allocatable      :: sr(:,:,:), srhex(:,:,:), Iz(:), thick(:)
complex(kind=dbl)               :: czero
logical                 :: usehex, switchmirror, verbose
! the following will need to be moved elsewhere at some point...
integer(kind=irg),parameter     :: LaueTest(11) = (/ 149, 151, 153, 156, 158, 160, 161, 164, 165, 166, 167 /)  ! space groups with 2 or mirror at 30 degrees

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

tpi = 2.D0*cPi
czero = dcmplx(0.D0,0.D0)

!=============================================
!=============================================
! crystallography section

nullify(cell)
allocate(cell)

 verbose = .TRUE.
 call Initialize_Cell(cell,Dyn,rlp, kmnl%xtalname, kmnl%dmin, kmnl%voltage, verbose)

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
! this is where we determine the values for the thicknesses 
allocate(thick(kmnl%numthick),Iz(kmnl%numthick))

do i=1,kmnl%numthick
  thick(i) = kmnl%startthick + float(i-1)*kmnl%thickinc
end do
!=============================================
!=============================================


!=============================================
!=============================================
! ---------- a couple of initializations
   npy = kmnl%npix
   gzero = 1  ! index of incident beam
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
allocate(sr(-kmnl%npix:kmnl%npix,-npy:npy,1:kmnl%numthick),stat=istat)

! in the trigonal/hexagonal case, we need intermediate storage arrays
  if (usehex) then
   npyhex = nint(2.0*float(npy)/sqrt(3.0))
   allocate(srhex(-kmnl%npix:kmnl%npix,-npyhex:npyhex,1:kmnl%numthick),stat=istat)
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
! print a message to indicate we're starting the computation
call Message('Starting computation', frm = "(/A)")

!=============================================
! ---------- create the incident beam directions list
! determine all independent incident beam directions (use a linked list starting at khead)
! numk is the total number of k-vectors to be included in this computation;
! note that this needs to be redone for each energy, since the wave vector changes with energy
   nullify(khead)
   if (usehex) then
    call Calckvectors(khead,cell, (/ 0.D0, 0.D0, 1.D0 /), (/ 0.D0, 0.D0, 0.D0 /),0.D0,kmnl%npix,npyhex,numk, &
                isym,ijmax,'RoscaLambert',usehex)
   else 
! Calckvectors(k,ga,ktmax,npx,npy,numk,isym,ijmax,mapmode,usehex)
    call Calckvectors(khead,cell, (/ 0.D0, 0.D0, 1.D0 /), (/ 0.D0, 0.D0, 0.D0 /),0.D0,kmnl%npix,npy,numk, &
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
  call OMP_SET_NUM_THREADS(kmnl%nthreads)
  io_int(1) = kmnl%nthreads
  call WriteValue(' Attempting to set number of threads to ',io_int, 1, frm = "(I4)")

! use OpenMP to run on multiple cores ... 
!!$OMP PARALLEL default(shared) COPYIN(rlp) &
!$OMP PARALLEL COPYIN(rlp) &
!$OMP& PRIVATE(DynMat,ik,FN,TID,kn,ipx,ipy,ix,iequiv,nequiv,reflist,firstw) &
!$OMP& PRIVATE(kkk,nns,nnw,nref,io_int,Iz)

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
     kkk = karray(1:3,ik)
     FN = kkk
     call Initialize_ReflectionList(cell, reflist, BetheParameters, FN, kkk, kmnl%dmin, nref, verbose)
! ---------- end of "create the master reflection list"
!=============================================

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

! for now, we're disabling the kinematical part
! solve the dynamical eigenvalue equation for this beam direction  Lgh,thick,kn,nn,gzero,kin,debug
     kn = karray(4,ik)
     call CalcKint(DynMat,kn,nns,kmnl%numthick,thick,Iz)
     deallocate(DynMat)

! and store the resulting values, applying point group symmetry where needed.
     ipx = kij(1,ik)
     ipy = kij(2,ik)
     call Apply2DLaueSymmetry(ipx,ipy,isym,iequiv,nequiv)
!$OMP CRITICAL
     if (usehex) then
       do ix=1,nequiv
         srhex(iequiv(1,ix),iequiv(2,ix),1:kmnl%numthick) = Iz(1:kmnl%numthick)
        end do
     else
         do ix=1,nequiv
           sr(iequiv(1,ix),iequiv(2,ix),1:kmnl%numthick) = Iz(1:kmnl%numthick)
          end do
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

!  we need to convert srhex to sr !!!


  open(unit=dataunit,file=trim(kmnl%outname),status='unknown',action='write',form = 'unformatted')
! write the program identifier
  write (dataunit) progname
! write the version number
  write (dataunit) scversion
! then the name of the crystal data file
  write (dataunit) kmnl%xtalname
! pattern size parameter
  write (dataunit) kmnl%npix
! thickness parameters
  write (dataunit) kmnl%numthick
  write (dataunit) kmnl%startthick, kmnl%thickinc
! is this a regular (square) or hexagonal projection ?
  if (usehex) then 
    write (dataunit) 'hexago'
  else
    write (dataunit) 'square'
  end if
! and finally the results array
  write (dataunit) sr
  close(unit=dataunit,status='keep')

  call Message('Final data stored in file '//trim(kmnl%outname), frm = "(A/)")

end subroutine ComputeKosselMasterPattern
