! ###################################################################
! Copyright (c) 2014, Marc De Graef/Carnegie Mellon University
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
! CTEMsoft2013:CTEMPED.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMPED 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Zone axis Precession Electron Diffraction
!
!> @todo implement full symmetry use; implement full Bloch wave output
!>
!
!> @date 11/29/01 MDG 1.0 original
!> @date 04/08/13 MDG 2.0 rewrite
!> @date 05/08/13 MDG 2.1 forked from mbcbed and adapted for large angle CBED patterns
!> @date 05/14/13 MDG 2.2 replaced all IO by namelist file and added command line argument handling
!> @date 09/04/13 MDG 2.3 all command line argument handling now via files.f90 routine
!> @date 03/05/14 MDG 3.0 new version to deal with precession electron diffraction
!> @date 09/07/14 MDG 4.0 removal of all global variables; new TypeDef and NamelistHandler; OpenMP version
!--------------------------------------------------------------------------
program CTEMPED


use local
use NameListTypedefs
use NameListHandlers
use files
use io

IMPLICIT NONE

character(fnlen)                        :: nmldeffile, progname, progdesc
type(EBSDMasterNameListType)            :: pednl

nmldeffile = 'CTEMped.nml'
progname = 'CTEMped.f90'
progdesc = 'Precession Electron Diffraction Pattern Simulation'

! deal with the command line arguments, if any
call Interpret_Program_Arguments(nmldeffile,1,(/ 13 /), progname)

! deal with the namelist stuff
call GetPEDNameList(nmldeffile,pednl)

! print some information
call CTEMsoft(progname, progdesc)

! generate a set of master EBSD patterns
 call PEDPattern(pednl,progname)

end program CTEMPED

!--------------------------------------------------------------------------
!
! SUBROUTINE:PEDpattern
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute a zone axis precession electron diffraction pattern
!
!> @param nmlfile namelist file name
!
!> @date 03/05/14 MDG 1.0 original, based on CTEMlacbed program
!> @date 09/07/14 MDG 2.0 complete rewrite using new CTEMSoftLib library
!--------------------------------------------------------------------------
subroutine PEDpattern(pednl,progname)

use local
use typedefs
use NameListTypedefs
use initializers
use constants
use crystal
use diffraction
use gvectors
use kvectors
use MBmodule
use postscript, ONLY: GetIndex
use symmetry
use math
use io
use error
use files
use omp_lib

IMPLICIT NONE

type(PEDNameListType),INTENT(IN)        :: pednl
character(fnlen),INTENT(IN)             :: progname

real(kind=sgl)                  :: ktmax, io_real(3), bragg, thetac, sc, minten, pxy(2), galen, DM(2,2), DD, X(2), &
                                   frac, startthick, thick(1), klaue(2), thetam, kk(3) 
integer(kind=irg)               :: ijmax,ga(3),gb(3),cnt, PX, numthick, ss, icnt, pgnum, ih, nunique, famnum, &
                                   newcount,count_rate,count_max, io_int(6), i, j, isym, ir, skip, ghkl(3), &
                                   npx, npy, numt, numk, ik, ip, jp, istat, dgn, nbeams, refcnt, &
                                   ifamily, famhkl(3), inum, maxHOLZ, numksame
character(3)                    :: method

real(kind=sgl),allocatable      ::  disk(:,:)
integer(kind=irg),allocatable   :: familymult(:), familyhkl(:,:), whichHOLZ(:), gequiv(:,:)
real(kind=sgl),allocatable      :: inten(:,:)
real(kind=dbl)                  :: s(3)
logical,allocatable             :: ksame(:)

type(unitcell),pointer          :: cell
type(DynType),save              :: Dyn
type(gnode),save                :: rlp
type(reflisttype),pointer       :: reflist,firstw, rltmpa
type(BetheParameterType)        :: BetheParameters
type(kvectorlist),pointer       :: khead, ktmp
real(kind=sgl),allocatable      :: karray(:,:)
integer(kind=irg),allocatable   :: kij(:,:)
complex(kind=dbl),allocatable   :: DynMat(:,:)


!=============================================
!=============================================
! crystallography section
nullify(cell)
allocate(cell)

verbose = .TRUE.
call Initialize_Cell(cell,Dyn,rlp,pednl%xtalname, pednl%dmin, pednl%voltage, verbose)

! set the foil normal 
 Dyn%FN = float(pednl%fn)
 call NormVec(cell, Dyn%FN, 'd')
 numt = 1
 thick(1) = pednl%thickness
 
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
 call ShortestG(pednl%k,ga,gb,isym)
 io_int(1:3)=ga(1:3)
 io_int(4:6)=gb(1:3)
 call WriteValue(' Reciprocal lattice vectors : ', io_int, 6,"('(',3I3,') and (',3I3,')',/)")

! construct the list of all possible reflections
! method = 'ALL'
! thetac = precangle/1000.0
! call Compute_ReflectionList(dmin,pednl%k,ga,gb,method,.FALSE.,maxHOLZ,thetac)
! galen = CalcLength(float(ga),'r')
!=============================================
!=============================================

!=============================================
!=============================================
! determine the list of contributing wave vectors
  call CalckvectorsPrecession(khead,cell,dble(pednl%k),dble(ga),pednl%precangle,pednl%prechalfwidth,pednl%precsample,&
        pednl%precazimuthal,numk)
!=============================================
!=============================================

!=============================================
!=============================================
! up to this point, everything is nearly identical to the mbcbed program,
! except that we do not use a camera length explicitly.  Now we
! need to do things a little differently.  First of all, we need a master 
! list for all the reflections that contribute in one way or another to the
! complete PED pattern.  We do this by going through the entire list of incident
! beam directions and flagging all reflections that contribute, either as weak
! or as strong reflections.
!  mess = ' Pruning reflection list (this takes a while ...) '
!  call Message("(A)")
!  call Prune_ReflectionList(numk,nbeams)
!  io_int(1) = nbeams
!  call WriteValue('Number of contributing beams  : ', io_int, 1, '(I8)')

! set the scale parameter for a default camera length of 1000 mm.
  sc = cell%mLambda * 1000.0 * 300.0 / 25.4  ! the absolute value does not matter and is derived from legacy Postscript code
! The original code used 300 dpi (hence 300/25.4) which was convenient for Postscript output; in the current case, we
! do not actually use the true value, but in the IDL visualization program, we scale the user defined camera length by
! 1000.0, and use this ratio to scale the diskoffset coordinates.  So, there's no absolute length scale, only a relative scale.


! next we create the output array (not sure yet how to do this)
  allocate(disk(-npix:npix,-npix:npix))
  disk=0.0

  io_int(1)=numk
  call WriteValue('Starting computation for # beam directions = ', io_int, 1, "(I8)")

! time the computation
  cnt = 0
  call system_clock(cnt,count_rate,count_max)

! go through the linked reflection list and set all intensities to zero; we'll use the 
! xg entry in the rltmpa linked list to accumulate the intensities

if (pednl%filemode.eq.'total') then
  rltmpa => reflist%next
  do i=1,DynNbeamsLinked
    rltmpa%xg = 0.D0
    rltmpa => rltmpa%next
  end do
end if


! the final bit of the program involves dumping all the results into a file,
! binary for now, but HDF5 in the future, for the IDL visualization program 
! to read.
  open(unit=dataunit,file=trim(outname),status='unknown',action='write',form='unformatted')
! write the program identifier
  write (dataunit) trim(progname)
! write the version number
  write (dataunit) scversion
! then the name of the crystal data file
  write (dataunit) pednl%xtalname
! the accelerating voltage [V]
  write (dataunit) pednl%voltage
! precangle angle [mrad]
  write (dataunit) pednl%precangle
! prechalfwidth angle [mrad]
  write (dataunit) pednl%prechalfwidth
! precsample
  write (dataunit) pednl%precsample
! precazimuthal
  write (dataunit) pednl%precazimuthal
! the zone axis indices
  write (dataunit) pednl%k
! the foil normal indices
  write (dataunit) pednl%fn
! number of k-values in disk
  write (dataunit) numk
! dmin value
  write (dataunit) pednl%dmin
! horizontal reciprocal lattice vector
  write (dataunit) ga  
! second reciprocal lattice vector
  write (dataunit) gb
! length horizontal reciprocal lattice vector (need for proper Laue center coordinate scaling)
  write (dataunit) galen
! length second reciprocal lattice vector (need for proper Laue center coordinate scaling)
  write (dataunit) CalcLength(cell,float(gb),'r')
! angle between these two vectors (in radians)
  write (dataunit) CalcAngle(cell,float(ga),float(gb),'r')
! thickness data
  write (dataunit) pednl%thickness
  
 DM(1,1) = CalcDot(cell, float(gb),float(gb),'c')
 DM(1,2) = -CalcDot(cell, float(ga),float(gb),'c')
 DM(2,1) = DM(1,2)
 DM(2,2) = CalcDot(cell, float(ga),float(ga),'c')
 DD = 1.0/(DM(1,1)*DM(2,2) - DM(1,2)*DM(2,1))

! first convert the linked wave vector list to a set of regular arrays, so that we can 
! use OpenMP for this computation
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

! set up OpenMP parameters
  io_int(1)=numk
  call WriteValue('Starting computation for # beam directions = ', io_int, 1, "(I8)")
  gzero = 1
  
! time the computation
  cnt = 0
  call system_clock(cnt,count_rate,count_max)

! create the list of potential reflections; this is actually the special PED routine that uses
! an analytical expression to determine which reflections will be part of the swept Ewald sphere 
! volume.  For each one, we set the xg parameter ot zero, and this parameter will accumulate the 
! intensities for all wave vectors on the precession cone.
  goffset = pmult * 2.0/pednl%thickness    ! reciprocal lattice point offset based on relrod expression 
  verbose = .TRUE.
  call Initialize_ReflectionList(cell, reflist, FN, pednl%k, nref, pednl%pedangle, goffset, verbose)

! set the number of OpenMP threads 
  call OMP_SET_NUM_THREADS(pednl%nthreads)
  io_int(1) = pednl%nthreads
  call WriteValue(' Attempting to set number of threads to ',io_int, 1, frm = "(I4)")

! use OpenMP to run on multiple cores ... 
!!$OMP PARALLEL default(shared) PRIVATE(

! set the foil normal 

  NUMTHREADS = OMP_GET_NUM_THREADS()
  TID = OMP_GET_THREAD_NUM()

  nullify(reflist)
  nullify(firstw)

  nns = 0
  nnw = 0

! loop over all beam orientations, selecting them from the linked list; but skipping the first one, which is 
! the center axis of the cone...
kvectorloop:  do ik = 1,numk

! generate the reflectionlist
        kk(1:3) = karray(1:3,ik)
        FN = kk
        call Initialize_ReflectionList(cell, reflist, BetheParameters, FN, kk, pednl%dmin, nref, verbose)

! determine strong and weak reflections
        call Apply_BethePotentials(cell, reflist, firstw, BetheParameters, nref, nns, nnw)

! generate the dynamical matrix
        allocate(DynMat(nns,nns))
        call GetDynMat(cell, reflist, firstw, rlp, DynMat, nns, nnw)

! allocate the intensity array to include both strong beams and weak beams (in that order)
        allocate(inten(numt,nns))
        inten = 0.0
 
! solve the dynamical eigenvalue equation for this beam direction
        kn = karray(4,ik)
        CalcPEDint(DynMat,cell,kn,nns,nt,thick,inten)
        deallocate(DynMat)

! now we need to write these intensities to an appropriate array ... 

! and remove the intensity array
     deallocate(inten)

! remove all the computed intensities for per pattern storage    
if (filemode.eq.'eachp') then

! first we count how many reflections have none-zero intensity  
  refcnt = 0
  rltmpa => reflist%next
  do i=1,DynNbeamsLinked
    if (rltmpa%xg.ne.0.D0) refcnt = refcnt + 1
    rltmpa => rltmpa%next
  end do

! write refcnt
  write (dataunit) refcnt
  rltmpa => reflist%next
  do i=1,DynNbeamsLinked
    if (rltmpa%xg.ne.0.D0) then
! decompose this point w.r.t ga and gb
      X(1) = CalcDot(float(rltmpa%hkl),float(ga),'c')
      X(2) = CalcDot(float(rltmpa%hkl),float(gb),'c')
      X = matmul(DM,X) * DD
      write (dataunit) rltmpa%hkl
      write (dataunit) rltmpa%xg
      write (dataunit) X
    end if
    rltmpa => rltmpa%next
  end do
  
! and reset the intensities for the next run
  rltmpa => reflist%next
  do i=1,DynNbeamsLinked
    rltmpa%xg = 0.D0
    rltmpa => rltmpa%next
  end do
end if

! select next beam direction
   if (ik.ne.numk) ktmp => ktmp%next

! update computation progress
   if (float(ik)/float(numk) .gt. frac) then
    io_int(1) = nint(100.0*frac) 
    call WriteValue('       ', io_int, 1, "(1x,I3,' percent completed')") 
    frac = frac + 0.05
   end if  

  end do kvectorloop



! stop the clock and report the total time     
  call system_clock(newcount,count_rate,count_max)
  io_real(1)=float(newcount-cnt)/float(count_rate)
  mess = ' Program run completed '; call Message("(/A/)")
  call WriteValue('Total computation time [s] ' , io_real, 1, "(F10.5)")


! remove all the computed intensities for per pattern storage    
if (filemode.eq.'total') then

! first we count how many reflections have none-zero intensity  
  refcnt = 0
  rltmpa => reflist%next
  do i=1,DynNbeamsLinked
    if (rltmpa%xg.ne.0.D0) refcnt = refcnt + 1
    rltmpa => rltmpa%next
  end do

! write refcnt
write (*,*) "LINE 402 DEBUG" !PGC
  write (dataunit) refcnt
  rltmpa => reflist%next
  do i=1,DynNbeamsLinked
    if (rltmpa%xg.ne.0.D0) then
! decompose this point w.r.t ga and gb
      X(1) = CalcDot(float(rltmpa%hkl),float(ga),'c')
      X(2) = CalcDot(float(rltmpa%hkl),float(gb),'c')
      X = matmul(DM,X) * DD
      write (dataunit) rltmpa%hkl
      write (dataunit) rltmpa%xg
      write (dataunit) X
    end if
    rltmpa => rltmpa%next
  end do
end if
  
  close(UNIT=dataunit,STATUS='keep')

  mess = ' Data stored in '//outname; call Message("(/A/)") 

end subroutine PEDpattern




!--------------------------------------------------------------------------
!
! SUBROUTINE: CalckvectorsPrecession
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief create a linked list of wave vectors for precession electron diffraction
!
!> @details This is a new version to test whether or not we can use the whole pattern
!> symmetry to determine the relevant list of incident wave vectors; this should be a 
!> general routine, so that we do not need to consider each symmetry case separately.
!> This will require a floating point version of the Apply2DPGSymmetry routine in symmetry.f90.
!
!> @param khead wave vector list pointer
!> @param cell unit cell pointer
!> @param k central wave vector
!> @param ga reciprocal lattice vector normal to k
!> @param precangle main precession angle [mrad]
!> @param prechalfwidth precession beam half width [mrad]
!> @param precsample number of samples in half width 
!> @param precazimuthal number of samples around each precession circle
!> @param numk total number of wave vectors in list
!
!> @date  03/05/14 MDG 1.0 original
!> @date  11/28/14 MDG 2.0 rewrite without global variables
!--------------------------------------------------------------------------
recursive subroutine CalckvectorsPrecession(khead,cell,k,ga,precangle,prechalfwidth,precsample,precazimuthal,numk)

use local
use typedefs
use io
use error
use constants
use diffraction
use crystal
use Lambert
use kvectors

IMPLICIT NONE

type(kvectorlist),pointer,INTENT(IN)    :: khead
type(unitcell),pointer,INTENT(IN)       :: cell
real(kind=dbl),INTENT(IN)               :: k(3)         !< initial wave vector
real(kind=dbl),INTENT(IN)               :: ga(3)        !< "horizontal" reciprocal lattice vector
real(kind=sgl),INTENT(IN)               :: precangle    !< precession angle in [mrad]
real(kind=sgl),INTENT(IN)               :: prechalfwidth        !< halfwidth of tilted beam [mrad]
integer(kind=irg),INTENT(IN)            :: precsample           !< number of kvectors along beam tilt
integer(kind=irg),INTENT(IN)            :: precazimuthal        !< number of kvectors along circumference
integer(kind=irg),INTENT(OUT)           :: numk         !< total number of kvectors in linked list

type(kvectorlist),pointer               :: ktail
integer(kind=irg)                       :: istat,i,j, iequiv(2,12), nequiv, jj, nx, ny, il, ith
real(kind=dbl)                          :: gp, dgp, glen, gan(3), gperp(3), kstar(3), dth
real(kind=dbl),allocatable              :: gw(:), ct(:), st(:), th(:)
logical                                 :: hexgrid = .FALSE.
real(kind=sgl)                          :: kt(3),kr(3)
real(kind=sgl)                          :: ktlen

write (*,*) precangle, prechalfwidth, precsample, precazimuthal, cell%mLambda

! compute geometrical factors 
 glen = CalcLength(cell,ga,'r')                 ! length of ga
 gan = ga/glen                                  ! normalized ga
 gp = 2.0*sin(precangle/1000.0)/cell%mLambda         ! precession angle converted to reciprocal length gp in units of glen
 dgp = 0.0
 if (precsample.gt.0) then
   dgp = 2.0*sin(0.001*(precangle-prechalfwidth))/cell%mLambda/glen/float(precsample)        ! half width step size converted to reciprocal length dgp in units of glen
 end if
 allocate(gw(2*precsample+1))                           ! sampling radii
 gw = gp + dgp * (/ (i,i=-precsample,precsample) /)     ! sampling radii

! pre-compute cosines and sines
 allocate(ct(precazimuthal),st(precazimuthal), th(precazimuthal))
 dth = 2.D0*cPi / dble(precazimuthal)
 th = (/ (i-1,i=1,precazimuthal) /) * dth
 ct = cos(th)
 st = sin(th)
 
 call TransSpace(cell,k,kstar,'d','r')               ! transform incident direction to reciprocal space
 call CalcCross(cell,ga,kstar,gperp,'r','r',0)       ! compute g_perp = ga x k
 call NormVec(cell,gperp,'r')                        ! normalize g_perp
 call NormVec(cell,kstar,'r')                        ! normalize reciprocal beam vector

! allocate the head and tail of the linked list
 allocate(khead,stat=istat)                             ! allocate new value
 if (istat.ne.0) call FatalError('CalckvectorsPrecession','unable to allocate khead pointer')
 ktail => khead                                         ! tail points to new value
 nullify(ktail%next)                                    ! nullify next in new value
 numk = 0                                               ! keep track of number of k-vectors so far

 
! next loop around each of the precession circles
 do il = 1,2*precsample+1                               ! number of concentric circles
  do ith = 1,precazimuthal                              ! number of points along each circle
! make a new one in the list, except for the first one
   if (numk.ne.0) then
     allocate(ktail%next,stat=istat)                    ! allocate new value
     if (istat.ne.0) call FatalError('Add_knode:',' unable to allocate pointer')
     ktail => ktail%next                        ! tail points to new value
     nullify(ktail%next)                                ! nullify next in new value
   end if
! and populate the fields   
   kt = - gw(il)*ct(ith)*gan - gw(il)*st(ith)*gperp     ! tangential component of k
   ktail%kt = kt                                        ! store tangential component of k
   ktlen = CalcLength(cell,kt,'r')**2                   ! squared length of tangential component

   kr = kt + sqrt(1.0/cell%mLambda**2 - ktlen)*kstar         ! complete wave vector
   ktail%k = kr                                         ! store in pointer list
   ktail%kn = CalcDot(cell,ktail%k,kstar,'r')           ! normal component of k
   numk = numk + 1
  end do
 end do


end subroutine CalckvectorsPrecession

