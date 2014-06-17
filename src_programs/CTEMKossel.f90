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
! CTEMsoft2013:CTEMKossel.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMKossel 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Kossel patterns, taken from master EBSD pattern simulation, with adjustment for 
!> Bloch wave part ...
!
!> @todo implement OpenMP multithreading for the actual computation part; requires modifications
!> in CTEMlib.a routines (mostly THREADPRIVATE commands in several modules)
!
!> @details So, in this new version, not only are we getting rid of all globals,
!> we're also attempting to split the program into two parts, so that the actual
!> computation routine can be called from elsewhere (i.e., no namelist handling
!> in the main computation part).
!
!> @date 04/08/11 MDG 1.0 early version
!> @date 01/07/14 MDG 2.0 new version
!> @date 06/13/14 MDG 3.0 rewrite without globals ... 
!> @date 06/13/14 MDG 3.1 separation of nml handling from computation part
!--------------------------------------------------------------------------
program CTEMKossel

use local
use files
use NameListTypedefs
use NameListHandlers
use io

IMPLICIT NONE

character(fnlen)                        :: nmldeffile, progname, progdesc
type(KosselNameListType)                :: knl

nmldeffile = 'CTEMKossel.nml'
progname = 'CTEMKossel.f90'
progdesc = 'Zone axis Kossel pattern simulation'

! deal with the command line arguments, if any
call Interpret_Program_Arguments(nmldeffile,1,(/ 12 /), progname)

! deal with the namelist stuff
call GetKosselNameList(nmldeffile,knl)

! print some information
call CTEMsoft(progname, progdesc)

! perform the zone axis computations for the knl input parameters
call ComputeKosselpattern(knl, progname)

end program CTEMKossel

!--------------------------------------------------------------------------
!
! SUBROUTINE:ComputeKosselpattern
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief computation of a Kossel diffraction pattern
!
!> @param knl Kossel name list structure
!> @param progname program name
!
!> @date 01/07/14  MDG 1.0 complete rewrite
!> @date 06/13/14  MDG 2.0 another complete rewrite, removal of namelist stuff
!> @date 06/16/14  MDG 2.1 debug and removal of unnecessary code; OpenMP implementation
!--------------------------------------------------------------------------
subroutine ComputeKosselpattern(knl,progname)

use local
use typedefs
use NameListTypedefs
use initializers
use symmetry
use Lambert, ONLY: Apply2DPGSymmetry
use crystal
use constants
use io
use error
use files
use gvectors
use kvectors
use diffraction
use multibeams
use postscript
use timing
use MBmodule
!use omp_lib

IMPLICIT NONE

type(KosselNameListType),INTENT(IN)     :: knl
character(fnlen),INTENT(IN)             :: progname

real(kind=sgl)                  :: ktmax, io_real(3), bragg, thetac,  galen, delta, kstar(3), gperp(3), &
                                   frac, klaue(2), thetam, kk(3), kn
integer(kind=irg)               :: ijmax,ga(3),gb(3),cnt, pgnum, NUMTHREADS, TID, ki, kj, &
                                   newcount,count_rate,count_max, io_int(6), i, j, isym, skip, &
                                   npx, npy, numt, numk, ik, ip, jp, istat, dgn, &
                                   numset, nn, gzero, ipx, ipy, ii, iequiv(2,12), nequiv
character(3)                    :: method


real(kind=sgl),allocatable      :: thickarray(:)
real(kind=sgl),allocatable      :: Iz(:), Izsum(:,:,:)
real(kind=sgl),allocatable      :: karray(:,:)
integer(kind=irg),allocatable   :: kij(:,:)
real(kind=dbl)                  :: pre, tpi
complex(kind=dbl)               :: czero
logical                         :: verbose

type(unitcell),pointer          :: cell
type(gnode)                     :: rlp
type(DynType)                   :: Dyn
type(kvectorlist),pointer       :: khead, ktmp
type(symdata2D)                 :: TDPG
type(BetheParameterType)        :: BetheParameters

interface
        subroutine GetDynMat(cell, Dyn)
        
        use local
        use typedefs
        use io
        use crystal
        use diffraction
        use kvectors
        use gvectors
        use constants
        
        IMPLICIT NONE
        
        type(unitcell),pointer           :: cell
        type(DynType),INTENT(INOUT)      :: Dyn
        end subroutine GetDynMat
end interface

  allocate(cell)

  verbose = .TRUE.
  call Initialize_Cell(cell,Dyn,knl%xtalname, knl%dmin, knl%voltage, verbose)

! set the foil normal 
  Dyn%FN = float(knl%fn)
 
! determine the point group number
  j=0
  do i=1,32
   if (SGPG(i).le.cell % SYM_SGnum) j=i
  end do

! use the new routine to get the whole pattern 2D symmetry group, since that
! is the one that determines the independent beam directions.
  dgn = GetPatternSymmetry(cell,knl%k,j,.TRUE.)
  pgnum = j
  isym = WPPG(dgn) ! WPPG lists the whole pattern point group numbers vs. diffraction group numbers

! determine the shortest reciprocal lattice points for this zone
  call ShortestG(cell,knl%k,ga,gb,isym)
  io_int(1:3)=ga(1:3)
  io_int(4:6)=gb(1:3)
  call WriteValue(' Reciprocal lattice vectors : ', io_int, 6,"('(',3I3,') and (',3I3,')',/)")
 
 ! for some of the 2D point groups, the standard orientation of the group according to ITC vol A
! may not be the orientation that we have here, so we need to determine by how much the 2D point
! group is rotated (CCW) with respect to the standard setting...
  call CheckPatternSymmetry(cell,knl%k,ga,isym,thetam)


! determine range of incident beam directions
  bragg = CalcDiffAngle(cell,ga(1),ga(2),ga(3))*0.5
  
! convert to ktmax along ga
  ktmax = 0.5*knl%convergence*0.001/bragg

! the number of pixels across the disk is equal to 2*npix + 1
  npx = knl%npix
  npy = npx
  
! set thicknesses for which to compute the intensities
  numt = knl%numthick
  allocate(thickarray(numt),stat=istat)
  thickarray = knl%startthick + knl%thickinc* (/ (float(i),i=0,numt-1) /)

! set parameters for wave vector computation
  klaue = (/ 0.0, 0.0 /)
  ijmax = float(npx)**2   ! truncation value for beam directions

! get the linked list of wave vectors
  call CalckvectorsSymmetry(khead,cell,TDPG,dble(knl%k),dble(ga),dble(ktmax),npx,npy,numk,isym,ijmax,klaue,.FALSE.)

! force dynamical matrix routine to read new Bethe parameters from file
  call Set_Bethe_Parameters(BetheParameters,.TRUE.)

!----------------------------MAIN COMPUTATIONAL LOOP-----------------------
  czero = cmplx(0.D0,0.D0)
  pre = cmplx(0.D0,1.D0) * cPi
  numset = cell % ATOM_ntype  ! number of special positions in the unit cell
  tpi = 2.D0*cPi
! allocate space for the results
  allocate(Iz(knl%numthick),Izsum(2*npx+1,2*npy+1,knl%numthick))
  
! BetheParameters%c2 = 40.0
! BetheParameters%c1 = 40.0 
  verbose = .FALSE.

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

  io_int(1)=numk
  call WriteValue('Starting computation for # beam directions = ', io_int, 1, "(I8)")
  gzero = 1
  
! time the computation
  cnt = 0
  call system_clock(cnt,count_rate,count_max)

! set the number of OpenMP threads 
!  call OMP_SET_NUM_THREADS(knl%nthreads)
!  io_int(1) = knl%nthreads
!  call WriteValue(' Setting number of threads to ',io_int, 1, frm = "(I4)")

!! use OpenMP to run on multiple cores ... 
!!$OMP PARALLEL default(shared) PRIVATE(ik,TID,kk,kn,ipx,ipy,ii,iequiv,nequiv,ip,jp,myreflist)

!  NUMTHREADS = OMP_GET_NUM_THREADS()
!  TID = OMP_GET_THREAD_NUM()

! each thread needs to have its own reflist linked list !!!
!  allocate(mycell)
! and copy some of the contents of cell into mycell
!  call CellCopy(cell,mycell)


!!$OMP DO SCHEDULE(STATIC,1)    

!  work through the beam direction list
  beamloop: do ik=1,numk

! generate the reflectionlist
        kk(1:3) = karray(1:3,ik)
        call Initialize_ReflectionList(cell, BetheParameters, Dyn, kk, knl%dmin, verbose)

! determine strong and wel reflections
        call Apply_BethePotentials(cell, BetheParameters)

! generate the dynamical matrix
        call GetDynMat(cell, Dyn)

! solve the dynamical eigenvalue equation for this beam direction
        kn = karray(4,ik)
        call CalcKint(Dyn, cell, kn, cell%nns, knl%numthick, thickarray, Iz)

        ipx = kij(1,ik)
        ipy = kij(2,ik)

! apply pattern symmetry 
        call Apply2DPGSymmetry(TDPG,ipx,ipy,isym,iequiv,nequiv)

!!$OMP CRITICAL
        do ii=1,nequiv
! is this point inside the viewing square ?
          ip = iequiv(1,ii) + npx + 1
          jp = iequiv(2,ii) + npy + 1
          if (((ip.ge.1).and.(ip.le.2*npx+1)).and.((jp.ge.1).and.(jp.le.2*npx+1))) then
           Izsum(ip,jp,1:knl%numthick) = Iz(1:knl%numthick)
          end if
        end do
!!$OMP END CRITICAL
      
! select next beam direction
!       if (ik.ne.numk) ktmp => ktmp%next

! update computation progress
        if (float(ik)/float(numk) .gt. frac) then
          io_int(1) = nint(100.0*frac) 
          call WriteValue('       ', io_int, 1, "(1x,I3,' percent completed')") 
          frac = frac + 0.05
        end if  
        
  end do beamloop
    

!!$OMP END PARALLEL


! stop the clock and report the total time     
  call system_clock(newcount,count_rate,count_max)
  io_real(1)=float(newcount-cnt)/float(count_rate)
  call Message(' Program run completed ', frm = "(/A/)")
  call WriteValue('Total computation time [s] ' , io_real, 1, "(F10.3)")
  
  
! store additional information for the IDL interface  
  open(unit=dataunit,file=trim(knl%outname),status='unknown',action='write',form='unformatted')
! write the program identifier
  write (dataunit) trim(progname)
! write the version number
  write (dataunit) scversion
! first write the array dimensions
  write (dataunit) 2*npx+1,2*npy+1,knl%numthick
! then the name of the crystal data file
  write (dataunit) knl%xtalname
! the accelerating voltage [V]
  write (dataunit) knl%voltage
! convergence angle [mrad]
  write (dataunit) knl%convergence
! max kt value in units of ga
  write (dataunit) ktmax
! the zone axis indices
  write (dataunit) knl%k
! the foil normal indices
  write (dataunit) knl%fn
! number of k-values in disk
  write (dataunit) numk
! dmin value
  write (dataunit) knl%dmin
! horizontal reciprocal lattice vector
  write (dataunit) ga  
! length horizontal reciprocal lattice vector (need for proper Laue center coordinate scaling)
  galen = CalcLength(cell,float(ga),'r')
  write (dataunit) galen
! we need to store the gperp vectors
  delta = 2.0*ktmax*galen/float(2*npx+1)             ! grid step size in nm-1 
  call TransSpace(cell,float(knl%k),kstar,'d','r')        ! transform incident direction to reciprocal space
  call CalcCross(cell,float(ga),kstar,gperp,'r','r',0)! compute g_perp = ga x k
  call NormVec(cell,gperp,'r')                        ! normalize g_perp
  write (dataunit) delta
  write (dataunit) gperp
! eight integers with the labels of various symmetry groups
  write (dataunit) (/ pgnum, PGLaue(pgnum), dgn, PDG(dgn), BFPG(dgn), WPPG(dgn), DFGN(dgn), DFSP(dgn) /)
! thickness data
  write (dataunit) knl%startthick, knl%thickinc
! and finaly the full pattern array
  write (dataunit) Izsum ! sr  
  close(UNIT=dataunit,STATUS='keep')
  
  call Message(' Data stored in '//trim(knl%outname), frm = "(/A/)") 

end subroutine ComputeKosselpattern

!--------------------------------------------------------------------------
!
! SUBROUTINE: GetDynMat
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute the dynamical matrix, including Bethe potentials
!
!> @param 
!
!> @date  04/22/14 MDG 1.0 new library version
!> @date  06/15/14 MDG 2.0 updated for removal of globals
!--------------------------------------------------------------------------
recursive subroutine GetDynMat(cell, Dyn)

use local
use typedefs
use io
use crystal
use diffraction
use kvectors
use gvectors
use constants

IMPLICIT NONE

type(unitcell),pointer           :: cell
type(DynType),INTENT(INOUT)      :: Dyn

complex(kind=dbl)                :: czero, ughp, uhph, weaksum 
real(kind=dbl)                   :: weaksgsum
integer(kind=sgl)                :: ir, ic, ll(3), istat
type(reflisttype),pointer        :: rlr, rlc, rlw
type(gnode)                      :: rlp

czero = cmplx(0.0,0.0,dbl)      ! complex zero

! create the dynamical matrix for his reference case
! get the absorption coefficient
        if (allocated(Dyn%DynMat)) deallocate(Dyn%DynMat)
        allocate(Dyn%DynMat(cell%nns,cell%nns),stat=istat)
        Dyn%DynMat = czero
        call CalcUcg(cell, rlp, (/0,0,0/) )
        Dyn%Upz = rlp%Vpmod

        rlr => cell%reflist%next
        ir = 1
        do
          if (.not.associated(rlr)) EXIT
          rlc => cell%reflist%next
          ic = 1
          do
          if (.not.associated(rlc)) EXIT
          if (ic.ne.ir) then  ! not a diagonal entry
! here we need to do the Bethe corrections if necessary
            if (cell%nnw.ne.0) then
              weaksum = czero
              rlw => cell%firstw
              do
               if (.not.associated(rlw)) EXIT
               ll = rlr%hkl - rlw%hkl
               ughp = cell%LUT(ll(1),ll(2),ll(3)) 
               ll = rlw%hkl - rlc%hkl
               uhph = cell%LUT(ll(1),ll(2),ll(3)) 
               weaksum = weaksum +  ughp * uhph *cmplx(1.D0/rlw%sg,0.0,dbl)
               rlw => rlw%nextw
              end do
!        ! and correct the dynamical matrix element to become a Bethe potential coefficient
              ll = rlr%hkl - rlc%hkl
              Dyn%DynMat(ir,ic) = cell%LUT(ll(1),ll(2),ll(3))  - cmplx(0.5D0*mLambda,0.0D0,dbl)*weaksum
             else
              ll = rlr%hkl - rlc%hkl
              Dyn%DynMat(ir,ic) = cell%LUT(ll(1),ll(2),ll(3))
            end if
          else  ! it is a diagonal entry, so we need the excitation error and the absorption length
! determine the total contribution of the weak beams
            if (cell%nnw.ne.0) then
              weaksgsum = 0.D0
              rlw => cell%firstw
              do
               if (.not.associated(rlw)) EXIT
                ll = rlr%hkl - rlw%hkl
                ughp = cell%LUT(ll(1),ll(2),ll(3)) 
                weaksgsum = weaksgsum +  cdabs(ughp)**2/rlw%sg
                rlw => rlw%nextw
              end do
              weaksgsum = weaksgsum * mLambda/2.D0
              Dyn%DynMat(ir,ir) = cmplx(2.D0*rlr%sg/mLambda-weaksgsum,Dyn%Upz,dbl)
            else
              Dyn%DynMat(ir,ir) = cmplx(2.D0*rlr%sg/mLambda,Dyn%Upz,dbl)
            end if           
        
           end if       
           rlc => rlc%nexts
           ic = ic + 1
          end do        
          rlr => rlr%nexts
          ir = ir+1
        end do

end subroutine GetDynMat
