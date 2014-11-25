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
! CTEMsoft2013:CTEMECPthinfilm.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMECPthinfilm 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Zone axis electron channeling patterns
!
!> @date 03/18/10 MDG 1.0 f90
!> @date 08/09/10 MDG 2.0 corrected weight factors and g-vector ordering problem
!> @date 11/18/13 MDG 3.0 major rewrite with new libraries 
!> @date 06/27/14 MDG 4.0 removal of all globals; separation of namelist handling from computation
!> @date 06/30/14 MDG 4.1 added OpenMP
!> @date 11/25/14 MDG 5.0 forked from old CTEMECP to start film+substrate implementation
!--------------------------------------------------------------------------
program CTEMECPthinfilm

use local
use typedefs
use NameListTypedefs
use NameListHandlers
use files
use io

IMPLICIT NONE

character(fnlen)                        :: nmldeffile, progname, progdesc
type(ECPNameListType)                   :: ecpnl

nmldeffile = 'CTEMECPthinfilm.nml'
progname = 'CTEMECPthinfilm.f90'
progdesc = 'Electron channeling pattern simulation for film on substrate'

! deal with the command line arguments, if any
call Interpret_Program_Arguments(nmldeffile,2,(/ 0, 40 /), progname)

! deal with the namelist stuff
call GetECPNameList(nmldeffile,ecpnl)

! print some information
call CTEMsoft(progname, progdesc)

! perform the zone axis computations (fos = film on substrate)
call ECpatternfos(ecpnl, progname)

end program CTEMECPthinfilm

!--------------------------------------------------------------------------
!
! SUBROUTINE:ECpatternfos
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute a large angle zone axis electron channeling pattern
!
!> @note This is really very similar to a LACBED computation, except that 
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
!> @date 11/25/14  MDG 3.0 new version for film on substrate
!--------------------------------------------------------------------------
subroutine ECpatternfos(ecpnl, progname)

use local
use typedefs
use NameListTypedefs
use crystal
use math
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

character(3)                    :: method
real(kind=sgl)                  :: galen, bragg, klaue(2), io_real(6), kstar(3), gperp(3), delta, thetac, FNstar(3), gg(3), sg, &
                                   kk(3), ktmax, FN(3), kn, fnat, fnatS, TTinv(3,3), r(3), p(3)
integer(kind=irg)               :: nt, skip, dgn, pgnum, io_int(6), maxHOLZ, ik, numk, ga(3), gb(3), TID, izfilm, &
                                   nn, npx, npy, isym, numset, numsetS, it, ijmax, jp, istat, iequiv(2,12), nequiv, NUMTHREADS, & 
                                   nns, nnw, nref, tots, totw, ikf
real(kind=dbl)                  :: ctmp(192,3),arg, lambdaF, lambdaS
integer                         :: i,j,ir, n,ipx,ipy,gzero,ic,ip,ikk
real(kind=sgl)                  :: pre, tpi,Znsq, kkl, DBWF, frac
real,allocatable                :: thick(:), sr(:,:,:) 
complex(kind=dbl),allocatable   :: Lgh(:,:,:),Sgh(:,:),Sghtmp(:,:,:)
complex(kind=dbl)               :: czero
real(kind=dbl)    		:: dE(3,3)
real(kind=sgl),allocatable      :: karray(:,:)
integer(kind=irg),allocatable   :: kij(:,:), nat(:), natS(:)
complex(kind=dbl),allocatable   :: DynMat(:,:)
logical                         :: verbose

type(orientation) 		:: orel
type(unitcell),pointer          :: cell, cellS
type(gnode)                     :: rlp, rlpS
type(DynType)                   :: Dyn, DynS
type(kvectorlist),pointer       :: khead, ktmp
type(symdata2D)                 :: TDPG, TDPGS
type(BetheParameterType)        :: BetheParameters
type(reflisttype),pointer       :: reflist, firstw, rltmp, reflistS, firstwS, rltmpS


! init some parameters
  gzero = 1
  frac = 0.05

  nullify(cell)
  nullify(khead)
  nullify(ktmp)
  nullify(cellS)

! ======================================
! ======================================
! read parameters for the film; note that there is no way to use symmetry in this case !

! load the crystal structure and compute the Fourier coefficient lookup table
  allocate(cell)
  verbose = .TRUE.
  call Initialize_Cell(cell,Dyn,rlp,ecpnl%xtalname, ecpnl%dmin, ecpnl%voltage,verbose)
  lambdaF = mLambda

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
  call WriteValue(' Film reciprocal lattice vectors : ', io_int, 6,"('(',3I3,') and (',3I3,')',/)")

! ======================================
! ======================================
! read parameters for the substrate

! load the substrate crystal structure and compute the Fourier coefficient lookup table
  allocate(cellS)
  call Initialize_Cell(cellS,DynS,rlpS,ecpnl%xtalname2, ecpnl%dmin, ecpnl%voltage, verbose)
  lambdaS = mLambda
  mLambda = lambdaF

! ======================================
! ======================================
! define orientation relation between two lattices and get transformation matrix
  orel%gA = float(ecpnl%gF)
  orel%tA = float(ecpnl%tF)
  orel%gB = float(ecpnl%gS)
  orel%tB = float(ecpnl%tS)

 ! check for orthonormality using zone equation
  if (sum(orel%tA*orel%gA).ne.0.0) call FatalError('ECPpatternfos','Plane does not contain direction (Film)')
  if (sum(orel%tB*orel%gB).ne.0.0) call FatalError('ECPpatternfos','Plane does not contain direction (Substrate)')

! and compute the (inverse) transformation matrix
  TTinv = ComputeOR(orel, cell, cellS, 'BA')
 ! note that matmul(TTinv, t) transforms a direct space vector t in the film to the 
 ! equivalent direct space vector in the substrate.
! ======================================
! ======================================

! here is how we can use the transformation matrix from film to substrate to 
! determine the vectors k_g in the substrate reciprocal reference frame

! reciprocal foil normal
  FN = ecpnl%fn
  call TransSpace(cell,FN,FNstar,'d','r')
  call NormVec(cell,FNstar,'r')
   


! ======================================
! ======================================
! next, we need to read the Monte Carlo file to get the lambdaE curves etc ... 



! ======================================
! ======================================

! diffraction geometry
  bragg = CalcDiffAngle(cell,ga(1),ga(2),ga(3))*0.5
  if (ecpnl%ktmax.ne.0.0) then 
    thetac = (ecpnl%ktmax * 2.0 * bragg)*1000.0
    ktmax = ecpnl%ktmax
    io_real(1) = thetac
  else
    ktmax = ecpnl%thetac / (2000.0 * bragg)
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

! set the thickness array (note that there is a special thickness for the film)
  nt = ecpnl%numthick
  allocate(thick(nt))
  thick = ecpnl%startthick + ecpnl%thickinc * (/ (i-1,i=1,nt) /)
! we also need to determine at which index the integration over the film ends.
  izfilm = int(ecpnl%filmthickness/depthstep) 
  
! ======================================
! ======================================

!----------------------------MAIN COMPUTATIONAL LOOP-----------------------
  czero = cmplx(0.D0,0.D0)
  pre = cmplx(0.D0,1.D0) * cPi
  tpi = 2.D0*cPi

! film
  numset = cell % ATOM_ntype  ! number of special positions in the unit cell
  allocate(nat(numset))
  nat = 0
  fnat = 1.0/float(sum(cell%numat(1:numset)))
! substrate
  numsetS = cellS % ATOM_ntype  ! number of special positions in the unit cell
  allocate(natS(numsetS))
  natS = 0
  fnatS = 1.0/float(sum(cellS%numat(1:numsetS)))

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

!  work through the beam direction list for the Film
  beamloop: do ik=1,numk

! generate the reflectionlist
        kk(1:3) = karray(1:3,ik)
        mLambda = lambdaF
!       FN = kk
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
        call CalcSghFilm(cell,reflist,nns,numset,Sghtmp,nat)

! sum Sghtmp over the sites
        Sgh = sum(Sghtmp,3)

! solve the dynamical eigenvalue equation; this is identical to what we have for the single layer
! version, but we need to end the integration/summation at the point where the substrate begins,
! instead of going to the largest depth.  So, we need a new CalcLghECP routine that takes the 
! lambda weight function, as in the EBSD case, and also carries out the secondary summation/integration
! over the substrate reflections.
        kn = karray(4,ik)
        call CalcLghFilm(DynMat,Lgh,ecpnl%filmthickness,kn,nns,gzero,depthstep,lambdaE,izfilm,amps)
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
        deallocate(Lgh, Sgh) 
          
! and here we need to do the loop over all the beams entering the substrate; we'll limit ourselves to
! the strong Film beams only, and for each one we'll compare the modulus of the amplitude to a threshold
! value to see if we should take it into account or not.  
        rltmpS => reflist%next  ! point to the first strong beam
        mLambda = lambdaS
        filmbeamloop: do ikf=1,nns
          bamp = cabs(amps(ikf))
          if (bamp.gt.amplitudethreshold) then  ! we do need to consider this incident beam
! get the incident beam direction in the substrate reciprocal reference frame
            kg = Convert_kgs_to_Substrate(cell, cellS, kk+float(rltmpS%hkl)+Calcsg(cell,float(rltmpS%hkl),kk,FNstar)*FNstar, &
                    TTinv, lambdaS) 
            call Initialize_ReflectionList(cellS, reflistS, BetheParameters, FNstar, kg, ecpnl%dmin, nrefS)
! determine strong and weak reflections
            call Apply_BethePotentials(cellS, reflistS, firstwS, BetheParameters, nrefS, nnsS, nnwS)
! generate the dynamical matrix
            allocate(DynMat(nnsS,nnsS))
            call GetDynMat(cellS, reflistS, firstwS, rlpS, DynMat, nnsS, nnwS)
! then we need to initialize the Sgh and Lgh arrays
            if (allocated(SghS)) deallocate(SghS)
            if (allocated(LghS)) deallocate(LghS)
            if (allocated(SghtmpS)) deallocate(SghtmpS)
            allocate(SghtmpS(nnsS,nnsS,numsetS),LghS(nnsS,nnsS,ntS),Sgh(nnsS,nnsS))
            SghS = czero
            SghtmpS = czero
            LghS = czero
            natS = 0
! compute Sgh matrix for the substrate
            call CalcSghSubstrate(cellS,reflistS,nnsS,numsetS,SghtmpS,natS)
! sum Sghtmp over the sites
            SghS = sum(SghtmpS,3)
! compute the Lgh matrix; fix normal component of k and integration thickness here !!!!!!!!!
            kn = karray(4,ik)
            call CalcLghSubstrate(DynMat,LghS,ecpnl%filmthickness,kn,nnsS,gzero,depthstep,lambdaE,izfilm,amp)
            deallocate(DynMat,SghtmpS)
! add to integrated signal from film
!$OMP CRITICAL
            if (isym.ne.1) then 
              do ip=1,nequiv
               do jp=1,ntS
                sr(iequiv(1,ip),iequiv(2,ip),jp) = sr(iequiv(1,ip),iequiv(2,ip),jp) + &
                  real(sum(LghS(1:nnsS,1:nnsS,jp)*SghS(1:nnsS,1:nnsS)))
               end do
              end do
            else
               do jp=1,nt
                sr(ipx+npx+1,ipy+npy+1,jp) = sr(ipx+npx+1,ipy+npy+1,jp)+real(sum(Lgh(1:nns,1:nns,jp)*Sgh(1:nns,1:nns)))
               end do
            end if
!$OMP END CRITICAL
            deallocate(LghS, SghS) 
          end if ! consider this beam ?
          rltmpS => rltmpS%nexts
        end do filmbeamloop



  ! these lines are to be used later on...
! incident beam direction into wave vector in reciprocal frame 
! p = float(ecpnl%k)
! call TransSpace(cell,p,r,'d','r')
! call NormVec(cell,r,'r')
! kstar = r/lambdaF

! transform k_g=k_0+g+s_g FN to Substrate reciprocal reference frame for the given orientation relation




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

end subroutine ECpatternfos



!--------------------------------------------------------------------------
!
! SUBROUTINE: CalcSghFilm
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute structure factor-like array for EBSD, ECCI and ECP simulations
!
!> @param cell unit cell pointer
!> @param nn dimension of array
!> @param Sgh output array
!> @param nat normalization array
!
!> @date 03/05/14  MDG 1.0 original (used to be in-line in ECP and ECCI programs)
!> @date 03/11/14  MDG 1.1 converted to diagonal Sgh array only
!> @date 06/19/14  MDG 2.0 no globals, taken out of CTEMECCI.f90
!> @date 11/25/14  MDG 3.0 forked from original for film-substrate case
!--------------------------------------------------------------------------
recursive subroutine CalcSghFilm(cell,reflist,nn,numset,Sgh,nat)

use local
use typedefs
use crystal
use gvectors
use constants
use symmetry

IMPLICIT NONE

type(unitcell),pointer                  :: cell
type(reflisttype),pointer               :: reflist
integer(kind=irg),INTENT(IN)            :: nn
integer(kind=irg),INTENT(IN)            :: numset
complex(kind=dbl),INTENT(INOUT)         :: Sgh(nn,nn,numset)
integer(kind=irg),INTENT(INOUT)         :: nat(numset)

integer(kind=irg)                       :: ip, ir, ic, kkk(3), ikk, n
real(kind=sgl)                          :: Znsq, DBWF, kkl
complex(kind=dbl)                       :: carg
real(kind=dbl)                          :: ctmp(192,3),arg, tpi
type(reflisttype),pointer               :: rltmpa, rltmpb

  tpi = 2.D0 * cPi

! for each special position we need to compute its contribution to the Sgh array
  do ip=1,cell % ATOM_ntype
!   call CalcOrbit(cell,ip,n,ctmp)
    nat(ip) = cell%numat(ip)
! get Zn-squared for this special position, and include the site occupation parameter as well
    Znsq = float(cell%ATOM_type(ip))**2 * cell%ATOM_pos(ip,4)
! loop over all contributing reflections
! ir is the row index
    rltmpa => reflist%next    ! point to the front of the list
    do ir=1,nn
! ic is the column index
      rltmpb => reflist%next    ! point to the front of the list
      do ic=1,nn
        kkk = rltmpb%hkl - rltmpa%hkl
! We'll assume isotropic Debye-Waller factors for now ...
! That means we need the square of the length of s=  kk^2/4
        kkl = 0.25 * CalcLength(cell,float(kkk),'r')**2
! Debye-Waller exponential times Z^2
        DBWF = Znsq * exp(-cell%ATOM_pos(ip,5)*kkl)
! here is where we insert the proper weight factor, Z^2 exp[-M_{h-g}]
! and also the detector geometry...   For now, we do nothing with the detector
! geometry; the Rossouw et al 1994 paper lists a factor A that does not depend
! on anything in particular, so we assume it is 1. 
        do ikk=1,n
! get the argument of the complex exponential
          arg = tpi*sum(kkk(1:3)*cell%apos(ip,ikk,1:3))
          carg = dcmplx(dcos(arg),dsin(arg))
! multiply with the prefactor and add
          Sgh(ir,ic,ip) = Sgh(ir,ic,ip) + carg * dcmplx(DBWF,0.D0)
        end do
        rltmpb => rltmpb%nexts  ! move to next column-entry
      end do
     rltmpa => rltmpa%nexts  ! move to next row-entry
   end do  
  end do
 
end subroutine CalcSghFilm



! ###################################################################
! 
! SUBROUTINE: CalcLghFilm
! 
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute the Film Lgh matrix for EBSD, ECCI, ECP, etc simulations
!
!> @details This routine is different from the regular CalcLgh in that it must
!> also return the wave amplitudes for all scattered beams; these are then 
!> used in the summation over the substrate contributions...
!
!> @param DMat dynamical matrix
!> @param Lgh output array
!> @param thick integration thickness of the Film only !
!> @param kn normal wave vector component
!> @param nn number of strong beams
!> @param gzero index of incident beam (should always be 1)
!> @param depthstep depth step size from MC simulation
!> @param lambdaE energy weight factors
!> @param izz number of energy weight factors to consider for the Film
!
!> @date 10/13/98  MDG 1.0 original
!> @date 07/04/01  MDG 2.0 f90
!> @date 06/19/14  MDG 3.0 no globals
!> @date 06/23/14  MDG 4.0 moved to MBmodule
!> @date 11/25/14  MDG 5.0 forked for film-substrate case
! ###################################################################
recursive subroutine CalcLghFilm(DMat,Lgh,thick,kn,nn,gzero,depthstep,lambdaE,izz,amps)

use local
use io
use files
use constants
use error

IMPLICIT NONE

complex(kind=dbl),INTENT(IN)        :: DMat(nn,nn)
complex(kind=dbl),INTENT(OUT)       :: Lgh(nn,nn)
real(kind=dbl),INTENT(IN)           :: thick
real(kind=dbl),INTENT(IN)           :: kn
integer(kind=irg),INTENT(IN)        :: nn
integer(kind=irg),INTENT(IN)        :: gzero
real(kind=dbl),INTENT(IN)           :: depthstep
real(kind=sgl),INTENT(IN)           :: lambdaE(izz)
integer(kind=irg),INTENT(IN)        :: izz
complex(kind=dbl),INTENT(OUT)       :: amps(nn)

integer                             :: i,j, iz
complex(kind=dbl)                   :: CGinv(nn,nn), Minp(nn,nn), tmp3(nn,nn)

real(kind=dbl)                      :: tpi, dzt
complex(kind=dbl)                   :: Ijk(nn,nn), q, getMIWORK, qold

integer(kind=irg)                   :: INFO, LDA, LDVR, LDVL,  JPIV(nn), MILWORK
complex(kind=dbl)                   :: CGG(nn,nn), W(nn)
complex(kind=dbl),allocatable       :: MIWORK(:)

integer(kind=irg),parameter         :: LWMAX = 5000 
complex(kind=dbl)                   :: VL(nn,nn),  WORK(LWMAX)
real(kind=dbl)                      :: RWORK(2*nn)
character                           :: JOBVL, JOBVR
integer(kind=sgl)                   :: LWORK

! compute the eigenvalues and eigenvectors
! using the LAPACK CGEEV, CGETRF, and CGETRI routines
! 
 Minp = DMat

! set some initial LAPACK variables 
 LDA = nn
 LDVL = nn
 LDVR = nn
 INFO = 0
 
! first initialize the parameters for the LAPACK ZGEEV, CGETRF, and CGETRI routines
 JOBVL = 'N'   ! do not compute the left eigenvectors
 JOBVR = 'V'   ! do compute the right eigenvectors
 LWORK = -1 ! so that we can ask the routine for the actually needed value

! call the routine to determine the optimal workspace size
  call zgeev(JOBVL,JOBVR,nn,Minp,LDA,W,VL,LDVL,CGG,LDVR,WORK,LWORK,RWORK,INFO)
  LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

! then call the eigenvalue solver
  call zgeev(JOBVL,JOBVR,nn,Minp,LDA,W,VL,LDVL,CGG,LDVR,WORK,LWORK,RWORK,INFO)
  if (INFO.ne.0) call FatalError('Error in CalcLghFilm: ','ZGEEV return not zero')

 CGinv = CGG

 call zgetrf(nn,nn,CGinv,LDA,JPIV,INFO)
 MILWORK = -1
 call zgetri(nn,CGinv,LDA,JPIV,getMIWORK,MILWORK,INFO)
 MILWORK =  INT(real(getMIWORK))
 if (.not.allocated(MIWORK)) allocate(MIWORK(MILWORK))
 MIWORK = dcmplx(0.D0,0.D0)
 call zgetri(nn,CGinv,LDA,JPIV,MIWORK,MILWORK,INFO)
 deallocate(MIWORK)

! then compute the integrated intensity matrix
! the following scaling step should really be done when we construct the dynamical matrix:
 W = W/cmplx(2.0*kn,0.0)

! recall that alpha(1:nn) = CGinv(1:nn,gzero)

! first the Ijk matrix (this is Winkelmann's B^{ij}(t) matrix)
! combined with numerical integration over [0, z0] interval,
! taking into account depth profiles from Monte Carlo simulations ...
! the depth profile lambdaE must be added to the absorption 
! components of the Bloch wave eigenvalues.

tpi = 2.D0*cPi*depthstep
dzt = depthstep/thick
 do i=1,nn
  do j=1,nn
     q =  cmplx(0.D0,0.D0)
     qold = tpi * dcmplx(aimag(W(i))+aimag(W(j)),real(W(i))-real(W(j)))
     do iz = 1,izz
       q = q + dble(lambdaE(iz)) * cdexp( - qold * dble(iz) ) !MNS changed cexp to cdexp to be compatible with gfortran
     end do
     Ijk(i,j) = conjg(CGinv(i,gzero)) * q * CGinv(j,gzero)
  end do
 end do

Ijk = Ijk * dzt

! then the summations for Lgh and kin
tmp3 = matmul(CGG,transpose(Ijk)) 
Lgh = matmul(tmp3,transpose(conjg(CGG)))

! and finally we need to compute the array of amplitudes for all scattered beams at the film-substrate interface



end subroutine CalcLghFilm


!--------------------------------------------------------------------------
