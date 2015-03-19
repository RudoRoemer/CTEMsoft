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
! CTEMsoft2013:CTEMEBSD.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMEBSD
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief CTEMEBSD computes energy-weighted EBSD patterns
!
!> @date  08/01/12  MDG 1.0 EBSD extraction program for fundamental zone patterns
!> @date  08/17/12  MDG 1.1 generalized fundamental zone to other symmetries
!> @date  09/20/12  MDG 1.2 adapted for Lambert projection
!> @date  09/25/12  MDG 1.3 prepared for multithreaded version by separating computation steps
!> @date  12/11/12  MDG 2.0 new branch with energy-dependent Lambert projections (cubic only for now)
!> @date  02/26/14  MDG 3.0 incorporation into git and adapted to new libraries
!> @date  03/26/14  MDG 3.1 modification of file formats; made compatible with IDL visualization interface
!> @date  06/24/14  MDG 4.0 removal of all global variables; separation of nml from computation; OpenMP
!> @date  03/10/15  MDG 4.1 added output format selector
! ###################################################################
! 

program CTEMEBSD

use local
use files
use NameListTypedefs
use NameListHandlers
use io
use EBSDmod

IMPLICIT NONE

character(fnlen)                       :: nmldeffile, progname, progdesc
type(EBSDNameListType)                 :: enl
type(MCNameListType)                   :: mcnl

type(EBSDAngleType),pointer            :: angles
type(EBSDLargeAccumType),pointer       :: acc
type(EBSDMasterType),pointer           :: master


integer(kind=irg)                      :: istat
logical                                :: verbose

interface
        subroutine ComputeEBSDPatterns(enl, angles, acc, master, progname)
        
        use local
        use typedefs
        use NameListTypedefs
        use symmetry
        use crystal
        use constants
        use io
        use files
        use diffraction
        use EBSDmod
        use Lambert
        use quaternions
        use rotations
        use noise
        
        IMPLICIT NONE
        
        type(EBSDNameListType),INTENT(IN)       :: enl
        type(EBSDAngleType),pointer             :: angles
        type(EBSDLargeAccumType),pointer        :: acc
        type(EBSDMasterType),pointer            :: master
        character(fnlen),INTENT(IN)             :: progname
        end subroutine ComputeEBSDPatterns
end interface


nmldeffile = 'CTEMEBSD.nml'
progname = 'CTEMEBSD.f90'
progdesc = 'Dynamical EBSD patterns, using precomputed MC and master Lambert projections'

! deal with the command line arguments, if any
call Interpret_Program_Arguments(nmldeffile,1,(/ 21 /), progname)

! deal with the namelist stuff
call GetEBSDNameList(nmldeffile,enl)

! print some information
call CTEMsoft(progname, progdesc)

! this program needs a lot of data, and it also should be integrated 
! with Dream.3D, so we need to make sure that all data is loaded outside
! of the main computational routine, and passed in as pointers/arguments
! either by the fortran program or by Dream.3D calls.  

! 1. read the angle array from file
verbose = .TRUE.
allocate(angles)
call EBSDreadangles(enl, angles, verbose)

! 2. read the Monte Carlo data file
allocate(acc)
call EBSDreadMCfile(enl, acc, verbose)

! 3. read EBSD master pattern file
allocate(master)
call EBSDreadMasterfile(enl, master, verbose)

! 4. generate detector arrays
allocate(master%rgx(enl%numsx,enl%numsy), master%rgy(enl%numsx,enl%numsy), master%rgz(enl%numsx,enl%numsy), stat=istat)
allocate(acc%accum_e_detector(enl%numEbins,enl%numsx,enl%numsy), stat=istat)
call EBSDGenerateDetector(enl, acc, master, verbose)
deallocate(acc%accum_e)

! perform the zone axis computations for the knl input parameters
call ComputeEBSDpatterns(enl, angles, acc, master, progname)

deallocate(master, acc, angles)

end program CTEMEBSD

!--------------------------------------------------------------------------
!
! SUBROUTINE:ComputeEBSDPatterns
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute an energy-weighted EBSD pattern
!
!> @param enl name list
!> @param angles angle structure
!> @param acc energy accumulator arrays
!> @param master structure with master and detector arrays
!> @param progname program name string
!
!> @date 11/29/01  MDG 1.0 original
!> @date 04/08/13  MDG 2.0 rewrite
!> @date 05/14/13  MDG 2.1 replaced IO by namelist file
!> @date 08/01/13  MDG 3.0 complete rewrite, eliminated old Lambert projection
!> @date 09/25/13  MDG 3.1 replaced k-vector code by kvectors module
!> @date 02/26/14  MDG 4.0 new version
!> @date 03/26/14  MDG 4.1 adapted to new input and out file formats
!> @date 05/22/14  MDG 4.2 slight modification of angle input file; update for new CTEMEBSDMaster file format
!> @date 06/24/14  MDG 5.0 removal of global variables; removal of namelist stuff; 
!> @date 03/09/15  MDG 5.1 added OpenMP functionality for final loop
!> @date 03/10/15  MDG 5.2 added 'bin' and 'gui' outputformat; added mask support for 'bin' outputformat
!> @date 03/14/15  MDG 5.3 attempt at speeding up the program by performing approximate energy sums (energyaverage = 1)
!--------------------------------------------------------------------------
subroutine ComputeEBSDPatterns(enl, angles, acc, master, progname)

use local
use typedefs
use NameListTypedefs
use symmetry
use crystal
use constants
use io
use files
use diffraction
use EBSDmod
use Lambert
use quaternions
use rotations
use noise
use omp_lib

IMPLICIT NONE

type(EBSDNameListType),INTENT(IN)       :: enl
type(EBSDAngleType),pointer             :: angles
type(EBSDLargeAccumType),pointer        :: acc
type(EBSDMasterType),pointer            :: master
character(fnlen),INTENT(IN)             :: progname


! all geometrical parameters and filenames
real(kind=dbl)                          :: prefactor, qz(3)

! allocatable arrays
real(kind=sgl),allocatable              :: EBSDpattern(:,:), binned(:,:)        ! array with EBSD patterns
real(kind=sgl),allocatable              :: z(:,:)               ! used to store the computed patterns before writing to disk

! quaternion variables
real(kind=dbl)                          :: qq(4), qq1(4), qq2(4), qq3(4)

! various items
integer(kind=irg)                       :: i, j, iang, jang, k, io_int(6), etotal          ! various counters
integer(kind=irg)                       :: istat                ! status for allocate operations
integer(kind=irg)                       :: nix, niy, binx, biny,num_el       ! various parameters
integer(kind=irg)                       :: NUMTHREADS, TID   ! number of allocated threads, thread ID
integer(kind=irg)                       :: istart, istop, ninbatch, nbatches, nremainder, ibatch, nthreads, maskradius

real(kind=sgl)                          :: bindx, sig, ma, mi
real(kind=sgl),parameter                :: dtor = 0.0174533  ! convert from degrees to radians
real(kind=dbl),parameter                :: nAmpere = 6.241D+18   ! Coulomb per second
integer(kind=irg),parameter             :: storemax = 20        ! number of EBSD patterns stored in one output block
integer(kind=irg)                       :: Emin, Emax      ! various parameters
real(kind=sgl)                          :: dc(3), scl           ! direction cosine array
real(kind=sgl)                          :: sx, dx, dxm, dy, dym, rhos, x         ! various parameters
real(kind=sgl)                          :: ixy(2)

real(kind=sgl),allocatable              :: mask(:,:), lx(:), ly(:)
character(len=1),allocatable            :: batchpatterns(:,:,:), bpat(:,:)
integer(kind=irg),allocatable           :: acc_array(:,:)
real(kind=sgl),allocatable              :: master_array(:,:), wf(:) 
character(len=3)                        :: outputformat

! parameter for random number generator
integer, parameter                      :: K4B=selected_int_kind(9)      ! used by ran function in math.f90
integer(K4B)                            :: idum

!====================================
! what is the output format?  GUI or BIN ?
outputformat = enl%outputformat

! define some energy-related parameters derived from MC input parameters
!====================================
etotal = enl%num_el 
sig = enl%MCsig

! get the indices of the minimum and maximum energy
Emin = nint((enl%energymin - enl%Ehistmin)/enl%Ebinsize) +1
if (Emin.lt.1)  Emin=1
if (Emin.gt.enl%numEbins)  Emin=enl%numEbins

Emax = nint((enl%energymax - enl%Ehistmin)/enl%Ebinsize) +1
if (Emax.lt.1)  Emax=1
if (Emax.gt.enl%numEbins)  Emax=enl%numEbins

num_el = nint(sum(acc%accum_e_detector))

!====================================

!====================================
! init a bunch of parameters
!====================================
! binned pattern array
  binx = enl%numsx/enl%binning
  biny = enl%numsy/enl%binning
  bindx = 1.0/float(enl%binning)**2


! determine the scale factor for the Lambert interpolation; the square has
! an edge length of 2 x sqrt(pi/2)
  scl = float(enl%npx) / LPs%sPio2

! intensity prefactor
  prefactor = 0.25D0 * nAmpere * enl%beamcurrent * enl%dwelltime * 1.0D-15/ dble(num_el)
!====================================

!====================================
! ------ and open the output file for IDL visualization (only thread 0 can write to this file)
!====================================
open(unit=dataunit,file=trim(enl%datafile),status='unknown',form='unformatted',action='write')
! we need to write the image dimensions, and also how many of those there are...
if (enl%binning.eq.1) then 
  write (dataunit) enl%numsx, enl%numsy, enl%numangles
else 
  write (dataunit) binx, biny, enl%numangles
end if

!====================================
! ------ start the actual image computation loop
!====================================

!====================================
! to speed things up, we'll split the computation into batches of 1,000 patterns per thread; once those 
! are computed, we leave the OpenMP part to write them to a file (will be replaced with HDF5 output
! at a later stage)
!====================================


! and allocate space to store each batch
if (outputformat.eq.'gui') then
  nthreads = 1
  ninbatch = 1000
  nbatches = 0
  nremainder = mod(enl%numangles,ninbatch)
else
  nthreads = enl%nthreads
  ninbatch = 1000
  nbatches = enl%numangles/(ninbatch*nthreads)
  nremainder = mod(enl%numangles,ninbatch*nthreads)
  allocate(batchpatterns(enl%numsx,enl%numsy,ninbatch*nthreads),stat=istat)
! here we also create a mask if necessary
  allocate(mask(binx,biny),stat=istat)
  mask = 1.0
  if (enl%maskpattern.eq.'y') then
! create the circular mask in a potentially rectangular array
    maskradius = (minval( (/ binx, biny /) ) / 2 )**2
    allocate(lx(binx), ly(biny), stat=istat)
    lx = (/ (float(i),i=1,binx) /) - float(binx/2)
    ly = (/ (float(i),i=1,biny) /) - float(biny/2)
    do i=1,binx
      do j=1,biny
        if ((lx(i)**2+ly(j)**2).gt.maskradius) mask(i,j) = 0.0
      end do
    end do
    deallocate(lx, ly)
  end if
end if

! for dictionary computations, the patterns are usually rather small, so perhaps the explicit 
! energy sums can be replaced by an averaged approximate approach, in which all the energy bins
! are added together from the start, and all the master patterns are totaled as well...
if (enl%energyaverage.eq.1) then
  allocate(acc_array(enl%numsx,enl%numsy))
  acc_array = sum(acc%accum_e_detector,1)
  allocate(wf(enl%numEbins))
  wf = sum(sum(acc%accum_e_detector,2),2)
  wf = wf/sum(wf)

! this is a straightforward sum; we should probably do a weighted sum instead
  allocate(master_array(-enl%npx:enl%npx,-enl%npy:enl%npy))
  do k=Emin,Emax
    master%sr(-enl%npx:enl%npx,-enl%npy:enl%npy,k) = master%sr(-enl%npx:enl%npx,-enl%npy:enl%npy,k) * wf(k)
  end do
  master_array = sum(master%sr,3)
end if


! set the number of OpenMP threads and allocate the corresponding number of random number streams
io_int(1) = nthreads
call WriteValue(' Attempting to set number of threads to ',io_int,1,"(I4)")
call OMP_SET_NUM_THREADS(nthreads)

io_int(1) = ninbatch
io_int(2) = nbatches
io_int(3) = nremainder
call WriteValue('  OpenMP loop variables : ',io_int,3,"(I6,I6,I6)")

do ibatch=1,nbatches+1

  istart = 1
  if (ibatch.eq.nbatches+1) then ! take care of the remainder patterns
    istop = nremainder
  else
    istop = ninbatch*nthreads
  end if

! use OpenMP to run on multiple cores ... 
!$OMP PARALLEL default(shared)  PRIVATE(TID,iang,jang,i,j,qq,qq1,qq2,qq3,dc,ixy,istat,nix,niy,dx,dy,dxm,dym,EBSDpattern) &
!$OMP& PRIVATE(k,binned,idum,bpat)

! allocate the array that will hold the computed pattern
  allocate(EBSDpattern(enl%numsx,enl%numsy),stat=istat)
  allocate(binned(binx,biny),stat=istat)
  allocate(bpat(binx,biny),stat=istat)
! if (istat.ne.0) then ...

  NUMTHREADS = OMP_GET_NUM_THREADS()
  TID = OMP_GET_THREAD_NUM()

! initialize the random number generator for the Poison noise
  idum = -1-TID               

!$OMP DO SCHEDULE(STATIC,1000)  
  do iang=istart,istop
! convert the direction cosines to quaternions, include the 
! sample quaternion orientation, and then back to direction cosines...
! then convert these individually to the correct EBSD pattern location
        jang = (ibatch-1)*ninbatch*nthreads + iang
        EBSDpattern = 0.0
        binned = 0.0
        bpat = ' '

        do i=1,enl%numsx
            do j=1,enl%numsy
!  do the active coordinate transformation for this euler angle
              dc = sngl(quat_Lp(conjg(angles%quatang(1:4,jang)),  (/ master%rgx(i,j),master%rgy(i,j),master%rgz(i,j) /) )) 
! make sure the third one is positive; if not, switch all 
              dc = dc/sqrt(sum(dc**2))
              if (dc(3).lt.0.0) dc = -dc
! convert these direction cosines to coordinates in the Rosca-Lambert projection
              ixy = scl * LambertSphereToSquare( dc, istat )
! four-point interpolation (bi-quadratic)
              nix = int(enl%npx+ixy(1))-enl%npx
              niy = int(enl%npy+ixy(2))-enl%npy
              dx = ixy(1)-nix
              dy = ixy(2)-niy
              dxm = 1.0-dx
              dym = 1.0-dy
 ! interpolate the intensity 
              if (enl%energyaverage.eq.1) then
                EBSDpattern(i,j) = EBSDpattern(i,j) + acc_array(i,j) * ( master_array(nix,niy) * dxm * dym + &
                                           master_array(nix+1,niy) * dx * dym + master_array(nix,niy+1) * dxm * dy + &
                                           master_array(nix+1,niy+1) * dx * dy )
              else
                do k=Emin,Emax 
                  EBSDpattern(i,j) = EBSDpattern(i,j) + acc%accum_e_detector(k,i,j) * ( master%sr(nix,niy,k) * dxm * dym + &
                                             master%sr(nix+1,niy,k) * dx * dym + master%sr(nix,niy+1,k) * dxm * dy + &
                                             master%sr(nix+1,niy+1,k) * dx * dy )
                end do
              end if
          end do
       end do

        EBSDpattern = prefactor * EBSDpattern

! add sampling noise (Poisson noise in this case, so multiplicative, sort of)
        do i=1,enl%numsx
          do j=1,enl%numsy
              EBSDpattern(i,j) = POIDEV(EBSDpattern(i,j),idum)
          end do
        end do      

! if this is a GUI-computation, then we can directly store the current pattern in the output file;
! otherwise, we process it and turn it into a byte array of the right binning level.
        if (outputformat.eq.'gui') then 
          write (dataunit) EBSDpattern
        else

! we may need to deal with the energy sensitivity of the scintillator as well...


! all the following things should really be done in the IDL visualization program:
!
! - point spread function of camera
! - binning
! - brightness/contrast scaling
! - storage in image files of large HDF5 ? file
!  
! that means that at this point, we really only need to store all the patterns in a single
! file, at full resolution, as observed at the scintillator stage.
       
        if (enl%scalingmode.ne.'not') then

! apply any camera binning
         if (enl%binning.ne.1) then 
          do i=1,enl%numsx,enl%binning
            do j=1,enl%numsy,enl%binning
                binned(i/enl%binning+1,j/enl%binning+1) = sum(EBSDpattern(i:i+enl%binning-1,j:j+enl%binning-1))
            end do
          end do  
! and divide by binning^2
          binned = binned * bindx
         else
          binned = EBSDpattern
         end if
        
! and finally, before saving the patterns, apply the contrast function (linear or gamma correction)
! for the linear case, we do not need to do anything here...
         if (enl%scalingmode.eq.'gam') then
           binned = binned**enl%gammavalue
         end if

        end if ! scaling mode .ne. 'not'

! write either the EBSDpattern array or the binned array to the batchpatterns array (convert to bytes first)
        if (enl%scalingmode.eq.'not') then
         ma = maxval(EBSDpattern)
         mi = minval(EBSDpattern)
         EBSDpattern = mask * ((EBSDpattern - mi)/ (ma-mi))
         bpat = char(nint(255.0*EBSDpattern))
        else
         ma = maxval(binned)
         mi = minval(binned)
         binned = mask * ((binned - mi)/ (ma-mi))
         bpat = char(nint(255.0*binned))
        end if
       
        if (enl%binning.ne.1) then 
         batchpatterns(1:enl%numsx/enl%binning+1,1:enl%numsy/enl%binning+1, iang) = bpat
       else
         batchpatterns(1:enl%numsx,1:enl%numsy, iang) = bpat
       end if

      end if
  end do
!$OMP END DO

!$OMP END PARALLEL

! here we write all the entries in the batchpatterns array to a file, one at a time for now...
 if (outputformat.eq.'bin') then
  if (ibatch.le.nbatches) then 
    do iang=1,ninbatch*enl%nthreads 
      bpat = batchpatterns(1:binx,1:biny,iang)
      write (dataunit) bpat
    end do
  else
    do iang=1,nremainder
      bpat = batchpatterns(1:binx,1:biny,iang)
      write (dataunit) bpat
    end do
  end if
 end if
  write (*,*) 'completed cycle ',ibatch,' of ',nbatches+1
end do


close(unit=dataunit,status='keep')

i=6
call timestamp(i)

end subroutine ComputeEBSDPatterns

