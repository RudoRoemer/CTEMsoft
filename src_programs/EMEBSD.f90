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
! EMsoft:EMEBSD.f90
!--------------------------------------------------------------------------
!
! PROGRAM: EMEBSD
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief EMEBSD computes energy-weighted EBSD patterns
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
!> @date  04/02/15  MDG 5.0 changed program input & output to HDF format
! ###################################################################

program EMEBSD

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
        subroutine ComputeEBSDPatterns(enl, angles, acc, master, progname, nmldeffile)
        
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
        
        type(EBSDNameListType),INTENT(INOUT)    :: enl
        type(EBSDAngleType),pointer             :: angles
        type(EBSDLargeAccumType),pointer        :: acc
        type(EBSDMasterType),pointer            :: master
        character(fnlen),INTENT(IN)             :: progname
        character(fnlen),INTENT(IN)             :: nmldeffile
        end subroutine ComputeEBSDPatterns
end interface


nmldeffile = 'EMEBSD.nml'
progname = 'EMEBSD.f90'
progdesc = 'Dynamical EBSD patterns, using precomputed MC and master Lambert projections'

! print some information
call EMsoft(progname, progdesc)

! deal with the command line arguments, if any
call Interpret_Program_Arguments(nmldeffile,1,(/ 22 /), progname)

! deal with the namelist stuff
call GetEBSDNameList(nmldeffile,enl)

! this program needs a lot of data, and it also should be integrated 
! with Dream.3D, so we need to make sure that all data is loaded outside
! of the main computational routine, and passed in as pointers/arguments
! either by the fortran program or by Dream.3D calls.  

! 1. read the angle array from file
verbose = .TRUE.
allocate(angles)
call EBSDreadangles(enl, angles, verbose=.TRUE.)

! 2. read the Monte Carlo data file (including HDF format)
allocate(acc)
call EBSDreadMCfile(enl, acc, verbose=.TRUE.)

! 3. read EBSD master pattern file (including HDF format)
allocate(master)
call EBSDreadMasterfile(enl, master, verbose=.TRUE.)

! 3.1 twin the master pattern with equal weight (FZ == pg 622)
!call TwinCubicMasterPattern(enl,master)

! 4. generate detector arrays
allocate(master%rgx(enl%numsx,enl%numsy), master%rgy(enl%numsx,enl%numsy), master%rgz(enl%numsx,enl%numsy), stat=istat)
allocate(acc%accum_e_detector(enl%numEbins,enl%numsx,enl%numsy), stat=istat)
call EBSDGenerateDetector(enl, acc, master, verbose)
deallocate(acc%accum_e)

! call TwinCubicMasterPattern(enl, master)

! perform the zone axis computations for the knl input parameters
call ComputeEBSDpatterns(enl, angles, acc, master, progname, nmldeffile)

deallocate(master, acc, angles)

end program EMEBSD

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
!> @param nmldeffile name of nml file
!
!> @date 11/29/01  MDG 1.0 original
!> @date 04/08/13  MDG 2.0 rewrite
!> @date 05/14/13  MDG 2.1 replaced IO by namelist file
!> @date 08/01/13  MDG 3.0 complete rewrite, eliminated old Lambert projection
!> @date 09/25/13  MDG 3.1 replaced k-vector code by kvectors module
!> @date 02/26/14  MDG 4.0 new version
!> @date 03/26/14  MDG 4.1 adapted to new input and out file formats
!> @date 05/22/14  MDG 4.2 slight modification of angle input file; update for new EMEBSDMaster file format
!> @date 06/24/14  MDG 5.0 removal of global variables; removal of namelist stuff; 
!> @date 03/09/15  MDG 5.1 added OpenMP functionality for final loop
!> @date 03/10/15  MDG 5.2 added 'bin' and 'gui' outputformat; added mask support for 'bin' outputformat
!> @date 03/14/15  MDG 5.3 attempt at speeding up the program by performing approximate energy sums (energyaverage = 1)
!> @date 03/20/15  MDG 5.4 corrected out-of-bounds error in EBSDpattern array
!> @date 04/07/15  MDG 5.5 added HDF-formatted output
!> @date 05/08/15  MDG 5.6 added support for hexagonal/trigonal master pattern interpolation
!--------------------------------------------------------------------------
subroutine ComputeEBSDPatterns(enl, angles, acc, master, progname, nmldeffile)

use local
use typedefs
use NameListTypedefs
use NameListHDFwriters
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
use HDF5
use HDFsupport
use ISO_C_BINDING
use omp_lib

IMPLICIT NONE

type(EBSDNameListType),INTENT(INOUT)    :: enl
type(EBSDAngleType),pointer             :: angles
type(EBSDLargeAccumType),pointer        :: acc
type(EBSDMasterType),pointer            :: master
character(fnlen),INTENT(IN)             :: progname
character(fnlen),INTENT(IN)             :: nmldeffile


! all geometrical parameters and filenames
real(kind=dbl)                          :: prefactor, qz(3)

! allocatable arrays
real(kind=sgl),allocatable              :: EBSDpattern(:,:), binned(:,:)        ! array with EBSD patterns
real(kind=sgl),allocatable              :: z(:,:)               ! used to store the computed patterns before writing to disk

! quaternion variables
real(kind=dbl)                          :: qq(4), qq1(4), qq2(4), qq3(4)

! various items
integer(kind=irg)                       :: i, j, iang, jang, k, io_int(6), etotal, hdferr          ! various counters
integer(kind=irg)                       :: istat                ! status for allocate operations
integer(kind=irg)                       :: nix, niy, binx, biny,num_el, nixp, niyp       ! various parameters
integer(kind=irg)                       :: NUMTHREADS, TID   ! number of allocated threads, thread ID
integer(kind=irg)                       :: istart, istop, ninbatch, nbatches, nremainder, ibatch, nthreads, maskradius

real(kind=sgl)                          :: bindx, sig, ma, mi
real(kind=sgl),parameter                :: dtor = 0.0174533  ! convert from degrees to radians
real(kind=dbl),parameter                :: nAmpere = 6.241D+18   ! Coulomb per second
integer(kind=irg),parameter             :: storemax = 20        ! number of EBSD patterns stored in one output block
integer(kind=irg)                       :: Emin, Emax      ! various parameters
real(kind=dbl)                          :: dc(3), scl           ! direction cosine array
real(kind=dbl)                          :: sx, dx, dxm, dy, dym, rhos, x         ! various parameters
real(kind=dbl)                          :: ixy(2), tmp

real(kind=sgl),allocatable              :: mask(:,:), lx(:), ly(:)
character(len=1),allocatable            :: batchpatterns(:,:,:), bpat(:,:)
integer(kind=irg),allocatable           :: acc_array(:,:)
real(kind=sgl),allocatable              :: master_array(:,:), wf(:) 
character(len=3)                        :: outputformat

! parameter for random number generator
integer, parameter                      :: K4B=selected_int_kind(9)      ! used by ran function in math.f90
integer(K4B)                            :: idum

type(HDFobjectStackType),pointer        :: HDF_head

integer(HSIZE_T), dimension(1:3)        :: hdims, offset 
integer(HSIZE_T)                        :: dim0, dim1, dim2
character(fnlen,kind=c_char)            :: line2(1)
character(fnlen)                        :: groupname, dataset
character(11)                           :: dstr
character(15)                           :: tstrb
character(15)                           :: tstre
character(fnlen)                        :: datafile
logical                                 :: overwrite = .TRUE., insert = .TRUE.

nullify(HDF_head)

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


! intensity prefactor
  prefactor = 0.25D0 * nAmpere * enl%beamcurrent * enl%dwelltime * 1.0D-15/ dble(num_el)
!====================================

!====================================
! ------ and open the output file for IDL visualization (only thread 0 can write to this file)
!====================================
! we need to write the image dimensions, and also how many of those there are...

! Initialize FORTRAN interface.
!
CALL h5open_f(hdferr)

call timestamp(datestring=dstr, timestring=tstrb)
tstre = tstrb

! Create a new file using the default properties.
datafile = trim(EMdatapathname)//trim(enl%datafile)
hdferr =  HDF_createFile(datafile, HDF_head)

! write the EMheader to the file
call HDF_writeEMheader(HDF_head, dstr, tstrb, tstre, progname)

! create a namelist group to write all the namelist files into
groupname = "NMLfiles"
hdferr = HDF_createGroup(groupname, HDF_head)

! read the text file and write the array to the file
dataset = 'EMEBSDNML'
hdferr = HDF_writeDatasetTextFile(dataset, nmldeffile, HDF_head)

call HDF_pop(HDF_head)

! create a NMLparameters group to write all the namelist entries into
groupname = "NMLparameters"
hdferr = HDF_createGroup(groupname, HDF_head)

call HDFwriteEBSDNameList(HDF_head, enl)

! and leave this group
call HDF_pop(HDF_head)

! then the remainder of the data in a EMData group
groupname = 'EMData'
hdferr = HDF_createGroup(groupname, HDF_head)

! we need to write the image dimensions
dataset = 'binx'
hdferr = HDF_writeDatasetInteger(dataset, binx, HDF_head) 

dataset = 'biny'
hdferr = HDF_writeDatasetInteger(dataset, biny, HDF_head) 

dataset = 'enl%numangles'
hdferr = HDF_writeDatasetInteger(dataset, enl%numangles, HDF_head) 

! and we leave this group open for further data output ... 

!====================================
! ------ start the actual image computation loop
!====================================

!====================================
! to speed things up, we'll split the computation into batches of 1,024 patterns per thread; once those 
! are computed, we leave the OpenMP part to write them to a file 
!====================================


! and allocate space to store each batch
if (outputformat.eq.'gui') then
  nthreads = 1
  ninbatch = 1024
  nbatches = 0
  nremainder = mod(enl%numangles,ninbatch)
else
  nthreads = enl%nthreads
  ninbatch = 1024
  nbatches = enl%numangles/(ninbatch*nthreads)
  nremainder = mod(enl%numangles,ninbatch*nthreads)
  allocate(batchpatterns(enl%numsx/enl%binning,enl%numsy/enl%binning,ninbatch*nthreads),stat=istat)
! this is also the size of the hyperslabs that we will write to the HDF output file

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

! determine the scale factor for the Lambert interpolation; the square has
! an edge length of 2 x sqrt(pi/2), and a different scale factor for the 
! hexagonal/trigonal case
  if (enl%sqorhe.eq.'square') then 
    scl = dble(enl%npx) / LPs%sPio2
  else
    scl = dble(enl%npx) / LPs%preg
  end if


! set the number of OpenMP threads and allocate the corresponding number of random number streams
io_int(1) = nthreads
call WriteValue(' Attempting to set number of threads to ',io_int,1,"(I4)")
call OMP_SET_NUM_THREADS(nthreads)


io_int(1) = enl%numangles 
io_int(2) = ninbatch
io_int(3) = nthreads
io_int(4) = nbatches
io_int(5) = nremainder
call WriteValue('  OpenMP loop variables : ',io_int,5,"(I10,' = ',I4,' * ',I2,' * ',I4,' + ',I6)")

do ibatch=1,nbatches+1

  istart = 1
  if (ibatch.eq.nbatches+1) then ! take care of the remainder patterns
    istop = nremainder
  else
    istop = ninbatch*nthreads
  end if

! use OpenMP to run on multiple cores ... 
!$OMP PARALLEL default(shared)  PRIVATE(TID,iang,jang,i,j,qq,qq1,qq2,qq3,dc,ixy,istat,nix,niy,dx,dy,dxm,dym,EBSDpattern) &
!$OMP& PRIVATE(k,binned,idum,bpat,nixp,niyp)

! allocate the array that will hold the computed pattern
  allocate(EBSDpattern(enl%numsx,enl%numsy),stat=istat)
  allocate(binned(binx,biny),stat=istat)
  allocate(bpat(binx,biny),stat=istat)
! if (istat.ne.0) then ...

  NUMTHREADS = OMP_GET_NUM_THREADS()
  TID = OMP_GET_THREAD_NUM()

! initialize the random number generator for the Poison noise
  idum = -1-TID               

!$OMP DO SCHEDULE(STATIC,1024)  
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
! do the active coordinate transformation for this euler angle
              dc = quat_Lp(angles%quatang(1:4,jang),  (/ master%rgx(i,j),master%rgy(i,j),master%rgz(i,j) /) )
! make sure the third one is positive; if not, switch all 
              dc = dc/dsqrt(sum(dc*dc))
              if (dc(3).lt.0.D0) dc = -dc
! convert these direction cosines to coordinates in the Rosca-Lambert projection
              if (enl%sqorhe.eq.'square') then
                ixy = scl * LambertSphereToSquare( dc, istat )
              else
                ixy = scl * LambertSphereToHex( dc, istat )
                tmp = ixy(1)
                ixy(1) = ixy(2)
                ixy(2) = tmp
              end if

              if (istat.eq.0) then 
! four-point interpolation (bi-quadratic)
                nix = int(enl%npx+ixy(1))-enl%npx
                niy = int(enl%npy+ixy(2))-enl%npy
                nixp = nix+1
                niyp = niy+1
                if (nixp.gt.enl%npx) nixp = nix
                if (niyp.gt.enl%npy) niyp = niy
                if (nix.lt.-enl%npx) nix = nixp
                if (niy.lt.-enl%npy) niy = niyp
                dx = ixy(1)-nix
                dy = ixy(2)-niy
                dxm = 1.0-dx
                dym = 1.0-dy
 ! interpolate the intensity 
                if (enl%energyaverage.eq.1) then
                  EBSDpattern(i,j) = EBSDpattern(i,j) + acc_array(i,j) * ( master_array(nix,niy) * dxm * dym + &
                                             master_array(nixp,niy) * dx * dym + master_array(nix,niyp) * dxm * dy + &
                                             master_array(nixp,niyp) * dx * dy )
                else
                  do k=Emin,Emax 
                    EBSDpattern(i,j) = EBSDpattern(i,j) + acc%accum_e_detector(k,i,j) * ( master%sr(nix,niy,k) * dxm * dym + &
                                               master%sr(nixp,niy,k) * dx * dym + master%sr(nix,niyp,k) * dxm * dy + &
                                               master%sr(nixp,niyp,k) * dx * dy )
                  end do
                end if
              end if
          end do
       end do

        EBSDpattern = prefactor * EBSDpattern

! add sampling noise (Poisson noise in this case, so multiplicative, sort of)
        if (outputformat.eq.'gui') then 
         do i=1,enl%numsx
          do j=1,enl%numsy
              EBSDpattern(i,j) = POIDEV(EBSDpattern(i,j),idum)
          end do
         end do
! if this is a GUI-computation, then we can directly store the current pattern in the output file;
! otherwise, we process it and turn it into a byte array of the right binning level.
          dataset = 'EBSDpatterns'
          offset = (/ 0, 0, ibatch-1 /)
          hdims = (/ enl%numsx, enl%numsy, enl%numangles /)
          dim0 = enl%numsx
          dim1 = enl%numsy
          dim2 = 1
          if (ibatch.eq.1) then
            hdferr = HDF_writeHyperslabFloatArray3D(dataset, EBSDpattern, hdims, offset, dim0, dim1, dim2, &
                                          HDF_head)
          else
            hdferr = HDF_writeHyperslabFloatArray3D(dataset, EBSDpattern, hdims, offset, dim0, dim1, dim2, &
                                          HDF_head, insert)
          end if
        else ! if output mode not equal to gui
       

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


! write the binned array to the batchpatterns array (convert to bytes first)

         ma = maxval(binned)
         mi = minval(binned)
         binned = mask * ((binned - mi)/ (ma-mi))
         bpat = char(nint(255.0*binned))

         batchpatterns(1:enl%numsx/enl%binning,1:enl%numsy/enl%binning, iang) = bpat

      end if
  end do ! end of iang loop
!$OMP END DO

!$OMP END PARALLEL

! here we write all the entries in the batchpatterns array to the HDF file as a hyperslab
dataset = 'EBSDpatterns'



 if (outputformat.eq.'bin') then
  if (ibatch.le.nbatches) then 
   offset = (/ 0, 0, (ibatch-1)*ninbatch*enl%nthreads /)
   hdims = (/ binx, biny, enl%numangles /)
   dim0 = binx
   dim1 = biny
   dim2 = ninbatch*enl%nthreads
! write (*,*) ibatch,offset,hdims,dim0,dim1,dim2
   if (ibatch.eq.1) then
     hdferr = HDF_writeHyperslabCharArray3D(dataset, batchpatterns, hdims, offset, dim0, dim1, dim2, &
                                          HDF_head)
   else
     hdferr = HDF_writeHyperslabCharArray3D(dataset, batchpatterns, hdims, offset, dim0, dim1, dim2, &
                                          HDF_head, insert)
   end if
  else
   offset = (/ 0, 0, (ibatch-1)*ninbatch*enl%nthreads /)
   hdims = (/ binx, biny, enl%numangles /)
   dim0 = binx
   dim1 = biny
   dim2 = nremainder
! write (*,*) ibatch,offset,hdims,dim0,dim1,dim2
   if (ibatch.eq.1) then
     hdferr = HDF_writeHyperslabCharArray3D(dataset, batchpatterns(1:binx,1:biny,1:nremainder), hdims, offset, dim0, dim1, dim2, &
                                          HDF_head)
   else
     hdferr = HDF_writeHyperslabCharArray3D(dataset, batchpatterns(1:binx,1:biny,1:nremainder), hdims, offset, dim0, dim1, dim2, &
                                          HDF_head, insert)
   end if
  end if
 end if
  write (*,*) 'completed cycle ',ibatch,' of ',nbatches+1
end do

call HDF_pop(HDF_head)

! and update the end time
call timestamp(datestring=dstr, timestring=tstre)
groupname = "EMheader"
hdferr = HDF_openGroup(groupname, HDF_head)

! stop time /EMheader/StopTime 'character'
dataset = 'StopTime'
line2(1) = tstre
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, overwrite)

! close the datafile
call HDF_pop(HDF_head,.TRUE.)

! close the Fortran interface
call h5close_f(hdferr)


end subroutine ComputeEBSDPatterns

