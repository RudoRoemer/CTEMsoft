! ###################################################################
! Copyright (c) 2013-2014, Marc De Graef/Saransh Singh/Carnegie Mellon University
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
! EMsoft:EMECP.f90
!--------------------------------------------------------------------------
!
! PROGRAM: EMECP
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief Electron channeling patterns from master pattern
!
!> @date 08/26/14 SS 1.0 f90
!> @date 13/10/15 SS 2.0 added detector model+new GetVectorCone routine+OpenMP
!--------------------------------------------------------------------------

program EMECP

use local
use typedefs
use NameListTypedefs
use NameListHandlers
use files
use io

IMPLICIT NONE

character(fnlen)                        :: nmldeffile, progname, progdesc
type(ECPNameListType)                   :: ecpnl

nmldeffile = 'EMECP.nml'
progname = 'EMECP.f90'
progdesc = 'Electron channeling patterns from master pattern'

! print some information
call EMsoft(progname, progdesc)

! deal with the command line arguments, if any
call Interpret_Program_Arguments(nmldeffile,2,(/ 0, 40 /), progname)

! deal with the namelist stuff
call GetECPNameList(nmldeffile,ecpnl)

! perform the zone axis computations
call ECpattern(ecpnl, progname, nmldeffile)

end program EMECP

!-----------------------------------------------------------------------------------
!
! SUBROUTINE: ECpattern
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief Electron channeling patterns from master pattern
!
!> @date 08/27/14 SS 1.0 f90
!> @date 13/10/15 SS 2.0 added detector model+new GetVectorCone routine+OpenMP+hdf5
!-------------------------------------------------------------------------------------
subroutine ECpattern(ecpnl, progname, nmldeffile)

use local
use typedefs
use NameListTypedefs
use crystal
use constants
use symmetry
use Lambert
use initializers
use constants
use error
use io
use files
use ECPmod
use rotations
use quaternions
use HDF5
use HDFsupport
use NameListHDFwriters
use ISO_C_BINDING
use omp_lib

IMPLICIT NONE

type(ECPNameListType),INTENT(INOUT)     :: ecpnl
character(fnlen),INTENT(IN)             :: progname
character(fnlen),INTENT(IN)             :: nmldeffile

type(ECPLargeAccumType),pointer         :: acc
type(ECPMasterType),pointer             :: master
type(IncidentListECP),pointer           :: khead, ktmp
type(ECPAngleType),pointer              :: angles
real(kind=dbl),allocatable              :: klist(:,:)

integer(kind=irg)                       :: npx,npy,numset,istat,val
integer(kind=irg),allocatable           :: ATOM_type(:)
real(kind=dbl)                          :: EkeV
real(kind=sgl)                          :: dmin, FN(3)
real(kind=sgl),allocatable              :: mask(:,:), lx(:), ly(:)
integer(kind=irg)                       :: maskradius, io_int(1), hdferr
logical                                 :: verbose
real(kind=sgl),allocatable              :: master_arrayNH(:,:), master_arraySH(:,:)
integer(kind=irg)                       :: numk, nix, niy, nixp, niyp, i, j, ierr, &
                                           ipx, ipy, iang, idir
real(kind=dbl)                          :: scl, x, dx, dy, dxm, dym
real(kind=dbl)                          :: dc(3), ixy(2)
integer(kind=irg),allocatable           :: kij(:,:)
real(kind=sgl),allocatable              :: ECPpattern(:,:), ECPpatternintd(:,:)
real(kind=sgl)                          :: time_start, time_end, ma, mi
character(len=1),allocatable            :: bpat(:,:)
integer(kind=irg)                       :: TID, nthreads

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

!=================================================================
! read Monte Carlo output file and extract necessary parameters
! first, we need to load the data from the output of EMMCOpenCL
! used in the bse1 mode
!=================================================================

call Message('opening '//trim(ecpnl%energyfile), frm = "(A)" )
call  ECPreadMCfile(ecpnl, acc, verbose=.TRUE.)
call Message(' -> completed reading '//trim(ecpnl%energyfile), frm = "(A)")

!=================================================================
! read Master pattern output file and extract necessary parameters
! first, we need to load the data from the ECP master program
!=================================================================

call Message('opening '//trim(ecpnl%masterfile), frm = "(A)" )
call ECPreadMasterfile(ecpnl, master, verbose=.TRUE.)
call Message(' -> completed reading '//trim(ecpnl%masterfile), frm = "(A)")

!=================================================================
! reading the angle file as euler angles or quaternions
!=================================================================

call Message('opening '//trim(ecpnl%anglefile), frm = "(A)" )
call ECPreadangles(ecpnl, angles, verbose=.TRUE.)
call Message(' -> completed reading '//trim(ecpnl%anglefile), frm = "(A)")

!=================================================================
! completed reading the file; generating list of incident vectors
!=================================================================

numk = 0
call GetVectorsCone(ecpnl, khead, numk)
allocate(kij(2,numk),klist(3,numk),stat=istat)

!=============================================================================================
! ------ and open the output file for IDL visualization (only thread 0 can write to this file)
!=============================================================================================
! we need to write the image dimensions, and also how many of those there are...

! Initialize FORTRAN interface.
!
CALL h5open_f(hdferr)

call timestamp(datestring=dstr, timestring=tstrb)
tstre = tstrb
call CPU_TIME(time_start)

! Create a new file using the default properties.
datafile = trim(EMdatapathname)//trim(ecpnl%datafile)
hdferr =  HDF_createFile(datafile, HDF_head)

! write the EMheader to the file
call HDF_writeEMheader(HDF_head, dstr, tstrb, tstre, progname)

! create a namelist group to write all the namelist files into
groupname = "NMLfiles"
hdferr = HDF_createGroup(groupname, HDF_head)

! read the text file and write the array to the file
dataset = 'EMECPNML'
hdferr = HDF_writeDatasetTextFile(dataset, nmldeffile, HDF_head)

call HDF_pop(HDF_head)

! create a NMLparameters group to write all the namelist entries into
groupname = "NMLparameters"
hdferr = HDF_createGroup(groupname, HDF_head)

call HDFwriteECPNameList(HDF_head, ecpnl, .FALSE.)

! and leave this group
call HDF_pop(HDF_head)

! then the remainder of the data in a EMData group
groupname = 'EMData'
hdferr = HDF_createGroup(groupname, HDF_head)

! we need to write the image dimensions
dataset = 'npix'
hdferr = HDF_writeDatasetInteger(dataset, ecpnl%npix, HDF_head) 

dataset = 'numangle_dictionary'
hdferr = HDF_writeDatasetInteger(dataset, ecpnl%numangle_anglefile, HDF_head) 

! and we leave this group open for further data output ... 

ktmp => khead
! converting to array for OpenMP parallelization
do i = 1,numk
   klist(1:3,i) = ktmp%k(1:3)
   kij(1:2,i) = (/ktmp%i,ktmp%j/)
   ktmp => ktmp%next
end do

allocate(mask(ecpnl%npix, ecpnl%npix),stat=istat)
if (istat .ne. 0) then
   call FatalError('ECpattern','could not allocate mask array')
end if

mask = 1.0
if (ecpnl%maskpattern.eq.'y') then
! create the circular mask in a potentially rectangular array
  maskradius = (minval( (/ ecpnl%npix, ecpnl%npix /) ) / 2 )**2
  allocate(lx(ecpnl%npix), ly(ecpnl%npix), stat=istat)
  lx = (/ (float(i),i=1,ecpnl%npix) /) - float(ecpnl%npix/2)
  ly = (/ (float(i),i=1,ecpnl%npix) /) - float(ecpnl%npix/2)
  do i=1,ecpnl%npix
    do j=1,ecpnl%npix
      if ((lx(i)**2+ly(j)**2).gt.maskradius) mask(i,j) = 0.0
    end do
  end do
  deallocate(lx, ly)
end if

allocate(master_arrayNH(-ecpnl%npx:ecpnl%npx,-ecpnl%npy:ecpnl%npy))
allocate(master_arraySH(-ecpnl%npx:ecpnl%npx,-ecpnl%npy:ecpnl%npy))

master_arrayNH = master%mLPNH
master_arraySH = master%mLPSH

! determine the scale factor for the Lambert interpolation; the square has
! an edge length of 2 x sqrt(pi/2)
! This has been changed on 09/01/15 to accommodate the new Lambert module]

!scl = float(npx)/LPs%sPio2
scl = dble(ecpnl%npx)

allocate(ECPpattern(-ecpnl%npix:ecpnl%npix,-ecpnl%npix:ecpnl%npix),&
ECPpatternintd(-ecpnl%npix:ecpnl%npix,-ecpnl%npix:ecpnl%npix),stat=istat)
ECPpattern = 0.0

! set the number of OpenMP threads and allocate the corresponding number of random number streams
io_int(1) = ecpnl%nthreads
call WriteValue(' Attempting to set number of threads to ',io_int,1,"(I4)")
call OMP_SET_NUM_THREADS(ecpnl%nthreads)

! use OpenMP to run on multiple cores
!$OMP PARALLEL PRIVATE(TID,nthreads)

TID = OMP_GET_THREAD_NUM()
nthreads = OMP_GET_NUM_THREADS()
if (TID .eq. 0) write(*,*)'Number of threads = ',nthreads 

!$OMP DO SCHEDULE(DYNAMIC)
angleloop: do iang = 1,ecpnl%numangle_anglefile

    imageloop: do idir = 1,numk

! do the active coordinate transformation for this euler angle

        dc = klist(1:3,idir)
        dc = quat_LP(angles%quatang(1:4,iang),dc)
        dc = dc/dsqrt(sum(dc*dc))

! make sure the third one is positive; if not, switch all
! convert these direction cosines to coordinates in the Rosca-Lambert projection

        ixy = scl * LambertSphereToSquare( dc, istat )
        if (istat .ne. 0) call FatalError('ECpattern','Cannot convert to square Lambert projection')
        nix = int(ecpnl%npx+ixy(1))-ecpnl%npx
        niy = int(ecpnl%npy+ixy(2))-ecpnl%npy
        nixp = nix+1
        niyp = niy+1
        if (nixp.gt.ecpnl%npx) nixp = nix
        if (niyp.gt.ecpnl%npy) niyp = niy
        if (nix.lt.-ecpnl%npx) nix = nixp
        if (niy.lt.-ecpnl%npy) niy = niyp
        dx = ixy(1)-nix
        dy = ixy(2)-niy
        dxm = 1.0-dx
        dym = 1.0-dy

! interpolate the intensity
        ipx = kij(1,idir)
        ipy = kij(2,idir)
! including the detector model with some sample tilt
        
        if (dc(3).gt.0.0) then 
            ECPpattern(ipx,ipy) =   master_arrayNH(nix,niy) * dxm * dym + &
                             master_arrayNH(nixp,niy) * dx * dym + &
                             master_arrayNH(nix,niyp) * dxm * dy + &
                             master_arrayNH(nixp,niyp) * dx * dy
        else
            ECPpattern(ipx,ipy) =   master_arrayNH(nix,niy) * dxm * dym + &
                             master_arrayNH(nixp,niy) * dx * dym + &
                             master_arrayNH(nix,niyp) * dxm * dy + &
                             master_arrayNH(nixp,niyp) * dx * dy

        end if

       
    end do imageloop 

!$OMP CRITICAL
    if (ecpnl%outputformat .eq. 'bin') then
        allocate(bpat(ecpnl%npix,ecpnl%npix),stat=istat)
        if (istat .ne. 0) call FatalError('ECpatter','cannot allocate bpat array')
        ma = maxval(ECPpattern)
        mi = minval(ECPpattern)
        ECPpatternintd = mask * ((ECPpattern - mi)/ (ma-mi))
        bpat = char(nint(255.0*ECPpatternintd))

! write dictionary pattern to h5 file
        offset = (/ 0, 0, iang-1 /)
        hdims = (/ ecpnl%npix, ecpnl%npix, ecpnl%numangle_anglefile /)
        dim0 = ecpnl%npix
        dim1 = ecpnl%npix
        dim2 = 1
        hdferr = HDF_writeHyperslabCharArray3D(dataset, bpat, hdims, offset, dim0, dim1, dim2, &
                                          HDF_head, insert)
    end if
!$OMP END CRITICAL

    if (ecpnl%maskpattern.eq.'y')  ECPpattern = ECPpattern * mask

!$OMP CRITICAL
    if (ecpnl%outputformat .eq. 'gui') then
          offset = (/ 0, 0, iang-1 /)
          hdims = (/ ecpnl%npix, ecpnl%npix, ecpnl%numangle_anglefile /)
          dim0 = ecpnl%npix
          dim1 = ecpnl%npix
          dim2 = 1
          if (iang.eq.1) then
            hdferr = HDF_writeHyperslabFloatArray3D(dataset, ECPpattern, hdims, offset, dim0, dim1, dim2, &
                                          HDF_head)
          else
            hdferr = HDF_writeHyperslabFloatArray3D(dataset, ECPpattern, hdims, offset, dim0, dim1, dim2, &
                                          HDF_head, insert)
          end if
    end if
!$OMP END CRITICAL

end do angleloop 

!$OMP END DO
!$OMP END PARALLEL

call HDF_pop(HDF_head)

! and update the end time
call timestamp(datestring=dstr, timestring=tstre)
groupname = "EMheader"
hdferr = HDF_openGroup(groupname, HDF_head)

! stop time /EMheader/StopTime 'character'
dataset = 'StopTime'
line2(1) = tstre
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, overwrite)

call CPU_TIME(time_end)
dataset = 'Duration'
time_end = time_end - time_start
hdferr = HDF_writeDatasetFloat(dataset, time_end, HDF_head)

! close the datafile
call HDF_pop(HDF_head,.TRUE.)

! close the Fortran interface
call h5close_f(hdferr)


end subroutine ECpattern
