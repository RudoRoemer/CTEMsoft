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
! EMsoft:EBSDmod.f90
!--------------------------------------------------------------------------
!
! MODULE: EBSDmod
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief EMEBSD helper routines
!
!> @date  06/24/14  MDG 1.0 original, lifted from EMEBSD.f90 to simplify code
!> @date  09/01/15  MDG 1.1 modified EBSDMasterType definition to accommodate multiple Lambert maps
!> @date  09/15/15  SS  added accum_z to EBSDLargeAccumType
!--------------------------------------------------------------------------
module EBSDmod

use local
use typedefs

IMPLICIT NONE

type EBSDAngleType
        real(kind=sgl),allocatable      :: quatang(:,:)
end type EBSDAngleType

type EBSDLargeAccumType
        integer(kind=irg),allocatable   :: accum_e(:,:,:),accum_z(:,:,:,:)
        real(kind=irg),allocatable      :: accum_e_detector(:,:,:)
end type EBSDLargeAccumType

type EBSDMasterType
        real(kind=sgl),allocatable      :: mLPNH(:,:,:) , mLPSH(:,:,:)
        real(kind=sgl),allocatable      :: rgx(:,:), rgy(:,:), rgz(:,:)          ! auxiliary detector arrays needed for interpolation
end type EBSDMasterType
        
contains

!--------------------------------------------------------------------------
!
! SUBROUTINE:EBSDreadangles
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read angles from an angle file
!
!> @param enl EBSD name list structure
!> @param quatang array of unit quaternions (output)
!
!> @date 06/24/14  MDG 1.0 original
!--------------------------------------------------------------------------
subroutine EBSDreadangles(enl,angles,verbose)

use NameListTypedefs
use io
use files
use quaternions
use rotations

IMPLICIT NONE


type(EBSDNameListType),INTENT(INOUT)    :: enl
type(EBSDAngleType),pointer             :: angles
logical,INTENT(IN),OPTIONAL             :: verbose

integer(kind=irg)                       :: io_int(1), i
character(2)                            :: angletype
real(kind=sgl),allocatable              :: eulang(:,:)   ! euler angle array
real(kind=sgl)                          :: qax(4)        ! axis-angle rotation quaternion

real(kind=sgl),parameter                :: dtor = 0.0174533  ! convert from degrees to radians
integer(kind=irg)                       :: istat


!====================================
! get the angular information, either in Euler angles or in quaternions, from a file
!====================================
! open the angle file 
open(unit=dataunit,file=trim(EMsoft_toNativePath(enl%anglefile)),status='old',action='read')

! get the type of angle first [ 'eu' or 'qu' ]
read(dataunit,*) angletype
if (angletype.eq.'eu') then 
  enl%anglemode = 'euler'
else
  enl%anglemode = 'quats'
end if

! then the number of angles in the file
read(dataunit,*) enl%numangles

if (present(verbose)) then 
  io_int(1) = enl%numangles
  call WriteValue('Number of angle entries = ',io_int,1)
end if

if (enl%anglemode.eq.'euler') then
! allocate the euler angle array
  allocate(eulang(3,enl%numangles),stat=istat)
! if istat.ne.0 then do some error handling ... 
  do i=1,enl%numangles
    read(dataunit,*) eulang(1:3,i)
  end do
  close(unit=dataunit,status='keep')

  if (enl%eulerconvention.eq.'hkl') then
    if (present(verbose)) call Message('  -> converting Euler angles to TSL representation', frm = "(A/)")
    eulang(1,1:enl%numangles) = eulang(1,1:enl%numangles) + 90.0
  end if

! convert the euler angle triplets to quaternions
  allocate(angles%quatang(4,enl%numangles),stat=istat)
! if (istat.ne.0) then ...

  if (present(verbose)) call Message('  -> converting Euler angles to quaternions', frm = "(A/)")
  
  do i=1,enl%numangles
    angles%quatang(1:4,i) = eu2qu(eulang(1:3,i)*dtor)
  end do

else
! the input file has quaternions, not Euler triplets
  allocate(angles%quatang(4,enl%numangles),stat=istat)
  do i=1,enl%numangles
    read(dataunit,*) angles%quatang(1:4,i)
  end do
end if

close(unit=dataunit,status='keep')

!====================================
! Do we need to apply an additional axis-angle pair rotation to all the quaternions ?
!
if (enl%axisangle(4).ne.0.0) then
  enl%axisangle(4) = enl%axisangle(4) * dtor
  qax = ax2qu( enl%axisangle )
  do i=1,enl%numangles
    angles%quatang(1:4,i) = quat_mult(qax,angles%quatang(1:4,i))
  end do 
end if

write (*,*) 'completed reading Euler angles'

end subroutine EBSDreadangles

!--------------------------------------------------------------------------
!
! SUBROUTINE:EBSDreadMCfile
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read angles from an angle file
!
!> @param enl EBSD name list structure
!> @param acc energy structure
!
!> @date 06/24/14  MDG 1.0 original
!> @date 11/18/14  MDG 1.1 removed enl%MCnthreads from file read
!> @date 04/02/15  MDG 2.0 changed program input & output to HDF format
!> @date 04/29/15  MDG 2.1 add optional parameter efile
!> @date 09/15/15  SS  2.2 added accum_z reading 
!--------------------------------------------------------------------------
subroutine EBSDreadMCfile(enl,acc,efile,verbose)

use NameListTypedefs
use files
use io
use HDF5
use HDFsupport
use error

IMPLICIT NONE

type(EBSDNameListType),INTENT(INOUT)    :: enl
type(EBSDLargeAccumType),pointer        :: acc
character(fnlen),INTENT(IN),OPTIONAL    :: efile
logical,INTENT(IN),OPTIONAL             :: verbose

integer(kind=irg)                       :: istat, hdferr, nlines, nx
logical                                 :: stat, readonly
integer(HSIZE_T)                        :: dims3(3),dims4(4)
character(fnlen)                        :: groupname, dataset, energyfile 
character(fnlen),allocatable            :: stringarray(:)

integer(kind=irg),allocatable           :: acc_e(:,:,:),acc_z(:,:,:,:)

type(HDFobjectStackType),pointer        :: HDF_head


! is the efile parameter present? If so, use it as the filename, otherwise use the enl%energyfile parameter
if (PRESENT(efile)) then
  energyfile = efile
else
  energyfile = trim(Emdatapathname)//trim(enl%energyfile)
end if

allocate(acc)

! first, we need to check whether or not the input file is of the HDF5 format type; if
! it is, we read it accordingly, otherwise we use the old binary format.
!
call h5fis_hdf5_f(energyfile, stat, hdferr)

if (stat) then
! open the fortran HDF interface
  call h5open_f(hdferr)

  nullify(HDF_head)

! open the MC file using the default properties.
  readonly = .TRUE.
  hdferr =  HDF_openFile(energyfile, HDF_head, readonly)

! open the namelist group
  groupname = 'NMLparameters'
  hdferr = HDF_openGroup(groupname, HDF_head)

  groupname = 'MCCLNameList'
  hdferr = HDF_openGroup(groupname, HDF_head)

! read all the necessary variables from the namelist group
  dataset = 'xtalname'
  stringarray = HDF_readDatasetStringArray(dataset, nlines, HDF_head)
  enl%MCxtalname = trim(stringarray(1))
  deallocate(stringarray)

  dataset = 'mode'
  stringarray = HDF_readDatasetStringArray(dataset, nlines, HDF_head)
  enl%MCmode = trim(stringarray(1))
  deallocate(stringarray)
  if (enl%MCmode .ne. 'full') call FatalError('EBSDreadMCfile','This file is not in full mode. Please input correct HDF5 file')
  dataset = 'numsx'
  enl%nsx = HDF_readDatasetInteger(dataset, HDF_head)
  enl%nsx = (enl%nsx - 1)/2
  enl%nsy = enl%nsx

  dataset = 'EkeV'
  enl%EkeV = HDF_readDatasetDouble(dataset, HDF_head)

  dataset = 'Ehistmin'
  enl%Ehistmin = HDF_readDatasetDouble(dataset, HDF_head)

  dataset = 'Ebinsize'
  enl%Ebinsize = HDF_readDatasetDouble(dataset, HDF_head)

  dataset = 'depthmax'
  enl%depthmax = HDF_readDatasetDouble(dataset, HDF_head)

  dataset = 'depthstep'
  enl%depthstep = HDF_readDatasetDouble(dataset, HDF_head)

  dataset = 'sig'
  enl%MCsig = HDF_readDatasetDouble(dataset, HDF_head)

  dataset = 'omega'
  enl%MComega = HDF_readDatasetDouble(dataset, HDF_head)

! close the name list group
  call HDF_pop(HDF_head)
  call HDF_pop(HDF_head)

! read from the EMheader
  groupname = 'EMheader'
  hdferr = HDF_openGroup(groupname, HDF_head)

  dataset = 'ProgramName'
  stringarray = HDF_readDatasetStringArray(dataset, nlines, HDF_head)
  enl%MCprogname = trim(stringarray(1))
  deallocate(stringarray)

  dataset = 'Version'
  stringarray = HDF_readDatasetStringArray(dataset, nlines, HDF_head)
  enl%MCscversion = trim(stringarray(1))
  deallocate(stringarray)

  call HDF_pop(HDF_head)

! open the Data group
  groupname = 'EMData'
  hdferr = HDF_openGroup(groupname, HDF_head)

! read data items 
  dataset = 'numEbins'
  enl%numEbins = HDF_readDatasetInteger(dataset, HDF_head)

  dataset = 'numzbins'
  enl%numzbins = HDF_readDatasetInteger(dataset, HDF_head)

  dataset = 'accum_e'
  acc_e = HDF_readDatasetIntegerArray3D(dataset, dims3, HDF_head)
  enl%num_el = sum(acc_e)
  nx = (dims3(2)-1)/2
  allocate(acc%accum_e(1:dims3(1),-nx:nx,-nx:nx))
  acc%accum_e = acc_e
  deallocate(acc_e)

  dataset = 'accum_z'
  acc_z = HDF_readDatasetIntegerArray4D(dataset, dims4, HDF_head)
  allocate(acc%accum_z(1:dims4(1),1:dims4(2),1:dims4(3),1:dims4(4)))
  acc%accum_z = acc_z
  deallocate(acc_z)

! and close everything
  call HDF_pop(HDF_head,.TRUE.)

! close the fortran HDF interface
  call h5close_f(hdferr)

else
!====================================
! ------ read old format Monte Carlo data file
!====================================
  if (present(verbose)) call Message('opening '//trim(enl%energyfile), frm = "(A)")

  open(dataunit,file=trim(EMsoft_toNativePath(energyfile)),status='unknown',form='unformatted')

! read the program identifier
   read (dataunit) enl%MCprogname
! read the version number
   read (dataunit) enl%MCscversion
! then the name of the crystal data file
   read (dataunit) enl%MCxtalname
! energy information etc...
   read(dataunit) enl%numEbins, enl%numzbins, enl%nsx, enl%nsy, enl%num_el ! , enl%MCnthreads
   enl%nsx = (enl%nsx - 1)/2
   enl%nsy = (enl%nsy - 1)/2
! more energy information
   read (dataunit) enl%EkeV, enl%Ehistmin, enl%Ebinsize, enl%depthmax, enl%depthstep
! angular information
   read (dataunit) enl%MCsig, enl%MComega
! Monte Carlo mode ('CSDA' or other)
   read (dataunit) enl%MCmode
! and finally the actual energy histogram (in square Lambert projection format)
   allocate(acc%accum_e(enl%numEbins,-enl%nsx:enl%nsx,-enl%nsy:enl%nsy),stat=istat)
   read(dataunit) acc%accum_e
   enl%num_el = sum(acc%accum_e)
! we do not need the other array in this energyfile
! read(dataunit) accum_z    ! we only need this array for the depth integrations in EMEBSDmaster.f90
  close(dataunit,status='keep')
end if


if (present(verbose)) call Message(' -> completed reading '//trim(enl%energyfile), frm = "(A)")

end subroutine EBSDreadMCfile


!--------------------------------------------------------------------------
!
! SUBROUTINE:EBSDreadMasterfile
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read EBSD master pattern from file
!
!> @param enl EBSD name list structure
!> @param 
!
!> @date 06/24/14  MDG 1.0 original
!> @date 04/02/15  MDG 2.0 changed program input & output to HDF format
!> @date 09/01/15  MDG 3.0 changed Lambert maps to Northern + Southern maps; lots of changes...
!> @date 09/03/15  MDG 3.1 removed support for old file format (too difficult to maintain after above changes)
!--------------------------------------------------------------------------
subroutine EBSDreadMasterfile(enl, master, mfile, verbose)

use NameListTypedefs
use files
use io
use error
use HDF5
use HDFsupport


IMPLICIT NONE

type(EBSDNameListType),INTENT(INOUT)    :: enl
type(EBSDMasterType),pointer            :: master
character(fnlen),INTENT(IN),OPTIONAL    :: mfile
logical,INTENT(IN),OPTIONAL             :: verbose

real(kind=sgl),allocatable              :: mLPNH(:,:,:) 
real(kind=sgl),allocatable              :: mLPSH(:,:,:) 
real(kind=sgl),allocatable              :: EkeVs(:) 
integer(kind=irg),allocatable           :: atomtype(:)

real(kind=sgl),allocatable              :: srtmp(:,:,:,:)
integer(kind=irg)                       :: istat

logical                                 :: stat, readonly
integer(kind=irg)                       :: hdferr, nlines
integer(HSIZE_T)                        :: dims(1), dims4(4)
character(fnlen)                        :: groupname, dataset, masterfile
character(fnlen),allocatable            :: stringarray(:)

type(HDFobjectStackType),pointer        :: HDF_head

! open the fortran HDF interface
call h5open_f(hdferr)

nullify(HDF_head, HDF_head)

! is the mfile parameter present? If so, use it as the filename, otherwise use the enl%masterfile parameter
if (PRESENT(mfile)) then
  masterfile = mfile
else
  masterfile = trim(EMdatapathname)//trim(enl%masterfile)
end if

! is this a proper HDF5 file ?
call h5fis_hdf5_f(trim(masterfile), stat, hdferr)

if (stat) then 
! open the master file 
  readonly = .TRUE.
  hdferr =  HDF_openFile(masterfile, HDF_head, readonly)

! open the namelist group
  groupname = 'NMLparameters'
  hdferr = HDF_openGroup(groupname, HDF_head)

  groupname = 'EBSDMasterNameList'
  hdferr = HDF_openGroup(groupname, HDF_head)

! read all the necessary variables from the namelist group
  dataset = 'energyfile'
  stringarray = HDF_readDatasetStringArray(dataset, nlines, HDF_head)
  enl%Masterenergyfile = trim(stringarray(1))
  deallocate(stringarray)

  dataset = 'npx'
  enl%npx = HDF_readDatasetInteger(dataset, HDF_head)
  enl%npy = enl%npx

  call HDF_pop(HDF_head)
  call HDF_pop(HDF_head)

  groupname = 'EMData'
  hdferr = HDF_openGroup(groupname, HDF_head)

  dataset = 'numEbins'
  enl%nE = HDF_readDatasetInteger(dataset, HDF_head)
! make sure that MC and Master results are compatible
  if ((enl%numEbins.ne.enl%nE).and.(.not.PRESENT(mfile))) then
    call Message('Energy histogram and Lambert stack have different energy dimension; aborting program', frm = "(A)")
    call HDF_pop(HDF_head,.TRUE.)
    stop
  end if

  dataset = 'numset'
  enl%numset = HDF_readDatasetInteger(dataset, HDF_head)

  dataset = 'squhex'
  stringarray = HDF_readDatasetStringArray(dataset, nlines, HDF_head)
  enl%sqorhe = trim(stringarray(1))
  deallocate(stringarray)

  dataset = 'mLPNH'
  srtmp = HDF_readDatasetFloatArray4D(dataset, dims4, HDF_head)
  allocate(master%mLPNH(-enl%npx:enl%npx,-enl%npy:enl%npy,enl%nE),stat=istat)
  master%mLPNH = sum(srtmp,4)
  deallocate(srtmp)

  dataset = 'mLPSH'
  srtmp = HDF_readDatasetFloatArray4D(dataset, dims4, HDF_head)
  allocate(master%mLPSH(-enl%npx:enl%npx,-enl%npy:enl%npy,enl%nE),stat=istat)
  master%mLPSH = sum(srtmp,4)
  deallocate(srtmp)

  dataset = 'xtalname'
  stringarray = HDF_readDatasetStringArray(dataset, nlines, HDF_head)
  enl%Masterxtalname = trim(stringarray(1))
  deallocate(stringarray)

  call HDF_pop(HDF_head)

  groupname = 'EMheader'
  hdferr = HDF_openGroup(groupname, HDF_head)

  dataset = 'ProgramName'
  stringarray = HDF_readDatasetStringArray(dataset, nlines, HDF_head)
  enl%Masterprogname = trim(stringarray(1))
  deallocate(stringarray)
  
  dataset = 'Version'
  stringarray = HDF_readDatasetStringArray(dataset, nlines, HDF_head)
  enl%Masterscversion = trim(stringarray(1))
  deallocate(stringarray)
  
  call HDF_pop(HDF_head,.TRUE.)

! close the fortran HDF interface
  call h5close_f(hdferr)

else
  masterfile = 'File '//trim(masterfile)//' is not an HDF5 file'
  call FatalError('EBSDreadMasterfile',masterfile)
end if
!====================================

if (present(verbose)) call Message(' -> completed reading '//trim(enl%masterfile), frm = "(A)")

end subroutine EBSDreadMasterfile

!--------------------------------------------------------------------------
!
! SUBROUTINE:EBSDreadMasterfile_overlap
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read EBSD master pattern from file
!
!> @param enl EBSDoverlap name list structure
!> @param 
!
!> @date 06/24/14  MDG 1.0 original
!> @date 04/02/15  MDG 2.0 changed program input & output to HDF format
!> @date 09/03/15  MDG 2.1 removed old file format support
!--------------------------------------------------------------------------
subroutine EBSDreadMasterfile_overlap(enl, master, mfile, verbose)

use NameListTypedefs
use files
use io
use error
use HDF5
use HDFsupport


IMPLICIT NONE

type(EBSDoverlapNameListType),INTENT(INOUT)    :: enl
type(EBSDMasterType),pointer            :: master
character(fnlen),INTENT(IN),OPTIONAL    :: mfile
logical,INTENT(IN),OPTIONAL             :: verbose

real(kind=sgl),allocatable              :: sr(:,:,:) 
real(kind=sgl),allocatable              :: EkeVs(:) 
integer(kind=irg),allocatable           :: atomtype(:)

real(kind=sgl),allocatable              :: srtmp(:,:,:,:)
integer(kind=irg)                       :: istat

logical                                 :: stat, readonly
integer(kind=irg)                       :: hdferr, nlines
integer(HSIZE_T)                        :: dims(1), dims4(4)
character(fnlen)                        :: groupname, dataset, masterfile
character(fnlen),allocatable            :: stringarray(:)

type(HDFobjectStackType),pointer        :: HDF_head

! open the fortran HDF interface
call h5open_f(hdferr)

nullify(HDF_head, HDF_head)

! is the mfile parameter present? If so, use it as the filename, otherwise use the enl%masterfile parameter
if (PRESENT(mfile)) then
  masterfile = trim(EMdatapathname)//trim(mfile)
else
  masterfile = trim(EMdatapathname)//trim(enl%masterfile)
end if

! first, we need to check whether or not the input file is of the HDF5 forat type
call h5fis_hdf5_f(trim(masterfile), stat, hdferr)

if (stat) then 
! open the master file 
  readonly = .TRUE.
  hdferr =  HDF_openFile(masterfile, HDF_head, readonly)

! open the namelist group
  groupname = 'NMLparameters'
  hdferr = HDF_openGroup(groupname, HDF_head)

  groupname = 'EBSDMasterNameList'
  hdferr = HDF_openGroup(groupname, HDF_head)

! read all the necessary variables from the namelist group
  dataset = 'energyfile'
  stringarray = HDF_readDatasetStringArray(dataset, nlines, HDF_head)
  enl%Masterenergyfile = trim(stringarray(1))
  deallocate(stringarray)

  dataset = 'npx'
  enl%npx = HDF_readDatasetInteger(dataset, HDF_head)
  enl%npy = enl%npx

  call HDF_pop(HDF_head)
  call HDF_pop(HDF_head)

  groupname = 'EMData'
  hdferr = HDF_openGroup(groupname, HDF_head)

  dataset = 'numEbins'
  enl%nE = HDF_readDatasetInteger(dataset, HDF_head)

  dataset = 'numset'
  enl%numset = HDF_readDatasetInteger(dataset, HDF_head)

  dataset = 'squhex'
  stringarray = HDF_readDatasetStringArray(dataset, nlines, HDF_head)
  enl%sqorhe = trim(stringarray(1))
  deallocate(stringarray)

  dataset = 'mLPNH'
  srtmp = HDF_readDatasetFloatArray4D(dataset, dims4, HDF_head)
  allocate(master%mLPNH(-enl%npx:enl%npx,-enl%npy:enl%npy,enl%nE),stat=istat)
  master%mLPNH = sum(srtmp,4)
  deallocate(srtmp)

  dataset = 'mLPSH'
  srtmp = HDF_readDatasetFloatArray4D(dataset, dims4, HDF_head)
  allocate(master%mLPSH(-enl%npx:enl%npx,-enl%npy:enl%npy,enl%nE),stat=istat)
  master%mLPSH = sum(srtmp,4)
  deallocate(srtmp)

  dataset = 'xtalname'
  stringarray = HDF_readDatasetStringArray(dataset, nlines, HDF_head)
  enl%Masterxtalname = trim(stringarray(1))
  deallocate(stringarray)

  call HDF_pop(HDF_head)

  groupname = 'EMheader'
  hdferr = HDF_openGroup(groupname, HDF_head)

  dataset = 'ProgramName'
  stringarray = HDF_readDatasetStringArray(dataset, nlines, HDF_head)
  enl%Masterprogname = trim(stringarray(1))
  deallocate(stringarray)
  
  dataset = 'Version'
  stringarray = HDF_readDatasetStringArray(dataset, nlines, HDF_head)
  enl%Masterscversion = trim(stringarray(1))
  deallocate(stringarray)
  
  call HDF_pop(HDF_head,.TRUE.)

! close the fortran HDF interface
  call h5close_f(hdferr)

else
  masterfile = 'File '//trim(masterfile)//' is not an HDF5 file'
  call FatalError('EBSDreadMasterfile_overlap',masterfile)
end if
!====================================

if (present(verbose)) then
  if (verbose) call Message(' -> completed reading '//trim(masterfile), frm = "(A)")
end if

end subroutine EBSDreadMasterfile_overlap




!--------------------------------------------------------------------------
!
! SUBROUTINE:EBSDGenerateDetector
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief generate the detector arrays
!
!> @param enl EBSD name list structure
!
!> @date 06/24/14  MDG 1.0 original
!> @date 07/01/15   SS  1.1 added omega as the second tilt angle
!> @date 07/07/15   SS  1.2 correction to the omega tilt parameter; old version in the comments

!--------------------------------------------------------------------------
subroutine EBSDGenerateDetector(enl, acc, master, verbose)

use NameListTypedefs
use files
use constants
use io
use Lambert

IMPLICIT NONE

type(EBSDNameListType),INTENT(INOUT)    :: enl
type(EBSDLargeAccumType),pointer        :: acc
type(EBSDMasterType),pointer            :: master
logical,INTENT(IN),OPTIONAL             :: verbose

real(kind=sgl),allocatable              :: scin_x(:), scin_y(:)                 ! scintillator coordinate ararays [microns]
real(kind=sgl),parameter                :: dtor = 0.0174533  ! convert from degrees to radians
real(kind=sgl)                          :: alp, ca, sa, cw, sw
real(kind=sgl)                          :: L2, Ls, Lc     ! distances
real(kind=sgl),allocatable              :: z(:,:)           
integer(kind=irg)                       :: nix, niy, binx, biny , i, j, Emin, Emax, istat, k      ! various parameters
real(kind=sgl)                          :: dc(3), scl           ! direction cosine array
real(kind=sgl)                          :: sx, dx, dxm, dy, dym, rhos, x, bindx         ! various parameters
real(kind=sgl)                          :: ixy(2)


!====================================
! ------ generate the detector arrays
!====================================
! This needs to be done only once for a given detector geometry
allocate(scin_x(enl%numsx),scin_y(enl%numsy),stat=istat)
! if (istat.ne.0) then ...
scin_x = - ( enl%xpc - ( 1.0 - enl%numsx ) * 0.5 - (/ (i-1, i=1,enl%numsx) /) ) * enl%delta
scin_y = ( enl%ypc - ( 1.0 - enl%numsy ) * 0.5 - (/ (i-1, i=1,enl%numsy) /) ) * enl%delta

! auxiliary angle to rotate between reference frames
alp = 0.5 * cPi - (enl%MCsig - enl%thetac) * dtor
ca = cos(alp)
sa = sin(alp)

cw = cos(enl%omega * dtor)
sw = sin(enl%omega * dtor)

! we will need to incorporate a series of possible distortions 
! here as well, as described in Gert nolze's paper; for now we 
! just leave this place holder comment instead

! compute auxilliary interpolation arrays
! if (istat.ne.0) then ...

L2 = enl%L * enl%L
do j=1,enl%numsx
  sx = L2 + scin_x(j) * scin_x(j)
  Ls = -sw * scin_x(j) + enl%L*cw
  Lc = cw * scin_x(j) + enl%L*sw
  do i=1,enl%numsy
   rhos = 1.0/sqrt(sx + scin_y(i)**2)
   master%rgx(j,i) = (scin_y(i) * ca + sa * Ls) * rhos!Ls * rhos
   master%rgy(j,i) = Lc * rhos!(scin_x(i) * cw + Lc * sw) * rhos
   master%rgz(j,i) = (-sa * scin_y(i) + ca * Ls) * rhos!(-sw * scin_x(i) + Lc * cw) * rhos
  end do
end do
deallocate(scin_x, scin_y)

! normalize the direction cosines.
allocate(z(enl%numsx,enl%numsy))
  z = 1.0/sqrt(master%rgx*master%rgx+master%rgy*master%rgy+master%rgz*master%rgz)
  master%rgx = master%rgx*z
  master%rgy = master%rgy*z
  master%rgz = master%rgz*z
deallocate(z)
!====================================

!====================================
! ------ create the equivalent detector energy array
!====================================
! from the Monte Carlo energy data, we need to extract the relevant
! entries for the detector geometry defined above.  Once that is 
! done, we can get rid of the larger energy array
!
! in the old version, we either computed the background model here, or 
! we would load a background pattern from file.  In this version, we are
! using the background that was computed by the MC program, and has 
! an energy histogram embedded in it, so we need to interpolate this 
! histogram to the pixels of the scintillator.  In other words, we need
! to initialize a new accum_e array for the detector by interpolating
! from the Lambert projection of the MC results.
!

! determine the scale factor for the Lambert interpolation; the square has
! an edge length of 2 x sqrt(pi/2)
  scl = float(enl%nsx) !  / LPs%sPio2  [removed on 09/01/15 by MDG for new Lambert routines]

! get the indices of the minimum and maximum energy
  Emin = nint((enl%energymin - enl%Ehistmin)/enl%Ebinsize) +1
  if (Emin.lt.1)  Emin=1
  if (Emin.gt.enl%numEbins)  Emin=enl%numEbins

  Emax = nint((enl%energymax - enl%Ehistmin)/enl%Ebinsize) +1
  if (Emax.lt.1)  Emax=1
  if (Emax.gt.enl%numEbins)  Emax=enl%numEbins

  do i=1,enl%numsx
    do j=1,enl%numsy
! do the coordinate transformation for this detector pixel
       dc = (/ master%rgx(i,j),master%rgy(i,j),master%rgz(i,j) /)
! make sure the third one is positive; if not, switch all 
       if (dc(3).lt.0.0) dc = -dc
! convert these direction cosines to coordinates in the Rosca-Lambert projection
        ixy = scl * LambertSphereToSquare( dc, istat )
        x = ixy(1)
        ixy(1) = ixy(2)
        ixy(2) = -x
! four-point interpolation (bi-quadratic)
        nix = int(enl%nsx+ixy(1))-enl%nsx
        niy = int(enl%nsy+ixy(2))-enl%nsy
        dx = ixy(1)-nix
        dy = ixy(2)-niy
        dxm = 1.0-dx
        dym = 1.0-dy
! interpolate the intensity 
        do k=Emin,Emax 
          acc%accum_e_detector(k,i,j) =   acc%accum_e(k,nix,niy) * dxm * dym + &
                                        acc%accum_e(k,nix+1,niy) * dx * dym + &
                                        acc%accum_e(k,nix,niy+1) * dxm * dy + &
                                        acc%accum_e(k,nix+1,niy+1) * dx * dy
        end do
    end do
  end do 
  acc%accum_e_detector = acc%accum_e_detector * 0.25



! and finally, get rid of the original accum_e array which is no longer needed
! [we'll do that in the calling program ]
!  deallocate(accum_e)

!====================================
end subroutine EBSDGenerateDetector

!--------------------------------------------------------------------------
!
! SUBROUTINE:TwinCubicMasterPattern
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief Generate a master pattern with regular and twin master pattern overlapped, both with 50% weights
!
!> @param enl EBSD name list structure
!> @param master  EBSDMasterType pointer
!
!> @date 04/16/15  SS 1.0 original
!> @date 04/20/15 MDG 1.1 minor edits
!> @date 09/03/15 MDG 1.2 added support for Northern and Southern Lambert hemispheres
!--------------------------------------------------------------------------
subroutine TwinCubicMasterPattern(enl,master)

use local
use io
use quaternions
use Lambert
use rotations
use NameListTypedefs
use NameListHandlers
use constants

IMPLICIT NONE

type(EBSDNameListType),INTENT(INOUT)                :: enl
type(EBSDMasterType),pointer                        :: master

real(kind=dbl),allocatable                          :: master_twinNH(:,:,:), master_twinSH(:,:,:)
type(EBSDLargeAccumType),pointer                    :: acc
logical                                             :: verbose
real(kind=dbl)                                      :: q(4),Lamproj(2),dc(3),dc_new(3),dx,dy,dxm,dym,ixy(2),scl
integer(kind=irg)                                   :: nix,niy,nixp,niyp
integer(kind=irg)                                   :: ii,jj,kk,ierr,istat,pp,qq

allocate(master_twinNH(-enl%npx:enl%npx,-enl%npy:enl%npy,1:enl%nE),stat=istat)
allocate(master_twinSH(-enl%npx:enl%npx,-enl%npy:enl%npy,1:enl%nE),stat=istat)

q = (/ dsqrt(3.D0)/2.D0,1/dsqrt(3.D0)/2.D0,1/dsqrt(3.D0)/2.D0,1/dsqrt(3.D0)/2.D0 /)

scl = float(enl%npx) ! / LPs%sPio2 [removed 09/01/15 by MDG for new Lambert module]

    master_twinNH = 0.0
    master_twinSH = 0.0
    do jj = -enl%npx,enl%npx
        do kk = -enl%npy,enl%npy

            Lamproj = (/ float(jj)/scl,float(kk)/scl /)
            dc = LambertSquareToSphere(Lamproj,ierr)
            dc_new = quat_Lp(conjg(q),dc)
            dc_new = dc_new/sqrt(sum(dc_new**2))
            if (dc_new(3) .lt. 0.0) dc_new = -dc_new

! convert direction cosines to lambert projections
            ixy = scl * LambertSphereToSquare( dc_new, istat )
! interpolate intensity from the neighboring points

            nix = floor(ixy(1))
            niy = floor(ixy(2))
            nixp = nix+1
            niyp = niy+1
            if (nixp.gt.enl%npx) nixp = nix
            if (niyp.gt.enl%npy) niyp = niy
            dx = ixy(1) - nix
            dy = ixy(2) - niy
            dxm = 1.0 - dx
            dym = 1.0 - dy

            master_twinNH(jj,kk,1:enl%nE) = master%mLPNH(nix,niy,1:enl%nE)*dxm*dym + master%mLPNH(nixp,niy,1:enl%nE)*dx*dym + &
                                    master%mLPNH(nix,niyp,1:enl%nE)*dxm*dy + master%mLPNH(nixp,niyp,1:enl%nE)*dx*dy
            master_twinSH(jj,kk,1:enl%nE) = master%mLPSH(nix,niy,1:enl%nE)*dxm*dym + master%mLPSH(nixp,niy,1:enl%nE)*dx*dym + &
                                    master%mLPSH(nix,niyp,1:enl%nE)*dxm*dy + master%mLPSH(nixp,niyp,1:enl%nE)*dx*dy
        end do
    end do
master%mLPNH = 0.5D0 * (master_twinNH + master%mLPNH)
master%mLPSH = 0.5D0 * (master_twinSH + master%mLPSH)

call Message(' -> completed superimposing twin and regular master patterns', frm = "(A)")

end subroutine TwinCubicMasterPattern

!--------------------------------------------------------------------------
!
! SUBROUTINE:OverlapMasterPattern
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Generate a master pattern with regular and rotated master pattern overlapped, both with 50% weights
!
!> @param enl EBSD name list structure
!> @param master  EBSDMasterType pointer
!> @param q unit quaternion providing the necessary rotation
!
!> @date 04/20/15 MDG 1.0 original, based on Saransh's twin routine above
!> @date 09/03/15 MDG 1.2 added support for Northern and Southern Lambert hemispheres
!> @date 10/16/15  SS 1.3 added alpha parameter for degree of mixing; 0<=alpha<=1; also
!> added master_inp and master_out variables in the subroutine for NtAg dataset
!--------------------------------------------------------------------------
recursive subroutine OverlapMasterPattern(enl,master_in,master_out,q,alpha)

use local
use io
use error
use quaternions
use Lambert
use rotations
use NameListTypedefs
use NameListHandlers
use constants

IMPLICIT NONE

type(EBSDNameListType),INTENT(INOUT)                :: enl
type(EBSDMasterType),pointer                        :: master_in, master_out
real(kind=dbl),INTENT(IN)                           :: q(4)
real(kind=dbl),INTENT(IN)                           :: alpha

real(kind=dbl),allocatable                          :: master_rotatedNH(:,:,:), master_rotatedSH(:,:,:)
type(EBSDLargeAccumType),pointer                    :: acc
logical                                             :: verbose
real(kind=dbl)                                      :: Lamproj(2),dc(3),dc_new(3),dx,dy,dxm,dym,ixy(2),scl
integer(kind=irg)                                   :: nix,niy,nixp,niyp
integer(kind=irg)                                   :: ii,jj,kk,ierr,istat,pp,qq

if (alpha .lt. 0.D0) then
    call FatalError('OverlapMasterPattern','value of mixing paramter is less than zero')
end if

if (alpha .gt. 1.D0) then
    call FatalError('OverlapMasterPattern','value of mixing paramter is greater than one')
end if

allocate(master_rotatedNH(-enl%npx:enl%npx,-enl%npy:enl%npy,1:enl%nE),stat=istat)
allocate(master_rotatedSH(-enl%npx:enl%npx,-enl%npy:enl%npy,1:enl%nE),stat=istat)

if (allocated(master_out%mLPNH)) deallocate(master_out%mLPNH)
if (allocated(master_out%mLPSH)) deallocate(master_out%mLPSH)

allocate(master_out%mLPNH(-enl%npx:enl%npx,-enl%npy:enl%npy,1:enl%nE),stat=istat)
allocate(master_out%mLPSH(-enl%npx:enl%npx,-enl%npy:enl%npy,1:enl%nE),stat=istat)

master_out%mLPNH = 0.0
master_out%mLPSH = 0.0

scl = float(enl%npx) ! / LPs%sPio2 [ removed on 09/01/15 by MDG for new Lambert module]

master_rotatedNH = 0.0
master_rotatedSH = 0.0
do jj = -enl%npx,enl%npx
    do kk = -enl%npy,enl%npy

        Lamproj = (/ float(jj)/scl,float(kk)/scl /)
        dc = LambertSquareToSphere(Lamproj,ierr)
        dc_new = quat_Lp(conjg(q), dc)
        dc_new = dc_new/sqrt(sum(dc_new**2))
        if (dc_new(3) .lt. 0.0) dc_new = -dc_new

! convert direction cosines to lambert projections
        ixy = scl * LambertSphereToSquare( dc_new, istat )

! interpolate intensity from the neighboring points
        nix = floor(ixy(1))
        niy = floor(ixy(2))
        nixp = nix+1
        niyp = niy+1
        if (nixp.gt.enl%npx) nixp = nix
        if (niyp.gt.enl%npy) niyp = niy
        dx = ixy(1) - nix
        dy = ixy(2) - niy
        dxm = 1.0 - dx
        dym = 1.0 - dy

        master_rotatedNH(jj,kk,1:enl%nE) = master_in%mLPNH(nix,niy,1:enl%nE)*dxm*dym + master_in%mLPNH(nixp,niy,1:enl%nE)&
                                    *dx*dym + master_in%mLPNH(nix,niyp,1:enl%nE)*dxm*dy + master_in%mLPNH(nixp,niyp,1:enl%nE)&
                                    *dx*dy
        master_rotatedSH(jj,kk,1:enl%nE) = master_in%mLPSH(nix,niy,1:enl%nE)*dxm*dym + master_in%mLPSH(nixp,niy,1:enl%nE)&
                                    *dx*dym + master_in%mLPSH(nix,niyp,1:enl%nE)*dxm*dy + master_in%mLPSH(nixp,niyp,1:enl%nE)*dx*dy
    end do
end do

master_out%mLPNH = (1 - alpha) * master_rotatedNH + alpha * master_in%mLPNH
master_out%mLPSH = (1 - alpha) * master_rotatedSH + alpha * master_in%mLPSH

call Message(' -> completed superimposing rotated and regular master patterns', frm = "(A)")

end subroutine OverlapMasterPattern

!--------------------------------------------------------------------------
!
! SUBROUTINE:GenerateBackground
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief Generate a binned and normalized background for the dictionary patterns using the monte carlo simulation
!
!> @param enl EBSD name list structure
!> @param master  EBSDMasterType pointer
!> @param q unit quaternion providing the necessary rotation
!
!> @date 04/20/15 MDG 1.0 original, based on Saransh's twin routine above
!--------------------------------------------------------------------------
subroutine GenerateBackground(enl,acc,EBSDBackground)

use local 
use typedefs
use NameListTypedefs

type(EBSDNameListType),INTENT(IN)       :: enl
type(EBSDLargeAccumType),pointer        :: acc
real(kind=sgl),INTENT(OUT)              :: EBSDBackground(enl%numsx/enl%binning,enl%numsy/enl%binning)

integer(kind=irg)                       :: ii, jj, kk, istat
real(kind=sgl),allocatable              :: EBSDtmp(:,:)
integer(kind=irg)                       :: Emin, Emax
real(kind=sgl)                          :: bindx


allocate(EBSDtmp(enl%numsx,enl%numsy),stat=istat)
EBSDtmp = 0.0
EBSDBackground = 0.0
! get the indices of the minimum and maximum energy
Emin = nint((enl%energymin - enl%Ehistmin)/enl%Ebinsize) +1
if (Emin.lt.1)  Emin=1
if (Emin.gt.enl%numEbins)  Emin=enl%numEbins

Emax = nint((enl%energymax - enl%Ehistmin)/enl%Ebinsize) +1
if (Emax.lt.1)  Emax=1
if (Emax.gt.enl%numEbins)  Emax=enl%numEbins

bindx = 1.0/float(enl%binning)**2

do ii = 1,enl%numsx
   do jj = 1,enl%numsy
      do kk = Emin,Emax
         EBSDtmp(ii,jj) = EBSDtmp(ii,jj) + acc%accum_e_detector(kk,ii,jj) 
      end do
   end do
end do

if(enl%binning .ne. 1) then
  do ii=1,enl%numsx/enl%binning
      do jj=1,enl%numsy/enl%binning
           EBSDBackground(ii,jj) = sum(EBSDtmp((ii-1)*enl%binning+1:ii*enl%binning,(jj-1)*enl%binning:jj*enl%binning))
           if(isnan(EBSDBackground(ii,jj))) then
               stop 'Background pattern encountered NaN during binning'
           end if
      end do
  end do  
! and divide by binning^2
  EBSDBackground = EBSDBackground * bindx
else
   EBSDBackground = EBSDtmp
end if

! apply gamma scaling
EBSDBackground = EBSDBackground**enl%gammavalue

! normalize the pattern
EBSDBackground = EBSDBackground/NORM2(EBSDBackground)

end subroutine GenerateBackground


end module EBSDmod

