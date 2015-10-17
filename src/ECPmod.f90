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
!--------------------------------------------------------------------------
! EMsoft:ECPmod.f90
!--------------------------------------------------------------------------
!
! MODULE: ECPmod
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief EMECP helper routines
!
!> @date 09/15/15 SS 1.0 original
!---------------------------------------------------------------------------
module ECPmod

use local 
use typedefs

IMPLICIT NONE

type ECPAngleType
        real(kind=dbl),allocatable      :: quatang(:,:)
end type ECPAngleType

type ECPLargeAccumType
        integer(kind=irg),allocatable   :: accum_z(:,:,:,:)
end type ECPLargeAccumType

type ECPMasterType
        real(kind=sgl),allocatable      :: mLPNH(:,:) , mLPSH(:,:)
end type ECPMasterType

type IncidentListECP
        integer(kind=irg)               :: i, j
        real(kind=dbl)                  :: k(3)
        type(IncidentListECP),pointer   :: next
end type IncidentListECP


contains

!--------------------------------------------------------------------------
!
! SUBROUTINE:ECPreadMCfile
!
!> @author Marc De Graef/Saransh Singh, Carnegie Mellon University
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
!> @date 09/15/15  SS  3.0 made part of ECPmod module
!> @date 10/12/15  SS  3.1 changes to handle new mc program; old version of mc file
!>                         not supported anymore
!--------------------------------------------------------------------------
subroutine ECPreadMCfile(enl,acc,efile,verbose)

use NameListTypedefs
use files
use io
use HDF5
use HDFsupport
use error

IMPLICIT NONE

type(ECPNameListType),INTENT(INOUT)     :: enl
type(ECPLargeAccumType),pointer         :: acc
character(fnlen),INTENT(IN),OPTIONAL    :: efile
logical,INTENT(IN),OPTIONAL             :: verbose

integer(kind=irg)                       :: istat, hdferr, nlines, nx
logical                                 :: stat, readonly
integer(HSIZE_T)                        :: dims3(3),dims4(4)
character(fnlen)                        :: groupname, dataset, energyfile 
character(fnlen),allocatable            :: stringarray(:)

integer(kind=irg),allocatable           :: acc_z(:,:,:,:)

type(HDFobjectStackType),pointer        :: HDF_head

! is the efile parameter present? If so, use it as the filename, otherwise use the enl%energyfile parameter
if (PRESENT(efile)) then
  energyfile = efile
else
  energyfile = trim(Emdatapathname)//trim(enl%energyfile)
end if

allocate(acc)

! first, we need to check whether or not the input file is of the HDF5 format type; if
! it is, we read it accordingly, otherwise we give error. Old format not supported anymore
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

  if(enl%MCmode .ne. 'bse1') then
     call FatalError('ECPreadMCfile','This file is not bse1 mode. Please input correct HDF5 file')
  end if

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

  dataset = 'sigstart'
  enl%MCsigstart = HDF_readDatasetDouble(dataset, HDF_head)

  dataset = 'sigend'
  enl%MCsigend = HDF_readDatasetDouble(dataset, HDF_head)

  dataset = 'sigstep'
  enl%MCsigstep = HDF_readDatasetDouble(dataset, HDF_head)

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
  dataset = 'numangle'
  enl%numangle = HDF_readDatasetInteger(dataset, HDF_head)

  dataset = 'numzbins'
  enl%numzbins = HDF_readDatasetInteger(dataset, HDF_head)

  dataset = 'accum_z'
  acc_z = HDF_readDatasetIntegerArray4D(dataset, dims4, HDF_head)
  allocate(acc%accum_z(1:dims4(1),1:dims4(2),1:dims4(3),1:dims4(4)))
  acc%accum_z = acc_z
  deallocate(acc_z)
 
  enl%num_el = sum(acc%accum_z)

! and close everything
  call HDF_pop(HDF_head,.TRUE.)

! close the fortran HDF interface
  call h5close_f(hdferr)

else
!==============================================
! ------ read old format Monte Carlo data file
!==============================================

!==============================================
! OLD VERSION OF MC FILE NOT SUPPORTED ANYMORE
! COMMENTING OUT THE FOLLOWING LINES
! REPLACING WITH FATALERROR COMMENT
!==============================================

  call FatalError('ECPmod (ECPreadMCfile)','The file is not a h5 file. Old version of MC file not supported anymore!')
  !if (present(verbose)) call Message('opening '//trim(enl%energyfile), frm = "(A)")

  !open(dataunit,file=trim(energyfile),status='unknown',form='unformatted')

! read the program identifier
   !read (dataunit) enl%MCprogname
! read the version number
   !read (dataunit) enl%MCscversion
! then the name of the crystal data file
   !read (dataunit) enl%MCxtalname
! energy information etc...
   !read(dataunit) enl%numEbins, enl%numzbins, enl%nsx, enl%nsy, enl%num_el ! , enl%MCnthreads
   !enl%nsx = (enl%nsx - 1)/2
   !enl%nsy = (enl%nsy - 1)/2
! more energy information
   !read (dataunit) enl%EkeV, enl%Ehistmin, enl%Ebinsize, enl%depthmax, enl%depthstep
! angular information
   !read (dataunit) enl%MCsig, enl%MComega
! Monte Carlo mode ('CSDA' or other)
   !read (dataunit) enl%MCmode
! and finally the actual energy histogram (in square Lambert projection format)
   !allocate(acc%accum_z(enl%numEbins,enl%numzbins,-enl%nsx:enl%nsx,-enl%nsy:enl%nsy),stat=istat)
   !read(dataunit) acc%accum_z
   !enl%num_el = sum(acc%accum_z)
! we do not need the other array in this energyfile
! read(dataunit) accum_z    ! we only need this array for the depth integrations in EMEBSDmaster.f90
  !close(dataunit,status='keep')
end if

if (present(verbose)) then

if (verbose) call Message(' -> completed reading '//trim(enl%energyfile), frm = "(A)")

end if

end subroutine ECPreadMCfile

!--------------------------------------------------------------------------
!
! SUBROUTINE:ECPreadMasterfile
!
!> @author Saransh Singh/Marc De Graef, Carnegie Mellon University
!
!> @brief read EBSD master pattern from file
!
!> @param enl EBSD name list structure
!
!> @date 06/24/14  MDG 1.0 original
!> @date 04/02/15  MDG 2.0 changed program input & output to HDF format
!> @date 09/01/15  MDG 3.0 changed Lambert maps to Northern + Southern maps; lots of changes...
!> @date 09/03/15  MDG 3.1 removed support for old file format (too difficult to maintain after above changes)
!> @date 09/15/15  SS  4.0 modified for ECP master program
!--------------------------------------------------------------------------
subroutine ECPreadMasterfile(enl, master, mfile, verbose)

use NameListTypedefs
use files
use io
use error
use HDF5
use HDFsupport


IMPLICIT NONE

type(ECPNameListType),INTENT(INOUT)     :: enl
type(ECPMasterType),pointer             :: master
character(fnlen),INTENT(IN),OPTIONAL    :: mfile
logical,INTENT(IN),OPTIONAL             :: verbose

real(kind=sgl),allocatable              :: mLPNH(:,:) 
real(kind=sgl),allocatable              :: mLPSH(:,:) 
real(kind=sgl),allocatable              :: EkeVs(:) 
integer(kind=irg),allocatable           :: atomtype(:)

real(kind=sgl),allocatable              :: srtmp(:,:,:)
integer(kind=irg)                       :: istat

logical                                 :: stat, readonly
integer(kind=irg)                       :: hdferr, nlines
integer(HSIZE_T)                        :: dims(1), dims3(3)
character(fnlen)                        :: groupname, dataset, masterfile
character(fnlen),allocatable            :: stringarray(:)

type(HDFobjectStackType),pointer        :: HDF_head

allocate(master)

! open the fortran HDF interface
call h5open_f(hdferr)

nullify(HDF_head, HDF_head)

! is the mfile parameter present? If so, use it as the filename, otherwise use the enl%masterfile parameter
if (PRESENT(mfile)) then
  masterfile = mfile
else
  masterfile = trim(EMdatapathname)//trim(enl%masterfile)
end if

! is this a propoer HDF5 file ?
call h5fis_hdf5_f(trim(masterfile), stat, hdferr)

if (stat) then 
! open the master file 
  readonly = .TRUE.
  hdferr =  HDF_openFile(masterfile, HDF_head, readonly)

! open the namelist group
  groupname = 'NMLparameters'
  hdferr = HDF_openGroup(groupname, HDF_head)

  groupname = 'ECPMasterNameList'
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

  dataset = 'EkeV'
  enl%EkeV = HDF_readDatasetfloat(dataset, HDF_head)
  
  dataset = 'numset'
  enl%numset = HDF_readDatasetInteger(dataset, HDF_head)

  dataset = 'squhex'
  stringarray = HDF_readDatasetStringArray(dataset, nlines, HDF_head)
  enl%sqorhe = trim(stringarray(1))
  deallocate(stringarray)

  dataset = 'mLPNH'
  srtmp = HDF_readDatasetFloatArray3D(dataset, dims3, HDF_head)
  allocate(master%mLPNH(-enl%npx:enl%npx,-enl%npy:enl%npy),stat=istat)
  master%mLPNH = sum(srtmp,3)
  deallocate(srtmp)

  dataset = 'mLPSH'
  srtmp = HDF_readDatasetFloatArray3D(dataset, dims3, HDF_head)
  allocate(master%mLPSH(-enl%npx:enl%npx,-enl%npy:enl%npy),stat=istat)
  master%mLPSH = sum(srtmp,3)
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

end subroutine ECPreadMasterfile

!--------------------------------------------------------------------------
!
! SUBROUTINE:GetVectorsCone
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief generate list of incident vectors for interpolation of ECP
!
!> @param ecpnl ECP namelist structure
!> @param klist IncidentListECP pointer
!> @param rotmat rotation matrix for the microscope to grain reference frame
!> @param numk number of incident vectors in the linked list
!
!> @date 10/12/15  SS 1.0 original
!--------------------------------------------------------------------------
subroutine GetVectorsCone(ecpnl, klist, numk)

use local
use io
use NameListTypedefs
use error

type(ECPNameListType),INTENT(IN)                 :: ecpnl
type(IncidentListECP),pointer                    :: klist, ktmp
integer(kind=irg),INTENT(OUT)                    :: numk

real(kind=dbl)                                   :: kk(3), thetacr, delta, ktmax
real(kind=dbl),parameter                         :: DtoR = 0.01745329251D0
integer(kind=irg)                                :: imin, imax, jmin, jmax
integer(kind=irg)                                :: ii, jj, istat

numk = 0
kk = (/0.D0,0.D0,1.D0/)
thetacr = DtoR*ecpnl%thetac
ktmax = tan(thetacr)
delta = 2.0*ktmax/(2.0*float(ecpnl%npix)+1.0)

imin = -ecpnl%npix
imax = ecpnl%npix
jmin = -ecpnl%npix
jmax = ecpnl%npix

allocate(klist,stat=istat)
if (istat .ne. 0) then
    call FatalError('GetVectorsCone','Failed to allocate klist pointer')
end if

ktmp => klist
nullify(ktmp%next)
ktmp%i = 0
ktmp%j = 0
numk = numk + 1

do ii = imin, imax
    do jj = jmin, jmax
        if (abs(ii) + abs(jj) .ne. 0) then
           allocate(ktmp%next,stat=istat)
           ktmp => ktmp%next
           nullify(ktmp%next)
           ktmp%k(1:3) = (/delta*ii,delta*jj,0.D0/) + kk
           ktmp%k = ktmp%k/sqrt(sum(ktmp%k**2))
           ktmp%i = ii
           ktmp%j = jj
           numk = numk + 1
        end if
    end do
end do

end subroutine GetVectorsCone

!--------------------------------------------------------------------------
!
! SUBROUTINE:EBSDreadangles
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief read angles from an angle file
!
!> @param enl EBSD name list structure
!> @param quatang array of unit quaternions (output)
!
!> @date 10/12/15  SS 1.0 original
!--------------------------------------------------------------------------
subroutine ECPreadangles(enl,angles,verbose)

use NameListTypedefs
use io
use files
use quaternions
use rotations

IMPLICIT NONE


type(ECPNameListType),INTENT(INOUT)     :: enl
type(ECPAngleType),pointer              :: angles
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
open(unit=dataunit,file=trim(enl%anglefile),status='old',action='read')

! get the type of angle first [ 'eu' or 'qu' ]
read(dataunit,*) angletype
if (angletype.eq.'eu') then 
  enl%anglemode = 'euler'
else
  enl%anglemode = 'quats'
end if

! then the number of angles in the file
read(dataunit,*) enl%numangle_anglefile

if (present(verbose)) then 
  io_int(1) = enl%numangle_anglefile
  call WriteValue('Number of angle entries = ',io_int,1)
end if

if (enl%anglemode.eq.'euler') then
! allocate the euler angle array
  allocate(eulang(3,enl%numangle_anglefile),stat=istat)
! if istat.ne.0 then do some error handling ... 
  do i=1,enl%numangle_anglefile
    read(dataunit,*) eulang(1:3,i)
  end do
  close(unit=dataunit,status='keep')

  if (enl%eulerconvention.eq.'hkl') then
    if (present(verbose)) call Message('  -> converting Euler angles to TSL representation', frm = "(A/)")
    eulang(1,1:enl%numangle_anglefile) = eulang(1,1:enl%numangle_anglefile) + 90.0
  end if

! convert the euler angle triplets to quaternions
  allocate(angles%quatang(4,enl%numangle_anglefile),stat=istat)
! if (istat.ne.0) then ...

  if (present(verbose)) call Message('  -> converting Euler angles to quaternions', frm = "(A/)")
  
  do i=1,enl%numangle_anglefile
    angles%quatang(1:4,i) = eu2qu(eulang(1:3,i)*dtor)
  end do

else
! the input file has quaternions, not Euler triplets
  allocate(angles%quatang(4,enl%numangle_anglefile),stat=istat)
  do i=1,enl%numangle_anglefile
    read(dataunit,*) angles%quatang(1:4,i)
  end do
end if

close(unit=dataunit,status='keep')

!====================================
! Do we need to apply an additional axis-angle pair rotation to all the quaternions ?

! commented out for now; needs to be verified


!if (enl%axisangle(4).ne.0.0) then
!  enl%axisangle(4) = enl%axisangle(4) * dtor
!  qax = ax2qu( enl%axisangle )
!  do i=1,enl%numangles_
!    angles%quatang(1:4,i) = quat_mult(qax,angles%quatang(1:4,i))
!  end do 
!end if

write (*,*) 'completed reading Euler angles'

end subroutine ECPreadangles

end module ECPmod

