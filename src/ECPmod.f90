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

type ECPLargeAccumType
        integer(kind=irg),allocatable   :: accum_z(:,:,:,:)
end type ECPLargeAccumType

type ECPMasterType
        real(kind=sgl),allocatable      :: mLPNH(:,:) , mLPSH(:,:)
end type ECPMasterType
 
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

  dataset = 'MCmode'
  stringarray = HDF_readDatasetStringArray(dataset, nlines, HDF_head)
  enl%MCmode = trim(stringarray(1))
  deallocate(stringarray)

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
  dataset = 'numEbins'
  enl%numEbins = HDF_readDatasetInteger(dataset, HDF_head)

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
! SUBROUTINE:EBSDreadMasterfile
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

end module ECPmod

