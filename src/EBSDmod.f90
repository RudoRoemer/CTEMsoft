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
!--------------------------------------------------------------------------
module EBSDmod

use local
use typedefs

IMPLICIT NONE

type EBSDAngleType
        real(kind=sgl),allocatable      :: quatang(:,:)
end type EBSDAngleType

type EBSDLargeAccumType
        integer(kind=irg),allocatable   :: accum_e(:,:,:)
        real(kind=irg),allocatable      :: accum_e_detector(:,:,:)
end type EBSDLargeAccumType

type EBSDMasterType
        real(kind=sgl),allocatable      :: sr(:,:,:) 
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
open(unit=dataunit,file=trim(enl%anglefile),status='old',action='read')

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
!> @param accum_e energy histogram (square Lambert projection)
!
!> @date 06/24/14  MDG 1.0 original
!> @date 11/18/14  MDG 1.1 removed enl%MCnthreads from file read
!> @date 04/02/15  MDG 2.0 changed program input & output to HDF format
!--------------------------------------------------------------------------
subroutine EBSDreadMCfile(enl,acc,verbose)

use NameListTypedefs
use files
use io
use HDF5
use HDFsupport

IMPLICIT NONE

type(EBSDNameListType),INTENT(INOUT)    :: enl
type(EBSDLargeAccumType),pointer        :: acc
logical,INTENT(IN),OPTIONAL             :: verbose

integer(kind=irg)                       :: istat, hdferr, nlines, nx
logical                                 :: stat, readonly
integer(HSIZE_T)                        :: dims3(3)
character(fnlen)                        :: groupname, dataset 
character(fnlen),allocatable            :: stringarray(:)

integer(kind=irg),allocatable           :: acc_e(:,:,:)

type(HDFobjectStackType),pointer        :: HDF_head


! first, we need to check whether or not the input file is of the HDF5 forat type; if
! it is, we read it accordingly, otherwise we use the old binary format.
!
call h5fis_hdf5_f(trim(enl%energyfile), stat, hdferr)

if (stat) then 
! open the fortran HDF interface
  call h5open_f(hdferr)

  nullify(HDF_head)

! open the MC file using the default properties.
  readonly = .TRUE.
  hdferr =  HDF_openFile(enl%energyfile, HDF_head, readonly)

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

! and close everything
  call HDF_pop(HDF_head,.TRUE.)

! close the fortran HDF interface
  call h5close_f(hdferr)

else
!====================================
! ------ read old format Monte Carlo data file
!====================================
  if (present(verbose)) call Message('opening '//trim(enl%energyfile), frm = "(A)")

  open(dataunit,file=trim(enl%energyfile),status='unknown',form='unformatted')

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
!> @param quatang array of unit quaternions (output)
!
!> @date 06/24/14  MDG 1.0 original
!> @date 04/02/15  MDG 2.0 changed program input & output to HDF format
!--------------------------------------------------------------------------
subroutine EBSDreadMasterfile(enl, master, verbose)

use NameListTypedefs
use files
use io
use HDF5
use HDFsupport


IMPLICIT NONE

type(EBSDNameListType),INTENT(INOUT)    :: enl
type(EBSDMasterType),pointer            :: master
logical,INTENT(IN),OPTIONAL             :: verbose

real(kind=sgl),allocatable              :: sr(:,:,:) 
real(kind=sgl),allocatable              :: EkeVs(:) 
integer(kind=irg),allocatable           :: atomtype(:)

real(kind=sgl),allocatable              :: srtmp(:,:,:,:)
integer(kind=irg)                       :: istat

logical                                 :: stat, readonly
integer(kind=irg)                       :: hdferr, nlines
integer(HSIZE_T)                        :: dims(1), dims4(4)
character(fnlen)                        :: groupname, dataset
character(fnlen),allocatable            :: stringarray(:)

type(HDFobjectStackType),pointer        :: HDF_head

! open the fortran HDF interface
call h5open_f(hdferr)

nullify(HDF_head, HDF_head)

! first, we need to check whether or not the input file is of the HDF5 forat type; if
! it is, we read it accordingly, otherwise we use the old binary format.
!
call h5fis_hdf5_f(trim(enl%masterfile), stat, hdferr)

if (stat) then 
! open the master file 
  readonly = .TRUE.
  hdferr =  HDF_openFile(enl%masterfile, HDF_head, readonly)

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
  if (enl%numEbins.ne.enl%nE) then
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

  dataset = 'sr'
  srtmp = HDF_readDatasetFloatArray4D(dataset, dims4, HDF_head)
  allocate(master%sr(-enl%npx:enl%npx,-enl%npy:enl%npy,enl%nE),stat=istat)
  master%sr = sum(srtmp,4)
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
!====================================
! ----- Read energy-dispersed Lambert projections (master pattern)
! this has been updated on 3/26/14 to accommodate the new EBSDmaster file format
!====================================
  open(unit=dataunit,file=trim(enl%masterfile),status='old',form='unformatted')
  read (dataunit) enl%Masterprogname
! read the version number
  read (dataunit) enl%Masterscversion
! then the name of the crystal data file
  read (dataunit) enl%Masterxtalname
! then the name of the corresponding Monte Carlo data file
  read (dataunit) enl%Masterenergyfile
! energy information and array size    
  read (dataunit) enl%npx,enl%npy,enl%nE,enl%numset
! make sure that MC and Master results are compatible
  if (enl%numEbins.ne.enl%nE) then
    call Message('Energy histogram and Lambert stack have different energy dimension; aborting program', frm = "(A)")
!   write (*,*) 'energy histogram = ',shape(accum_e)
!   write (*,*) 'Lambert stack = ', nE, npx, npy
    stop
  end if
  allocate(master%sr(-enl%npx:enl%npx,-enl%npy:enl%npy,enl%nE),srtmp(-enl%npx:enl%npx,-enl%npy:enl%npy,enl%nE,enl%numset), &
           EkeVs(enl%nE),atomtype(enl%numset),stat=istat)
  read (dataunit) EkeVs
  read (dataunit) atomtype
  deallocate(EkeVs, atomtype)   ! arrays are only needed by IDL visualization routine
! is this a regular (square) or hexagonal projection ?
  read (dataunit) enl%sqorhe
! and finally the results array
  read (dataunit) srtmp
! convert to a smaller array by summing over all atom types 
! [in a later version of the program we might allow for the 
! user to request an element specific EBSD pattern calculation]
  master%sr = sum(srtmp,4)
  deallocate(srtmp)
  close(unit=dataunit,status='keep')
end if
!====================================

if (present(verbose)) call Message(' -> completed reading '//trim(enl%masterfile), frm = "(A)")

end subroutine EBSDreadMasterfile




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
real(kind=sgl)                          :: alp, ca, sa
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

! we will need to incorporate a series of possible distortions 
! here as well, as described in Gert nolze's paper; for now we 
! just leave this place holder comment instead

! compute auxilliary interpolation arrays
! if (istat.ne.0) then ...

L2 = enl%L * enl%L
do j=1,enl%numsy
  sx = L2 + scin_y(j) * scin_y(j)
  Ls = ca * scin_y(j) + enl%L*sa
  Lc = -sa * scin_y(j) + enl%L*ca
  do i=1,enl%numsx
   rhos = 1.0/sqrt(sx + scin_x(i)**2)
   master%rgx(i,j) = Ls * rhos
   master%rgy(i,j) = scin_x(i) * rhos
   master%rgz(i,j) = Lc * rhos
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
  scl = float(enl%nsx) / LPs%sPio2

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

real(kind=dbl),allocatable                          :: master_twin(:,:,:)
type(EBSDLargeAccumType),pointer                    :: acc
logical                                             :: verbose
real(kind=dbl)                                      :: q(4),Lamproj(2),dc(3),dc_new(3),dx,dy,dxm,dym,ixy(2),scl
integer(kind=irg)                                   :: nix,niy,nixp,niyp
integer(kind=irg)                                   :: ii,jj,kk,ierr,istat,pp,qq

allocate(master_twin(-enl%npx:enl%npx,-enl%npy:enl%npy,1:enl%nE),stat=istat)

q = (/ dsqrt(3.D0)/2.D0,1/dsqrt(3.D0)/2.D0,1/dsqrt(3.D0)/2.D0,1/dsqrt(3.D0)/2.D0 /)
scl = float(enl%npx) / LPs%sPio2

do ii = 1,enl%nE
    master_twin = 0.0
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

            master_twin(jj,kk,ii) = master%sr(nix,niy,ii)*dxm*dym + master%sr(nixp,niy,ii)*dx*dym + &
                                    master%sr(nix,niyp,ii)*dxm*dy + master%sr(nixp,niyp,ii)*dx*dy
        end do
    end do
end do
master%sr = 0.5D0 * (master_twin + master%sr)

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
!--------------------------------------------------------------------------
subroutine OverlapMasterPattern(enl,master,q)

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
real(kind=dbl),INTENT(IN)                           :: q(4)

real(kind=dbl),allocatable                          :: master_rotated(:,:,:)
type(EBSDLargeAccumType),pointer                    :: acc
logical                                             :: verbose
real(kind=dbl)                                      :: Lamproj(2),dc(3),dc_new(3),dx,dy,dxm,dym,ixy(2),scl
integer(kind=irg)                                   :: nix,niy,nxip,niyp
integer(kind=irg)                                   :: ii,jj,kk,ierr,istat,pp,qq

allocate(master_rotated(-enl%npx:enl%npx,-enl%npy:enl%npy,1:enl%nE),stat=istat)

scl = float(enl%npx) / LPs%sPio2

do ii = 1,enl%nE
    master_rotated = 0.0
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

            master_rotated(jj,kk,ii) = master%sr(nix,niy,ii)*dxm*dym + master%sr(nixp,niy,ii)*dx*dym + &
                                       master%sr(nix,niyp,ii)*dxm*dy + master%sr(nixp,niyp,ii)*dx*dy
        end do
    end do
end do
master%sr = 0.5D0 * (master_rotated + master%sr)

call Message(' -> completed superimposing rotated and regular master patterns', frm = "(A)")

end subroutine OverlapMasterPattern


end module EBSDmod

