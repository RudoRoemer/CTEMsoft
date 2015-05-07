! ###################################################################
! Copyright (c) 2013-2015, Marc De Graef/Carnegie Mellon University
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
! EMsoft:EMEBSDoverlap.f90
!--------------------------------------------------------------------------
!
! PROGRAM: EMEBSDoverlap
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief EMEBSDoverlap merges master pattern data for a given orientation relation
!>
!> @note This program takes a pair of master pattern files, and 
!> merges them together for a given orientation relation between the two phases; they
!> could be the same phase for twins, for instance, or two different phases, for instance
!> a mixture of hcp and fcc with a given OR.  There are some conditions on the array sizes
!> that have to be met in order for this to work. The output of this program should only
!> be used for visualization purposes, not for actual EBSD pattern simulations; the output
!> file is HDF5, with only the overlap master pattern in it for the highest electron energy.
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
!> @date  05/03/15  MDG 6.0 branched EMEBSD.f90 into this program; tested square-square case; others to be implemented
!> @date  05/04/15  MDG 6.1 added fracA volume fraction parameter; added stereographic projection and circular Lambert output
! ###################################################################

program EMEBSDoverlap

use local
use files
use crystal
use NameListTypedefs
use NameListHandlers
use io
use EBSDmod
use error
use gvectors
use symmetry
use rotations
use quaternions
use constants
use Lambert
use HDFsupport
use initializers


IMPLICIT NONE

character(fnlen)                       :: nmldeffile, progname, progdesc, dataset
type(EBSDoverlapNameListType)          :: enl

type(EBSDLargeAccumType),pointer       :: acc
type(EBSDMasterType),pointer           :: masterA, masterB

integer(kind=irg)                      :: istat, sA(3), sB(3), ierr, i, j , hdferr
logical                                :: verbose
character(6)                           :: sqorheA, sqorheB
character(fnlen)                       :: outstr, datafile
type(unitcell),pointer                 :: cellA, cellB
type(DynType),save                     :: DynA, DynB
type(gnode),save                       :: rlpA, rlpB
real(kind=sgl)                         :: dmin, voltage, TTAB(3,3), TT(3,3), io_real(3), fracB, scl
real(kind=dbl)                         :: edge, xy(2), xyz(3), txyz(3), txy(2), Radius, dc(3)
type(orientation)                      :: orel    
real(kind=sgl),allocatable             :: master(:,:), masterLC(:,:), masterSP(:,:)
type(HDFobjectStackType),pointer       :: HDF_head

interface
  function InterpolateMaster(dc, master, s, sqorhe) result(res)
  
  use local
  use Lambert
  use EBSDmod
  use constants
  
  IMPLICIT NONE
  
  real(kind=dbl),INTENT(INOUT)            :: dc(3)
  type(EBSDMasterType),INTENT(IN),pointer :: master
  integer(kind=irg),INTENT(IN)            :: s(3)
  character(6),INTENT(IN)                 :: sqorhe
  real(kind=sgl)                          :: res
  
  end function InterpolateMaster

  function InterpolateLambert(dc, master, npx) result(res)

  use local
  use Lambert
  use EBSDmod
  use constants
  
  IMPLICIT NONE
  
  real(kind=dbl),INTENT(INOUT)            :: dc(3)
  real(kind=sgl),INTENT(IN)               :: master(-npx:npx,-npx:npx)
  integer(kind=irg),INTENT(IN)            :: npx 
  real(kind=sgl)                          :: res
  end function InterpolateLambert

end interface



verbose = .FALSE.

nmldeffile = 'EMEBSDoverlap.nml'
progname = 'EMEBSDoverlap.f90'
progdesc = 'Merge EBSD master patterns for a particular orientation relation'

! print some information
call EMsoft(progname, progdesc)

! deal with the command line arguments, if any
call Interpret_Program_Arguments(nmldeffile,1,(/ 23 /), progname)

! deal with the namelist stuff
call GetEBSDoverlapNameList(nmldeffile,enl)

! read EBSD master pattern files 
allocate(masterA, masterB)
call EBSDreadMasterfile_overlap(enl, masterA, enl%masterfileA, verbose)
sqorheA = enl%sqorhe
enl%xtalnameA = enl%Masterxtalname
call EBSDreadMasterfile_overlap(enl, masterB, enl%masterfileB, verbose)
sqorheB = enl%sqorhe
enl%xtalnameB = enl%Masterxtalname

! make sure that the master pattern arrays have the same dimensions
sA = shape(masterA%sr)
sB = shape(masterB%sr)

if (sum(abs(sA-sB)).ne.0) then
  call FatalError('EMEBSDoverlap','master patterns have different dimensions')
end if

!=============================
!=============================
! ok, we're in business... let's initialize the crystal structures so that we
! can compute the orientation relation matrix
nullify(cellA,cellB)
allocate(cellA,cellB)
dmin=0.05
voltage = 30000.0

call Initialize_Cell(cellA,DynA,rlpA,enl%xtalnameA, dmin, voltage, verbose)
call Initialize_Cell(cellB,DynB,rlpB,enl%xtalnameB, dmin, voltage, verbose)

!=============================
!=============================
orel%gA = enl%gA
orel%gB = enl%gB
orel%tA = enl%tA
orel%tB = enl%tB

! check the OR for orthogonality
if (sum(orel%gA*orel%tA).ne.0.0) then
  call FatalError('EMEBSDoverlap','gA and tA must be orthogonal !!!')
end if

if (sum(orel%gB*orel%tB).ne.0.0) then
  call FatalError('EMEBSDoverlap','gB and tB must be orthogonal !!!')
end if

! compute the rotation matrix for this OR
TTAB = ComputeOR(orel, cellA, cellB, 'AB')
TT = transpose(matmul( cellA%rsm, matmul( TTAB, transpose(cellB%dsm))))

! output
outstr = ' '//trim(enl%xtalnameA)//' --> '//trim(enl%xtalnameB)
call WriteValue('Transformation Matrix : ',outstr)
do i=1,3
  io_real(1:3) = TT(i,1:3)
  call WriteValue('',io_real,3,"(3f10.6)")
end do

! we no longer need the crystal structure data
 deallocate(cellA, cellB)

!=============================
!=============================
! next, allocate a new master array into which we'll write the superimposed patterns
! we discard all energies except for the highest energy
allocate(master(-enl%npx:enl%npx,-enl%npy:enl%npy))

! the only difficulty here is that one or both of the structures could be hexagonal/trigonal,
! and the output master pattern will be a square Lambert projection, regardless of the input
! formats; we'll need to distinguish between all possible cases...
fracB = 1.0-enl%fracA

call WriteValue('','Each master pattern has its own intensity range.',"(/A)")
call WriteValue('','This means that one pattern may dominate over another')
call WriteValue('','even when the volume fractions of A and B are equal. ',"(A/)")
io_real(1) = maxval(masterA%sr)
call WriteValue('maximum intensity in master A: ',io_real, 1)
io_real(1) = maxval(masterB%sr)
call WriteValue('maximum intensity in master B: ',io_real, 1)

!=============================
!=============================
! A is square, B is either square or hexagonal
if (sqorheA.eq.'square') then 
  edge = LPs%sPio2 / dble(enl%npx)
  do i=-enl%npx,enl%npx
    do j=-enl%npy,enl%npy
! determine the spherical direction for this point
      xy = (/ dble(i), dble(j) /) * edge
      xyz = LambertSquareToSphere(xy, ierr)
! since A is already square Lambert, all we need to do is compute the 
! beam orientation in crystal B, and sample the master pattern for that
! location. 
      txyz = matmul(TT, xyz)
! normalize these direction cosines (they are already in a cartesian reference frame!)
      txyz = txyz/sqrt(sum(txyz*txyz))
! and interpolate the masterB pattern
      master(i,j) = enl%fracA*masterA%sr(i,j,sA(3)) + fracB*InterpolateMaster(txyz, masterB, sB, sqorheB)
    end do
  end do
end if

!=============================
!=============================
! A is hexagonal, B is square [to be implemented]


!=============================
!=============================
! both A and B are hexagonal patterns [to be implemented]


! free up the master pattern memory
deallocate(masterA, masterB)


!=============================
!=============================
! convert the square Lambert projection to a stereographic projection
! with the PatternAxis of structure A at the center
  allocate(masterSP(-enl%npx:enl%npx,-enl%npy:enl%npy))
  Radius = 1.0
  do i=-enl%npx,enl%npx 
    do j=-enl%npy,enl%npy
      xy = (/ float(i), float(j) /) / float(enl%npx)
      xyz = StereoGraphicInverse( xy, ierr, Radius )
      if (ierr.ne.0) then 
        masterSP(i,j) = 0.0
      else
        masterSP(i,j) = InterpolateLambert(xyz, master, enl%npx)
      end if
    end do
  end do

! convert the square Lambert projection to a circular Lambert projection
! with the PatternAxis of structure A at the center
  allocate(masterLC(-enl%npx:enl%npx,-enl%npy:enl%npy))
  Radius = 1.0
  do i=-enl%npx,enl%npx 
    do j=-enl%npy,enl%npy
      xy = sqrt(2.0) * (/ float(i), float(j) /) / float(enl%npx)
      if (sum(xy*xy).gt.2.0) then 
        masterLC(i,j) = 0.0
      else
        xyz = LambertInverse( xy, ierr, Radius )
        masterLC(i,j) = InterpolateLambert(xyz, master, enl%npx)
      end if
    end do
  end do


! finally, create simple HDF5 file with only the overlap master array in it
nullify(HDF_head)
! Initialize FORTRAN interface.
call h5open_f(hdferr)

! Create a new file using the default properties.
datafile = trim(EMdatapathname)//trim(enl%datafile)
hdferr =  HDF_createFile(datafile, HDF_head)

! create datasets 
dataset = 'MasterLambertSquare'
hdferr = HDF_writeDatasetFloatArray2D(dataset, master, 2*enl%npx+1, 2*enl%npx+1, HDF_head)
 
dataset = 'MasterLambertCircle'
hdferr = HDF_writeDatasetFloatArray2D(dataset, masterLC, 2*enl%npx+1, 2*enl%npx+1, HDF_head)
 
dataset = 'MasterStereographic'
hdferr = HDF_writeDatasetFloatArray2D(dataset, masterSP, 2*enl%npx+1, 2*enl%npx+1, HDF_head)
 
call HDF_pop(HDF_head,.TRUE.)

! and close the fortran hdf interface
call h5close_f(hdferr)

call WriteValue('','Output data stored in '//trim(datafile),"(//A/)")

end program EMEBSDoverlap




function InterpolateLambert(dc, master, npx) result(res)

use local
use Lambert
use EBSDmod
use constants

IMPLICIT NONE

real(kind=dbl),INTENT(INOUT)            :: dc(3)
real(kind=sgl),INTENT(IN)               :: master(-npx:npx,-npx:npx)
integer(kind=irg),INTENT(IN)            :: npx 
real(kind=sgl)                          :: res

integer(kind=irg)                       :: nix, niy, nixp, niyp, istat
real(kind=sgl)                          :: xy(2), dx, dy, dxm, dym, scl

scl = float(npx) / LPs%sPio2

if (dc(3).lt.0.0) dc = -dc

! convert direction cosines to lambert projections
xy = scl * LambertSphereToSquare( dc, istat )
res = 0.0

if (istat.eq.0) then 
! interpolate intensity from the neighboring points
  nix = floor(xy(1))
  niy = floor(xy(2))
  nixp = nix+1
  niyp = niy+1
  if (nixp.gt.npx) nixp = nix
  if (niyp.gt.npx) niyp = niy
  dx = xy(1) - nix
  dy = xy(2) - niy
  dxm = 1.0 - dx
  dym = 1.0 - dy
  
  res = master(nix,niy)*dxm*dym + master(nixp,niy)*dx*dym + &
        master(nix,niyp)*dxm*dy + master(nixp,niyp)*dx*dy
end if

end function InterpolateLambert


function InterpolateMaster(dc, master, s, sqorhe) result(res)

use local
use Lambert
use EBSDmod
use constants

IMPLICIT NONE

real(kind=dbl),INTENT(INOUT)            :: dc(3)
type(EBSDMasterType),INTENT(IN),pointer :: master
integer(kind=irg),INTENT(IN)            :: s(3)
character(6),INTENT(IN)                 :: sqorhe
real(kind=sgl)                          :: res

integer(kind=irg)                       :: nix, niy, nixp, niyp, istat, npx
real(kind=sgl)                          :: xy(2), dx, dy, dxm, dym, scl, tmp

npx = (s(1)-1)/2
if (dc(3).lt.0.0) dc = -dc

! convert direction cosines to lambert projections
if (sqorhe.eq.'square') then 
  scl = float(npx) / LPs%sPio2
  xy = scl * LambertSphereToSquare( dc, istat )
else
  scl = float(npx)/LPs%preg
  xy = scl * LambertSphereToHex( dc, istat )
  tmp = xy(1)
  xy(1) = xy(2)
  xy(2) = tmp
end if
res = 0.0

if (istat.eq.0) then 
! interpolate intensity from the neighboring points
  nix = floor(xy(1))
  niy = floor(xy(2))
  nixp = nix+1
  niyp = niy+1
  if (nixp.gt.npx) nixp = nix
  if (niyp.gt.npx) niyp = niy
  dx = xy(1) - nix
  dy = xy(2) - niy
  dxm = 1.0 - dx
  dym = 1.0 - dy
  
  res = master%sr(nix,niy,s(3))*dxm*dym + master%sr(nixp,niy,s(3))*dx*dym + &
        master%sr(nix,niyp,s(3))*dxm*dy + master%sr(nixp,niyp,s(3))*dx*dy
end if

end function InterpolateMaster
