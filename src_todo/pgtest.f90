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
! EMsoft:pgtest.f90
!--------------------------------------------------------------------------
!
! PROGRAM: pgtest
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief test the Lambert projections for all point groups
!
!> @details This program runs through all the different point groups and
!> their settings to make sure that the k-space sampling algorithm in 
!> the kvectors module is correct for all possible symmetries.
! 
!> @date   08/26/15 MDG 1.0 original
!--------------------------------------------------------------------------
program pgtest

use local
use symmetry
use crystal
use kvectors
use gvectors
use io
use diffraction
use HDF5
use HDFsupport
use Lambert

IMPLICIT NONE

type(unitcell),pointer          :: cell
character(fnlen)                :: outname, groupname, dataset
character(fnlen)                :: progname, progdesc
 
type(HDFobjectStackType),pointer:: HDF_head
character(11)                   :: dstr
character(15)                   :: tstrb
character(15)                   :: tstre
character(2)                    :: groupstring
integer(kind=irg)               :: pgnums(42), isymval(42), csvals(42), sgnums(42), sgset(42), istrig(42)
integer(kind=irg)               :: npx, npy, npyhex, skip, numk, ijmax, io_int(3), hdferr, i, j, l, nx, ny, SamplingType, &
                                   ix, iy, ierr, jj, nequiv, iequiv(3,48), jjj
real(kind=dbl)                  :: aval(42), bval(42), cval(42), alval(42), beval(42), gaval(42), xyz(3), xy(2), delta, srt, &
                                   xyz2(3), xx(6), yy(6), x, y, kstar(3)
integer(kind=irg),parameter     :: cnums = 42! number of tests to carry out
type(kvectorlist),pointer       :: khead, ktmp, ktmp2
logical                         :: usehex
integer(kind=irg),allocatable   :: mLPNH(:,:), mLPSH(:,:), spNH(:,:), spSH(:,:)
type(gnode),save                :: rlp

real(kind=dbl)                  :: stmp(48,3)           !< output array with equivalent vectors
integer(kind=irg)               :: n                    !< number of entries in equivalent vector array
character(1)                    :: space                !< 'd' or 'r'


progname = 'pgtest.f90'
progdesc = 'Point group k-space sampling tests'

call timestamp(datestring=dstr, timestring=tstrb)

! define the cell pointer and assign default values
nullify(cell)
allocate(cell)
call ResetCell(cell)

!  test the hexagonal case

cell%a = 0.4D0
cell%b = 0.4D0
cell%c = 0.6D0
cell%alpha = 90.0D0
cell%beta  = 90.0D0
cell%gamma = 120.0D0
cell%xtal_system = 4
cell%SYM_SGnum = 162
cell%hexset = .TRUE.
usehex = .TRUE.
cell%SG%SYM_trigonal = .FALSE.
call CalcMatrices(cell)
!
!
npx = 200
delta = 1.D0 / dble(npx)
srt = 0.5D0*dsqrt(3.D0)
cell%voltage = 20000.D0
call CalcWaveLength(cell, rlp, skip)
write (*,*) 'wave number : ', 1.D0/cell%mLambda

xx = dble( (/ npx, npx, 0, -npx, -npx, 0 /) )
yy = dble( (/ 0, npx, npx, 0, -npx, -npx /) ) 

do i=1,6
  xy = (/ xx(i), yy(i) /) * delta
  xyz = LambertHexToSphere(xy,ierr)
  write (*,*) xy, ' --> ',xyz

  call NormVec(cell,xyz,'c')                       ! normalize incident direction in cartesian space
  xyz = xyz/cell%mLambda                         ! divide by wavelength
! and transform to reciprocal crystal space using the direct structure matrix
  call TransSpace(cell, xyz, xyz2, 'c', 'r')
  write (*,*) 'xyz2 = ',xyz2, CalcLength(cell,xyz2,'r')

! then do the stereographic projections
  call TransSpace(cell,xyz2,xyz,'r','c')
  call NormVec(cell, xyz, 'c')
  write (*,*) 'normalized vector = ', xyz, sum(xyz*xyz)
  xy = LambertSphereToHex(xyz,ierr)*dble(npx)
  write (*,*) 'inverse hexmap : ',xy
     
  x = xyz(1) + 0.5D0*xyz(2)/srt
  y = xyz(2)/srt
  ix = nint(float(npx)*x/(1.D0+xyz(3)))
  iy = nint(float(npx)*y/(1.D0+xyz(3)))
 write (*,*) 'NH : ',ix, iy, xx(i)-yy(i)*0.5D0, yy(i)*srt 
write (*,*) '---'

end do


! initialize the necessary parameters for all 42 tests
!point group numbers
pgnums = (/ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 14, 15, &
           16, 16, 17, 17, 18, 18, 18, 19, 19, 19, 20, 20, 20, 21, 22, 23, &
           24, 25, 26, 26, 27, 28, 29, 30, 31, 32 /)

! correct SamplingType numbers (must be compared to the returned numbers)
isymval = (/ 1,  2,  3,  4,  5,  5,  5,  6,  5,  5,  6,  6,  7,  8,  6,  9, &
            10, 11, 12, 13, 12, 12, 13, 14, 14, 15, 16, 16, 17, 15, 12, 17, &
            16, 18, 16, 17, 19,  3,  6,  6,  8,  9 /)

! trigonal case? if 1 then rhombohedral setting, if 2 then hexagonal setting
istrig = (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            2, 1, 2, 1, 2, 1, 2, 2, 1, 2, 2, 1, 2, 0, 0, 0, &
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)

! crystal system values for all the tests
!  1. Cubic
!  2. Tetragonal
!  3. Orthorhombic
!  4. Hexagonal
!  5. Trigonal
!  6. Monoclinic
!  7. Triclinic
csvals = (/ 7, 7, 6, 6, 6, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, &
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, &
            4, 4, 4, 4, 4, 1, 1, 1, 1, 1 /)

! define some nominal lattice parameters for each case
aval = (/ (0.4D0,i=1,42) /)
bval = (/ 0.3D0, 0.3D0, 0.3D0, 0.3D0, 0.3D0, 0.3D0, 0.3D0, 0.3D0,  &
          0.4D0, 0.4D0, 0.4D0, 0.4D0, 0.4D0, 0.4D0, 0.4D0, 0.4D0,  &
          0.4D0, 0.4D0, 0.4D0, 0.4D0, 0.4D0, 0.4D0, 0.4D0, 0.4D0,  &
          0.4D0, 0.4D0, 0.4D0, 0.4D0, 0.4D0, 0.4D0, 0.4D0, 0.4D0,  &
          0.4D0, 0.4D0, 0.4D0, 0.4D0, 0.4D0, 0.4D0, 0.4D0, 0.4D0, 0.4D0, 0.4D0 /)
cval = (/ 0.5D0, 0.5D0, 0.5D0, 0.5D0, 0.5D0, 0.5D0, 0.5D0, 0.5D0,  &
          0.5D0, 0.5D0, 0.5D0, 0.5D0, 0.5D0, 0.5D0, 0.5D0, 0.5D0,  &
          0.5D0, 0.4D0, 0.5D0, 0.4D0, 0.5D0, 0.4D0, 0.5D0, 0.5D0,  &
          0.4D0, 0.5D0, 0.5D0, 0.4D0, 0.5D0, 0.5D0, 0.5D0, 0.5D0,  &
          0.5D0, 0.5D0, 0.5D0, 0.5D0, 0.5D0, 0.4D0, 0.4D0, 0.4D0, 0.4D0, 0.4D0 /)
alval = (/70.D0, 70.D0, 90.D0, 90.D0, 90.D0, 90.D0, 90.D0, 90.D0,  &
          90.D0, 90.D0, 90.D0, 90.D0, 90.D0, 90.D0, 90.D0, 90.D0,  &
          90.D0, 50.D0, 90.D0, 50.D0, 90.D0, 50.D0, 90.D0, 90.D0,  &
          50.D0, 90.D0, 90.D0, 50.D0, 90.D0, 90.D0, 90.D0, 90.D0,  &
          90.D0, 90.D0, 90.D0, 90.D0, 90.D0, 90.D0, 90.D0, 90.D0, 90.D0, 90.D0 /)
beval = (/80.D0, 80.D0, 80.D0, 80.D0, 80.D0, 90.D0, 90.D0, 90.D0,  &
          90.D0, 90.D0, 90.D0, 90.D0, 90.D0, 90.D0, 90.D0, 90.D0,  &
          90.D0, 50.D0, 90.D0, 50.D0, 90.D0, 50.D0, 90.D0, 90.D0,  &
          50.D0, 90.D0, 90.D0, 50.D0, 90.D0, 90.D0, 90.D0, 90.D0,  &
          90.D0, 90.D0, 90.D0, 90.D0, 90.D0, 90.D0, 90.D0, 90.D0, 90.D0, 90.D0 /)
gaval = (/75.D0, 75.D0, 90.D0, 90.D0, 90.D0, 90.D0, 90.D0, 90.D0,  &
          90.D0, 90.D0, 90.D0, 90.D0, 90.D0, 90.D0, 90.D0, 90.D0,  &
         120.D0, 50.D0,120.D0, 50.D0,120.D0, 50.D0,120.D0,120.D0,  &
          50.D0,120.D0,120.D0, 50.D0,120.D0,120.D0,120.D0,120.D0,  &
         120.D0,120.D0,120.D0,120.D0,120.D0, 90.D0, 90.D0, 90.D0, 90.D0, 90.D0 /)
 
! for all cases, we need to set a representative space group number
sgnums = (/   1,   2,   3,   6,  10,  16,  25,  47,  75,  81,  83,  89,  99, 111, 115, 123, &
            143, 146, 147, 148, 150, 155, 149, 156, 160, 157, 164, 166, 162, 168, 174, 175, &
            177, 183, 187, 189, 191, 195, 200, 207, 215, 221 /)

! open an HDF5 file for program output; the output will be a series of 
! stereographic projections, one pair (Northern and Southern hemispheres)
! for each of the 42 point group cases (this includes rhombohedral and 
! hexagonal settings for the trigonal groups, and special cases for tetragonal
! -42m/-4m2 and hexagonal -6m2/-62m).
!
! Initialize FORTRAN interface.
!
call h5open_f(hdferr)

! Create a new file using the default properties.
nullify(HDF_head)
outname = "pgtest.h5"
hdferr =  HDF_createFile(outname, HDF_head)

! write the EMheader to the file
call HDF_writeEMheader(HDF_head, dstr, tstrb, tstre, progname)

! create a group to write all the stereographic projections into
groupname = "PGprojections"
hdferr = HDF_createGroup(groupname, HDF_head)

! various other parameters
skip = 3
npx = 200
npy = 200
npyhex = nint(2.0*float(npy)/sqrt(3.0))
ijmax = 0 
nx = 2*npx+1



! loop through the point groups and all their settings (a total of 42 cases)
do i=1,cnums
 if (istrig(i).ne.1) then
! initialize the components of cell that will be needed for this test
  call ResetCell(cell)
  cell%a = aval(i)
  cell%b = bval(i)
  cell%c = cval(i)
  cell%alpha = alval(i)
  cell%beta = beval(i)
  cell%gamma = gaval(i)
  cell%xtal_system = csvals(i)
  cell%SYM_SGnum = sgnums(i)
  cell%hexset = .FALSE.
  cell%SG%SYM_trigonal = .FALSE.
  cell%SG%SYM_second = .FALSE.
  if (istrig(i).gt.0) cell%SG%SYM_trigonal = .TRUE.
  if (istrig(i).eq.1) cell%SG%SYM_second = .TRUE.
  usehex = .FALSE.
  if ((cell%xtal_system.eq.4).or.(cell%xtal_system.eq.5)) then
    usehex = .TRUE.
    cell%hexset = .TRUE.
  end if

  nullify(khead)

  SamplingType = PGSamplingType(pgnums(i))

! next, intercept the special cases (hexagonal vs. rhombohedral cases that require special treatment)
  if (SamplingType.eq.-1) then 
    SamplingType = getHexvsRho(cell,pgnums(i))
  end if

  write (*,*) 'Starting point group ',i,pgnums(i), csvals(i), sgnums(i), PGTHD(pgnums(i)), SamplingType, &
              cell%SYM_SGset, cell%SG%SYM_trigonal
  write (*,*) ' --> latparm = ',cell%a, cell%b, cell%c, cell%alpha, cell%beta, cell%gamma


! compute the metric matrices
  call CalcMatrices(cell)

! generate the symmetry elements (all groups in sgnum are symmorphic)
  call GenerateSymmetry(cell,.TRUE.)

! ok, so we have a unit cell and all the symmetry elements; we do not need to populate the 
! cell with atoms, so we can directly proceed to generating a list of independent k-vectors
! which we will then scan and plot on a pair of modified Lambert projections; we'll convert
! them to stereographic projections and write those to the HDF5 file with the case number 
! as a groupname.

! set the accelerating voltage
   cell%voltage = 20000.D0
   call CalcWaveLength(cell, rlp, skip)

! ---------- create the incident beam directions list
! determine all independent incident beam directions (use a linked list starting at khead)
! numk is the total number of k-vectors to be included in this computation;
! note that this needs to be redone for each energy, since the wave vector changes with energy
   nullify(khead)
!  if (usehex) then
!   call Calckvectors(khead,cell, (/ 0.D0, 0.D0, 1.D0 /), (/ 0.D0, 0.D0, 0.D0 /),0.D0,npx,npyhex,numk, &
!               SamplingType,ijmax,'RoscaLambert',usehex)
!  else 
    call Calckvectors(khead,cell, (/ 0.D0, 0.D0, 1.D0 /), (/ 0.D0, 0.D0, 0.D0 /),0.D0,npx,npy,numk, &
                SamplingType,ijmax,'RoscaLambert',usehex)
!  end if
   io_int(1)=numk
   call WriteValue('# independent beam directions to be considered = ', io_int, 1, "(I8)")

! next, fill the Northern and Southern hemisphere arrays
!  if (usehex) then 
!    allocate(mLPNH(-npx:npx,-npyhex:npyhex),mLPSH(-npx:npx,-npyhex:npyhex))
!  else
     allocate(mLPNH(-npx:npx,-npy:npy),mLPSH(-npx:npx,-npy:npy))
!  end if
   allocate(spNH(-npx:npx,-npy:npy),spSH(-npx:npx,-npy:npy))
   mLPNH = 0
   mLPSH = 0
   spNH = 0
   spSH = 0
   
write (*,*) 'starting loop over k-vectors'

ktmp => khead
! loop through the list and delete it at the same time
   do j=1,numk
! fill in the point and its symmetrically equivalent versions
     if (usehex) then 
       call Apply3DPGsymmetry(cell,ktmp%i,ktmp%j,ktmp%hs,npx,iequiv,nequiv,usehex)
     else
       if ((cell%SYM_SGnum.ge.195).and.(cell%SYM_SGnum.le.230)) then
         call Apply3DPGsymmetry(cell,ktmp%i,ktmp%j,ktmp%hs,npx,iequiv,nequiv,cubictype=SamplingType)
       else
         call Apply3DPGsymmetry(cell,ktmp%i,ktmp%j,ktmp%hs,npx,iequiv,nequiv)
       end if
     end if
     do jj=1,nequiv
       if (iequiv(3,jj).eq.-1) mLPSH(iequiv(1,jj),iequiv(2,jj)) = jj
       if (iequiv(3,jj).eq.1) mLPNH(iequiv(1,jj),iequiv(2,jj)) = jj
     end do
! stereographic projections
     if (usehex) then 
       call Apply3DPGsymmetry(cell,ktmp%i,ktmp%j,ktmp%hs,npx,iequiv,nequiv,usehex,stereographic=.TRUE.)
     else
       if ((cell%SYM_SGnum.ge.195).and.(cell%SYM_SGnum.le.230)) then
         call Apply3DPGsymmetry(cell,ktmp%i,ktmp%j,ktmp%hs,npx,iequiv,nequiv,cubictype=SamplingType,stereographic=.TRUE.)
       else
         call Apply3DPGsymmetry(cell,ktmp%i,ktmp%j,ktmp%hs,npx,iequiv,nequiv,stereographic=.TRUE.)
       end if
     end if
     do jj=1,nequiv
       if (iequiv(3,jj).eq.-1) SPSH(iequiv(1,jj),iequiv(2,jj)) = jj
       if (iequiv(3,jj).eq.1)  SPNH(iequiv(1,jj),iequiv(2,jj)) = jj
     end do

! delete the linked list entry
    ktmp2 => ktmp%next
    deallocate(ktmp)
    ktmp => ktmp2 
   end do

write (*,*) 'ready to write to file'

! write the output arrays to the HDF5 file
   write (groupstring,"(I2.2)") i
!  if (usehex) then 
!    ny = 2*npyhex+1
!    dataset = groupstring//"NH"
!    hdferr = HDF_writeDatasetIntegerArray2D(dataset, mLPNH, nx, ny, HDF_head)
!    dataset = groupstring//"SH"
!    hdferr = HDF_writeDatasetIntegerArray2D(dataset, mLPSH, nx, ny, HDF_head)
!  else 
     ny = 2*npy+1
     dataset = groupstring//"NH"
     hdferr = HDF_writeDatasetIntegerArray2D(dataset, mLPNH, nx, ny, HDF_head)
     dataset = groupstring//"SH"
     hdferr = HDF_writeDatasetIntegerArray2D(dataset, mLPSH, nx, ny, HDF_head)
!  end if
   ny = 2*npy+1
   dataset = groupstring//"NHsp"
   hdferr = HDF_writeDatasetIntegerArray2D(dataset, spNH, nx, ny, HDF_head)
   dataset = groupstring//"SHsp"
   hdferr = HDF_writeDatasetIntegerArray2D(dataset, spSH, nx, ny, HDF_head)
write (*,*) '  --> done'
! and get rid of the arrays
   deallocate(mLPNH, mLPSH, spNH, spSH)
 end if
end do


! and close the HDF5 file
call HDF_pop(HDF_head,.TRUE.)

! and close the fortran hdf interface
call h5close_f(hdferr)



end program pgtest
