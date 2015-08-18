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
! EMsoft:EMEBSDmaster.f90
!--------------------------------------------------------------------------
!
! PROGRAM: EMEBSDmaster
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief EMEBSDmaster computes the energy-dependent master EBSD pattern for a given structure
!
!> @todo implement full symmetry use; implement multiple reflection output
!> or, easier perhaps, selection of one reflection;
!>
!> implement OpenMP multithreading for the actual computation part; requires modifications
!> in EMlib.a routines (mostly THREADPRIVATE commands in several modules)
!
!> @date  03/08/12  MDG 1.0 EBSD program for fundamental zone patterns
!> @date  08/17/12  MDG 1.1 added generalized fundamental zone for all crystal symmetries
!> @date  08/20/12  MDG 1.2 modifid FZ to use the Lambert projection of a square sampling grid
!> @date  09/06/12  MDG 1.3 added support for second setting of Laue group -3m
!> @date  11/21/12  MDG 2.0 added full Lambert projection support for both square and hexagonal grids
!>				   the older code is still available for now, but will be removed after validation
!>				   of the newer code.
!> @date  12/04/12  MDG 2.1 added support for equilateral triangle mapping; needs to be validated.
!>				   also modified the structure of the output file, so that EBSD.f90 will
!>				   know which of the inverse mapping methods it should use.
!> @date  12/10/12  MDG 3.0 expanded EBSDFZ program to include energy-dependencies from Monte Carlo
!> @date  12/12/12  MDG 3.1 test to do an actual numerical integration for the I_jk integrals, using MC profiles
!> @date  08/01/13  MDG 4.0 complete rewrite with Lambert format for MC output and new ctemlib.a routines 
!>                          also, the earlier versions would do only one energy value, whereas this new 
!>                          implementation does the complete energy-dependent master pattern
!> @date  01/27/14  MDG 4.1 continued rewrite, fixed problem with kvector list, replaced gvector routines
!>		   	     with updated routines; changed program name to EMEBSDmaster.
!> @date  05/03/14  MDG 4.2 test version to resolve bug in the Sgh matrix part (solved)
!> @date  06/19/14  MDG 4.3 rewrite, removal of all globals, split of namelist handling from computation; add OpenMP
!> @date  03/26/15  MDG 5.0 all output now in HDF5 format (this is the first converted program)
!> @date  05/09/15  MDG 5.1 completed conversion from hexagonal Lambert to square Labert
!--------------------------------------------------------------------------
program EMEBSDmaster

use local
use NameListTypedefs
use NameListHandlers
use files
use io

IMPLICIT NONE

character(fnlen)                        :: nmldeffile, progname, progdesc
type(EBSDMasterNameListType)            :: emnl

nmldeffile = 'EMEBSDmaster.nml'
progname = 'EMEBSDmaster.f90'
progdesc = 'EBSD Energy-dependent Master Pattern Simulation'

! print some information
call EMsoft(progname, progdesc)

! deal with the command line arguments, if any
call Interpret_Program_Arguments(nmldeffile,1,(/ 21 /), progname)

! deal with the namelist stuff
call GetEBSDMasterNameList(nmldeffile,emnl)

! generate a set of master EBSD patterns
 call ComputeMasterPattern(emnl, progname, nmldeffile)

end program EMEBSDmaster

!--------------------------------------------------------------------------
!
! SUBROUTINE:ComputeMasterPattern
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute an EBSD master pattern as a function of energy
!
!> @param emnl namelist 
!> @param progname program name
!> @param nmldeffile namelist file name (so that the entire file can be stored inside the HDF5 file)
!
!> @date 11/29/01  MDG 1.0 original
!> @date 04/08/13  MDG 2.0 rewrite
!> @date 05/14/13  MDG 2.1 replaced IO by namelist file
!> @date 08/01/13  MDG 3.0 complete rewrite, eliminated old Lambert projection
!> @date 09/25/13  MDG 3.1 replaced k-vector code by kvectors module
!> @date 06/19/14  MDG 4.0 no more globals, nml split, added OpenMP
!> @date 03/26/15  MDG 5.0 all output now in HDF5 format (this is the first converted program)
!> @date 03/28/15  MDG 5.1 converted Monte Carlo input to HDF5 format and tested
!> @date 04/14/15  MDG 5.2 modified HDF5 output using hyperslabs
!> @date 04/15/15  MDG 5.3 debugged file closing problem, mostly in HDFsupport.f90
!> @date 04/15/15  MDG 5.4 corrected offset reading of accum_z array
!> @date 04/27/15  MDG 5.5 reactivate the hexagonal code; needs to be debugged
!> @date 05/08/15  MDG 5.6 added automated conversion from hexagonal to square Lambert; debugged
!--------------------------------------------------------------------------
subroutine ComputeMasterPattern(emnl, progname, nmldeffile)

use typedefs
use NameListTypedefs
use initializers
use MBmodule
use symmetry
use crystal
use constants
use error
use gvectors
use kvectors
use io
use local
use files
use diffraction
use multibeams
use timing
use Lambert
use HDF5
use NameListHDFwriters
use HDFsupport
use ISO_C_BINDING
use omp_lib
 
IMPLICIT NONE

type(EBSDMasterNameListType),INTENT(INOUT) :: emnl
character(fnlen),INTENT(IN)                :: progname
character(fnlen),INTENT(IN)                :: nmldeffile

real(kind=dbl)          :: ctmp(192,3), arg
integer(HSIZE_T)        :: dims4(4), cnt4(4), offset4(4)
integer(kind=irg)      :: isym,i,j,ik,npy,ipx,ipy,debug,iE,izz, izzmax, iequiv(2,12), nequiv, num_el, MCnthreads, & ! counters
                        numk, & ! number of independent incident beam directions
                        ir,nat(100),kk(3), npyhex, skip, ijmax, one, NUMTHREADS, TID, &
                        numset,n,ix,iy,iz, io_int(6), nns, nnw, nref,  &
                        istat,gzero,ic,ip,ikk, totstrong, totweak, jh, ierr, nix, niy, nixp, niyp     ! counters
real(kind=dbl)         :: tpi,Znsq, kkl, DBWF, kin, delta, h, lambda, omtl, srt, dc(3), xy(2), edge, scl, tmp, dx, dxm, dy, dym !
real(kind=sgl)          :: io_real(5), selE, kn, FN(3), kkk(3)
real(kind=sgl),allocatable      :: sr(:,:,:,:), srhex(:,:,:,:), srh(:,:,:,:), EkeVs(:), svals(:) ! results
complex(kind=dbl)               :: czero
complex(kind=dbl),allocatable   :: Lgh(:,:), Sgh(:,:,:)
logical                 :: usehex, switchmirror, verbose
character(fnlen)        :: xtalname
! the following will need to be moved elsewhere at some point...
integer(kind=irg),parameter     :: LaueTest(11) = (/ 149, 151, 153, 156, 158, 160, 161, 164, 165, 166, 167 /)  ! space groups with 2 or mirror at 30 degrees

! Monte Carlo derived quantities
integer(kind=irg)       :: numEbins, numzbins, nsx, nsy, hdferr, nlines    ! variables used in MC energy file
real(kind=dbl)          :: EkeV, Ehistmin, Ebinsize, depthmax, depthstep, etotal ! enery variables from MC program
integer(kind=irg),allocatable :: accum_e(:,:,:), accum_z(:,:,:,:), thick(:), acc_z(:,:,:,:)
real(kind=sgl),allocatable :: lambdaE(:,:)
character(fnlen)        :: oldprogname, groupname, energyfile, outname
character(8)            :: MCscversion
character(11)           :: dstr
character(15)           :: tstrb
character(15)           :: tstre
logical                 :: f_exists, readonly, overwrite=.TRUE., insert=.TRUE.
character(fnlen, KIND=c_char),allocatable,TARGET :: stringarray(:)
character(fnlen,kind=c_char)                     :: line2(1)

type(unitcell),pointer          :: cell
type(DynType),save              :: Dyn
type(gnode),save                :: rlp
type(reflisttype),pointer       :: reflist,firstw, rltmp
type(BetheParameterType)        :: BetheParameters
type(kvectorlist),pointer       :: khead, ktmp
real(kind=sgl),allocatable      :: karray(:,:)
integer(kind=irg),allocatable   :: kij(:,:)
complex(kind=dbl),allocatable   :: DynMat(:,:)
character(fnlen)                :: dataset, instring

type(HDFobjectStackType),pointer  :: HDF_head

!$OMP THREADPRIVATE(rlp) 

nullify(HDF_head)

call timestamp(datestring=dstr, timestring=tstrb)

! note that the EMMC.f90 program creates a statistical output file that 
! must be read by the present program, so that things like energy etc are 
! known to the program...  The content of the MC output file is as follows:
! write(dataunit) numEbins, numzbins, numsx, numsy
! write (dataunit) EkeV, Ehistmin, Ebinsize, depthmax, depthstep
! write(dataunit) accum_e
! write (dataunit) accum_z
! So that means that the beam voltage range is already known; we do not need to
! ask for this again...

tpi = 2.D0*cPi
czero = dcmplx(0.D0,0.D0)

!=============================================
!=============================================
! ---------- read Monte Carlo .h5 output file and extract necessary parameters
call Message('opening '//trim(emnl%energyfile), frm = "(A)" )

! this is now an HDF5 file as of version 5.1
! [since this is no longer a sequential access file, we do not need to read everything, just
! the quantities that we need...]

!
! Initialize FORTRAN interface.
!
call h5open_f(hdferr)

! first of all, if the file exists, then delete it and rewrite it on each energyloop
energyfile = trim(EMdatapathname)//trim(emnl%energyfile)
inquire(file=energyfile, exist=f_exists)

if (.not.f_exists) then
  call FatalError('ComputeMasterPattern','Monte Carlo input file does not exist')
end if

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
xtalname = trim(stringarray(1))

dataset = 'numsx'
nsx = HDF_readDatasetInteger(dataset, HDF_head)
nsx = (nsx - 1)/2
nsy = nsx

dataset = 'EkeV'
EkeV = HDF_readDatasetDouble(dataset, HDF_head)

dataset = 'Ehistmin'
Ehistmin = HDF_readDatasetDouble(dataset, HDF_head)

dataset = 'Ebinsize'
Ebinsize = HDF_readDatasetDouble(dataset, HDF_head)

dataset = 'depthmax'
depthmax = HDF_readDatasetDouble(dataset, HDF_head)

dataset = 'depthstep'
depthstep = HDF_readDatasetDouble(dataset, HDF_head)

! close the name list group
call HDF_pop(HDF_head)
call HDF_pop(HDF_head)

! open the Data group
groupname = 'EMData'
hdferr = HDF_openGroup(groupname, HDF_head)

! read data items 
dataset = 'numEbins'
numEbins = HDF_readDatasetInteger(dataset, HDF_head)

dataset = 'numzbins'
numzbins = HDF_readDatasetInteger(dataset, HDF_head)

dataset = 'totnum_el'
num_el = HDF_readDatasetInteger(dataset, HDF_head)

dataset = 'accum_z'
! dims4 =  (/ numEbins, numzbins, 2*(nsx/10)+1,2*(nsy/10)+1 /)
acc_z = HDF_readDatasetIntegerArray4D(dataset, dims4, HDF_head)
allocate(accum_z(numEbins,numzbins,-nsx/10:nsx/10,-nsy/10:nsy/10),stat=istat)
accum_z = acc_z
deallocate(acc_z)

! and close everything
call HDF_pop(HDF_head,.TRUE.)

! close the fortran interface
call h5close_f(hdferr)

! print some information to the console
io_int(1:5) = (/ numEbins, numzbins, nsx, nsy, num_el /)
call WriteValue(' NumEbins, numzbins, nsx, nsy, num_el',io_int,5,"(4I8,', ',I14)")
etotal = num_el 

io_real(1:5) = (/ EkeV, Ehistmin, Ebinsize, depthmax, depthstep /)
call WriteValue(' EkeV, Ehistmin, Ebinsize, depthmax, depthstep ',io_real,5,"(4F10.5,',',F10.5)")

call Message(' -> completed reading '//trim(emnl%energyfile), frm = "(A//)")

!=============================================
!=============================================


!=============================================
!=============================================
! crystallography section

nullify(cell)
allocate(cell)

 verbose = .TRUE.
 call Initialize_Cell(cell,Dyn,rlp,xtalname, emnl%dmin, sngl(1000.0*EkeV), verbose)

! the following line needs to be verified ... 
!  hexset = .TRUE.    ! if hexagonal structure this switch selects between three and four index notation (4 if true)

! determine the point group number
 j=0
 do i=1,32
  if (SGPG(i).le.cell%SYM_SGnum) j=i
 end do
 isym = j

! and convert this to the corresponding  Laue point group number since diffraction 
! patterns are always centrosymmetric (hence, there are only 11 different cases for
! the symmetry of the incident beam).   
  isym = PGLaueinv(isym)  

! If the Laue group is # 7, then we need to determine the orientation of the mirror plane.
! The second orientation of the mirror plane is represented by "Laue group" # 12 in this program.
 switchmirror = .FALSE.
 if (isym.eq.7) then
  do i=1,11
    if (cell%SYM_SGnum.eq.LaueTest(i)) switchmirror = .TRUE.
  end do
 end if
 if (switchmirror) then
  isym = 12
  call Message(' Switching computational wedge to second setting for this space group', frm = "(A)")
 end if
 write (*,*) ' Laue group # ',isym, PGTHD(j)

! if this point group is trigonal or hexagonal, we need to switch usehex to .TRUE. so that
! the program will use the hexagonal sampling method
usehex = .FALSE.
if (((isym.ge.6).and.(isym.le.9)).or.(isym.eq.12)) usehex = .TRUE.
! ---------- end of symmetry and crystallography section
!=============================================
!=============================================

!=============================================
!=============================================
! this is where we determine the value for the thickness integration limit for the CalcLgh3 routine...
allocate(EkeVs(numEbins),thick(numEbins))

do i=1,numEbins
  EkeVs(i) = Ehistmin + float(i-1)*Ebinsize
end do

! then, for each energy determine the 95% histogram thickness
izzmax = 0
do iE = 1,numEbins
 do ix=-nsx/10,nsx/10
  do iy=-nsy/10,nsy/10
   istat = sum(accum_z(iE,:,ix,iy))
   izz = 1
   do while (sum(accum_z(iE,1:izz,ix,iy)).lt.(0.99*istat)) 
    izz = izz+1
   end do
   if (izz.gt.izzmax) izzmax = izz
  end do
 end do
 thick(iE) = dble(izzmax) * depthstep
end do

izz = nint(maxval(thick)/depthstep)
allocate(lambdaE(1:numEbins,1:izz),stat=istat)
do iE=1,numEbins
 do iz=1,izz
  lambdaE(iE,iz) = float(sum(accum_z(iE,iz,-nsx/10:nsx/10,-nsy/10:nsy/10)))/etotal
 end do
end do

! and get rid of the accum_z array
deallocate(accum_z)
! ---------- end of 'read Monte Carlo output file and extract necessary parameters' section
!=============================================
!=============================================


!=============================================
!=============================================
! ---------- a couple of initializations
   numset = cell % ATOM_ntype  
   npy = emnl%npx
   allocate(svals(numset),stat=istat)
   gzero = 1  ! index of incident beam
   debug = 0  ! no longer used
! ----------
!=============================================
!=============================================

!=============================================
!=============================================
! ---------- allocate memory for the master pattern
! we need to sample the stereographic projection Northern hemisphere or a portion
! thereoff, depending on the order of the Laue group.  There are 11 Laue groups, 
! which leads to 9 different shapes for the stereographic asymmetric unit for the 
! independent incident beam directions.  
! allocate space for the results (needs to be altered for general symmetry case)
!if (emnl%Esel.eq.-1) then
!  allocate(sr(-emnl%npx:emnl%npx,-npy:npy,1:numEbins,1:numset),stat=istat)
!else
  allocate(sr(-emnl%npx:emnl%npx,-npy:npy,1,1:numset),stat=istat)
!end if 

! in the trigonal/hexagonal case, we need intermediate storage arrays
  if (usehex) then
   npyhex = nint(2.0*float(npy)/sqrt(3.0))
   allocate(srhex(-emnl%npx:emnl%npx,-npyhex:npyhex,1,1:numset),stat=istat)
   allocate(srh(-emnl%npx:emnl%npx,-npy:npy,1,1:numset),stat=istat)
   allocate(sr(-emnl%npx:emnl%npx,-npy:npy,1,1:numset),stat=istat)
  end if

! set various arrays to zero
   sr = 0.0
   if (usehex) then
     srhex = 0.0
   end if
! ---------- end allocate memory for the master pattern
!=============================================
!=============================================

! force dynamical matrix routine to read new Bethe parameters from file
! this will all be changed with the new version of the Bethe potentials
!  call Set_Bethe_Parameters(BetheParameters)


!=============================================
! create the HDF5 output file
!=============================================
! delete the file if it already exists
! inquire(file=trim(emnl%outname), exist=f_exists)

! if (f_exists) then
!   open(unit=dataunit, file=trim(emnl%outname), status='old',form='unformatted')
!   close(unit=dataunit, status='delete')
!   call Message('Old HDF file deleted')
! end if

  nullify(HDF_head)
! Initialize FORTRAN interface.
  call h5open_f(hdferr)

! Create a new file using the default properties.
  outname = trim(EMdatapathname)//trim(emnl%outname)
  hdferr =  HDF_createFile(outname, HDF_head)

! write the EMheader to the file
  call HDF_writeEMheader(HDF_head, dstr, tstrb, tstre, progname)

! create a namelist group to write all the namelist files into
  groupname = "NMLfiles"
  hdferr = HDF_createGroup(groupname, HDF_head)

! read the text file and write the array to the file
  dataset = 'EBSDmasterNML'
  hdferr = HDF_writeDatasetTextFile(dataset, nmldeffile, HDF_head)

! leave this group
  call HDF_pop(HDF_head)
  
! create a namelist group to write all the namelist files into
  groupname = "NMLparameters"
  hdferr = HDF_createGroup(groupname, HDF_head)
  call HDFwriteEBSDMasterNameList(HDF_head, emnl)

! leave this group
  call HDF_pop(HDF_head)

! then the remainder of the data in a EMData group
  groupname = 'EMData'
  hdferr = HDF_createGroup(groupname, HDF_head)

  dataset = 'xtalname'
  stringarray(1)= trim(xtalname)
  hdferr = HDF_writeDatasetStringArray(dataset, stringarray, 1, HDF_head)

  dataset = 'numset'
  hdferr = HDF_writeDatasetInteger(dataset, numset, HDF_head)

  if (emnl%Esel.eq.-1) then
    dataset = 'numEbins'
    hdferr = HDF_writeDatasetInteger(dataset, numEbins, HDF_head)
  
    dataset = 'EkeVs'
    hdferr = HDF_writeDatasetFloatArray1D(dataset, EkeVs, numEbins, HDF_head)
  else
    dataset = 'numEbins'
    hdferr = HDF_writeDatasetInteger(dataset, one, HDF_head)
  
    dataset = 'selE'
    hdferr = HDF_writeDatasetFloat(dataset, selE, HDF_head)
  end if 

  dataset = 'cell%ATOM_type'
  hdferr = HDF_writeDatasetIntegerArray1D(dataset, cell%ATOM_type(1:numset), numset, HDF_head)

  dataset = 'squhex'
  if (usehex) then 
    stringarray(1)= 'hexago'
    hdferr = HDF_writeDatasetStringArray(dataset, stringarray, 1, HDF_head)
  else
    stringarray(1)= 'square'
    hdferr = HDF_writeDatasetStringArray(dataset, stringarray, 1, HDF_head)
  end if  

! create the hyperslab and write zeroes to it for now
  dataset = 'sr'
  dims4 = (/  2*emnl%npx+1, 2*emnl%npx+1, numEbins, numset /)
  cnt4 = (/ 2*emnl%npx+1, 2*emnl%npx+1, 1, numset /)
  offset4 = (/ 0, 0, 0, 0 /)
  hdferr = HDF_writeHyperslabFloatArray4D(dataset, sr, dims4, offset4, cnt4(1), cnt4(2), cnt4(3), cnt4(4), HDF_head)

! create the hexagonal hyperslab hyperslab and write zeroes to it for now
  if (usehex) then 
    dataset = 'srhex'
    dims4 = (/  2*emnl%npx+1, 2*emnl%npx+1, numEbins, numset /)
    cnt4 = (/ 2*emnl%npx+1, 2*emnl%npx+1, 1, numset /)
    offset4 = (/ 0, 0, 0, 0 /)
    hdferr = HDF_writeHyperslabFloatArray4D(dataset, srh, dims4, offset4, cnt4(1), cnt4(2), cnt4(3), cnt4(4), HDF_head)
  end if

  call HDF_pop(HDF_head,.TRUE.)

! and close the fortran hdf interface
 call h5close_f(hdferr)

!=============================================
!=============================================
! ---------- from here on, we need to repeat the entire computation for each energy value
! so this is where we could in principle implement an OpenMP approach; alternatively, 
! we could do the inner loop over the incident beam directions in OpenMP (probably simpler)

energyloop: do iE=numEbins,1,-1
! is this a single-energy run ?
   if (emnl%Esel.ne.-1) then
     if (emnl%Esel.ne.iE) CYCLE energyloop
   end if
   
! print a message to indicate where we are in the computation
   io_int(1)=iE
   call Message('Starting computation for energy bin', frm = "(/A$)")
   call WriteValue(' ',io_int,1,"(I4$)")
   io_real(1) = EkeVs(iE)
   call WriteValue('; energy [keV] = ',io_real,1,"(F6.2/)")
   selE = EkeVs(iE)

! set the accelerating voltage
   skip = 3
   cell%voltage = dble(EkeVs(iE)*1000.0)
   call CalcWaveLength(cell, rlp, skip)

!=============================================
! ---------- create the incident beam directions list
! determine all independent incident beam directions (use a linked list starting at khead)
! numk is the total number of k-vectors to be included in this computation;
! note that this needs to be redone for each energy, since the wave vector changes with energy
   nullify(khead)
   if (usehex) then
    call Calckvectors(khead,cell, (/ 0.D0, 0.D0, 1.D0 /), (/ 0.D0, 0.D0, 0.D0 /),0.D0,emnl%npx,npyhex,numk, &
                isym,ijmax,'RoscaLambert',usehex)
   else 
! Calckvectors(k,ga,ktmax,npx,npy,numk,isym,ijmax,mapmode,usehex)
    call Calckvectors(khead,cell, (/ 0.D0, 0.D0, 1.D0 /), (/ 0.D0, 0.D0, 0.D0 /),0.D0,emnl%npx,npy,numk, &
                isym,ijmax,'RoscaLambert',usehex)
   end if
   io_int(1)=numk
   call WriteValue('# independent beam directions to be considered = ', io_int, 1, "(I8)")

! convert the kvector linked list into arrays for OpenMP
  allocate(karray(4,numk), kij(2,numk),stat=istat)
! point to the first beam direction
  ktmp => khead
! and loop through the list, keeping k, kn, and i,j
  karray(1:3,1) = sngl(ktmp%k(1:3))
  karray(4,1) = sngl(ktmp%kn)
  kij(1:2,1) = (/ ktmp%i, ktmp%j /)
open(unit=dataunit,file='kvecs.txt',status='unknown',form='formatted')
ik = 1
write(dataunit,"(I5,'->',2(F14.5,','),F14.5)") ik, karray(1:3,ik)
  do ik=2,numk
    ktmp => ktmp%next
    karray(1:3,ik) = sngl(ktmp%k(1:3))
    karray(4,ik) = sngl(ktmp%kn)
    kij(1:2,ik) = (/ ktmp%i, ktmp%j /)
write(dataunit,"(I5,'->',2(F14.5,','),F14.5)") ik, karray(1:3,ik)
  end do
close(unit=dataunit,status='keep')
! and remove the linked list
  call Delete_kvectorlist(khead)

  verbose = .FALSE.
  totstrong = 0
  totweak = 0

! ---------- end of "create the incident beam directions list"
!=============================================

! here's where we introduce the OpenMP calls, to spead up the overall calculations...

! set the number of OpenMP threads 
  call OMP_SET_NUM_THREADS(emnl%nthreads)
  io_int(1) = emnl%nthreads
  call WriteValue(' Attempting to set number of threads to ',io_int, 1, frm = "(I4)")

! use OpenMP to run on multiple cores ... 
!!$OMP PARALLEL default(shared) COPYIN(rlp) &
!$OMP PARALLEL COPYIN(rlp) &
!$OMP& PRIVATE(DynMat,Sgh,Lgh,ik,FN,TID,kn,ipx,ipy,ix,iequiv,nequiv,reflist,firstw) &
!$OMP& PRIVATE(kkk,nns,nnw,nref,svals,nat,io_int)

  NUMTHREADS = OMP_GET_NUM_THREADS()
  TID = OMP_GET_THREAD_NUM()

!$OMP DO SCHEDULE(DYNAMIC,100)    
! ---------- and here we start the beam direction loop
   beamloop:do ik = 1,numk

!=============================================
! ---------- create the master reflection list for this beam direction
! Then we must determine the masterlist of reflections (also a linked list);
! This list basically samples a large reciprocal space volume; it does not 
! distinguish between zero and higher order Laue zones, since that 
! distinction becomes meaningless when we consider the complete 
! reciprocal lattice.  
     nullify(reflist)
     kkk = karray(1:3,ik)
     FN = kkk

! [110] beam direction
!kkk =  (/ 44.73017,      44.73017,       0.00000 /) 
!FN = kkk
!    call TransSpace(cell,kkk,FN,'r','d')
!    call NormVec(cell,FN,'d')
!    FN = karray(1:3,ik)
!verbose = .TRUE.
     call Initialize_ReflectionList(cell, reflist, BetheParameters, FN, kkk, emnl%dmin, nref, verbose)
!write (*,*) 'original number of reflections in list ',nref
!write (*,*) 'zone axis = ',kkk



! ---------- end of "create the master reflection list"
!=============================================

! lines copied from EMKossel, may need to be modified ... 

! determine strong and weak reflections
     nullify(firstw)
     nns = 0
     nnw = 0
     call Apply_BethePotentials(cell, reflist, firstw, BetheParameters, nref, nns, nnw)

!rltmp=>reflist%next
!write (*,*) 'number of reflections in list ',nref, nns, nnw
!write(*,*) ' double diffraction for 0, 0, 2 ? ',cell%dbdiff(0,0,2),  IsGAllowed(cell, (/ 0, 0, 2 /) )
!do i=1,nref
!  write (*,*) i, rltmp%hkl, rltmp%sg, cell%LUT( rltmp%hkl(1), rltmp%hkl(2), rltmp%hkl(3) ), rltmp%strong, rltmp%weak, &
!              cell%dbdiff( rltmp%hkl(1), rltmp%hkl(2), rltmp%hkl(3) ), rltmp%hkl(1)+rltmp%hkl(2)
!  rltmp => rltmp%next
!end do
!
!
!stop


! generate the dynamical matrix
     allocate(DynMat(nns,nns))
     call GetDynMat(cell, reflist, firstw, rlp, DynMat, nns, nnw)
     totstrong = totstrong + nns
     totweak = totweak + nnw

! if (TID.eq.0) write (*,*) ik, FN, nref, nns, nnw

! then we need to initialize the Sgh and Lgh arrays
     if (allocated(Sgh)) deallocate(Sgh)
     if (allocated(Lgh)) deallocate(Lgh)
     allocate(Sgh(nns,nns,numset),Lgh(nns,nns))
     Sgh = czero
     Lgh = czero
     nat = 0
     call CalcSgh(cell,reflist,nns,numset,Sgh,nat)

! for now, we're disabling the kinematical part
! solve the dynamical eigenvalue equation for this beam direction  Lgh,thick,kn,nn,gzero,kin,debug
     kn = karray(4,ik)
     call CalcLgh(DynMat,Lgh,dble(thick(iE)),dble(kn),nns,gzero,depthstep,lambdaE(iE,1:izzmax),izzmax)
     deallocate(DynMat)

! dynamical contribution
     svals = 0.0
     do ix=1,numset
       svals(ix) = real(sum(Lgh(1:nns,1:nns)*Sgh(1:nns,1:nns,ix)))
     end do
     svals = svals/float(sum(nat(1:numset)))

! and store the resulting values, applying point group symmetry where needed.
     ipx = kij(1,ik)
     ipy = kij(2,ik)
     call Apply2DLaueSymmetry(ipx,ipy,isym,iequiv,nequiv)
!$OMP CRITICAL
     if (usehex) then
       do ix=1,nequiv
!        srhex(iequiv(1,ix),iequiv(2,ix),iE,1:numset) = svals(1:numset)
         srhex(iequiv(1,ix),iequiv(2,ix),1,1:numset) = svals(1:numset)
        end do
     else
        if (emnl%Esel.eq.-1) then
         do ix=1,nequiv
!          sr(iequiv(1,ix),iequiv(2,ix),iE,1:numset) = svals(1:numset)
           sr(iequiv(1,ix),iequiv(2,ix),1,1:numset) = svals(1:numset)
          end do
        else
         do ix=1,nequiv
           sr(iequiv(1,ix),iequiv(2,ix),1,1:numset) = svals(1:numset)
          end do
        endif
     end if
!$OMP END CRITICAL
  
     if (mod(ik,2500).eq.0) then
       io_int(1) = ik
       call WriteValue('  completed beam direction ',io_int, 1, "(I8)")
     end if

     call Delete_gvectorlist(reflist)

    end do beamloop

! end of OpenMP portion
!$OMP END PARALLEL

  deallocate(karray, kij)

! stop the clock and report the total time     
!   call Time_stop(numk)

! write (*,*) 'Some statistics :'
! write (*,*) 'Average number of strong beams : ',float(BetheParameter%totstrong)/float(numk)
! write (*,*) '          (min,max) : ',BetheParameter%minstrong,BetheParameter%maxstrong
! write (*,*) 'Average number of weak beams : ',float(BetheParameter%totweak)/float(numk)
! write (*,*) '          (min,max) : ',BetheParameter%minweak,BetheParameter%maxweak


! Finally, if this was sampled on a hexagonal array, we need to do barycentric interpolation
! to the standard square array for final output and use of the subsequent program.
! [this interpolation scheme must be verified; it is possible that there is an off-by-one error somewhere ...]
if (usehex) then
  delta = dsqrt(2.D0)/dble(emnl%npx)
  srt = 2.D0/dsqrt(3.D0)
! copy the central row without modifications
  srh(-emnl%npx:emnl%npx,0,1,1:numset) = srhex(-emnl%npx:emnl%npx,0,1,1:numset)
! we'll go through the array with pairs of horizontal rows at a time
  do j=1,npy-1
! determine which way the triangle is oriented for this row of the square array
    jh = floor(j*srt)
    if (mod(jh,2).eq.0) then ! even numbers mean triangle points down
      h = delta/srt - (j*delta - float(jh)*delta/srt)
      lambda = 0.5D0 - h/delta/dsqrt(3.D0)
      omtl = 1.D0-2.D0*lambda
      do i=-emnl%npx+1,emnl%npx-1  ! perform the barycentric interpolation
! positive row, pay attention to hexagonal coordinate transformation !
        srh(i,j,1,1:numset) = ( srhex(i-1,jh+1,1,1:numset) + srhex(i,jh+1,1,1:numset) )*lambda + omtl * srhex(i,jh,1,1:numset)
! negative row
        srh(i,-j,1,1:numset) = ( srhex(i-1,-jh-1,1,1:numset) + srhex(i,-jh-1,1,1:numset) )*lambda + omtl * srhex(i,-jh,1,1:numset)
      end do
    else
      h = j*delta - float(jh)*delta/srt
      lambda = 0.5D0 - h/delta/dsqrt(3.D0)
      omtl = 1.D0-2.D0*lambda
      do i=-emnl%npx+1,emnl%npx-1  ! perform the barycentric interpolation
!! positive row, pay attention to hexagonal coordinate transformation !
        srh(i,j,1,1:numset) = ( srhex(i-1,jh,1,1:numset) + srhex(i,jh,1,1:numset) )*lambda + omtl * srhex(i,jh+1,1,1:numset)
! negative row
        srh(i,-j,1,1:numset) = ( srhex(i-1,-jh,1,1:numset) + srhex(i,-jh,1,1:numset) )*lambda + omtl * srhex(i,-jh-1,1,1:numset)
      end do
    end if
  end do
  
! and finally, we convert this to a square Lambert projection which will be used 
! for all EBSD pattern interpolations; we do store the srh array in the HDF5 file as well for 
! visualization purposes.
  edge = LPs%sPio2 / dble(emnl%npx)
  scl = float(emnl%npx) / LPs%preg
  do i=-emnl%npx,emnl%npx
    do j=-npy,npy
! determine the spherical direction for this point
      xy = (/ dble(i), dble(j) /) * edge
      dc = LambertSquareToSphere(xy, ierr)
! and interpolate the hexagonal master pattern
      if (dc(3).lt.0.D0) dc = -dc
! convert direction cosines to hexagonal Lambert projections
      xy = scl * LambertSphereToHex( dc, ierr )
      tmp = xy(1)
      xy(1) = xy(2)
      xy(2) = tmp
! interpolate intensity from the neighboring points
      if (ierr.eq.0) then 
        nix = floor(xy(1))
        niy = floor(xy(2))
        nixp = nix+1
        niyp = niy+1
        if (nixp.gt.emnl%npx) nixp = nix
        if (niyp.gt.emnl%npx) niyp = niy
        dx = xy(1) - nix
        dy = xy(2) - niy
        dxm = 1.D0 - dx
        dym = 1.D0 - dy
        sr(i,j,1,1:numset) = srh(nix,niy,1,1:numset)*dxm*dym + srh(nixp,niy,1,1:numset)*dx*dym + &
                             srh(nix,niyp,1,1:numset)*dxm*dy + srh(nixp,niyp,1,1:numset)*dx*dy
      end if
    end do
  end do
end if

! since these computations can take a long time, here we store 
! all the output at the end of each pass through the energyloop.


  io_int(1) = nint(float(totstrong)/float(numk))
  call WriteValue(' -> Average number of strong reflections = ',io_int, 1, "(I5)")
  io_int(1) = nint(float(totweak)/float(numk))
  call WriteValue(' -> Average number of weak reflections   = ',io_int, 1, "(I5)")

! and here is where the major changes are for this version 5.0: all output now in HDF5 format
  call timestamp(timestring=tstre)

  nullify(HDF_head)
! Initialize FORTRAN HDF interface.
  call h5open_f(hdferr)

! open the existing file using the default properties.
  hdferr =  HDF_openFile(outname, HDF_head)

! update the time string
  groupname = 'EMheader'
  hdferr = HDF_openGroup(groupname, HDF_head)

  dataset = 'StopTime'
  line2(1) = tstre
  hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, overwrite)

  call HDF_pop(HDF_head)
  
  groupname = 'EMData'
  hdferr = HDF_openGroup(groupname, HDF_head)

! add data to the hyperslab
  dataset = 'sr'
  dims4 = (/  2*emnl%npx+1, 2*emnl%npx+1, numEbins, numset /)
  cnt4 = (/ 2*emnl%npx+1, 2*emnl%npx+1, 1, numset /)
  offset4 = (/ 0, 0, iE-1, 0 /)
  hdferr = HDF_writeHyperslabFloatArray4D(dataset, sr, dims4, offset4, cnt4(1), cnt4(2), cnt4(3), cnt4(4), HDF_head, insert)

  if (usehex) then
    dataset = 'srhex'
    dims4 = (/  2*emnl%npx+1, 2*emnl%npx+1, numEbins, numset /)
    cnt4 = (/ 2*emnl%npx+1, 2*emnl%npx+1, 1, numset /)
    offset4 = (/ 0, 0, iE-1, 0 /)
    hdferr = HDF_writeHyperslabFloatArray4D(dataset, srh, dims4, offset4, cnt4(1), cnt4(2), cnt4(3), cnt4(4), HDF_head, insert)
  end if

  call HDF_pop(HDF_head,.TRUE.)

! and close the fortran hdf interface
  call h5close_f(hdferr)

 if ((emnl%Esel.eq.-1).and.(iE.ne.1)) then 
  call Message('Intermediate data stored in file '//trim(emnl%outname), frm = "(A/)")
 end if

 if ((emnl%Esel.eq.-1).and.(iE.eq.1)) then 
  call Message('Final data stored in file '//trim(emnl%outname), frm = "(A/)")
 end if

end do energyloop

if (emnl%Esel.ne.-1) then
  call Message('Final data stored in file '//trim(emnl%outname), frm = "(A/)")
end if

end subroutine ComputeMasterPattern
