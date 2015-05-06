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
! EMsoft:EMEBSDmaster.f90
!--------------------------------------------------------------------------
!
! PROGRAM: EMEBSDmaster
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief ebsd master pattern on using scattering matrix GPU
!
!> @date 03/11/15 SS  1.0 original
!> @date 04/03/15 MDG 2.0 converted input/output to HDF format
!> @date 04/20/15 SS  3.0 made changes for approximate master pattern calculation
!--------------------------------------------------------------------------

program EMEBSDmaster

use local
use typedefs
use NameListTypedefs
use NameListHandlers
use files
use io

IMPLICIT NONE

character(fnlen)                                :: nmldeffile, progname, progdesc
type(EBSDMasterNameListType)                    :: emnl

nmldeffile = 'EMEBSDmaster.nml'
progname = 'EMEBSDmasterOpenCL.f90'
progdesc = 'Master pattern generation for Electron backscatter diffraction using scattering matrix on a GPU'

! print some information
call EMsoft(progname, progdesc)

! deal with the command line arguments, if any
call Interpret_Program_Arguments(nmldeffile,1,(/ 21 /), progname)

! deal with the namelist stuff
call GetEBSDMasterNameList(nmldeffile,emnl)

! perform the zone axis computations
call EBSDmasterpatternOpenCL(emnl, progname, nmldeffile)

end program EMEBSDmaster

!--------------------------------------------------------------------------
!
! SUBROUTINE:EBSDmasterpatternOpenCL
!
!> @author Marc De Graef/Saransh Singh, Carnegie Mellon University
!
!> @brief compute a master electron channeling pattern using scattering matrix on a GPU
!
!> @note Rewrite of the older ECPmaster program to perform the calculations on
!> a GPU. The bethe potential approximation combined with the GPU should really
!> things up a lot!
!
!> @param emnl name list structure
!> @param progname program name
!
!> @date 04/20/15  SS 1.0 original
!> @date 05/15/15 MDG 1.1 removed getenv() call; replaced by global path string
!--------------------------------------------------------------------------
subroutine EBSDmasterpatternOpenCL(emnl, progname, nmldeffile)

use local
use typedefs
use NameListTypedefs
use crystal
use symmetry
use Lambert
use initializers
use constants
use gvectors
use kvectors
use error
use io
use NameListHDFwriters
use files
use diffraction
use MBModule
use HDF5
use HDFsupport
use ISO_C_BINDING
use cl

IMPLICIT NONE

type(EBSDMasterNameListType),INTENT(INOUT)      :: emnl
character(fnlen),INTENT(IN)                     :: progname
character(fnlen),INTENT(IN)                     :: nmldeffile

real(kind=dbl)                                  :: frac
integer(kind=irg)                               :: gzero, istat

integer(kind=irg)                               :: numEbins, numzbins, nx, ny, totnum_el ! reading from MC file
real(kind=dbl)                                  :: EkeV, Ehistmin, Ebinsize, depthmax, depthstep, sig, omega  ! reading from MC file
integer(kind=irg), allocatable                  :: accum_e(:,:,:), accum_z(:,:,:,:) ! reading from MC file

integer(kind=irg)                               :: io_int_sgl(1), io_int(6) ! integer output variable
real(kind=dbl)                                  :: io_real(5) ! real output variable

integer(kind=irg)                               :: i, j, isym, pgnum,ix,iy ! variables for point group and Laue group
integer(kind=irg),parameter                     :: LaueTest(11) = (/ 149, 151, 153, 156, 158, 160, 161, 164, 165, 166, 167 /)  ! space groups with 2 or mirror at 30 degrees
integer(kind=irg)                               :: npyhex, ijmax, numk, skip ! parameters for calckvectors and calcwavelength subroutine

integer(kind=irg)                               :: ga(3), gb(3) ! shortest reciprocal lattice vector for zone axis
real(kind=sgl), allocatable                     :: thick(:), sr(:,:,:,:), lambdaE(:,:)
real(kind=dbl)                                  :: intthick
complex(kind=dbl),allocatable                   :: Lgh(:,:),Sgh(:,:),Sghtmp(:,:,:)
complex(kind=dbl),allocatable                   :: DynMat(:,:)
complex(kind=dbl)                               :: czero

integer(kind=irg)                               :: nt, nns, nnw, tots, totw ! thickness array and BetheParameters strong and weak beams
real(kind=sgl)                                  :: FN(3), k(3), fnat, kn
integer(kind=irg)                               :: numset, nref, ipx, ipy, iequiv(2,12), nequiv, ip, jp, izz, iE, iz, one,ierr,ss
integer(kind=irg)                               :: izzmax
integer(kind=irg),allocatable                   :: kij(:,:), nat(:)
real(kind=dbl)                                  :: res(2),selE

character(fnlen)                                :: oldprogname
character(fnlen)                                :: xtalname
character(8)                                    :: MCscversion
character(4)                                    :: MCmode
character(6)                                    :: projtype

logical                                         :: verbose, usehex, switchmirror
type(unitcell), pointer                         :: cell
type(gnode)                                     :: rlp
type(DynType)                                   :: Dyn
type(kvectorlist), pointer                      :: khead,ktmp,kheadtmp ! linked list for incident wave vectors for master list
type(kvectorlist), pointer                      :: kheadcone,ktmpcone ! linked list for incident wave vectors for individual pattern
real(kind=dbl),allocatable                      :: ecpattern(:,:)
type(BetheParameterType)                        :: BetheParameters
type(reflisttype),pointer                       :: reflist, firstw,rltmp

integer(kind=8)                                 :: size_in_bytes,size_in_bytes_coeff,size_in_bytes_lambda,&
size_in_bytes_dynmatsz,size_in_bytes_wavefn,size_in_bytes_offdiagonal
integer, parameter                              :: iunit = 10
real(kind=4)                                    :: absmax,e,pf,eps
character(len = 100)                            :: info ! info about the GPU
integer(kind=8)                                 :: globalsize(2),localsize(2)
integer, parameter                              :: source_length = 100000
character(len = source_length)                  :: source
real(kind=4),allocatable                        :: lambdas(:),EkeVs(:)
integer(kind=4)                                 :: num,irec,ii,jj,kk,pp,qq,ll,mm,npx,npy,npiximgx,npiximgy

complex(kind=4)                                 :: coef(3,3),cvals(9)
integer(kind=4)                                 :: numdepth
integer(kind=4),allocatable                     :: locpix(:,:)
integer(kind=4),allocatable                     :: arrsize(:),arrsizesum(:),offset(:),ns(:)
complex(kind=4),allocatable                     :: LghCumulative(:),SghCumulative(:,:),A(:)
complex(kind=4),allocatable                     :: diagonal(:),off_diagonal(:)

integer(kind=4)                                 :: size1,size2

type(cl_platform_id)                            :: platform
type(cl_device_id)                              :: device
type(cl_context)                                :: context
type(cl_command_queue)                          :: command_queue
type(cl_program)                                :: prog
type(cl_kernel)                                 :: kernel
type(cl_event)                                  :: event
type(cl_mem)                                    :: cl_expA,cl_AA,cl_AAA,cl_coeff,cl_T1,cl_T2,&
                                                   cl_wavefncoeff,cl_wavefncoeffintd,&
                                                   cl_sqrsize,cl_lambdas,cl_diagonal,cl_offdiagonal
integer(kind=irg)                               :: istart,iend,jstart,jend,lamcenter(2)

character(11)                                   :: dstr
character(15)                                   :: tstrb
character(15)                                   :: tstre
integer(HSIZE_T)                                :: dims4(4)
logical                                         :: f_exists, readonly
character(fnlen, KIND=c_char),allocatable,TARGET:: stringarray(:)
character(fnlen)                                :: dataset, groupname
integer                                         :: hdferr, nlines, nsx, nsy, num_el

type(HDFobjectStackType),pointer                :: HDF_head

nullify(HDF_head)

call timestamp(datestring=dstr, timestring=tstrb)


gzero = 1
frac = 0.05
eps = 1.0
!dataunit = 10

globalsize = (/ 32, 32 /)
localsize = (/ 4, 4 /)

npiximgx = emnl%npx
npiximgy = emnl%npx
npx = globalsize(1)
npy = globalsize(2)


size_in_bytes_dynmatsz = npx*npy*sizeof(npx)
size_in_bytes_coeff = 9*sizeof(coef(1,1))


cvals(1) = cmplx(-3.3335514852690488032942739163345055,0.0)
cvals(2) = cmplx(-3.0386480729366970892124687564926859,-1.5868011957588383288038677051222921)
cvals(3) = cmplx(-3.0386480729366970892124687564926859,+1.5868011957588383288038677051222921)
cvals(4) = cmplx(-2.1108398003026547374987047865183922,-3.0899109287255009227777015426228801)
cvals(5) = cmplx(-2.1108398003026547374987047865183922,+3.0899109287255009227777015426228801)
cvals(6) = cmplx(-0.38106984566311299903129424501333242,-4.3846445331453979503692027283066828)
cvals(7) = cmplx(-0.38106984566311299903129424501333242,+4.3846445331453979503692027283066828)
cvals(8) = cmplx(2.6973334615369892273896047461916633,-5.1841620626494141778340870727109629)
cvals(9) = cmplx(2.6973334615369892273896047461916633,+5.1841620626494141778340870727109629)

coef(1,1) = cvals(1)+cvals(2)+cvals(3)
coef(2,1) = cvals(1)*cvals(2)+cvals(2)*cvals(3)+cvals(3)*cvals(1)
coef(3,1) = cvals(1)*cvals(2)*cvals(3)

coef(1,2) = cvals(4)+cvals(5)+cvals(6)
coef(2,2) = cvals(4)*cvals(5)+cvals(5)*cvals(6)+cvals(6)*cvals(4)
coef(3,2) = cvals(4)*cvals(5)*cvals(6)

coef(1,3) = cvals(7)+cvals(8)+cvals(9)
coef(2,3) = cvals(7)*cvals(8)+cvals(8)*cvals(9)+cvals(9)*cvals(7)
coef(3,3) = cvals(7)*cvals(8)*cvals(9)

!=============================================================
!read Monte Carlo output file and extract necessary parameters
! first, we need to load the data from the MC program.
!=============================================================

call Message('opening '//trim(emnl%energyfile), frm = "(A)" )

! [since this is no longer a sequential access file, we do not need to read everything, just
! the quantities that we need...]

!
! Initialize FORTRAN interface.
!
call h5open_f(hdferr)

! first of all, if the file exists, then delete it and rewrite it on each energyloop
inquire(file=trim(emnl%energyfile), exist=f_exists)

if (.not.f_exists) then
  call FatalError('ComputeMasterPattern','Monte Carlo input file does not exist')
end if

! open the MC file using the default properties.
readonly = .TRUE.
hdferr =  HDF_openFile(emnl%energyfile, HDF_head, readonly)

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

!allocate(accum_z(numEbins,numzbins,-nsx/10:nsx/10,-nsy/10:nsy/10),stat=istat)

dataset = 'accum_z'
!dims4 =  (/ numEbins, numzbins, 2*(nsx/10)+1,2*(nsy/10)+1 /)
accum_z = HDF_readDatasetIntegerArray4D(dataset, dims4, HDF_head)

! and close everything
call HDF_pop(HDF_head,.TRUE.)

! close the fortran interface
call h5close_f(hdferr)

call Message(' -> completed reading '//trim(emnl%energyfile), frm = "(A)")
!=============================================
! completed reading monte carlo file
!=============================================

nullify(cell)
allocate(cell)

! load the crystal structure and compute the Fourier coefficient lookup table
verbose = .TRUE.
call Initialize_Cell(cell,Dyn,rlp, xtalname, emnl%dmin, sngl(1000.0*EkeV),verbose)

! determine the point group and Laue group number
j=0
do i=1,32
    if (SGPG(i).le.cell%SYM_SGnum) j=i
end do
isym = j
pgnum = j

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

if(usehex) then
    npyhex = nint(2.0*float(emnl%npx)/sqrt(3.0))
    npiximgy = npyhex
end if
ijmax = float(emnl%npx)**2   ! truncation value for beam directions

czero = cmplx(0.D0, 0.D0)
! force dynamical matrix routine to read new Bethe parameters from file
call Set_Bethe_Parameters(BetheParameters,.TRUE.)

!nt = nint((depthmax - emnl%startthick)/depthstep)
!allocate(thick(nt))
!thick = emnl%startthick + depthstep * (/ (i-1,i=1,nt) /)
nullify(reflist)
nullify(firstw)

nns = 0
nnw = 0
tots = 0
totw = 0

numset = cell % ATOM_ntype  ! number of special positions in the unit cell
izz = numzbins

!allocate(lambdaZ(1:izz),stat=istat)
allocate(nat(numset),stat=istat)
!allocate(kij(2,numk),stat=istat)

numdepth = 0

allocate(EkeVs(numEbins),thick(numEbins))

do ii=1,numEbins
    EkeVs(ii) = Ehistmin + float(ii-1)*Ebinsize
end do

! then, for each energy determine the 95% histogram thickness
izzmax = 0
do iE = 1,numEbins
    do ix=-nx/10,nx/10
        do iy=-ny/10,ny/10
            istat = sum(accum_z(iE,:,ix,iy))
            izz = 1
            do while (sum(accum_z(iE,1:izz,ix,iy)).lt.(0.95*istat))
                izz = izz+1
            end do
            if (izz.gt.izzmax) izzmax = izz
        end do
    end do
    thick(iE) = dble(izzmax) * depthstep
end do

izz = nint(maxval(thick)/depthstep)

allocate(lambdaE(1:numEbins,1:izz),stat=istat)
allocate(lambdas(1:izz),stat=istat)

do iE=1,numEbins
    do iz=1,izz
        lambdaE(iE,iz) = float(sum(accum_z(iE,iz,:,:)))/float(sum(accum_z(:,:,:,:)))
    end do
end do

size_in_bytes_lambda = izz*sizeof(lambdas(1))

nat = 0
fnat = 1.0/float(sum(cell%numat(1:numset)))
intthick = dble(depthmax)

!==========================
! INITIALIZATION OF GPU
!==========================

! get the platform ID
call clGetPlatformIDs(platform, num, ierr)
if(ierr /= CL_SUCCESS) stop "Cannot get CL platform."

! get the device ID
call clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, device, num, ierr)
if(ierr /= CL_SUCCESS) stop "Cannot get CL device."

! get the device name and print it
call clGetDeviceInfo(device, CL_DEVICE_NAME, info, ierr)
write(6,*) "CL device: ", info

! create the context and the command queue
context = clCreateContext(platform, device, ierr)
if(ierr /= CL_SUCCESS) stop "Cannot create context"
command_queue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, ierr)
if(ierr /= CL_SUCCESS) stop "Cannot create command queue"

!=====================
! BUILD THE KERNEL
!=====================

! read the source file
open(unit = iunit, file = trim(openclpathname)//'MBmoduleOpenCL.cl', access='direct', status = 'old', &
action = 'read', iostat = ierr, recl = 1)
if (ierr /= 0) stop 'Cannot open file MBmoduleOpenCL.cl'

source = ''
irec = 1
do
read(unit = iunit, rec = irec, iostat = ierr) source(irec:irec)
if (ierr /= 0) exit
if(irec == source_length) stop 'Error: CL source file is too big'
irec = irec + 1
end do
close(unit=iunit)

! create the program
prog = clCreateProgramWithSource(context, source, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot create program from source.'

! build
call clBuildProgram(prog, '-cl-no-signed-zeros', ierr)

! get the compilation log
call clGetProgramBuildInfo(prog, device, CL_PROGRAM_BUILD_LOG, source, irec)
if(len(trim(source)) > 0) print*, trim(source)

if(ierr /= CL_SUCCESS) stop 'Error: program build failed.'

write(6,*) "Kernel Build Successful...."

! finally get the kernel and release the program
kernel = clCreateKernel(prog, 'CalcLghMaster', ierr)
call clReleaseProgram(prog, ierr)

allocate(sr(2*emnl%npx+1,2*emnl%npx+1,1:numEbins,1:numset),stat=istat)
sr = 0.0

allocate(ns(npx*npy),locpix(npx*npy,2),stat=istat)

select case (isym)

    case (1)  ! triclinic symmetry
        istart = -npiximgx
        iend = npiximgx
        jstart = -npiximgy
        jend = npiximgy

    case (2)   !  monoclinic symmetry
        istart = -npiximgx
        iend = npiximgx
        jstart = 0
        jend = npiximgy

    case (3,4,10,5,11)  ! orthorhombic mmm, tetragonal 4/m, cubic m-3
        istart = 0
        iend = npiximgx
        jstart = 0
        jend = npiximgy

!======================================================================
! NOT READY TO BE USED STILL
!======================================================================
    case (6,7,8,9,12)   ! npy is now npyhex ! These cases use the hexagonal lambert maps and have to be dealt with seperately. The CalckvectorsGPU is only for square lambert projections
        istart = 0
        iend = npy
        jstart = 0
        jend = npy

end select

energyloop: do iE=numEbins,1,-1

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
! generating list of incident wave vectors
!=============================================

! determine all independent incident beam directions (use a linked list starting at khead)
! numk is the total number of k-vectors to be included in this computation;


!if (usehex) then
!    call Calckvectors(khead,cell, (/ 0.D0, 0.D0, 1.D0 /), (/ 0.D0, 0.D0, 0.D0 /),0.D0,emnl%npx,npyhex,numk, &
!    isym,ijmax,'RoscaLambert',usehex)
!else
! Calckvectors(k,ga,ktmax,npx,npy,numk,isym,ijmax,mapmode,usehex)
!call Calckvectors(khead,cell, (/ 0.D0, 0.D0, 1.D0 /), (/ 0.D0, 0.D0, 0.D0 /),0.D0,emnl%npx,emnl%npx,numk, &
!    isym,ijmax,'RoscaLambert',usehex)
!end if

!io_int_sgl(1)=numk

!call WriteValue('# independent beam directions to be considered = ', io_int_sgl, 1, "(I8)")

    lambdas = 0.0
    lambdas(:) = lambdaE(iE,:)
    numdepth = thick(iE)

    beamloopx : do ii = -floor(float(istart)/float(npx)),ceiling(float(iend)/float(npx))
        beamloopy: do jj = -floor(float(jstart)/float(npy)),ceiling(float(jend)/float(npy))



            ns = 0
            locpix = 0

                nullify(khead)

                lamcenter = (/istart+(ii+floor(float(istart)/float(npx))+1)*npx/2,&
                              jstart+(jj+floor(float(jstart)/float(npy))+1)*npy/2 /)
                call CalckvectorsGPU(khead,cell,npiximgx,npx/2,lamcenter,numk)

                ktmp => khead

                k = ktmp%k
                FN = k
                call Initialize_ReflectionList(cell, reflist, BetheParameters, FN,  k, emnl%dmin, nref)

! we have to ignore the bethe approximation in this case. hopefully, the GPU kernel is fast enough to make up for the offset

                if (allocated(Sghtmp)) deallocate(Sghtmp)
                if (allocated(Sgh)) deallocate(Sgh)
                if (allocated(DynMat)) deallocate(DynMat)
                if (allocated(diagonal)) deallocate(diagonal)
                if (allocated(off_diagonal)) deallocate(off_diagonal)


                allocate(Sghtmp(nref,nref,numset),stat=istat)
                allocate(Sgh(nref*nref,numset),stat=istat)
                allocate(diagonal(nref*npx*npy),stat=istat)
                allocate(off_diagonal(nref*nref),stat=istat)
                allocate(DynMat(nref,nref),stat=istat)

                call CalcSghMaster(cell,reflist,nref,numset,Sghtmp,nat)

                do pp = 1,numset
                    do kk = 1,nref
                        do ll = 1,nref
                            Sgh((kk-1)*nref+ll,pp) = Sghtmp(kk,ll,pp)
                        end do
                    end do
                end do

            call GetDynMatmaster(cell,reflist,DynMat,nref)

            do kk = 1,nref
                do ll = 1,nref
                    off_diagonal((kk-1)*nref+ll) = DynMat(kk,ll)
                end do
            end do

            ktmp => khead
            do kk = 1,npx
                do ll = 1,npy
                    k = ktmp%k
                    FN = k
                    locpix((kk-1)*npy+ll,:) = (/ktmp%i,ktmp%j/)
                    rltmp => reflist%next

                    do pp = 1,nref
                        diagonal(((kk-1)*npy+ll-1)*nref+pp) = Calcsg(cell,float(rltmp%hkl),k,FN)
                        rltmp => rltmp%next
                    end do

                    absmax = maxval(abs(diagonal(((kk-1)*npy+ll-1)*nref+1:((kk-1)*npy+ll-1)*nref+nref)))
                    e = ceiling(log(absmax)/log(2.0))
                    ns((kk-1)*npy+ll) = e+1
                    if (ns((kk-1)*npy+ll) .le. 0) then
                        ns((kk-1)*npy+ll) = 1
                    else

                        diagonal = diagonal/2**ns((kk-1)*npy+ll)
                    end if

                    ktmp => ktmp%next

                end do
            end do

            diagonal = diagonal*cmplx(0.0,cPi*cell%mLambda)
            diagonal = diagonal*cmplx(eps,0.0)

            off_diagonal = off_diagonal*cmplx(0.0,cPi*cell%mLambda)
            off_diagonal = off_diagonal*cmplx(eps,0.0)

! allocate device memory
            size_in_bytes = nref*nref*npx*npy*sizeof(diagonal(1))
            size_in_bytes_wavefn = nref*npx*npy*sizeof(diagonal(1))
            size_in_bytes_offdiagonal = nref*nref*sizeof(off_diagonal(1))

            if (allocated(LghCumulative)) deallocate(LghCumulative)
            allocate(LghCumulative(1:nref*nref*npx*npy),stat=istat)

            cl_expA = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory for exponential.'

            cl_offdiagonal = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_offdiagonal, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory for input offdiagonal matrix.'

            cl_diagonal = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_wavefn, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory for input diagonal matrix.'

            cl_AA = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory for square of A.'

            cl_AAA = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory for cube of A.'

            cl_T1 = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory for T1.'

            cl_T2 = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory for T2.'

            cl_wavefncoeff = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_wavefn, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory for fourier coefficients .'

            cl_wavefncoeffintd = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_wavefn, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory for intermediate array.'

            cl_coeff = clCreateBuffer(context, CL_MEM_READ_ONLY, size_in_bytes_coeff, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory for ceofficients.'

            cl_sqrsize = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_dynmatsz, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory for wavefunction coefficient offset information.'

            cl_lambdas = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_dynmatsz, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory for lambda values.'

! write list of initial arrays to buffer
            call clEnqueueWriteBuffer(command_queue, cl_diagonal, cl_bool(.true.), 0_8, size_in_bytes_wavefn, diagonal(1), ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot write to buffer.'

            call clEnqueueWriteBuffer(command_queue, cl_offdiagonal, cl_bool(.true.), 0_8,&
            size_in_bytes_offdiagonal, off_diagonal(1), ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot write to buffer.'

            call clEnqueueWriteBuffer(command_queue, cl_sqrsize, cl_bool(.true.), 0_8, size_in_bytes_dynmatsz, ns(1), ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot write to buffer.'

            call clEnqueueWriteBuffer(command_queue, cl_coeff, cl_bool(.true.), 0_8, size_in_bytes_coeff, coef(1,1), ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot write to buffer.'

            call clEnqueueWriteBuffer(command_queue, cl_lambdas, cl_bool(.true.), 0_8, size_in_bytes_lambda, lambdas(1), ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot write to buffer.'

! set kernel arguments

            call clSetKernelArg(kernel, 0, cl_expA, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot set first kernel argument.'

            call clSetKernelArg(kernel, 1, cl_offdiagonal, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot set second kernel argument.'

            call clSetKernelArg(kernel, 2, cl_diagonal, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot set third kernel argument.'

            call clSetKernelArg(kernel, 3, cl_AA, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot set fourth kernel argument.'

            call clSetKernelArg(kernel, 4, cl_AAA, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot set fifth kernel argument.'

            call clSetKernelArg(kernel, 5, cl_coeff, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot set sixth kernel argument.'

            call clSetKernelArg(kernel, 6, cl_T1, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot set seventhth kernel argument.'

            call clSetKernelArg(kernel, 7, cl_T2, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot set eighth kernel argument.'

            call clSetKernelArg(kernel, 8, cl_wavefncoeff, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot set nineth kernel argument.'

            call clSetKernelArg(kernel, 9, cl_wavefncoeffintd, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot set tenth kernel argument.'

            call clSetKernelArg(kernel, 10, nref, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot set eleventh kernel argument.'

            call clSetKernelArg(kernel, 11, cl_sqrsize, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot set twelveth kernel argument.'

            call clSetKernelArg(kernel, 12, numdepth, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot set thirteenth kernel argument.'

            call clSetKernelArg(kernel, 13, cl_lambdas, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot set fourteenth kernel argument.'

!execute the kernel

            call clEnqueueNDRangeKernel(command_queue, kernel, globalsize, localsize, event, ierr)
! wait for the commands to finish

            call clFinish(command_queue, ierr)
! read the resulting vector from device memory

            call clEnqueueReadBuffer(command_queue, cl_T1, cl_bool(.true.), 0_8, size_in_bytes, LghCumulative(1), ierr)

! divide by integration depth explicitly (OpenCL kernel giving some problems)
            LghCumulative = LghCumulative/float(izz-1)
            do pp = 1,npx
                do qq = 1,npy
                    ipx = locpix((pp-1)*npy+qq,1)
                    ipy = locpix((pp-1)*npy+qq,2)
                    size1 = nref
                    size2 = ((pp-1)*npy+(qq-1))*nref*nref
                    if (isym.ne.1) then
                        call Apply2DLaueSymmetry(ipx,ipy,isym,iequiv,nequiv)
                        iequiv(1,1:nequiv) = iequiv(1,1:nequiv) + npiximgx + 1
                        iequiv(2,1:nequiv) = iequiv(2,1:nequiv) + npiximgy + 1
                        do ss = 1,numset
                            do ip=1,nequiv
                                sr(iequiv(1,ip),iequiv(2,ip),iE,ss) = real(sum(Sgh(1:nref*nref,ss)*&
                                LghCumulative(size2+1:size2+size1**2)))!/float(sum(nat))
                            end do
                        end do
                    else
                        do ss = 1,numset
                            sr(ipx+npiximgx+1,ipy+npiximgy+1,iE,ss) = real(sum(Sgh(1:nref*nref,ss)*&
                            LghCumulative(size2+1:size2+size1**2)))!/float(sum(nat))
                        end do
                    end if
                end do
            end do

      ! if (mod(ii,4) .eq. 0) then
      !      write(6,'(A,I8,A)') 'Completed ',(ii*npx*npy),' beams.'
      ! end if
            call clReleaseMemObject(cl_expA, ierr)
            call clReleaseMemObject(cl_diagonal, ierr)
            call clReleaseMemObject(cl_offdiagonal, ierr)
            call clReleaseMemObject(cl_AA, ierr)
            call clReleaseMemObject(cl_AAA, ierr)
            call clReleaseMemObject(cl_T1, ierr)
            call clReleaseMemObject(cl_T2, ierr)
            call clReleaseMemObject(cl_wavefncoeff, ierr)
            call clReleaseMemObject(cl_wavefncoeffintd, ierr)
            call clReleaseMemObject(cl_coeff, ierr)
            call clReleaseMemObject(cl_lambdas, ierr)

            !call Delete_kvectorlist(khead)
        end do beamloopy
write(6,*) 'Ended one batch'
    end do beamloopx
open(unit=11,file='test.txt',action='write')
do i=1,2*npiximgx+1
do j=1,2*npiximgy+1
write(11, '(F15.6)', advance='no') sr(i,j,3,1)
end do
write(11, *) ''  ! this gives you the line break
end do
close(11)
end do energyloop

! Initialize FORTRAN interface.
!
call h5open_f(hdferr)
call timestamp(timestring=tstre)

! first of all, if the file exists, then delete it and rewrite it on each energyloop
inquire(file=trim(emnl%outname), exist=f_exists)

if (f_exists) then
    open(unit=dataunit, file=trim(emnl%outname), status='old',form='unformatted')
    close(unit=dataunit, status='delete')
end if

! Create a new file using the default properties.
hdferr =  HDF_createFile(emnl%outname, HDF_head)

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
    hdferr = HDF_writeDatasetFloat(dataset, sngl(selE), HDF_head)
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
dataset = 'sr'
hdferr = HDF_writeDatasetFloatArray4D(dataset, sr, 2*emnl%npx+1, 2*emnl%npx+1, numEbins, numset, HDF_head)

call HDF_pop(HDF_head,.TRUE.)
! and close the fortran hdf interface
call h5close_f(hdferr)

call Message('Final data stored in file '//trim(emnl%outname), frm = "(A/)")


end subroutine EBSDmasterpatternOpenCL
