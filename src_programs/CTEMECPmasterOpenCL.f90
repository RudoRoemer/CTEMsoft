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
! CTEMsoft2013:CTEMECPmaster.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMECPmaster
!
!> @author Marc De Graef/ Saransh Singh, Carnegie Mellon University
!
!> @brief Zone axis electron channeling master patterns on using scattering matrix GPU
!
!> @date 03/06/14 SS 1.0 original
!--------------------------------------------------------------------------

program CTEMECPmaster

use local
use typedefs
use NameListTypedefs
use NameListHandlers
use files
use io

IMPLICIT NONE

character(fnlen)                                :: nmldeffile, progname, progdesc
type(ECPMasterNameListType)                     :: ecpnl

nmldeffile = 'CTEMECPmaster.nml'
progname = 'CTEMECPmasterOpenCL.f90'
progdesc = 'Master pattern generation for Electron channeling pattern using scattering matrix on a GPU'

! deal with the command line arguments, if any
call Interpret_Program_Arguments(nmldeffile,2,(/ 0, 40 /), progname)

! deal with the namelist stuff
call GetECPMasterNameList(nmldeffile,ecpnl)

! print some information
call CTEMsoft(progname, progdesc)

! perform the zone axis computations
call ECmasterpatternOpenCL(ecpnl, progname)

end program CTEMECPmaster

!--------------------------------------------------------------------------
!
! SUBROUTINE:ECmasterpatternOpenCL
!
!> @author Marc De Graef/Saransh Singh, Carnegie Mellon University
!
!> @brief compute a master electron channeling pattern using scattering matrix on a GPU
!
!> @note Rewrite of the older ECPmaster program to perform the calculations on
!> a GPU. The bethe potential approximation combined with the GPU should really
!> things up a lot!
!
!> @param ecpnl name list structure
!> @param progname program name
!
!> @date 03/06/14  SS 1.0 original
!--------------------------------------------------------------------------
subroutine ECmasterpatternOpenCL(ecpnl, progname)

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
use files
use diffraction
use MBModule
use cl

IMPLICIT NONE

type(ECPMasterNameListType),INTENT(IN)          :: ecpnl
character(fnlen),INTENT(IN)                     :: progname

real(kind=dbl)                                  :: frac
integer(kind=irg)                               :: gzero, istat

integer(kind=irg)                               :: numEbins, numzbins, nx, ny, totnum_el ! reading from MC file
real(kind=dbl)                                  :: EkeV, Ehistmin, Ebinsize, depthmax, depthstep, sig, omega  ! reading from MC file
integer(kind=irg), allocatable                  :: accum_e(:,:,:), accum_z(:,:,:,:) ! reading from MC file

integer(kind=irg)                               :: io_int_sgl(1), io_int(6) ! integer output variable
real(kind=dbl)                                  :: io_real(5) ! real output variable

integer(kind=irg)                               :: i, j, isym, pgnum ! variables for point group and Laue group
integer(kind=irg),parameter                     :: LaueTest(11) = (/ 149, 151, 153, 156, 158, 160, 161, 164, 165, 166, 167 /)  ! space groups with 2 or mirror at 30 degrees
integer(kind=irg)                               :: npyhex, ijmax, numk, skip ! parameters for calckvectors and calcwavelength subroutine

integer(kind=irg)                               :: ga(3), gb(3) ! shortest reciprocal lattice vector for zone axis
real(kind=sgl), allocatable                     :: thick(:), sr(:,:), lambdaZ(:)
real(kind=dbl)                                  :: intthick
complex(kind=dbl),allocatable                   :: Lgh(:,:),Sgh(:,:),Sghtmp(:,:,:)
complex(kind=dbl),allocatable                   :: DynMat(:,:)
complex(kind=dbl)                               :: czero

integer(kind=irg)                               :: nt, nns, nnw, tots, totw ! thickness array and BetheParameters strong and weak beams
real(kind=sgl)                                  :: FN(3), k(3), fnat, kn
integer(kind=irg)                               :: numset, nref, ipx, ipy, iequiv(2,12), nequiv, ip, jp, izz, IE, iz, one,ierr
integer(kind=irg),allocatable                   :: kij(:,:), nat(:)
real(kind=dbl)                                  :: res(2)

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
size_in_bytes_dynmatsz,size_in_bytes_wavefn
integer, parameter                              :: iunit = 10
real(kind=4)                                    :: absmax,e,pf,eps
character(len = 100)                            :: info ! info about the GPU
integer(kind=8)                                 :: globalsize(2),localsize(2)
integer, parameter                              :: source_length = 100000
character(len = source_length)                  :: source
real(kind=4),allocatable                        :: lambdas(:)
integer(kind=4)                                 :: num,irec,ii,jj,kk,pp,qq,ll,mm,npx,npy,npiximgx,npiximgy

complex(kind=4)                                 :: coef(3,3),cvals(9)
integer(kind=4)                                 :: numdepth
integer(kind=4),allocatable                     :: locpix(:,:)
integer(kind=4),allocatable                     :: arrsize(:),arrsizesum(:),offset(:),ns(:)
complex(kind=4),allocatable                     :: LghCumulative(:),SghCumulative(:),A(:)

integer(kind=4)                                 :: size1,size2

type(cl_platform_id)                            :: platform
type(cl_device_id)                              :: device
type(cl_context)                                :: context
type(cl_command_queue)                          :: command_queue
type(cl_program)                                :: prog
type(cl_kernel)                                 :: kernel
type(cl_event)                                  :: event
type(cl_mem)                                    :: cl_expA,cl_A,cl_AA,cl_AAA,cl_coeff,cl_T1,cl_T2,&
cl_wavefncoeff,cl_wavefncoeffintd,&
cl_offset,cl_arrsize,cl_arrsizesum,cl_sqrsize,cl_lambdas


gzero = 1
frac = 0.05
eps = 1.0
!dataunit = 10

globalsize = (/ 32, 32 /)
localsize = (/ 4, 4 /)

npiximgx = ecpnl%npx
npiximgy = ecpnl%npx
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


allocate(cell)

!=============================================================
!read Monte Carlo output file and extract necessary parameters
! first, we need to load the data from the MC program.
!=============================================================

call Message('opening '//trim(ecpnl%energyfile), frm = "(A)" )

open(dataunit,file=trim(ecpnl%energyfile),status='unknown',form='unformatted')

! lines from CTEMMCCL.f90... these are the things we need to read in...
! write (dataunit) progname
!! write the version number
! write (dataunit) MCscversion
!! then the name of the crystal data file
! write (dataunit) xtalname
!! energy information etc...
! write (dataunit) numEbins, numzbins, numsx, numsy, , totnum_el
! write (dataunit) EkeV, Ehistmin, Ebinsize, depthmax, depthstep
! write (dataunit) sig, omega
! write (dataunit) MCmode
!! and here are the actual results
! write (dataunit) accum_e
! write (dataunit) accum_z

read (dataunit) oldprogname
read (dataunit) MCscversion
read (dataunit) xtalname

read(dataunit) numEbins, numzbins, nx, ny, totnum_el
nx = (nx - 1)/2
ny = (ny - 1)/2

read (dataunit) EkeV, Ehistmin, Ebinsize, depthmax, depthstep
io_real(1:5) = (/ EkeV, Ehistmin, Ebinsize, depthmax, depthstep /)
call WriteValue(' EkeV, Ehistmin, Ebinsize, depthmax, depthstep ',io_real,5,"(4F10.5,',',F10.5)")

read (dataunit) sig, omega
read (dataunit) MCmode

allocate(accum_e(numEbins,-nx:nx,-nx:nx),accum_z(numEbins,numzbins,-nx/10:nx/10,-nx/10:nx/10),stat=istat)

read(dataunit) accum_e
! actually, we do not yet need the accum_e array for ECP. This will be removed with an updated version of the MC code
! but we need to skip it in this unformatted file so that we can read the accum_z array ...
!deallocate(accum_e)

read(dataunit) accum_z    ! we only need this array for the depth integrations
close(dataunit,status='keep')
call Message(' -> completed reading '//trim(ecpnl%energyfile), frm = "(A)")
!=============================================
! completed reading monte carlo file
!=============================================

nullify(cell)
allocate(cell)

! load the crystal structure and compute the Fourier coefficient lookup table
verbose = .TRUE.
call Initialize_Cell(cell,Dyn,rlp, xtalname, ecpnl%dmin, sngl(1000.0*EkeV),verbose)

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

if(usehex)  npyhex = nint(2.0*float(ecpnl%npx)/sqrt(3.0))
ijmax = float(ecpnl%npx)**2   ! truncation value for beam directions


!=============================================
! generating list of incident wave vectors
!=============================================

! determine all independent incident beam directions (use a linked list starting at khead)
! numk is the total number of k-vectors to be included in this computation;
nullify(khead)
if (usehex) then
call Calckvectors(khead,cell, (/ 0.D0, 0.D0, 1.D0 /), (/ 0.D0, 0.D0, 0.D0 /),0.D0,ecpnl%npx,npyhex,numk, &
isym,ijmax,'RoscaLambert',usehex)
else
! Calckvectors(k,ga,ktmax,npx,npy,numk,isym,ijmax,mapmode,usehex)
call Calckvectors(khead,cell, (/ 0.D0, 0.D0, 1.D0 /), (/ 0.D0, 0.D0, 0.D0 /),0.D0,ecpnl%npx,ecpnl%npx,numk, &
isym,ijmax,'RoscaLambert',usehex)
end if
io_int_sgl(1)=numk

call WriteValue('# independent beam directions to be considered = ', io_int_sgl, 1, "(I8)")

ktmp => khead
czero = cmplx(0.D0, 0.D0)
! force dynamical matrix routine to read new Bethe parameters from file
call Set_Bethe_Parameters(BetheParameters,.TRUE.)

nt = nint((depthmax - ecpnl%startthick)/depthstep)
allocate(thick(nt))
thick = ecpnl%startthick + depthstep * (/ (i-1,i=1,nt) /)

nullify(reflist)
nullify(firstw)

nns = 0
nnw = 0
tots = 0
totw = 0

numset = cell % ATOM_ntype  ! number of special positions in the unit cell
izz = numzbins

allocate(lambdaZ(1:izz),stat=istat)
allocate(nat(numset),stat=istat)
allocate(kij(2,numk),stat=istat)

numdepth = 0

do iz=1,izz
    lambdaZ(iz) = float(sum(accum_z(:,iz,:,:)))/float(sum(accum_z(:,:,:,:)))
    if ((lambdaZ(iz)) .ne. 0.0) numdepth = numdepth+1
end do

deallocate(lambdaZ)
allocate(lambdaZ(1:numdepth),stat=istat)

do iz=1,numdepth
    lambdaZ(iz) = float(sum(accum_z(:,iz,:,:)))/float(sum(accum_z(:,:,:,:)))
end do

size_in_bytes_lambda = numdepth*sizeof(lambdaZ(1))

!kij(1:2,1) = (/ ktmp%i, ktmp%j /)

!do i = 2,numk
!    ktmp => ktmp%next
!    kij(1:2,i) = (/ ktmp%i, ktmp%j /)
!end do

open(unit=11,file="test.txt",action="write")

nat = 0
fnat = 1.0/float(sum(cell%numat(1:numset)))
intthick = dble(depthmax)


open(unit=dataunit,file=trim(ecpnl%outname),status='unknown',action='write',form = 'unformatted')
! write the program identifier
write (dataunit) progname
! write the version number
write (dataunit) scversion
! then the name of the crystal data file
write (dataunit) xtalname
! then the name of the corresponding Monte Carlo data file
write (dataunit) ecpnl%energyfile
! energy information and array size
write (dataunit) ecpnl%npx,ecpnl%npx,numset
write (dataunit) EkeV
write (dataunit) ecpnl%dmin
! atom type array for asymmetric unit
write (dataunit) cell%ATOM_type(1:numset)
! is this a regular (square) or hexagonal projection ?
!if (usehex) then
 !   projtype = 'hexago'
 !  write (dataunit) projtype
!else
 !   projtype = 'square'
 !   write (dataunit) projtype
!end if


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
call getenv("CTEMsoft2013opencl",openclpathname)
open(unit = iunit, file = trim(openclpathname)//'/MBmoduleOpenCL.cl', access='direct', status = 'old', &
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
kernel = clCreateKernel(prog, 'CalcLgh', ierr)
call clReleaseProgram(prog, ierr)

allocate(sr(2*ecpnl%npx+1,2*ecpnl%npx+1),stat=istat)
sr = 0.0

allocate(arrsize(npx*npy),arrsizesum(npx*npy),offset(npx*npy),ns(npx*npy),locpix(npx*npy,2),stat=istat)

ktmp => khead
kheadtmp => ktmp

do ii = 1,floor(float(numk)/float(npx*npy))
    arrsize = 0
    arrsizesum = 0
    offset = 0
    ns = 0
    locpix = 0
    do kk = 1,npx
        do ll = 1,npy
            k = ktmp%k
            FN = kk
            locpix((kk-1)*npy+ll,1) = ktmp%i
            locpix((kk-1)*npy+ll,2) = ktmp%j

            call Initialize_ReflectionList(cell, reflist, BetheParameters, FN, k, ecpnl%dmin, nref)

! determine strong and weak reflections
            call Apply_BethePotentials(cell, reflist, firstw, BetheParameters, nref, nns, nnw)
            arrsize((kk-1)*npy+ll) = nns

            if ((kk-1)*npy+ll .lt. npx*npy) then
                arrsizesum((kk-1)*npy+ll+1) = sum(arrsize)
                offset((kk-1)*npy+ll+1) = sum(arrsize**2)
            end if

            ktmp => ktmp%next
        end do
    end do

    ktmp => kheadtmp

    if (allocated(SghCumulative)) deallocate(SghCumulative)
    if (allocated(A)) deallocate(A)

    allocate(SghCumulative(1:offset(npx*npy)+arrsize(npx*npy)**2),A(1:offset(npx*npy)+arrsize(npx*npy)**2),stat=istat)

    do kk = 1,npx
        do ll = 1,npy
                k = ktmp%k
                FN = kk

                call Initialize_ReflectionList(cell, reflist, BetheParameters, FN, k, ecpnl%dmin, nref)

! determine strong and weak reflections
                call Apply_BethePotentials(cell, reflist, firstw, BetheParameters, nref, nns, nnw)

                if (allocated(DynMat)) deallocate(DynMat)

                allocate(DynMat(nns,nns),stat=istat)

                call GetDynMat(cell, reflist, firstw, rlp, DynMat, nns, nnw)

                if (allocated(Sgh)) deallocate(Sgh)
                if (allocated(Sghtmp)) deallocate(Sghtmp)

                allocate(Sghtmp(nns,nns,numset),Sgh(nns,nns))

                Sgh = dcmplx(0.D0,0.D0)
                Sghtmp = dcmplx(0.D0,0.D0)
                nat = 0

                call CalcSgh(cell,reflist,nns,numset,Sghtmp,nat)
! sum Sghtmp over the sites
                Sgh = sum(Sghtmp,3)

                do pp = 1,nns
                    do qq = 1,nns
                        SghCumulative(offset((kk-1)*npy+ll)+(pp-1)*nns+qq) = Sgh(pp,qq)
                    end do
                end do

                DynMat = DynMat*cmplx(0.0,cPi*cell%mLambda)
                DynMat = DynMat*cmplx(eps,0.0)

                absmax = maxval(abs(DynMat))
                e = ceiling(log(absmax)/log(2.0))
                ns((kk-1)*npy+ll) = e+1

                if (ns((kk-1)*npy+ll) .le. 0) then
                    ns((kk-1)*npy+ll) = 1
                else

                DynMat = DynMat/2**ns((kk-1)*npy+ll)

                end if

                do pp = 1,nns
                    do qq = 1,nns
                        A(offset((kk-1)*npy+ll)+(pp-1)*nns+qq) = DynMat(pp,qq)
                    end do
                end do
                ktmp => ktmp%next
            end do
        end do
        kheadtmp => ktmp

! allocate device memory
        size_in_bytes = (offset(npx*npy)+arrsize(npx*npy)**2)*sizeof(A(1))
        size_in_bytes_wavefn = (arrsizesum(npx*npy)+arrsize(npx*npy))*sizeof(A(1))

        if (allocated(LghCumulative)) deallocate(LghCumulative)
        allocate(LghCumulative(1:offset(npx*npy)+arrsize(npx*npy)**2),stat=istat)

        cl_expA = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory for exponential.'

        cl_A = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory for input matrix.'

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

        cl_arrsize = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_dynmatsz, ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory for array size information.'

        cl_offset = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_dynmatsz, ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory for offset information.'

        cl_arrsizesum = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_dynmatsz, ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory for wavefunction coefficient offset information.'

        cl_sqrsize = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_dynmatsz, ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory for wavefunction coefficient offset information.'

        cl_lambdas = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_dynmatsz, ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory for thickness of each pixel.'

! write list of initial arrays to buffer
        call clEnqueueWriteBuffer(command_queue, cl_A, cl_bool(.true.), 0_8, size_in_bytes, A(1), ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot write to buffer.'

        call clEnqueueWriteBuffer(command_queue, cl_arrsize, cl_bool(.true.), 0_8, size_in_bytes_dynmatsz, arrsize(1), ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot write to buffer.'

        call clEnqueueWriteBuffer(command_queue, cl_offset, cl_bool(.true.), 0_8, size_in_bytes_dynmatsz, offset(1), ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot write to buffer.'

        call clEnqueueWriteBuffer(command_queue, cl_arrsizesum, cl_bool(.true.), 0_8, size_in_bytes_dynmatsz, arrsizesum(1), ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot write to buffer.'

        call clEnqueueWriteBuffer(command_queue, cl_sqrsize, cl_bool(.true.), 0_8, size_in_bytes_dynmatsz, ns(1), ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot write to buffer.'

        call clEnqueueWriteBuffer(command_queue, cl_coeff, cl_bool(.true.), 0_8, size_in_bytes_coeff, coef(1,1), ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot write to buffer.'

        call clEnqueueWriteBuffer(command_queue, cl_lambdas, cl_bool(.true.), 0_8, size_in_bytes_lambda, lambdaZ(1), ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot write to buffer.'

! set kernel arguments

        call clSetKernelArg(kernel, 0, cl_expA, ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot set first kernel argument.'

        call clSetKernelArg(kernel, 1, cl_A, ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot set second kernel argument.'

        call clSetKernelArg(kernel, 2, cl_AA, ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot set third kernel argument.'

        call clSetKernelArg(kernel, 3, cl_AAA, ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot set fourth kernel argument.'

        call clSetKernelArg(kernel, 4, cl_arrsize, ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot set fifth kernel argument.'

        call clSetKernelArg(kernel, 5, cl_coeff, ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot set sixth kernel argument.'

        call clSetKernelArg(kernel, 6, cl_T1, ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot set seventh kernel argument.'

        call clSetKernelArg(kernel, 7, cl_T2, ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot set eighth kernel argument.'

        call clSetKernelArg(kernel, 8, cl_wavefncoeff, ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot set nineth kernel argument.'

        call clSetKernelArg(kernel, 9, cl_wavefncoeffintd, ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot set eleventh kernel argument.'

        call clSetKernelArg(kernel, 10, cl_sqrsize, ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot set tenth kernel argument.'

        call clSetKernelArg(kernel, 11, cl_offset, ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot set eleventh kernel argument.'

        call clSetKernelArg(kernel, 12, cl_arrsizesum, ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot set twelveth kernel argument.'

        call clSetKernelArg(kernel, 13, numdepth, ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot set thirteenth kernel argument.'

        call clSetKernelArg(kernel, 14, cl_lambdas, ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot set thirteenth kernel argument.'

!execute the kernel

        call clEnqueueNDRangeKernel(command_queue, kernel, globalsize, localsize, event, ierr)
! wait for the commands to finish

        call clFinish(command_queue, ierr)
! read the resulting vector from device memory

        call clEnqueueReadBuffer(command_queue, cl_T1, cl_bool(.true.), 0_8, size_in_bytes, LghCumulative(1), ierr)

        do pp = 1,npx
            do qq = 1,npy
                ipx = locpix((pp-1)*npy+qq,1)
                ipy = locpix((pp-1)*npy+qq,2)
                size1 = arrsize((pp-1)*npy+qq)
                size2 = offset((pp-1)*npy+qq)

                if (isym.ne.1) then
                    call Apply2DLaueSymmetry(ipx,ipy,isym,iequiv,nequiv)
                    iequiv(1,1:nequiv) = iequiv(1,1:nequiv) + npiximgx + 1
                    iequiv(2,1:nequiv) = iequiv(2,1:nequiv) + npiximgy + 1
                do ip=1,nequiv
                    sr(iequiv(1,ip),iequiv(2,ip)) = real(sum(SghCumulative(size2+1:size2+size1**2)*&
                        LghCumulative(size2+1:size2+size1**2)))/float(sum(nat))
                end do
                else
                    sr(ipx+npiximgx+1,ipy+npiximgy+1) = real(sum(SghCumulative(size2+1:size2+size1**2)*&
                        LghCumulative(size2+1:size2+size1**2)))/float(sum(nat))
                end if

            end do
        end do

        if (mod(ii,4) .eq. 0) then
            write(6,'(A,I8,A)') 'Completed ',(ii*npx*npy),' beams.'
        end if

        call clReleaseMemObject(cl_expA, ierr)
        call clReleaseMemObject(cl_A, ierr)
        call clReleaseMemObject(cl_AA, ierr)
        call clReleaseMemObject(cl_AAA, ierr)
        call clReleaseMemObject(cl_T1, ierr)
        call clReleaseMemObject(cl_T2, ierr)
        call clReleaseMemObject(cl_wavefncoeff, ierr)
        call clReleaseMemObject(cl_wavefncoeffintd, ierr)
        call clReleaseMemObject(cl_coeff, ierr)
        call clReleaseMemObject(cl_arrsize, ierr)
        call clReleaseMemObject(cl_arrsizesum, ierr)
        call clReleaseMemObject(cl_offset, ierr)
        call clReleaseMemObject(cl_lambdas, ierr)

    !end do
end do

write(6,'(A)')'Finished computation...quitting now'

do i=1,2*ecpnl%npx+1
   do j=1,2*ecpnl%npx+1
        write(11, '(F15.6)', advance='no') sr(i,j)
    end do
    write(11, *) ''  ! this gives you the line break
end do

write (dataunit) sr
!deallocate(sr)

close(unit=dataunit,status='keep')
close(unit=11,status='keep')

end subroutine ECmasterpatternOpenCL