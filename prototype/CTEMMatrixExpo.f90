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
! PROGRAM: CTEMMatrixExpo
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief Host code for matrix exponential on a GPU
!
!> @date 10/08/14 MDG 1.0 original
!--------------------------------------------------------------------------

program CTEMMatrixExpo

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

real(kind=sgl)                  :: dmin,EkeV


type(unitcell), pointer         :: cell
type(gnode)                     :: rlp
type(DynType)                   :: Dyn
integer(kind=4)                 :: num,ierr,irec,iunit
integer(kind=4),allocatable     :: DynMatSz(:),ns(:),result1(:)
integer(kind=8)                 :: size_in_bytes,globalsize(2),localsize(2),&
                                   size_in_bytes_coeff,size_in_bytes_param
character(fnlen)                :: xtalname
logical                         :: verbose,usehex,switchmirror

complex(kind=dbl),allocatable   :: DynMat(:,:),result(:),A(:),DynMat1D(:)
complex(kind=sgl)               :: coefficient(9),cvals(9)
type(reflisttype),pointer       :: reflist, firstw
type(BetheParameterType)        :: BetheParameters
real(kind=sgl)                  :: FN(3),kk(3)
integer(kind=irg)               :: nref,nns,nnw
integer(kind=irg),parameter     :: maxDynMatsz = 55

type(kvectorlist),pointer       :: khead,ktmp
integer(kind=irg)               :: numk,npx,ijmax,isym,pgnum,i,j,npyhex
integer(kind=irg),parameter     :: LaueTest(11) = (/ 149, 151, 153, 156, 158, 160, 161, 164, 165, 166, 167 /)  ! space groups with 2 or mirror at 30 degrees

! OpenCL variables
type(cl_platform_id)            :: platform
type(cl_device_id)              :: device
type(cl_context)                :: context
type(cl_command_queue)          :: command_queue
type(cl_program)                :: prog
type(cl_kernel)                 :: kernel
type(cl_mem)                    :: cl_expA,cl_A,cl_AA,cl_AAA,cl_coeff,cl_T1,cl_T2,&
                                   cl_T3, cl_T1T2,cl_T1T2T3,cl_matsz,cl_ns
type(cl_event)                  :: event


character(fnlen)                :: info ! info about the GPU
integer, parameter              :: source_length = 10000000
character(len = source_length)  :: source


dmin = 0.05
xtalname = 'Ni.xtal'

globalsize = (/32,32/)
localsize = (/4,4/)

EkeV = 30.0
npx = 512
ijmax = npx**2
size_in_bytes = globalsize(1)*globalsize(2)*sizeof(DynMat(1,1))*(maxDynMatsz**2)
size_in_bytes_coeff = 9*sizeof(DynMat(1,1))
size_in_bytes_param = globalsize(1)*globalsize(2)*sizeof(ns(1))
allocate(result1(globalsize(1)*globalsize(2)))

allocate(result(globalsize(1)*globalsize(2)*maxDynMatsz**2),A(globalsize(1)*globalsize(2)*maxDynMatsz**2))
result = cmplx(0.0,0.0)
A = cmplx(0.0,0.0)
allocate(ns(globalsize(1)*globalsize(2)),DynMatSz(globalsize(1)*globalsize(2)))

!=============================================
! initialize the cell
!=============================================
nullify(cell)
allocate(cell)

! load the crystal structure and compute the Fourier coefficient lookup table
verbose = .TRUE.

call Initialize_Cell(cell,Dyn,rlp, xtalname, dmin, sngl(1000.0*EkeV),verbose)

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

if(usehex)  npyhex = nint(2.0*float(npx)/sqrt(3.0))



nullify(khead)
call Calckvectors(khead,cell, (/ 0.D0, 0.D0, 1.D0 /), (/ 0.D0, 0.D0, 0.D0 /),0.D0,npx,npx,numk, &
isym,ijmax,'RoscaLambert',usehex)

call Set_Bethe_Parameters(BetheParameters,.TRUE.)
nullify(reflist)
nullify(firstw)
ktmp => khead


beamloop: do i = 1, 1024
    kk = ktmp%k(1:3)
    FN = kk

    call Initialize_ReflectionList(cell, reflist, BetheParameters, FN, kk, dmin, nref)

! determine strong and weak reflections
    call Apply_BethePotentials(cell, reflist, firstw, BetheParameters, nref, nns, nnw)

    allocate(DynMat(nns,nns))
    allocate(DynMat1D(nns**2))
    call GetDynMat(cell, reflist, firstw, rlp, DynMat, nns, nnw)
    DynMatSz(i) = nns
    ns(i) = ceiling(log(maxval(abs(DynMat)))/log(2.0)) + 1
    DynMat1D = Reshape(DynMat,(/nns**2/))
    DynMat1D = DynMat1D/(2**ns(i))
    !if (i .eq. 1) print*,DynMat1D
    A((i-1)*maxDynMatsz+1:(i-1)*maxDynMatsz+nns**2) = DynMat1D(1:nns**2)
    ktmp => ktmp%next
    deallocate(DynMat)
    deallocate(DynMat1D)
end do beamloop


cvals(1) = cmplx(-3.3335514852690488032942739163345055,0.d0)
cvals(2) = cmplx(-3.0386480729366970892124687564926859,-1.5868011957588383288038677051222921)
cvals(3) = cmplx(-3.0386480729366970892124687564926859,+1.5868011957588383288038677051222921)
cvals(4) = cmplx(-2.1108398003026547374987047865183922,-3.0899109287255009227777015426228801)
cvals(5) = cmplx(-2.1108398003026547374987047865183922,+3.0899109287255009227777015426228801)
cvals(6) = cmplx(-0.38106984566311299903129424501333242,-4.3846445331453979503692027283066828)
cvals(7) = cmplx(-0.38106984566311299903129424501333242,+4.3846445331453979503692027283066828)
cvals(8) = cmplx(2.6973334615369892273896047461916633,-5.1841620626494141778340870727109629)
cvals(9) = cmplx(2.6973334615369892273896047461916633,+5.1841620626494141778340870727109629)

coefficient(1) = cvals(1)+cvals(2)+cvals(3)
coefficient(2) = cvals(1)*cvals(2)+cvals(2)*cvals(3)+cvals(3)*cvals(1)
coefficient(3) = cvals(1)*cvals(2)*cvals(3)

coefficient(4) = cvals(4)+cvals(5)+cvals(6)
coefficient(5) = cvals(4)*cvals(5)+cvals(5)*cvals(6)+cvals(6)*cvals(4)
coefficient(6) = cvals(4)*cvals(5)*cvals(6)

coefficient(7) = cvals(7)+cvals(8)+cvals(9)
coefficient(8) = cvals(7)*cvals(8)+cvals(8)*cvals(9)+cvals(9)*cvals(7)
coefficient(9) = cvals(7)*cvals(8)*cvals(9)

!=====================
! INITIALIZATION
!=====================

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
open(unit = iunit, file = trim(openclpathname)//'/CTEMMatExp.cl', access='direct', status = 'old', &
action = 'read', iostat = ierr, recl = 1)
if (ierr /= 0) stop 'Cannot open file CTEMMatExp.cl'

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

write(*,*) "Kernel Build Successful...."

! finally get the kernel and release the program
kernel = clCreateKernel(prog, 'MatExpImg', ierr)
!call clReleaseProgram(prog, ierr)

! allocate device memory
cl_expA = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

cl_A = clCreateBuffer(context, CL_MEM_READ_ONLY, size_in_bytes, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

cl_AA = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

cl_AAA = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

cl_coeff = clCreateBuffer(context, CL_MEM_READ_ONLY, size_in_bytes, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

cl_T1 = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

cl_T2 = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

cl_T3 = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

cl_T1T2 = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

cl_matsz = clCreateBuffer(context, CL_MEM_READ_ONLY, size_in_bytes_param, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

cl_ns = clCreateBuffer(context, CL_MEM_READ_ONLY, size_in_bytes_param, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'


call clEnqueueWriteBuffer(command_queue, cl_coeff, cl_bool(.true.), 0_8, size_in_bytes_coeff, coefficient(1), ierr)

call clEnqueueWriteBuffer(command_queue, cl_A, cl_bool(.true.), 0_8, size_in_bytes, A(1), ierr)

call clEnqueueWriteBuffer(command_queue, cl_matsz, cl_bool(.true.), 0_8, size_in_bytes_param, DynMatSz(1), ierr)

call clEnqueueWriteBuffer(command_queue, cl_ns, cl_bool(.true.), 0_8, size_in_bytes_param, ns(1), ierr)


call clSetKernelArg(kernel, 0, cl_expA, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument 0.'

call clSetKernelArg(kernel, 1, cl_A, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument 1.'

call clSetKernelArg(kernel, 2, cl_AA, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument 2.'

call clSetKernelArg(kernel, 3, cl_AAA, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument 3.'

call clSetKernelArg(kernel, 4, cl_coeff, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument 5.'

call clSetKernelArg(kernel, 5, cl_T1, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument 6.'

call clSetKernelArg(kernel, 6, cl_T2, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument 7.'

call clSetKernelArg(kernel, 7, cl_T3, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument 8.'

call clSetKernelArg(kernel, 8, cl_T1T2, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument 9.'

call clSetKernelArg(kernel, 9, cl_ns, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument 10.'

call clSetKernelArg(kernel, 10,cl_matsz , ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument 11.'

call clEnqueueNDRangeKernel(command_queue, kernel, globalsize, localsize, event, ierr)

call clFinish(command_queue, ierr)

call clEnqueueReadBuffer(command_queue, cl_expA, cl_bool(.true.), 0_8, size_in_bytes, result(1), ierr)

print*,result(1:100)

call clReleaseMemObject(cl_expA,ierr)
call clReleaseMemObject(cl_A,ierr)
call clReleaseMemObject(cl_AA,ierr)
call clReleaseMemObject(cl_AAA,ierr)
call clReleaseMemObject(cl_T1,ierr)
call clReleaseMemObject(cl_T2,ierr)
call clReleaseMemObject(cl_T3,ierr)
call clReleaseMemObject(cl_coeff,ierr)
call clReleaseKernel(kernel,ierr)
call clReleaseCommandQueue(command_queue,ierr)
call clReleaseContext(context,ierr)

end program