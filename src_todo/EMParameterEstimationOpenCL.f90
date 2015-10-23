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

program EMDictIndxOpenCL

use math
use typedefs
use dictmod
use cl

IMPLICIT NONE


! OpenCL variables
type(cl_platform_id)            :: platform
type(cl_device_id)              :: device
type(cl_context)                :: context
type(cl_command_queue)          :: command_queue
type(cl_program)                :: prog
type(cl_kernel)                 :: kernel
type(cl_mem)                    :: cl_quaternion,cl_mean,cl_ConcParam,cl_symmetry,cl_lookup,cl_seeds
type(cl_event)                  :: event

character(len = 100)            :: info ! info about the GPU
integer, parameter              :: iunit = 10
integer, parameter              :: source_length = 10000000
character(len = source_length)  :: source

integer(kind=4)                 :: num,ierr,istat,irec,numk,numsym,count,numinit
integer(kind=4)                 :: npx,npy
integer(kind=8)                 :: size_in_bytes_quat,size_in_bytes_mean,size_in_bytes_ConcParam,size_in_bytes_sym,&
    size_in_bytes_lookup,size_in_bytes_seeds
integer(kind=8)                 :: globalsize(2),localsize(2)

real(kind=sgl),allocatable      :: meanres(:),kappares(:),quat_rand(:),symmetry(:),quat_sample_lin(:),lookup(:)
real(kind=dbl),allocatable      :: quat_sample(:,:)
integer(kind=irg)               :: nseeds,ii,jj,kk,seed
integer(kind=irg),allocatable   :: init_seeds(:),rnseeds(:)

numk = 40
numsym = 48
npx = 32
npy = 32
seed = 12345
numinit = 3

size_in_bytes_quat = 4*npx*npy*numk*sizeof(npx)
size_in_bytes_mean = 4*npx*npy*sizeof(npx)
size_in_bytes_ConcParam = npx*npy*sizeof(npx)
size_in_bytes_sym = 4*numsym*sizeof(npx)
size_in_bytes_lookup = 590*sizeof(npx)
size_in_bytes_seeds = 4*npx*npy*sizeof(nseeds)

globalsize = (/npx,npy/)
localsize = (/16,16/)

allocate(meanres(4*npx*npy),kappares(npx*npy),stat=istat)
allocate(quat_rand(4*npx*npy*numk),symmetry(4*numsym),stat=istat)
allocate(quat_sample(4,numk),quat_sample_lin(4*numk),stat=istat)
allocate(lookup(590),stat=istat)

do ii = 0,numsym/2 - 1
    symmetry(4*ii+1:4*ii+4) = SYM_Qsymop(1:4,ii+1)
    symmetry(4*ii+97:4*ii+100) = -SYM_Qsymop(1:4,ii+1)
end do

do ii = 1,590
    lookup(ii) = BesselIn(0.05D0*ii,2)/BesselI1(0.05D0*ii)
end do

!===================================
! INITIALIZE QUATERNION ARRAY
!===================================

quat_sample = DI_SampleVMF(numk,seed,(/1.D0,0.D0,0.D0,0.D0/),200.D0)

do ii = 0,numk-1
    quat_sample_lin(4*ii+1:4*ii+4) = quat_sample(:,ii+1)
end do

!call RANDOM_NUMBER(quat_rand)

do ii = 0,npx*npy-1
    quat_rand(4*numk*ii+1:4*numk*ii+4*numk) = quat_sample_lin(:)
end do
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
open(unit = iunit, file = trim(openclpathname)//'DictIndx.cl', access='direct', status = 'old', &
action = 'read', iostat = ierr, recl = 1)
if (ierr /= 0) stop 'Cannot open file DictIndx.cl'

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
kernel = clCreateKernel(prog, 'ParamEstm', ierr)
call clReleaseProgram(prog, ierr)

open(unit = iunit, file = 'RandomSeeds.data', form='unformatted', status='old')
read(iunit) nseeds
allocate(rnseeds(nseeds))
read(iunit) rnseeds
close(unit=iunit,status='keep')

if (4*npx*npy .gt. nseeds) Stop 'Insufficient number of primes in the list'

allocate(init_seeds(4*npx*npy),stat=istat)
init_seeds = 0
do ii = 1,npy
    do jj = 1,npx
        do kk = 1,4
            init_seeds(4*((ii-1)*npx+(jj-1))+kk) = rnseeds(4*((ii-1)*npx+(jj-1))+kk)
        end do
    end do
end do


! create device memory buffers
cl_quaternion = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_quat, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

cl_mean = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_mean, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

cl_symmetry = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_sym, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

cl_ConcParam = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_ConcParam, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

cl_lookup = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_lookup, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

cl_seeds = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_seeds, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

! write the set of quaternions for all the pixels in the image
call clEnqueueWriteBuffer(command_queue, cl_quaternion, cl_bool(.true.), 0_8, size_in_bytes_quat, quat_rand(1), ierr)

call clEnqueueWriteBuffer(command_queue, cl_symmetry, cl_bool(.true.), 0_8, size_in_bytes_sym, symmetry(1), ierr)

call clEnqueueWriteBuffer(command_queue, cl_lookup, cl_bool(.true.), 0_8, size_in_bytes_lookup, lookup(1), ierr)

call clEnqueueWriteBuffer(command_queue, cl_seeds, cl_bool(.true.), 0_8, size_in_bytes_seeds, init_seeds(1), ierr)


! set the kernel arguments
call clSetKernelArg(kernel, 0, cl_quaternion, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot set first kernel argument.'

call clSetKernelArg(kernel, 1, cl_mean, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot set second kernel argument.'

call clSetKernelArg(kernel, 2, cl_ConcParam, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot set third kernel argument.'

call clSetKernelArg(kernel, 3, cl_symmetry, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot set fourth kernel argument.'

call clSetKernelArg(kernel, 4, numk, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot set fifth kernel argument.'

call clSetKernelArg(kernel, 5, numsym, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot set sixth kernel argument.'

call clSetKernelArg(kernel, 6, cl_lookup, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot set sixth kernel argument.'

call clSetKernelArg(kernel, 7, cl_seeds, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot set sixth kernel argument.'

call clSetKernelArg(kernel, 8, numinit, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot set sixth kernel argument.'

! execute the kernel
call clEnqueueNDRangeKernel(command_queue, kernel, globalsize, localsize, event, ierr)

! wait for the commands to finish
call clFinish(command_queue, ierr)

call clEnqueueReadBuffer(command_queue, cl_mean, cl_bool(.true.), 0_8, size_in_bytes_mean, meanres(1), ierr)

call clEnqueueReadBuffer(command_queue, cl_ConcParam, cl_bool(.true.), 0_8, size_in_bytes_ConcParam, kappares(1), ierr)

count = 0
do ii = 1,npx*npy
    if (kappares(ii) .eq. 0.0) count = count + 1
end do

!print*,count,npx*npy,float(npx*npy)/float(count)
print*,kappares
!print*,meanres
!print*,quat_rand(1:4)
call clReleaseKernel(kernel, ierr)
call clReleaseCommandQueue(command_queue, ierr)
call clReleaseContext(context, ierr)
call clReleaseMemObject(cl_quaternion, ierr)
call clReleaseMemObject(cl_mean, ierr)
call clReleaseMemObject(cl_ConcParam, ierr)
call clReleaseMemObject(cl_symmetry, ierr)
call clReleaseMemObject(cl_lookup, ierr)


end program EMDictIndxOpenCL
