program MatrixExpo

use cl
use local

IMPLICIT NONE

complex   :: a
complex,allocatable     :: result(:)
! OpenCL variables
type(cl_platform_id)    :: platform
type(cl_device_id)      :: device
type(cl_context)        :: context
type(cl_command_queue)  :: command_queue
type(cl_program)        :: prog
type(cl_kernel)         :: kernel
type(cl_mem)            :: ret
type(cl_event)          :: event


character(len = 100)    :: info ! info about the GPU
integer(kind=4)         :: num,ierr,irec,iunit
integer(kind=8)         :: size_in_bytes,globalsize(1),localsize(1)
integer, parameter              :: source_length = 10000000
character(len = source_length)  :: source

globalsize = 10
localsize = 1

allocate(result(10))
a = cmplx(1,2)
size_in_bytes = 10*sizeof(a)
write(*,*) a,size_in_bytes

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
open(unit = iunit, file = trim(openclpathname)//'/test.cl', access='direct', status = 'old', &
action = 'read', iostat = ierr, recl = 1)
if (ierr /= 0) stop 'Cannot open file test.cl'

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
kernel = clCreateKernel(prog, 'test', ierr)
!call clReleaseProgram(prog, ierr)

! allocate device memory
ret = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size_in_bytes, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'


call clSetKernelArg(kernel, 0, a, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument.'

call clSetKernelArg(kernel, 1, ret, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument.'

call clEnqueueNDRangeKernel(command_queue, kernel, globalsize, localsize, event, ierr)

call clFinish(command_queue, ierr)

call clEnqueueReadBuffer(command_queue, ret, cl_bool(.true.), 0_8, size_in_bytes, result(1), ierr)

print*,result

end program