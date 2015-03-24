! ###################################################################
! Copyright (c) 2013-2014, Saransh Singh/Carnegie Mellon University
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

program InnerProdOpenCL

use local
use typedefs
use NameListTypedefs
use NameListHandlers
use files
use io

IMPLICIT NONE

character(fnlen)                          :: nmldeffile, progname, progdesc
type(DictIndxOpenCLListType)              :: dictindxnl

nmldeffile = 'EMInnerProdOpenCL.nml'
progname = 'EMInnerProdOpenCL.f90'
progdesc = 'Program to compute inner product of observed and calculated patterns on a GPU'

! deal with the command line arguments, if any
call Interpret_Program_Arguments(nmldeffile,2,(/ 0, 40 /), progname)

! deal with the namelist stuff
call GetDictIndxOpenCLNameList(nmldeffile,dictindxnl)

! print some information
call EMsoft(progname, progdesc)

! perform the zone axis computations
call InnerProdGPU(dictindxnl, progname)


end program InnerProdOpenCL



!--------------------------------------------------------------------------
!
! SUBROUTINE:InnerProdGPU
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief Perform the inner product computations for the dictionary approach
!
!> @param dictindxnl dictionary indexing namelist pointer
!> @param progname name of the program
!
!> @date 12/09/14  SS 1.0 original
!--------------------------------------------------------------------------

subroutine InnerProdGPU(dictindxnl,progname)!(result,exp,dict,Ne,Nd,L)

use local
use NameListTypedefs
use NameListHandlers
use typedefs
use files
use io
use cl

IMPLICIT NONE

type(DictIndxOpenCLListType),INTENT(IN)             :: dictindxnl
character(fnlen),INTENT(IN)                         :: progname

real(kind=4),allocatable                            :: result(:)
real(kind=4),allocatable                            :: expt(:)
real(kind=4),allocatable                            :: dict(:)
real(kind=4),allocatable                            :: eulerangles(:,:),image_lin_flt(:)
integer(kind=1),allocatable                         :: image(:,:),image_lin(:),header(:)
integer(kind=4)                                     :: Ne
integer(kind=4)                                     :: Nd
integer(kind=4)                                     :: L
integer(kind=irg)                                   :: totnumdict,numdictsingle,numexptsingle,imght,imgwd,totnumexpt
integer(kind=4)                                     :: recordsize,filesize

type(cl_platform_id)                                :: platform
type(cl_device_id)                                  :: device
type(cl_context)                                    :: context
type(cl_command_queue)                              :: command_queue
type(cl_program)                                    :: prog
type(cl_kernel)                                     :: kernel
type(cl_mem)                                        :: cl_expt,cl_dict,cl_result,cl_As,cl_Bs
type(cl_event)                                      :: event

integer(kind=4),parameter                           :: iunit = 10
integer(kind=4),parameter                           :: iunitdict = 11
integer(kind=4),parameter                           :: iunitexpt = 12
character(len = 100)                                :: info ! info about the GPU
integer(kind=8)                                     :: globalsize(2),localsize(2)
integer, parameter                                  :: source_length = 10000000
character(len = source_length)                      :: source
integer(kind=4)                                     :: num,ierr,istat,irec,Wexp,Wdict,i,j,ii,jj,kk,ll,mm,nn,pp,qq
integer(kind=8)                                     :: size_in_bytes_expt,size_in_bytes_dict,size_in_bytes_result
real(kind=sgl),allocatable                          :: exptsorted(:,:),topk(:,:,:)


Ne = dictindxnl%numexptsingle
Nd = dictindxnl%numdictsingle
L = (dictindxnl%imght*dictindxnl%imgwd)
totnumdict = dictindxnl%totnumdict
totnumexpt = dictindxnl%totnumexpt
numdictsingle = dictindxnl%numdictsingle
numexptsingle = dictindxnl%numexptsingle
imght = dictindxnl%imght
imgwd = dictindxnl%imgwd
recordsize = 25 + L
filesize = 948633600

size_in_bytes_expt = L*Ne*sizeof(L)
size_in_bytes_dict = L*Nd*sizeof(L)
size_in_bytes_result = Ne*Nd*sizeof(result(1))!Ne*Nd*sizeof(L)
Wexp = L
Wdict = Nd
localsize = (/16,16/)
globalsize = (/Ne,Nd/)

allocate(image(imgwd,imght),image_lin(1:L),image_lin_flt(1:L),header(25),stat=istat)
image = 0
image_lin = 0
image_lin_flt = 0.0

allocate(eulerangles(totnumdict,3),stat=istat)
eulerangles = 0.0

!=====================================
! I/O FOR EULER ANGLE FILE
!=====================================

call Message('opening '//trim(dictindxnl%eulerfile), frm = "(A)" )
open(dataunit,file=trim(dictindxnl%eulerfile),status='old',form='formatted',action='read',iostat=ierr)

read(dataunit,'(I15)') totnumdict
do ii = 1,totnumdict
    read(dataunit,*) eulerangles(ii,1:3)
end do
close(dataunit,status='keep')
call Message(' -> completed reading '//trim(dictindxnl%eulerfile), frm = "(A)")

!=====================================
! OPENING DICT AND EXPT FOR I/O
!=====================================


call Message('opening '//trim(dictindxnl%dictfile), frm = "(A)" )
open(unit=iunitdict,file=trim(dictindxnl%dictfile),status='old',form='unformatted',access='stream',iostat=ierr)

call Message('opening '//trim(dictindxnl%exptfile), frm = "(A)" )
open(unit=iunitexpt,file=trim(dictindxnl%exptfile),status='old',form='unformatted',access='stream',iostat=ierr)


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
call getenv("EMsoftopencl",openclpathname)
open(unit = iunit, file = trim(openclpathname)//'/DictIndx.cl', access='direct', status = 'old', &
action = 'read', iostat = ierr, recl = 1)
if (ierr /= 0) stop 'Cannot open file InnerProd.cl'

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
call clBuildProgram(prog, '-cl-fast-relaxed-math', ierr)

! get the compilation log
call clGetProgramBuildInfo(prog, device, CL_PROGRAM_BUILD_LOG, source, irec)
if(len(trim(source)) > 0) print*, trim(source)

if(ierr /= CL_SUCCESS) stop 'Error: program build failed.'

write(6,*) "Kernel Build Successful...."

! finally get the kernel and release the program
kernel = clCreateKernel(prog, 'InnerProd', ierr)
!call clReleaseProgram(prog, ierr)

! allocate device memory
cl_expt = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_expt, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

cl_dict = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_dict, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

cl_result = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_result, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

!cl_As = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_local, ierr)
!if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

!cl_Bs = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_local, ierr)
!if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

!=====================================
! MAIN LOOP OVER ALL PATTERNS
!=====================================

experimentalloop: do ll = 1,nint(float(totnumexpt)/float(numexptsingle))

    if (allocated(expt)) deallocate(expt)
    allocate(expt(Ne*L),stat=istat)
    expt = 0.0

    if (allocated(exptsorted)) deallocate(exptsorted)
    allocate(exptsorted(numexptsingle,totnumdict),stat=istat)
    exptsorted = 0.0

    if (allocated(topk)) deallocate(topk)
    allocate(topk(numexptsingle,100,2),stat=istat)
    do ii = 1,numexptsingle
        read(iunitexpt) header
        read(iunitexpt) image
        do jj = 1,imght
            image_lin((jj-1)*imght+1:jj*imght) = image(:,jj)
        end do
        image_lin_flt = float(image_lin)
        do mm = 1,L
            if (image_lin_flt(mm) .lt. 0) image_lin_flt(mm) = image_lin_flt(mm) + 255.0
        end do
        expt((ii-1)*L+1:ii*L) = image_lin_flt(1:L)/sqrt(sum(image_lin_flt**2))
    end do

    call clEnqueueWriteBuffer(command_queue, cl_expt, cl_bool(.true.), 0_8, size_in_bytes_expt, expt(1), ierr)
    if(ierr /= CL_SUCCESS) stop 'Error: cannot write to buffer.'

    dictionaryloop: do kk = 1,2!floor(float(totnumdict)/float(numdictsingle))

        if (allocated(result)) deallocate(result)
        if (allocated(dict)) deallocate(dict)

        allocate(result(Ne*Nd),dict(Nd*L),stat=istat)

        result = 0.0
        dict = 0.0

        do ii = 1,numdictsingle
            read(iunitdict) image
            do jj = 1,imght
                image_lin((jj-1)*imght+1:jj*imght) = image(:,jj)
            end do
            image_lin_flt = float(image_lin)
            do mm = 1,L
                if (image_lin_flt(mm) .lt. 0) image_lin_flt(mm) = image_lin_flt(mm) + 255.0
            end do
            dict((ii-1)*L+1:ii*L) = image_lin_flt(1:L)/sqrt(sum(image_lin_flt**2))
        end do
!if (kk .eq. 2) then
!open(unit=13,file='test.txt',action='write')
!do i=1,imgwd
!do j=1,imght
!write(13, '(I4)', advance='no') image(j,i)
!end do
!write(13, *) ''  ! this gives you the line break
!end do
!close(13)
!end if

        call clEnqueueWriteBuffer(command_queue, cl_dict, cl_bool(.true.), 0_8, size_in_bytes_dict, dict(1), ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot write to buffer.'

! set kernel argument
        call clSetKernelArg(kernel, 0, cl_expt, ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot set first kernel argument.'

        call clSetKernelArg(kernel, 1, cl_dict, ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot set second kernel argument.'

        call clSetKernelArg(kernel, 2, Wexp, ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot set third kernel argument.'

        call clSetKernelArg(kernel, 3, Wdict, ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot set fourth kernel argument.'

        call clSetKernelArg(kernel, 4, cl_result, ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot set fifth kernel argument.'
!execute the kernel
        call clEnqueueNDRangeKernel(command_queue, kernel, globalsize, localsize, event, ierr)

! wait for the commands to finish
        call clFinish(command_queue, ierr)
! read the resulting vector from device memory
        call clEnqueueReadBuffer(command_queue, cl_result, cl_bool(.true.), 0_8, size_in_bytes_result, result(1), ierr)
! copying the dot products to an array for sorting later
        do i = 1,numexptsingle
            exptsorted(i,numdictsingle*(kk-1)+1:numdictsingle*kk) = result((i-1)*numdictsingle+1:i*numdictsingle)
        end do

        if (mod(kk,1) .eq. 0) then
            write(6,'(A26,I7,A26,I7)'),'Completed dot product of',kk*1024,'dictionary patterns with',ll*1024,'experimental patterns'
        end if

    end do dictionaryloop

    if (allocated(result)) deallocate(result)
    if (allocated(dict)) deallocate(dict)

    allocate(result(Ne*Nd),dict(Nd*L),stat=istat)

    result = 0.0
    dict = 0.0

! do the leftover dot products
if (1 .eq. 0) then
    do ii = 1,MODULO(totnumdict,numdictsingle)
        read(iunitdict) image
        do jj = 1,imght
            image_lin((jj-1)*imght+1:jj*imght) = image(:,jj)
        end do
        image_lin_flt = float(image_lin)
        do mm = 1,L
            if (image_lin_flt(mm) .lt. 0) image_lin_flt(mm) = image_lin_flt(mm) + 255.0
        end do
        dict((ii-1)*L+1:ii*L) = image_lin_flt(1:L)/sqrt(sum(image_lin_flt**2))
    end do

    call clEnqueueWriteBuffer(command_queue, cl_dict, cl_bool(.true.), 0_8, size_in_bytes_dict, dict(1), ierr)
    if(ierr /= CL_SUCCESS) stop 'Error: cannot write to buffer.'

! set kernel argument
    call clSetKernelArg(kernel, 0, cl_expt, ierr)
    if(ierr /= CL_SUCCESS) stop 'Error: cannot set first kernel argument.'

    call clSetKernelArg(kernel, 1, cl_dict, ierr)
    if(ierr /= CL_SUCCESS) stop 'Error: cannot set second kernel argument.'

    call clSetKernelArg(kernel, 2, Wexp, ierr)
    if(ierr /= CL_SUCCESS) stop 'Error: cannot set third kernel argument.'

    call clSetKernelArg(kernel, 3, Wdict, ierr)
    if(ierr /= CL_SUCCESS) stop 'Error: cannot set fourth kernel argument.'

    call clSetKernelArg(kernel, 4, cl_result, ierr)
    if(ierr /= CL_SUCCESS) stop 'Error: cannot set fifth kernel argument.'
!execute the kernel
    call clEnqueueNDRangeKernel(command_queue, kernel, globalsize, localsize, event, ierr)

! wait for the commands to finish
    call clFinish(command_queue, ierr)
! read the resulting vector from device memory
    call clEnqueueReadBuffer(command_queue, cl_result, cl_bool(.true.), 0_8, size_in_bytes_result, result(1), ierr)
! adding dot products of leftover dictionary patterns
    do i = 1,numexptsingle
        exptsorted(i,numdictsingle*numexptsingle+1:numdictsingle*(numexptsingle+1)) = result((i-1)*numdictsingle+1:i*numdictsingle)
    end do
end if

    !sort the exptsorted array and return the top 100 or so values
    ! code goes here
        do pp = 1,100
            topk(1:numexptsingle,pp,1) = maxloc(exptsorted,2)
            topk(1:numexptsingle,pp,2) = maxval(exptsorted,2)
            do nn = 1,numexptsingle
                exptsorted(nn,int(topk(nn,pp,1))) = 0.0
            end do
        end do
print*,topk(1,:,1)

end do experimentalloop

print*,result(1)
call clReleaseKernel(kernel, ierr)
call clReleaseCommandQueue(command_queue, ierr)
call clReleaseContext(context, ierr)
call clReleaseMemObject(cl_expt, ierr)
call clReleaseMemObject(cl_dict, ierr)
call clReleaseMemObject(cl_result, ierr)

call Message('Finished reading '//trim(dictindxnl%dictfile), frm = "(A)" )
call Message('Finished '//trim(dictindxnl%exptfile), frm = "(A)" )

close(iunitdict,status='keep')
close(iunitexpt,status='keep')



end subroutine InnerProdGPU

!--------------------------------------------------------------------------
!
! SUBROUTINE:Sort
!
!> @author Taken from rosettacode http://rosettacode.org/wiki/Sorting_algorithms/Insertion_sort#Fortran
!
!> @brief sort the array
!
!> @param array input array
!> @param size size of array
!
!--------------------------------------------------------------------------

