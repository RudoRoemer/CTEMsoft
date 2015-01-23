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

nmldeffile = 'CTEMInnerProdOpenCL.nml'
progname = 'CTEMInnerProdOpenCL.f90'
progdesc = 'Program to compute inner product of observed and calculated patterns on a GPU'

! deal with the command line arguments, if any
call Interpret_Program_Arguments(nmldeffile,2,(/ 0, 40 /), progname)

! deal with the namelist stuff
call GetDictIndxOpenCLNameList(nmldeffile,dictindxnl)

! print some information
call CTEMsoft(progname, progdesc)

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
real(kind=4),allocatable                            :: dict(:),dicttranspose(:)
real(kind=4),allocatable                            :: eulerangles(:,:),imagedictflt(:),imageexptflt(:),meandict(:),meanexpt(:)
integer(kind=1),allocatable                         :: imagedict(:),imageexpt(:)
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

integer(kind=4),parameter                           :: iunit = 40
integer(kind=4),parameter                           :: iunitdict = 41
integer(kind=4),parameter                           :: iunitexpt = 42
character(len = 100)                                :: info ! info about the GPU
integer(kind=8)                                     :: globalsize(2),localsize(2)
integer, parameter                                  :: source_length = 10000000
character(len = source_length)                      :: source
integer(kind=4)                                     :: num,ierr,istat,irec,Wexp,Wdict,i,j,ii,jj,kk,ll,mm,nn,pp,qq
integer(kind=8)                                     :: size_in_bytes_expt,size_in_bytes_dict,size_in_bytes_result
real(kind=sgl),allocatable                          :: exptsorted(:,:),topk(:,:,:)
integer(kind=irg),allocatable                       :: exptindex(:,:),indexlist(:)
!real(kind=4),allocatable                            :: exptsr(:,:),dictsr(:,:),res(:,:)
integer(kind=irg)                                   :: t1,t2,clock_rate,clock_max


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

allocate(imagedict(imgwd*imght),imageexpt(recordsize),stat=istat)
imagedict = 0
imageexpt = 0

allocate(imagedictflt(imgwd*imght),imageexptflt(imgwd*imght),stat=istat)
imagedict = 0.0
imageexpt = 0.0

allocate(meandict(imgwd*imght),meanexpt(imgwd*imght),stat=istat)
meanexpt = 0.0
meandict = 0.0

allocate(eulerangles(totnumdict,3),stat=istat)
eulerangles = 0.0

allocate(indexlist(1:totnumdict),stat=istat)
indexlist = 0

do ii = 1,totnumdict
    indexlist(ii) = ii
end do

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

print*,trim(dictindxnl%exptfile)
call Message('opening '//trim(dictindxnl%dictfile), frm = "(A)" )
open(unit=iunitdict,file=trim(dictindxnl%dictfile),status='old',form='unformatted',access='direct',recl=imght*imgwd,iostat=ierr)

call Message('opening '//trim(dictindxnl%exptfile), frm = "(A)" )
open(unit=iunitexpt,file=trim(dictindxnl%exptfile),status='old',form='unformatted',access='direct',recl=recordsize,iostat=ierr)

call Message('Start calculating mean pattern for dictionary '//trim(dictindxnl%dictfile), frm = "(A)" )
!if (1 .eq. 0) then
do ii = 1,totnumdict
    read(iunitdict,rec=ii) imagedict
    imagedictflt = float(imagedict)
    do jj = 1,L
        if (imagedictflt(jj) .lt. 0) imagedictflt(jj) = imagedictflt(jj)+255.0
    end do
    meandict = meandict+imagedictflt
end do
meandict = meandict/totnumdict


call Message('Finished calculating mean pattern for dictionary '//trim(dictindxnl%dictfile), frm = "(A)" )

call Message('Start calculating mean pattern for experiments '//trim(dictindxnl%exptfile), frm = "(A)" )

do ii = 1,totnumexpt
    read(iunitexpt,rec=ii) imageexpt
    imageexptflt = float(imageexpt(26:recordsize))
    do jj = 1,L
        if (imageexptflt(jj) .lt. 0) imageexptflt(jj) = imageexptflt(jj)+255.0
    end do
    meanexpt = meanexpt+imageexptflt
end do
meanexpt = meanexpt/totnumexpt


call Message('Finished calculating mean pattern for experiments '//trim(dictindxnl%exptfile), frm = "(A)" )
!end if
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
open(unit = iunit, file = trim(openclpathname)//'/DictIndx.cl', access='direct', status = 'old', &
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

!allocate(exptsr(Ne,L),dictsr(Nd,L),res(Ne,Nd),stat=istat)
!exptsr = 0.0
!dictsr = 0.0
!res = 0.0


experimentalloop: do ll = 1,nint(float(totnumexpt)/float(numexptsingle))

    if (allocated(expt)) deallocate(expt)
    allocate(expt(Ne*L),stat=istat)
    expt = 0.0

    if (allocated(exptsorted)) deallocate(exptsorted)
    allocate(exptsorted(numexptsingle,totnumdict),stat=istat)
    exptsorted = 0.0

    if (allocated(exptindex)) deallocate(exptindex)
    allocate(exptindex(numexptsingle,totnumdict),stat=istat)
    exptindex = 0

    if (allocated(topk)) deallocate(topk)
    allocate(topk(numexptsingle,100,2),stat=istat)
    do ii = 1,numexptsingle
        exptindex(ii,:) = indexlist(:)
        read(iunitexpt,rec=(ll-1)*numexptsingle+ii) imageexpt
        imageexptflt = float(imageexpt(26:recordsize))
        do mm = 1,L
            if (imageexptflt(mm) .lt. 0) imageexptflt(mm) = imageexptflt(mm) + 255.0
        end do
        expt((ii-1)*L+1:ii*L) = (imageexptflt(1:L)-meanexpt(1:L))/sqrt(sum((imageexptflt(1:L)-meanexpt(1:L))**2))
    end do

    call clEnqueueWriteBuffer(command_queue, cl_expt, cl_bool(.true.), 0_8, size_in_bytes_expt, expt(1), ierr)
    if(ierr /= CL_SUCCESS) stop 'Error: cannot write to buffer.'

!if (ll .eq. 1) then
!do ii = 1,numexptsingle
!exptsr(ii,1:L) = expt((ii-1)*L+1:ii*L)
!end do
!end if
!if (ll .eq. 1) then
!imageexptflt = expt(616*L+1:617*L)
!open(unit=13,file='testexpt.txt',action='write')
!do i=1,imgwd
!do j=1,imght
!write(13, '(F15.6)', advance='no') imageexptflt((i-1)*imght+j)
!end do
!write(13, *) ''  ! this gives you the line break
!end do
!close(13)
!end if
    dictionaryloop: do kk = 1,floor(float(totnumdict)/float(numdictsingle))

        if (allocated(result)) deallocate(result)
        if (allocated(dict)) deallocate(dict)
        if (allocated(dicttranspose)) deallocate(dicttranspose)

        allocate(result(Ne*Nd),dict(Nd*L),dicttranspose(Nd*L),stat=istat)

        result = 0.0
        dict = 0.0
        dicttranspose = 0.0

        do ii = 1,numdictsingle
            read(iunitdict,rec=(kk-1)*numdictsingle+ii) imagedict
            imagedictflt = float(imagedict)
            do mm = 1,L
                if (imagedictflt(mm) .lt. 0) imagedictflt(mm) = imagedictflt(mm) + 255.0
            end do
            dict((ii-1)*L+1:ii*L) = (imagedictflt(1:L)-meandict(1:L))/sqrt(sum((imagedictflt(1:L)-meandict(1:L))**2))
        end do

        do ii = 1,L
            do mm = 1,numdictsingle
                dicttranspose((ii-1)*numdictsingle+mm) = dict((mm-1)*L+ii)
            end do
        end do
!if (kk .eq. 1) then
!imagedictflt = dict(91*L+1:92*L)
!open(unit=13,file='testdict.txt',action='write')
!do i=1,imgwd
!do j=1,imght
!write(13, '(F15.6)', advance='no') imagedictflt((i-1)*imght+j)
!end do
!write(13, *) ''  ! this gives you the line break
!end do
!close(13)
!end if

        call clEnqueueWriteBuffer(command_queue, cl_dict, cl_bool(.true.), 0_8, size_in_bytes_dict, dicttranspose(1), ierr)
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
!if (kk .eq. 2) then
!do ii = 1,numdictsingle
!dictsr(ii,1:L) = dict((ii-1)*L+1:ii*L)
!end do
!res = matmul(exptsr,transpose(dictsr))
!print*,res(2,1:20)
!print*,''
!print*,result(1025:1044)
!end if
! copying the dot products to an array for sorting later
        do i = 1,numexptsingle
            exptsorted(i,numdictsingle*(kk-1)+1:numdictsingle*kk) = result((i-1)*numdictsingle+1:i*numdictsingle)
        end do

        if (mod(kk,30) .eq. 0) then
            write(6,'(A26,I7,A26,I7)'),'Completed dot product of',kk*1024,'dictionary patterns with',ll*1024,'experimental patterns'
        end if
    !print*,maxval(result),maxloc(result)
    end do dictionaryloop

    if (allocated(result)) deallocate(result)
    if (allocated(dict)) deallocate(dict)
    if (allocated(dicttranspose)) deallocate(dicttranspose)

    allocate(result(Ne*Nd),dict(Nd*L),dicttranspose(Nd*L),stat=istat)

    result = 0.0
    dict = 0.0
    dicttranspose = 0.0
! do the leftover dot products

    leftoverloop: do ii = 1,MODULO(totnumdict,numdictsingle)
        read(iunitdict,rec=(kk-1)*numdictsingle+ii) imagedict
        imagedictflt = float(imagedict)
        do mm = 1,L
            if (imagedictflt(mm) .lt. 0) imagedictflt(mm) = imagedictflt(mm) + 255.0
        end do
        dict((ii-1)*L+1:ii*L) = (imagedictflt(1:L)-meandict(1:L))/sqrt(sum((imagedictflt(1:L)-meandict(1:L))**2))
    end do leftoverloop

    do ii = 1,L
        do mm = 1,numdictsingle
            dicttranspose((ii-1)*numdictsingle+mm) = dict((mm-1)*L+ii)
        end do
    end do

    call clEnqueueWriteBuffer(command_queue, cl_dict, cl_bool(.true.), 0_8, size_in_bytes_dict, dicttranspose(1), ierr)
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
        exptsorted(i,numdictsingle*(kk-1)+1:numdictsingle*kk) = result((i-1)*numdictsingle+1:i*numdictsingle)
    end do
    !sort the exptsorted array and return the top 100 or so values
    do ii = 1,numexptsingle
    !call system_clock(t1,clock_rate,clock_max)
        call SSORT(exptsorted(ii,:),exptindex(ii,:),totnumdict,-2)
    !call system_clock(t2,clock_rate,clock_max)
!print*,'Time for sorting one image array =',real(t2-t1)/real(clock_rate),'seconds'
    end do
!print*,exptindex(1,1:100)
print*,''
!print*,exptsorted(1,1:100)!maxval(result),maxloc(result)
end do experimentalloop

!print*,result(1)
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
!> @author Taken from http://www.personal.psu.edu/jhm/f90/examples/sort/sorthalf.f
!> @brief sort the array
!
!> Parameters described in the program
!
!--------------------------------------------------------------------------
SUBROUTINE SSORT (X, Y, N, KFLAG)
!***BEGIN PROLOGUE  SSORT
!***PURPOSE  Sort an array and optionally make the same interchanges in
!            an auxiliary array.  The array may be sorted in increasing
!            or decreasing order.  A slightly modified QUICKSORT
!            algorithm is used.
!***LIBRARY   SLATEC
!***CATEGORY  N6A2B
!***TYPE      SINGLE PRECISION (SSORT-S, DSORT-D, ISORT-I)
!***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
!***AUTHOR  Jones, R. E., (SNLA)
!           Wisniewski, J. A., (SNLA)
!***DESCRIPTION
!
!   SSORT sorts array X and optionally makes the same interchanges in
!   array Y.  The array X may be sorted in increasing order or
!   decreasing order.  A slightly modified quicksort algorithm is used.
!
!   Description of Parameters
!      X - array of values to be sorted   (usually abscissas)
!      Y - array to be (optionally) carried along
!      N - number of values in array X to be sorted
!      KFLAG - control parameter
!            =  2  means sort X in increasing order and carry Y along.
!            =  1  means sort X in increasing order (ignoring Y)
!            = -1  means sort X in decreasing order (ignoring Y)
!            = -2  means sort X in decreasing order and carry Y along.
!
!***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
!                 for sorting with minimal storage, Communications of
!                 the ACM, 12, 3 (1969), pp. 185-187.
!***REVISION HISTORY  (YYMMDD)
!   761101  DATE WRITTEN
!   761118  Modified to use the Singleton quicksort algorithm.  (JAW)
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891009  Removed unreferenced statement labels.  (WRB)
!   891024  Changed category.  (WRB)
!   891024  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   901012  Declared all variables; changed X,Y to SX,SY. (M. McClain)
!   920501  Reformatted the REFERENCES section.  (DWL, WRB)
!   920519  Clarified error messages.  (DWL)
!   920801  Declarations section rebuilt and code restructured to use
!           IF-THEN-ELSE-ENDIF.  (RWC, WRB)
!***END PROLOGUE  SSORT
!     .. Scalar Arguments ..
INTEGER KFLAG, N
!     .. Array Arguments ..
REAL X(*)
INTEGER Y(*)
!     .. Local Scalars ..
REAL R, T, TT, TTY, TY
INTEGER I, IJ, J, K, KK, L, M, NN
!     .. Local Arrays ..
INTEGER IL(21), IU(21)
!     .. External Subroutines ..
!     None
!     .. Intrinsic Functions ..
INTRINSIC ABS, INT
!***FIRST EXECUTABLE STATEMENT  SSORT
NN = N
IF (NN .LT. 1) THEN
PRINT *,      'The number of values to be sorted is not positive.'
RETURN
ENDIF
!
KK = ABS(KFLAG)
IF (KK.NE.1 .AND. KK.NE.2) THEN
PRINT *,      'The sort control parameter, K, is not 2, 1, -1, or -2.'
RETURN
ENDIF
!
!     Alter array X to get decreasing order if needed
!
IF (KFLAG .LE. -1) THEN
DO 10 I=1,NN
X(I) = -X(I)
10    CONTINUE
ENDIF
!
IF (KK .EQ. 2) GO TO 100
!
!     Sort X only
!
M = 1
I = 1
J = NN
R = 0.375E0
!
20 IF (I .EQ. J) GO TO 60
IF (R .LE. 0.5898437E0) THEN
R = R+3.90625E-2
ELSE
R = R-0.21875E0
ENDIF
!
30 K = I
!
!     Select a central element of the array and save it in location T
!
IJ = I + INT((J-I)*R)
T = X(IJ)
!
!     If first element of array is greater than T, interchange with T
!
IF (X(I) .GT. T) THEN
X(IJ) = X(I)
X(I) = T
T = X(IJ)
ENDIF
L = J
!
!     If last element of array is less than than T, interchange with T
!
IF (X(J) .LT. T) THEN
X(IJ) = X(J)
X(J) = T
T = X(IJ)
!
!        If first element of array is greater than T, interchange with T
!
IF (X(I) .GT. T) THEN
X(IJ) = X(I)
X(I) = T
T = X(IJ)
ENDIF
ENDIF
!
!     Find an element in the second half of the array which is smaller
!     than T
!
40 L = L-1
IF (X(L) .GT. T) GO TO 40
!
!     Find an element in the first half of the array which is greater
!     than T
!
50 K = K+1
IF (X(K) .LT. T) GO TO 50
!
!     Interchange these elements
!
IF (K .LE. L) THEN
TT = X(L)
X(L) = X(K)
X(K) = TT
GO TO 40
ENDIF
!
!     Save upper and lower subscripts of the array yet to be sorted
!
IF (L-I .GT. J-K) THEN
IL(M) = I
IU(M) = L
I = K
M = M+1
ELSE
IL(M) = K
IU(M) = J
J = L
M = M+1
ENDIF
GO TO 70
!
!     Begin again on another portion of the unsorted array
!
60 M = M-1
IF (M .EQ. 0) GO TO 190
I = IL(M)
J = IU(M)
!
70 IF (J-I .GE. 1) GO TO 30
IF (I .EQ. 1) GO TO 20
I = I-1
!
80 I = I+1
IF (I .EQ. J) GO TO 60
T = X(I+1)
IF (X(I) .LE. T) GO TO 80
K = I
!
90 X(K+1) = X(K)
K = K-1
IF (T .LT. X(K)) GO TO 90
X(K+1) = T
GO TO 80
!
!     Sort X and carry Y along
!
100 M = 1
I = 1
J = NN
R = 0.375E0
!
110 IF (I .EQ. J) GO TO 150
IF (R .LE. 0.5898437E0) THEN
R = R+3.90625E-2
ELSE
R = R-0.21875E0
ENDIF
!
120 K = I
!
!     Select a central element of the array and save it in location T
!
IJ = I + INT((J-I)*R)
T = X(IJ)
TY = Y(IJ)
!
!     If first element of array is greater than T, interchange with T
!
IF (X(I) .GT. T) THEN
X(IJ) = X(I)
X(I) = T
T = X(IJ)
Y(IJ) = Y(I)
Y(I) = TY
TY = Y(IJ)
ENDIF
L = J
!
!     If last element of array is less than T, interchange with T
!
IF (X(J) .LT. T) THEN
X(IJ) = X(J)
X(J) = T
T = X(IJ)
Y(IJ) = Y(J)
Y(J) = TY
TY = Y(IJ)
!
!        If first element of array is greater than T, interchange with T
!
IF (X(I) .GT. T) THEN
X(IJ) = X(I)
X(I) = T
T = X(IJ)
Y(IJ) = Y(I)
Y(I) = TY
TY = Y(IJ)
ENDIF
ENDIF
!
!     Find an element in the second half of the array which is smaller
!     than T
!
130 L = L-1
IF (X(L) .GT. T) GO TO 130
!
!     Find an element in the first half of the array which is greater
!     than T
!
140 K = K+1
IF (X(K) .LT. T) GO TO 140
!
!     Interchange these elements
!
IF (K .LE. L) THEN
TT = X(L)
X(L) = X(K)
X(K) = TT
TTY = Y(L)
Y(L) = Y(K)
Y(K) = TTY
GO TO 130
ENDIF
!
!     Save upper and lower subscripts of the array yet to be sorted
!
IF (L-I .GT. J-K) THEN
IL(M) = I
IU(M) = L
I = K
M = M+1
ELSE
IL(M) = K
IU(M) = J
J = L
M = M+1
ENDIF
GO TO 160
!
!     Begin again on another portion of the unsorted array
!
150 M = M-1
IF (M .EQ. 0) GO TO 190
I = IL(M)
J = IU(M)
!
160 IF (J-I .GE. 1) GO TO 120
IF (I .EQ. 1) GO TO 110
I = I-1
!
170 I = I+1
IF (I .EQ. J) GO TO 150
T = X(I+1)
TY = Y(I+1)
IF (X(I) .LE. T) GO TO 170
K = I
!
180 X(K+1) = X(K)
Y(K+1) = Y(K)
K = K-1
IF (T .LT. X(K)) GO TO 180
X(K+1) = T
Y(K+1) = TY
GO TO 170
!
!     Clean up
!
190 IF (KFLAG .LE. -1) THEN
DO 200 I=1,NN
X(I) = -X(I)
200    CONTINUE
ENDIF
RETURN
END
