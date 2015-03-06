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

program DictionaryIndexing

use local
use typedefs
use NameListTypedefs
use NameListHandlers
use files
use io

IMPLICIT NONE

character(fnlen)                          :: nmldeffile, progname, progdesc
type(DictIndxOpenCLListType)              :: dictindxnl

nmldeffile = 'CTEMDictionaryIndexing.nml'
progname = 'CTEMInnerProdOpenCL.f90'
progdesc = 'Program to compute inner product of observed and calculated patterns on a GPU'

! deal with the command line arguments, if any
call Interpret_Program_Arguments(nmldeffile,2,(/ 0, 40 /), progname)

! deal with the namelist stuff
call GetDictIndxOpenCLNameList(nmldeffile,dictindxnl)

! print some information
call CTEMsoft(progname, progdesc)

! perform the dictionary indexing computations
call MasterSubroutine(dictindxnl, progname)


end program DictionaryIndexing

!--------------------------------------------------------------------------
!
! SUBROUTINE:MasterSubroutine
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief Master subroutine to control Inner products computations, sorting values
!> and indexing of points, all in parallel using openMP
!
!> @param dictindxnl dictionary indexing namelist pointer
!> @param progname name of the program
!
!> @date 01/27/15  SS 1.0 original
!--------------------------------------------------------------------------
subroutine MasterSubroutine(dictindxnl,progname)

use local
use NameListTypedefs
use NameListHandlers
use typedefs
use files
use io
use dictmod
use others, only: SSORT
use rotations, only: eu2qu,qu2eu
use constants
use cl
use omp_lib

IMPLICIT NONE

type(DictIndxOpenCLListType),INTENT(IN)             :: dictindxnl
character(fnlen),INTENT(IN)                         :: progname

real(kind=4),allocatable                            :: result(:),resultcpy(:)
real(kind=4),allocatable                            :: expt(:)
real(kind=4),allocatable                            :: dict(:),dicttranspose(:)
real(kind=4),allocatable                            :: eulerangles(:,:),imagedictflt(:),imageexptflt(:),meandict(:),meanexpt(:)
integer(kind=1),allocatable                         :: imagedict(:),imageexpt(:)
integer(kind=4)                                     :: Ne
integer(kind=4)                                     :: Nd
integer(kind=4)                                     :: L
logical                                             :: MeanSub
integer(kind=irg)                                   :: totnumdict,numdictsingle,numexptsingle,imght,imgwd,totnumexpt,nnk
integer(kind=4)                                     :: recordsize,filesize

integer(kind=4),parameter                           :: iunit = 40
integer(kind=4),parameter                           :: iunitdict = 41
integer(kind=4),parameter                           :: iunitexpt = 42
character(len = 100)                                :: info ! info about the GPU
integer(kind=8)                                     :: globalsize(2),localsize(2)
integer, parameter                                  :: source_length = 1000000
character(len = source_length)                      :: source
integer(kind=4)                                     :: num,ierr,istat,irec,Wexp,Wdict,i,j,ii,jj,kk,ll,mm,nn,pp,qq
integer(kind=8)                                     :: size_in_bytes_expt,size_in_bytes_dict,size_in_bytes_result
!real(kind=sgl),allocatable                          :: topk(:,:),topkintd(:,:),exptsorted(:,:)
!integer(kind=irg),allocatable                       :: indextopk(:,:),indextopkintd(:,:),indexlist(:),exptindex(:,:)
real(kind=sgl),allocatable                          :: arr(:,:),prevarr(:,:)
integer(kind=4),allocatable                         :: auxarr(:,:),prevauxarr(:,:),indexlist(:)
type(cl_platform_id)                                :: platform
type(cl_device_id)                                  :: device
type(cl_context)                                    :: context
type(cl_command_queue)                              :: command_queue
type(cl_mem)                                        :: cl_expt,cl_dict
integer(kind=irg)                                   :: TID,state1,state2,nthreads,index
integer(kind=irg)                                   :: seed
real(kind=dbl),allocatable                          :: samples(:,:)
type(dicttype),pointer                              :: dictlist
real(kind=dbl)                                      :: q0, q1, q2, q3, muhat(4), kappahat, qu(4)
integer(kind=irg)                                   :: t1,t2,rate,max
real(kind=4),allocatable                            :: resultarray(:,:)
integer(kind=4),allocatable                         :: indexarray(:,:)
real(kind=sgl),allocatable                          :: res(:,:),exptsr(:,:),dictsr(:,:)
character                                           :: TAB
character(fnlen)                                    :: str1,str2,str3,str4,str5,str6,str7,str8,str9,str10
real(kind=sgl)                                      :: euler(3)


TAB = CHAR(9)
seed = 432514
allocate(dictlist,stat=istat)
dictlist%Num_of_init = 3
dictlist%Num_of_iterations = 30
dictlist%pgnum = 32

call DI_Init(dictlist,'WAT')

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

Ne = dictindxnl%numexptsingle
Nd = dictindxnl%numdictsingle
L = (dictindxnl%imght*dictindxnl%imgwd)
totnumdict = dictindxnl%totnumdict
totnumexpt = dictindxnl%totnumexpt
numdictsingle = dictindxnl%numdictsingle
numexptsingle = dictindxnl%numexptsingle
imght = dictindxnl%imght
imgwd = dictindxnl%imgwd
nnk = dictindxnl%nnk
MeanSub = dictindxnl%MeanSubtraction
recordsize = 25 + L
filesize = 948633600

size_in_bytes_expt = L*Ne*sizeof(L)
size_in_bytes_dict = L*Nd*sizeof(L)
!size_in_bytes_result = Ne*Nd*sizeof(result(1))!Ne*Nd*sizeof(L)
Wexp = L
Wdict = Nd
localsize = (/16,16/)
globalsize = (/Ne,Nd/)

allocate(exptsr(Ne,L),dictsr(Nd,L),res(Ne,Nd))

allocate(resultarray(numdictsingle*ceiling(float(totnumdict)/float(numdictsingle)),numexptsingle),stat=istat)
resultarray = 0.0

if (istat .ne. 0) stop 'Could not allocate big result array'

allocate(indexarray(numdictsingle*ceiling(float(totnumdict)/float(numdictsingle)),numexptsingle),stat=istat)
indexarray = 0

if (istat .ne. 0) stop 'Could not allocate big result array'


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

!allocate(arr(numexptsingle,numdictsingle),auxarr(numexptsingle,numdictsingle),stat=istat)
!arr = 0.0
!auxarr = 0

!allocate(prevarr(numexptsingle,nnk),prevauxarr(numexptsingle,nnk),stat=istat)
!prevarr = 0.0
!prevauxarr = 0

allocate(indexlist(1:numdictsingle*(floor(float(totnumdict)/float(numdictsingle))+1)),stat=istat)
indexlist = 0

do ii = 1,numdictsingle*ceiling(float(totnumdict)/float(numdictsingle))
    indexlist(ii) = ii
end do


allocate(result(Ne*Nd),resultcpy(Ne*Nd),stat=istat)
result = 0.0
resultcpy = 0.0

if (allocated(dict)) deallocate(dict)
allocate(dict(Nd*L),stat=istat)
dict = 0.0

if (allocated(dicttranspose)) deallocate(dicttranspose)
allocate(dicttranspose(Nd*L),stat=istat)
dicttranspose = 0.0

if (allocated(expt)) deallocate(expt)
allocate(expt(Ne*L),stat=istat)
expt = 0.0

allocate(samples(4,nnk),stat=istat)
samples = 0

! allocate device memory
cl_expt = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_expt, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory for experimental data.'

cl_dict = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_dict, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory for dictionary data.'

!=====================================
! I/O FOR EULER ANGLE FILE
!=====================================

call Message(' -> opening '//trim(dictindxnl%eulerfile), frm = "(A)" )
open(dataunit,file=trim(dictindxnl%eulerfile),status='old',form='formatted',action='read',iostat=ierr)

read(dataunit,'(I15)') totnumdict
do ii = 1,totnumdict
read(dataunit,*) eulerangles(ii,1:3)
end do
close(dataunit,status='keep')
call Message(' -> completed reading '//trim(dictindxnl%eulerfile), frm = "(A)")

!=====================================================
! I/O FOR DICTIONARY AND EXPERIMENTAL DATA SET
!=====================================================

call Message(' -> opening '//trim(dictindxnl%dictfile), frm = "(A)" )
open(unit=iunitdict,file=trim(dictindxnl%dictfile),status='old',form='unformatted',access='direct',recl=imght*imgwd,iostat=ierr)

call Message(' -> opening '//trim(dictindxnl%exptfile), frm = "(A)" )
open(unit=iunitexpt,file=trim(dictindxnl%exptfile),status='old',form='unformatted',access='direct',recl=recordsize,iostat=ierr)

!==================================================================
! CALCULATING MEAN FOR BOTH DICTIONARY AND EXPERIMENTAL DATA SET
!==================================================================

if (MeanSub .eqv. .TRUE.) then

!$OMP PARALLEL PRIVATE(TID,imagedict,imagedictflt) &
!$OMP& SHARED(meandict,meanexpt)
    TID = OMP_GET_THREAD_NUM()

    if (TID .eq. 0) then
        call Message(' -> Start calculating mean pattern for dictionary '//trim(dictindxnl%dictfile), frm = "(A)" )
        do ii = 1,totnumdict
            read(iunitdict,rec=ii) imagedict
            imagedictflt = float(imagedict)
            do jj = 1,L
                if (imagedictflt(jj) .lt. 0) imagedictflt(jj) = imagedictflt(jj)+255.0
            end do
            meandict = meandict+imagedictflt/sqrt(sum(imagedictflt**2))
        end do
        meandict = meandict/totnumdict
        call Message(' -> Finished calculating mean pattern for dictionary pattern '//trim(dictindxnl%dictfile), frm = "(A)" )
    end if


    if (TID .eq. 1) then
        call Message(' -> Start calculating mean pattern for observed patterns '//trim(dictindxnl%exptfile), frm = "(A)" )
        do ii = 1,totnumexpt
            read(iunitexpt,rec=ii) imageexpt
            imageexptflt = float(imageexpt(26:recordsize))
            do jj = 1,L
                if (imageexptflt(jj) .lt. 0) imageexptflt(jj) = imageexptflt(jj)+255.0
            end do
            meanexpt = meanexpt+imageexptflt/sqrt(sum(imageexptflt**2))
        end do
        meanexpt = meanexpt/totnumexpt
        call Message(' -> Finished calculating mean pattern for observed patterns '//trim(dictindxnl%exptfile), frm = "(A)" )
    end if
!$OMP END PARALLEL
end if

!=====================================================
! MAIN LOOP FOR EVALUATING DOT PRODUCTS AND INDEXING
!=====================================================
open(unit=iunit,file='Results.ctf',action='write',status='unknown')

call WriteHeader(iunit,'ctf')

experimentalloop: do ll = 1,ceiling(float(totnumexpt)/float(numexptsingle))

    expt = 0.0
    resultarray = 0
    indexarray = 0

    if (ll .le. floor(float(totnumexpt)/float(numexptsingle))) then
        do ii = 1,numexptsingle
            read(iunitexpt,rec=(ll-1)*numexptsingle+ii) imageexpt
            imageexptflt = float(imageexpt(26:recordsize))
            do mm = 1,L
                if (imageexptflt(mm) .lt. 0) imageexptflt(mm) = imageexptflt(mm) + 255.0
            end do
            expt((ii-1)*L+1:ii*L) = (imageexptflt(1:L)/sqrt(sum(imageexptflt(1:L)**2)))-meanexpt
            expt((ii-1)*L+1:ii*L) = expt((ii-1)*L+1:ii*L)/sqrt(sum(expt((ii-1)*L+1:ii*L)**2))
        end do

    else if (ll .eq. ceiling(float(totnumexpt)/float(numexptsingle))) then
        do ii = 1,MODULO(totnumexpt,numexptsingle)
            read(iunitexpt,rec=(ll-1)*numexptsingle+ii) imageexpt
            imageexptflt = float(imageexpt(26:recordsize))
            do mm = 1,L
                if (imageexptflt(mm) .lt. 0) imageexptflt(mm) = imageexptflt(mm) + 255.0
            end do
            expt((ii-1)*L+1:ii*L) = (imageexptflt(1:L)/sqrt(sum(imageexptflt(1:L)**2)))-meanexpt
            expt((ii-1)*L+1:ii*L) = expt((ii-1)*L+1:ii*L)/sqrt(sum(expt((ii-1)*L+1:ii*L)**2))
        end do
    end if


    call clEnqueueWriteBuffer(command_queue, cl_expt, cl_bool(.true.), 0_8, size_in_bytes_expt, expt(1), ierr)
    if(ierr /= CL_SUCCESS) stop 'Error: cannot write to buffer.'
    state1 = 0
    state2 = 0

    dictionaryloop: do kk = 1,ceiling(float(totnumdict)/float(numdictsingle))

        dict = 0.0

        if (kk .le. floor(float(totnumdict)/float(numdictsingle))) then
            do ii = 1,numdictsingle
                read(iunitdict,rec=(kk-1)*numdictsingle+ii) imagedict
                imagedictflt = float(imagedict)
                do mm = 1,L
                    if (imagedictflt(mm) .lt. 0) imagedictflt(mm) = imagedictflt(mm) + 255.0
                end do
                dict((ii-1)*L+1:ii*L) = (imagedictflt(1:L)/sqrt(sum(imagedictflt(1:L)**2))) - meandict

                dict((ii-1)*L+1:ii*L) = dict((ii-1)*L+1:ii*L)/sqrt(sum(dict((ii-1)*L+1:ii*L)**2))

            end do
        else if (kk .eq. ceiling(float(totnumdict)/float(numdictsingle))) then
            do ii = 1,MODULO(totnumdict,numdictsingle)

                read(iunitdict,rec=(kk-1)*numdictsingle+ii) imagedict
                imagedictflt = float(imagedict)
                do mm = 1,L
                    if (imagedictflt(mm) .lt. 0) imagedictflt(mm) = imagedictflt(mm) + 255.0
                end do
                dict((ii-1)*L+1:ii*L) = (imagedictflt(1:L)/sqrt(sum(imagedictflt(1:L)**2))) - meandict
                dict((ii-1)*L+1:ii*L) = dict((ii-1)*L+1:ii*L)/sqrt(sum(dict((ii-1)*L+1:ii*L)**2))

            end do
        end if

        dicttranspose = 0.0

        do ii = 1,L
            do jj = 1,Nd
                dicttranspose((ii-1)*Nd+jj) = dict((jj-1)*L+ii)
            end do
        end do

        call clEnqueueWriteBuffer(command_queue, cl_dict, cl_bool(.true.), 0_8, size_in_bytes_dict, dicttranspose(1), ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot write to buffer.'

        call InnerProdGPU(cl_expt,cl_dict,Ne,Nd,L,result,source,source_length,platform,device,context,command_queue)

        if (mod(kk,floor(float(totnumdict)/float(numdictsingle)/10.0)) .eq. 0) then
            if (ll .le. floor(float(totnumexpt)/float(numexptsingle))) then
                write(6,'(A26,I7,A26,I7)'),'Completed dot product of',kk*numdictsingle,'dictionary patterns with',ll*numexptsingle,&
                    'experimental patterns'
            else if (ll .eq. ceiling(float(totnumexpt)/float(numexptsingle))) then
                write(6,'(A26,I7,A30,I7)'),'Completed dot product of',kk*numdictsingle,'dictionary patterns with all',totnumexpt,&
                'experimental patterns'
            end if
        end if

        do ii = 1,numexptsingle
            resultarray((kk-1)*numdictsingle+1:kk*numdictsingle,ii) = result((ii-1)*numdictsingle+1:ii*numdictsingle)
        end do

    end do dictionaryloop

    if (ll .le. floor(float(totnumexpt)/float(numexptsingle))) then

        write(6,'(A28,I7,A26,I7)'),'Completed dot product of all',totnumdict,&
        'dictionary patterns with',ll*numexptsingle,'experimental patterns'

    else if (ll .eq. ceiling(float(totnumexpt)/float(numexptsingle))) then

        write(6,'(A28,I7,A30,I7)'),'Completed dot product of all',totnumdict,'dictionary patterns with all',totnumexpt,&
        'experimental patterns'

    end if

    write(6,'(A41,I8,A10)'),' -> Starting sorting and indexing of the',numexptsingle,'  patterns'

    if (ll .le. floor(float(totnumexpt)/float(numexptsingle))) then
!$OMP PARALLEL PRIVATE(TID,arr,auxarr,prevarr,prevauxarr,index,samples,muhat,kappahat) &
!$OMP& PRIVATE(str1,str2,str3,str4,str5,str6,str7,str8,str9,str10,euler) &
!$OMP& SHARED(indexlist,nthreads)
    TID = OMP_GET_THREAD_NUM()
    nthreads = OMP_GET_NUM_THREADS()
    if (TID .eq. 0) write(6,'(A20,I3)'),'Number of threads =',nthreads

!$OMP DO SCHEDULE(DYNAMIC)
        do ii = 1,numexptsingle
            indexarray(:,ii) = indexlist(:)
            call SSORT(resultarray(:,ii),indexarray(:,ii),numdictsingle*ceiling(float(totnumdict)/float(numdictsingle)),-2)
            do jj = 1,nnk
                index = indexarray(jj,ii)
                samples(1:4,jj) = eu2qu((cPi/180.D0)*eulerangles(index,1:3))
            end do
            call DI_EMforDD(samples, dictlist, nnk, seed, muhat, kappahat,'WAT')
            euler = 180.0/cPi*qu2eu(muhat)
            write(str1,'(F12.3)') float(MODULO(ii-1,512))
            write(str2,'(F12.3)') float(floor(float((ll-1)*numexptsingle+ii-1)/512.0))
            write(str3,'(I2)') 10
            write(str4,'(I2)') 0
            write(str5,'(F12.3)') euler(1)
            write(str6,'(F12.3)') euler(2)
            write(str7,'(F12.3)') euler(3)
            write(str8,'(F12.3)') 180.0*acos(1.0-1.0/kappahat)/cPi
            write(str9,'(I3)') 255
            write(str10,'(I3)') 255
!$OMP CRITICAL
            write(iunit,'(A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A)')'1',TAB,trim(adjustl(str1)),TAB,&
            trim(adjustl(str2)),TAB,trim(adjustl(str3)),TAB,trim(adjustl(str4)),TAB,trim(adjustl(str5)),&
            ' ',trim(adjustl(str6)),' ',trim(adjustl(str7)),TAB,trim(adjustl(str8)),TAB,trim(adjustl(str9)),&
            TAB,trim(adjustl(str10))
!$OMP END CRITICAL
        !write (*,*) 'mu    = ',muhat
        !write (*,*) 'kappa = ',kappahat
        !write (*,*) 'equivalent angular precision : ',180.D0*dacos(1.D0-1.D0/kappahat)/cPi
        end do
!$OMP END DO
!$OMP END PARALLEL

    else if (ll .eq. ceiling(float(totnumexpt)/float(numexptsingle))) then
!$OMP PARALLEL PRIVATE(TID,arr,auxarr,prevarr,prevauxarr,index,samples,muhat,kappahat) &
!$OMP& SHARED(indexlist,nthreads)
    TID = OMP_GET_THREAD_NUM()
    nthreads = OMP_GET_NUM_THREADS()
    if (TID .eq. 0) write(6,'(A20,I3)'),'Number of threads =',nthreads

!$OMP DO SCHEDULE(DYNAMIC)
        do ii = 1,MODULO(totnumexpt,numexptsingle)
            indexarray(:,ii) = indexlist(:)
            call SSORT(resultarray(:,ii),indexarray(:,ii),numdictsingle*ceiling(float(totnumdict)/float(numdictsingle)),-2)
            do jj = 1,nnk
                index = indexarray(jj,ii)
                samples(1:4,jj) = eu2qu((cPi/180.D0)*eulerangles(index,1:3))
            end do
            call DI_EMforDD(samples, dictlist, nnk, seed, muhat, kappahat,'WAT')
!$OMP CRITICAL
            write(iunit,'(A,A,2F12.4,A,I2,A,I1,A,3F12.4,A,F12.4,A,I3,A,I3)')'1',TAB,float(MODULO(ii-1,512)),&
            float(floor(float((ll-1)*numexptsingle+ii-1)/512.0)),&
            TAB,10,TAB,0,TAB,180.0/cPi*qu2eu(muhat),TAB,180.0*acos(1.0-1.0/kappahat)/cPi,TAB,255,TAB,255
!$OMP END CRITICAL
        end do
!$OMP END DO
!$OMP END PARALLEL

    end if

    write(6,'(A)'),' -> Finished sorting, indexing and I/O of the patterns'

end do experimentalloop
write(6,'(A)'),' -> Finished indexing....quitting code now'
close(iunit)
call clReleaseCommandQueue(command_queue, ierr)
call clReleaseContext(context, ierr)
call clReleaseMemObject(cl_expt, ierr)
call clReleaseMemObject(cl_dict, ierr)


end subroutine MasterSubroutine

!--------------------------------------------------------------------------
!
! SUBROUTINE:WriteHeader
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief Write the header for *.ctf or *.ang file file format
!
!> @param iunit unit to write to
!> @param fileformat string to tell which type of header to write
!
!> @date 02/07/15  SS 1.0 original
!--------------------------------------------------------------------------
subroutine WriteHeader(iunit,fileformat)

use local

IMPLICIT NONE

integer(kind=irg),INTENT(IN)                        :: iunit
character(len=3),INTENT(IN)                         :: fileformat

logical                                             :: isopen
integer(kind=irg)                                   :: ierr
character(fnlen)                                    :: filename
character                                           :: TAB

TAB = CHAR(9)
inquire(unit=iunit,OPENED=isopen)
if (isopen .eqv. .FALSE.) then
    filename = 'test.'//fileformat
    write(6,*),'Warning: No such file exists...creating file with name test.'//fileformat
    open(unit=iunit,file=trim(filename),status='unknown',action='write',iostat=ierr)
end if

if (fileformat .eq. 'ctf') then
    write(iunit,'(A)'),'Channel Text File'
    write(iunit,'(A)'),'Prj Test'
    write(iunit,'(3A)'),'Author	[Unknown]'
    write(iunit,'(A)'),'JobMode	Grid'
    write(iunit,'(3A)'),'XCells',Tab,'512'
    write(iunit,'(3A)'),'YCells',TAB,'384'
    write(iunit,'(3A)'),'XStep',TAB,'1'
    write(iunit,'(3A)'),'YStep',TAB,'1'
    write(iunit,'(A)'),'AcqE1	0'
    write(iunit,'(A)'),'AcqE2	0'
    write(iunit,'(A)'),'AcqE3	0'
    write(iunit,'(A)',advance='no'),'Euler angles refer to Sample Coordinate system (CS0)!  '
    write(iunit,'(A)')'Mag	30	Coverage	100	Device	0	KV	288.9	TiltAngle	-1	TiltAxis	0'
    write(iunit,'(A)'),'Phases	1'
    write(iunit,'(A)'),'3.524;3.524;3.524	90;90;90	Nickel	11	225'
    write(iunit,'(A)'),'Phase	X	Y	Bands	Error	Euler1	Euler2	Euler3	MAD	BC	BS'

else if (fileformat .eq. 'ang') then
    write(iunit,'(A)'),'# TEM_PIXperUM          1.000000'
    write(iunit,'(A)'),'# x-star                0.372300'
    write(iunit,'(A)'),'# y-star                0.689300'
    write(iunit,'(A)'),'# z-star                0.970100'
    write(iunit,'(A)'),'# WorkingDistance       5.000000'
    write(iunit,'(A)'),'#'
    write(iunit,'(A)'),'# Phase 1'
    write(iunit,'(A)'),'# MaterialName  	Nickel'
    write(iunit,'(A)'),'# Formula     	Ni'
    write(iunit,'(A)'),'# Info'
    write(iunit,'(A)'),'# Symmetry              43'
    write(iunit,'(A)'),'# LatticeConstants      3.520 3.520 3.520  90.000  90.000  90.000'
    write(iunit,'(A)'),'# NumberFamilies        4'
    write(iunit,'(A)'),'# hklFamilies   	 1  1  1 1 0.000000'
    write(iunit,'(A)'),'# hklFamilies   	 2  0  0 1 0.000000'
    write(iunit,'(A)'),'# hklFamilies   	 2  2  0 1 0.000000'
    write(iunit,'(A)'),'# hklFamilies   	 3  1  1 1 0.000000'
    write(iunit,'(A)'),'# Categories 0 0 0 0 0'
    write(iunit,'(A)'),'#'
    write(iunit,'(A)'),'# GRID: SqrGrid'
    write(iunit,'(A)'),'# XSTEP: 1.000000'
    write(iunit,'(A)'),'# YSTEP: 1.000000'
    write(iunit,'(A)'),'# NCOLS_ODD: 189'
    write(iunit,'(A)'),'# NCOLS_EVEN: 189'
    write(iunit,'(A)'),'# NROWS: 201'
    write(iunit,'(A)'),'#'
    write(iunit,'(A)'),'# OPERATOR: 	Administrator'
    write(iunit,'(A)'),'#'
    write(iunit,'(A)'),'# SAMPLEID:'
    write(iunit,'(A)'),'#'
    write(iunit,'(A)'),'# SCANID:'
    write(iunit,'(A)'),'#'
else

    stop 'Error: Can not recognize specified file format'

end if

end subroutine WriteHeader

!--------------------------------------------------------------------------
!
! SUBROUTINE:InnerProdGPU
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief Perform the inner product computations for the dictionary approach
!
!> @param expt vector with list of observed patterns
!> @param dict vector with list of calculated patterns
!> @param Ne number of patterns in the expt vector
!> @param Nd number of patterns in the dict vector
!> @param L size of one single pattern
!> @param result result of the matrix multiplication
!
!> @date 12/09/14  SS 1.0 original
!> @date 27/01/15  SS 1.1 modified to call the subroutine from mastersubroutine
!--------------------------------------------------------------------------
subroutine InnerProdGPU(cl_expt,cl_dict,Ne,Nd,L,result,source,source_length,platform,device,context,command_queue)

use local
use cl

IMPLICIT NONE

real(kind=4),INTENT(OUT)                            :: result(Ne*Nd)
type(cl_mem),INTENT(INOUT)                          :: cl_expt
type(cl_mem),INTENT(INOUT)                          :: cl_dict

integer(kind=4),INTENT(IN)                          :: Ne
integer(kind=4),INTENT(IN)                          :: Nd
integer(kind=4),INTENT(IN)                          :: L
character(len=source_length),INTENT(IN)             :: source
integer(kind=4),INTENT(IN)                          :: source_length
type(cl_platform_id),INTENT(IN)                     :: platform
type(cl_device_id),INTENT(INOUT)                    :: device
type(cl_context),INTENT(INOUT)                      :: context
type(cl_command_queue),INTENT(INOUT)                :: command_queue

type(cl_program)                                    :: prog
type(cl_kernel)                                     :: kernel
type(cl_mem)                                        :: cl_result
type(cl_event)                                      :: event

real(kind=4)                                        :: dicttranspose(Nd*L)
integer(kind=4),parameter                           :: iunit = 40
character(len = 100)                                :: info ! info about the GPU
integer(kind=8)                                     :: globalsize(2),localsize(2)
integer, parameter                                  :: source_length_build_info = 10000
character(len = source_length)                      :: source_build_info
integer(kind=4)                                     :: num,ierr,istat,irec,Wexp,Wdict,i,j,ii,jj,kk
integer(kind=8)                                     :: size_in_bytes_expt,size_in_bytes_dict,size_in_bytes_result

!size_in_bytes_expt = L*Ne*sizeof(L)
!size_in_bytes_dict = L*Nd*sizeof(L)
size_in_bytes_result = Ne*Nd*sizeof(result(1))
Wexp = L
Wdict = Nd
localsize = (/16,16/)
globalsize = (/Ne,Nd/)
!allocate(indexlist(1:numdictsingle),stat=istat)
!indexlist = 0

!=====================
! INITIALIZATION
!=====================
! get the platform ID
!call clGetPlatformIDs(platform, num, ierr)
!if(ierr /= CL_SUCCESS) stop "Cannot get CL platform."
! get the device ID
!call clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, device, num, ierr)
!if(ierr /= CL_SUCCESS) stop "Cannot get CL device."
! create the context and the command queue
!context = clCreateContext(platform, device, ierr)
!if(ierr /= CL_SUCCESS) stop "Cannot create context"

!command_queue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, ierr)
!if(ierr /= CL_SUCCESS) stop "Cannot create command queue"


!=====================
! BUILD THE KERNEL
!=====================

! create the program
prog = clCreateProgramWithSource(context, source, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot create program from source.'

! build
call clBuildProgram(prog, '-cl-fast-relaxed-math', ierr)

! get the compilation log
call clGetProgramBuildInfo(prog, device, CL_PROGRAM_BUILD_LOG, source_build_info, irec)
if(len(trim(source_build_info)) > 0) print*, trim(source_build_info)

if(ierr /= CL_SUCCESS) stop 'Error: program build failed.'

!write(6,*) "Kernel Build Successful...."

! finally get the kernel and release the program
kernel = clCreateKernel(prog, 'InnerProd', ierr)
call clReleaseProgram(prog, ierr)
!print*,'First three values in expt is ',expt(1:3)
! allocate device memory
!cl_expt = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_expt, ierr)
!if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

!cl_dict = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_dict, ierr)
!if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

cl_result = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_result, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory for result.'

!call clEnqueueWriteBuffer(command_queue, cl_expt, cl_bool(.true.), 0_8, size_in_bytes_expt, expt(1), ierr)
!if(ierr /= CL_SUCCESS) stop 'Error: cannot write to buffer.'

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
!dicttranspose = 0.0

!do ii = 1,L
!    do jj = 1,Nd
!        dicttranspose((ii-1)*Nd+jj) = dict((jj-1)*L+ii)
!    end do
!end do

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

!call clEnqueueWriteBuffer(command_queue, cl_dict, cl_bool(.true.), 0_8, size_in_bytes_dict, dicttranspose(1), ierr)
!if(ierr /= CL_SUCCESS) stop 'Error: cannot write to buffer.'

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
!do ii = 1,Ne
!exptsr(ii,1:L) = expt((ii-1)*L+1:ii*L)
!end do
!do ii = 1,Nd
!dictsr(ii,1:L) = dict((ii-1)*L+1:ii*L)
!end do
!res = matmul(exptsr,transpose(dictsr))
!print*,res(353,356:456)-result(352*numdictsingle+356:352*numdictsingle+456)
!print*,''
!print*,maxval(res),maxloc(res),maxval(result),maxloc(result)
!end if

!print*,result(1)
call clReleaseKernel(kernel, ierr)
!call clReleaseCommandQueue(command_queue, ierr)
!call clReleaseContext(context, ierr)
!call clReleaseMemObject(cl_expt, ierr)
!call clReleaseMemObject(cl_dict, ierr)
call clReleaseMemObject(cl_result, ierr)

end subroutine InnerProdGPU


