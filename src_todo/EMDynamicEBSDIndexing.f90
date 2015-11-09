! ###################################################################
! Copyright (c) 2015, Marc De Graef/Carnegie Mellon University
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
! EMsoft:CTEMPEDIndexing.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMPEDIndexing
!
!> @author Saransh Singh/Marc De Graef, Carnegie Mellon University
!
!> @brief Indexing of EBSD patterns using the dictionary approach. Later, this program
!> will be used for dynamic pattern center correction to hopefully make
!> the dictionary indexing even more robust.
!
!> @date 03/25/15 SS 1.0 original
!--------------------------------------------------------------------------

program EBSDIndexing

use local
use typedefs
use NameListTypedefs
use NameListHandlers
use files
use io
use initializers
use EBSDmod

IMPLICIT NONE

character(fnlen)                            :: nmldeffile, progname, progdesc
type(EBSDNameListType)                      :: ebsdnl
type(EBSDLargeAccumType),pointer            :: acc
type(EBSDMasterType),pointer                :: master
logical                                     :: verbose
integer(kind=irg)                           :: istat

interface
        subroutine MasterSubroutine(ebsdnl,acc,master,progname)

        use local
        use typedefs
        use NameListTypedefs
        use NameListHandlers
        use files
        use dictmod
        use others, only: SSORT
        use crystal
        use initializers
        use gvectors
        use io
        use diffraction
        use symmetry
        use quaternions
        use constants
        use rotations
        use so3
        use math
        use EBSDmod
        use cl
        use omp_lib

        IMPLICIT NONE

        type(EBSDNameListType),INTENT(IN)                   :: ebsdnl
        type(EBSDLargeAccumType),pointer,INTENT(IN)         :: acc
        type(EBSDMasterType),pointer,Intent(IN)             :: master
        character(fnlen),INTENT(IN)                         :: progname

        end subroutine MasterSubroutine
end interface

nmldeffile = 'EMEBSD.nml'
progname = 'EMEBSDIndexing.f90'
progdesc = 'Program to index EBSD patterns using the dynamically calculated dictionary'
verbose = .TRUE.

! print some information
call EMsoft(progname, progdesc)

! deal with the command line arguments, if any
call Interpret_Program_Arguments(nmldeffile,2,(/ 0, 40 /), progname)

! deal with the namelist stuff
!call GetEBSDIndxNameList(nmldeffile,ebsdnl)
call GetEBSDNameList(nmldeffile,ebsdnl)

! 1. read the Monte Carlo data file
allocate(acc)
call EBSDreadMCfile(ebsdnl, acc)

! 2. read EBSD master pattern file
allocate(master)
call EBSDreadMasterfile(ebsdnl, master)

! 3. generate detector arrays
allocate(master%rgx(ebsdnl%numsx,ebsdnl%numsy), master%rgy(ebsdnl%numsx,ebsdnl%numsy), &
    master%rgz(ebsdnl%numsx,ebsdnl%numsy), stat=istat)
allocate(acc%accum_e_detector(ebsdnl%numEbins,ebsdnl%numsx,ebsdnl%numsy), stat=istat)
call EBSDGenerateDetector(ebsdnl, acc, master, verbose)
deallocate(acc%accum_e)

! perform the dictionary indexing computations
call MasterSubroutine(ebsdnl,acc,master,progname)

deallocate(master, acc)

end program EBSDIndexing

!--------------------------------------------------------------------------
!
! SUBROUTINE:MasterSubroutine
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief Master subroutine to control dictionary generation, inner products computations, sorting values
!> and indexing of points, all in parallel using OpenCL/openMP
!
!> @param ebsdnl ped indexing namelist pointer
!> @param acc accumulator pointer containing MC results
!> @param master master pattern is read into this pointer
!> @param progname name of the program
!
!> @date 03/30/15  SS 1.0 original
!> @date 05/05/15 MDG 1.1 removed getenv() call; replaced by global path strings
!--------------------------------------------------------------------------

subroutine MasterSubroutine(ebsdnl,acc,master,progname)

use local
use typedefs
use NameListTypedefs
use NameListHandlers
use files
use dictmod
use Lambert, only: LambertSphereToSquare
use others, only: SSORT
use crystal
use initializers
use gvectors
use io
use diffraction
use symmetry
use quaternions
use constants
use rotations
use so3
use math
use EBSDmod
use cl
use omp_lib
use HDF5
use HDFsupport
use NameListHDFwriters
use ISO_C_BINDING

IMPLICIT NONE

type(EBSDNameListType),INTENT(IN)                   :: ebsdnl
type(EBSDLargeAccumType),pointer,INTENT(IN)         :: acc
type(EBSDMasterType),pointer,Intent(IN)             :: master
character(fnlen),INTENT(IN)                         :: progname

type(unitcell),pointer                              :: cell
type(DynType)                                       :: Dyn
type(gnode)                                         :: rlp
logical                                             :: verbose

type(cl_platform_id)                                :: platform
type(cl_device_id)                                  :: device
type(cl_context)                                    :: context
type(cl_command_queue)                              :: command_queue
type(cl_mem)                                        :: cl_expt,cl_dict

integer(kind=irg)                                   :: num,ierr,irec,istat
integer(kind=irg),parameter                         :: iunit = 40
integer(kind=irg),parameter                         :: iunitexpt = 41
character(len = 100)                                :: info ! info about the GPU
integer, parameter                                  :: source_length = 50000
real(kind=dbl),parameter                            :: nAmpere = 6.241D+18   ! Coulomb per second

character(len = source_length)                      :: source

integer(kind=irg)                                   :: Ne,Nd,L,totnumexpt,numdictsingle,numexptsingle,imght,imgwd,nnk,recordsize
integer(kind=8)                                     :: size_in_bytes_dict,size_in_bytes_expt
integer(kind=1),allocatable                         :: imageexpt(:)
real(kind=sgl),allocatable                          :: imageexptflt(:),binned(:,:),imagedictflt(:),imagedictfltflip(:)
real(kind=sgl),allocatable                          :: result(:),expt(:),dict(:),dicttranspose(:),resultarray(:),&
eulerarray(:,:),resultmain(:,:),resulttmp(:,:)
integer(kind=irg),allocatable                       :: accum_e_MC(:,:)
real(kind=4),allocatable                            :: meandict(:),meanexpt(:),wf(:),master_array(:,:)
real(kind=sgl),allocatable                          :: EBSDpattern(:,:),FZarray(:,:)
integer(kind=irg)                                   :: i,j,ii,jj,kk,ll,mm,pp,qq
integer(kind=irg)                                   :: FZcnt, pgnum, io_int(3), ncubochoric
type(FZpointd),pointer                              :: FZlist, FZtmp
integer(kind=irg),allocatable                       :: indexlist(:),indexarray(:),indexmain(:,:),indextmp(:,:)
real(kind=sgl)                                      :: dmin,voltage,scl,prefactor
character(fnlen)                                    :: xtalname
integer(kind=irg)                                   :: binx,biny,TID,totnum_el,nthreads,Emin,Emax
real(kind=sgl)                                      :: sx,dx,dxm,dy,dym,rhos,x,projweight
real(kind=sgl)                                      :: dc(3),angle(4),ixy(2),bindx
integer(kind=irg)                                   :: nix,niy,nixp,niyp
character(fnlen)                                    :: str1,str2,str3,str4,str5,str6,str7,str8,str9,str10
real(kind=sgl)                                      :: euler(3)
integer(kind=irg)                                   :: index
character                                           :: TAB


TAB = CHAR(9)
verbose = .TRUE.
Ne = 512
Nd = 512
L = 80*60
totnumexpt = 196608
numdictsingle = 512
numexptsingle = 512
imght =80
imgwd = 60
nnk = 40
xtalname = 'Ni.xtal'
dmin = 0.04
voltage = 30000.0
ncubochoric = 50
recordsize = imght*imgwd+25

size_in_bytes_dict = Nd*L*sizeof(dict(1))
size_in_bytes_expt = Ne*L*sizeof(expt(1))

totnum_el = sum(acc%accum_e_detector)

! get the indices of the minimum and maximum energy
Emin = nint((ebsdnl%energymin - ebsdnl%Ehistmin)/ebsdnl%Ebinsize) +1
if (Emin.lt.1)  Emin=1
if (Emin.gt.ebsdnl%numEbins)  Emin=ebsdnl%numEbins

!====================================
! init a bunch of parameters
!====================================
! binned pattern array
binx = ebsdnl%numsx/ebsdnl%binning
biny = ebsdnl%numsy/ebsdnl%binning
bindx = 1.0/float(ebsdnl%binning)**2

! determine the scale factor for the Lambert interpolation; the square has
! an edge length of 2 x sqrt(pi/2)
scl = float(ebsdnl%npx) / LPs%sPio2

! intensity prefactor
prefactor = 0.25D0 * nAmpere * ebsdnl%beamcurrent * ebsdnl%dwelltime * 1.0D-15/ dble(totnum_el)

allocate(accum_e_MC(ebsdnl%numsx,ebsdnl%numsy),stat=istat)
accum_e_MC = sum(acc%accum_e_detector,1)
allocate(wf(ebsdnl%numEbins))
wf = sum(sum(acc%accum_e_detector,2),2)
wf = wf/sum(wf)

! for dictionary computations, the patterns are usually rather small, so perhaps the explicit
! energy sums can be replaced by an averaged approximate approach, in which all the energy bins
! are added together from the start, and all the master patterns are totaled as well...
! this is a straightforward sum; we should probably do a weighted sum instead
allocate(master_array(-ebsdnl%npx:ebsdnl%npx,-ebsdnl%npy:ebsdnl%npy))

do ii=Emin,Emax
    master%sr(-ebsdnl%npx:ebsdnl%npx,-ebsdnl%npy:ebsdnl%npy,ii) = &
    master%sr(-ebsdnl%npx:ebsdnl%npx,-ebsdnl%npy:ebsdnl%npy,ii) * wf(ii)
end do

master_array = sum(master%sr,3)

!=========================================
! ALLOCATION AND INITIALIZATION OF ARRAYS
!=========================================

allocate(expt(Ne*L),stat=istat)
if (istat .ne. 0) stop 'Could not allocate array for experimental patterns'
expt = 0.0

allocate(dict(Nd*L),dicttranspose(Nd*L),stat=istat)
if (istat .ne. 0) stop 'Could not allocate array for dictionary patterns'
dict = 0.0
dicttranspose = 0.0

allocate(result(Ne*Nd),stat=istat)
if (istat .ne. 0) stop 'Could not allocate array for results'
result = 0.0

allocate(imageexpt(recordsize),imageexptflt(imght*imgwd),imagedictflt(imght*imgwd),imagedictfltflip(imght*imgwd),stat=istat)
if (istat .ne. 0) stop 'Could not allocate array for reading experimental image patterns'
imageexpt = 0.0
imageexptflt = 0.0

allocate(meandict(imght*imgwd),meanexpt(imght*imgwd),&
    binned(ebsdnl%numsx/ebsdnl%binning,ebsdnl%numsy/ebsdnl%binning),stat=istat)
if (istat .ne. 0) stop 'Could not allocate array for mean dictionary and experimental patterns'
meandict = 0
meanexpt = 0

allocate(EBSDpattern(ebsdnl%numsx,ebsdnl%numsy),stat=istat)
if (istat .ne. 0) stop 'Could not allocate array for EBSD pattern'
EBSDpattern = 0.0


!================================
! INITIALIZATION OF CELL POINTER
!================================

nullify(cell)
allocate(cell)

verbose = .TRUE.
call Initialize_Cell(cell,Dyn,rlp,xtalname, dmin, voltage, verbose)

! determine the point group number
jj=0

do ii=1,32
    if (SGPG(ii).le.cell % SYM_SGnum) jj=ii
end do

pgnum = jj
write (*,*) 'point group number = ', pgnum

!==========================
! INITIALIZATION OF DEVICE
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

! allocate device memory
cl_expt = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_expt, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory for experimental data.'

cl_dict = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_dict, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory for dictionary data.'


!=====================================================
! SAMPLING OF RODRIGUES FUNDAMENTAL ZONE
!=====================================================

nullify(FZlist)
FZcnt = 0
write (*,*) 'pgnum = ', pgnum
write (*,*) 'N = ',ncubochoric

call sampleRFZ(ncubochoric, pgnum, FZcnt, FZlist)

! allocate and fill FZarray for OpenMP parallelization
allocate(FZarray(3,FZcnt),stat=istat)
FZarray = 0.0

FZtmp => FZlist

do ii = 1,FZcnt
    FZarray(1:3,ii) = FZtmp%rod(1:3)
    FZtmp => FZtmp%next
end do

io_int(1) = FZcnt
call WriteValue(' Number of incident beam directions       : ', io_int, 1, "(I8)")

! allocate some arrays

allocate(resultarray(numdictsingle),stat=istat)
if (istat .ne. 0) stop 'could not allocate result arrays'

resultarray = 0.0

allocate(indexarray(1:numdictsingle),stat=istat)
if (istat .ne. 0) stop 'could not allocate index arrays'

indexarray = 0

allocate(indexlist(1:numdictsingle*(ceiling(float(FZcnt)/float(numdictsingle)))),stat=istat)
if (istat .ne. 0) stop 'could not allocate indexlist arrays'

indexlist = 0

do ii = 1,numdictsingle*ceiling(float(FZcnt)/float(numdictsingle))
    indexlist(ii) = ii
end do

allocate(resultmain(nnk,numexptsingle*ceiling(float(totnumexpt)/float(numexptsingle))),stat=istat)
if (istat .ne. 0) stop 'could not allocate main result array'

resultmain = 0.0

allocate(indexmain(nnk,numexptsingle*ceiling(float(totnumexpt)/float(numexptsingle))),stat=istat)
if (istat .ne. 0) stop 'could not allocate main index array'

indexmain = 0

allocate(resulttmp(2*nnk,numexptsingle*ceiling(float(totnumexpt)/float(numexptsingle))),stat=istat)
if (istat .ne. 0) stop 'could not allocate temporary result array'

resulttmp = 0.0

allocate(indextmp(2*nnk,numexptsingle*ceiling(float(totnumexpt)/float(numexptsingle))),stat=istat)
if (istat .ne. 0) stop 'could not allocate temporary index array'

indextmp = 0

allocate(eulerarray(1:3,numdictsingle*ceiling(float(FZcnt)/float(numdictsingle))),stat=istat)
if (istat .ne. 0) stop 'could not allocate euler array'

!=========================================================
! I/O FOR EXPERIMENTAL DATA SET AND CALCULATION OF MEAN
!=========================================================

call Message(' -> opened experimental file for I/O', frm = "(A)" )

open(unit=iunitexpt,file=trim(EMdatapathname)//'Experiments/FrameData',&
    status='old',form='unformatted',access='direct',recl=recordsize,iostat=ierr)

!!!!$OMP PARALLEL PRIVATE(TID,ii,jj,imageexpt,imageexptflt) &
!!!!$OMP& SHARED(meandict,meanexpt)

!TID = OMP_GET_THREAD_NUM()

!if (TID .eq. 0) then
call Message(' -> Start calculating mean pattern for observed patterns', frm = "(A)" )

do ii = 1,totnumexpt
    read(iunitexpt,rec=ii) imageexpt
    imageexptflt = float(imageexpt(26:recordsize))

    do jj = 1,L
        if (imageexptflt(jj) .lt. 0) imageexptflt(jj) = imageexptflt(jj)+256.0
    end do

    meanexpt = meanexpt+imageexptflt

    end do

meanexpt = meanexpt/totnumexpt
meanexpt = meanexpt/NORM2(meanexpt)

! the mean of the dictionary is assumed to be the same as the experiments
meandict = meanexpt

!open(unit=13,file='testmean.txt',action='write')
!do ll = 1,80
!do mm = 1,60
!write(13, '(F15.6)', advance='no') meanexpt((mm-1)*80+ll)
!end do
!write(13, *) ''  ! this gives you the line break
!end do
!close(13)

call Message(' -> Finished calculating mean pattern for observed patterns', frm = "(A)" )
!end if

!if (TID .eq. 1) then
!    call Message(' -> Start estimating mean pattern for dictionnary from Monte carlo data file', frm = "(A)" )

!    do pp = 1,ebsdnl%numsx,ebsdnl%binning
!        do qq = 1,ebsdnl%numsy,ebsdnl%binning
!            binned(pp/ebsdnl%binning+1,qq/ebsdnl%binning+1) = &
!            float(sum(accum_e_MC(pp:pp+ebsdnl%binning-1,qq:qq+ebsdnl%binning-1)))
!        end do
!    end do

!    binned = binned/NORM2(binned)

!    do pp = 1,ebsdnl%numsx/ebsdnl%binning
!        do qq = 1,ebsdnl%numsy/ebsdnl%binning
!            meandict((qq-1)*ebsdnl%numsx/ebsdnl%binning+pp) = binned(pp,qq)
!        end do
!    end do

!    do pp = 1,imgwd
!        imagedictfltflip((pp-1)*imght+1:pp*imght) = meandict((imgwd-pp)*imght+1:(imgwd-pp+1)*imght)
!    end do

!    meandict(1:L) = imagedictfltflip(1:L)

!open(unit=13,file='testmean.txt',action='write')
!do pp = 1,80
!do qq = 1,60
!write(13, '(F15.6)', advance='no') meandict((qq-1)*80+pp)
!end do
!write(13, *) ''  ! this gives you the line break
!end do
!close(13)

!call Message(' -> Finished calculating mean pattern for dictionary', frm = "(A)" )

!end if
!!!!$OMP END PARALLEL


!FZtmp => FZlist

dictionaryloop: do ii = 1,ceiling(float(FZcnt)/float(numdictsingle))
    dict = 0.0
    result = 0.0

    if (ii .le. floor(float(FZcnt)/float(numdictsingle))) then

!$OMP PARALLEL PRIVATE(TID,binned,angle,ll,mm,dc,ixy,istat,nix,niy) &
!$OMP& PRIVATE(nixp,niyp,dx,dy,dxm,dym,EBSDpattern,imagedictflt,projweight) &
!$OMP& SHARED(pp,dict,master_array,accum_e_MC,meandict,eulerarray,ii,numdictsingle,FZcnt,prefactor)

        TID = OMP_GET_THREAD_NUM()
!$OMP BARRIER
!$OMP DO SCHEDULE(DYNAMIC)

        do pp = 1,numdictsingle

            binned = 0.0
            angle = ro2qu(FZarray(1:3,(ii-1)*numdictsingle+pp))
            do ll=1,ebsdnl%numsx
                do mm=1,ebsdnl%numsy
!  do the active coordinate transformation for this euler angle
                    dc = sngl(quat_Lp(conjg(angle(1:4)),  (/ master%rgx(ll,mm),master%rgy(ll,mm),master%rgz(ll,mm) /) ))
! make sure the third one is positive; if not, switch all
                    dc = dc/sqrt(sum(dc**2))
                    if (dc(3).lt.0.0) dc = -dc
! convert these direction cosines to coordinates in the Rosca-Lambert projection
                    ixy = scl * LambertSphereToSquare( dc, istat )
! four-point interpolation (bi-quadratic)
                    nix = int(ebsdnl%npx+ixy(1))-ebsdnl%npx
                    niy = int(ebsdnl%npy+ixy(2))-ebsdnl%npy
                    nixp = nix+1
                    niyp = niy+1
                    if (nixp.gt.ebsdnl%npx) nixp = nix
                    if (niyp.gt.ebsdnl%npy) niyp = niy
                    dx = ixy(1)-nix
                    dy = ixy(2)-niy
                    dxm = 1.0-dx
                    dym = 1.0-dy
! interpolate the intensity
                    EBSDpattern(ll,mm) = EBSDpattern(ll,mm) + accum_e_MC(ll,mm) * ( master_array(nix,niy) * dxm * dym + &
                    master_array(nixp,niy) * dx * dym + master_array(nix,niyp) * dxm * dy + &
                    master_array(nixp,niyp) * dx * dy )

                end do
            end do

            EBSDpattern = prefactor * EBSDpattern

            if (ebsdnl%binning .ne. 1) then
                do ll=1,ebsdnl%numsx,ebsdnl%binning
                    do mm=1,ebsdnl%numsy,ebsdnl%binning
                        binned(ll/ebsdnl%binning+1,mm/ebsdnl%binning+1) = &
                        sum(EBSDpattern(ll:ll+ebsdnl%binning-1,mm:mm+ebsdnl%binning-1))
                    end do
                end do
! and divide by binning^2
                binned = binned * bindx
            else
                binned = EBSDpattern
            end if

            do ll = 1,ebsdnl%numsx/ebsdnl%binning
                do mm = 1,ebsdnl%numsy/ebsdnl%binning
                    imagedictflt((mm-1)*ebsdnl%numsx/ebsdnl%binning+ll) = binned(ll,mm)
                end do
            end do

            do ll = 1,imgwd
                imagedictfltflip((ll-1)*imght+1:ll*imght) = imagedictflt((imgwd-ll)*imght+1:(imgwd-ll+1)*imght)
            end do
            imagedictflt = imagedictfltflip

            dict((pp-1)*L+1:pp*L) = imagedictflt(1:L)/NORM2(imagedictflt(1:L))
            projweight = DOT_PRODUCT(meandict,dict((pp-1)*L+1:pp*L))
            dict((pp-1)*L+1:pp*L) = dict((pp-1)*L+1:pp*L) - projweight*meandict
            dict((pp-1)*L+1:pp*L) = dict((pp-1)*L+1:pp*L)/NORM2(dict((pp-1)*L+1:pp*L))
!if ((ii-1)*numdictsingle+pp .eq. 500) then
!print*,180.0/cPi*ro2eu(FZarray(1:3,(ii-1)*numdictsingle+pp))
!open(unit=13,file='testdict.txt',action='write')
!do ll = 1,80
!do mm = 1,60
!write(13, '(F15.6)', advance='no') dict((mm-1)*80+ll)
!end do
!write(13, *) ''  ! this gives you the line break
!end do
!close(13)
!end if
            eulerarray(1:3,(ii-1)*numdictsingle+pp) = 180.0/cPi*ro2eu(FZarray(1:3,(ii-1)*numdictsingle+pp))
        end do
!$OMP END DO
!$OMP END PARALLEL


    else if (ii .eq. ceiling(float(FZcnt)/float(numdictsingle))) then

!$OMP PARALLEL PRIVATE(TID,pp,binned,angle,ll,mm,dc,ixy,istat,nix,niy) &
!$OMP& PRIVATE(nixp,niyp,dx,dy,dxm,dym,EBSDpattern,imagedictflt,projweight) &
!$OMP& SHARED(dict,master_array,accum_e_MC,meandict,eulerarray,ii,numdictsingle,FZcnt,prefactor)

TID = OMP_GET_THREAD_NUM()
!$OMP BARRIER
!$OMP DO SCHEDULE(DYNAMIC)

        do pp = 1,MODULO(FZcnt,numdictsingle)

            binned = 0.0
            angle = ro2qu(FZarray(1:3,(ii-1)*numdictsingle+pp))

            do ll=1,ebsdnl%numsx
                do mm=1,ebsdnl%numsy
!  do the active coordinate transformation for this euler angle
                    dc = sngl(quat_Lp(conjg(angle(1:4)),  (/ master%rgx(ll,mm),master%rgy(ll,mm),master%rgz(ll,mm) /) ))
! make sure the third one is positive; if not, switch all
                    dc = dc/sqrt(sum(dc**2))
                    if (dc(3).lt.0.0) dc = -dc
! convert these direction cosines to coordinates in the Rosca-Lambert projection
                    ixy = scl * LambertSphereToSquare( dc, istat )
! four-point interpolation (bi-quadratic)
                    nix = int(ebsdnl%npx+ixy(1))-ebsdnl%npx
                    niy = int(ebsdnl%npy+ixy(2))-ebsdnl%npy
                    nixp = nix+1
                    niyp = niy+1
                    if (nixp.gt.ebsdnl%npx) nixp = nix
                    if (niyp.gt.ebsdnl%npy) niyp = niy
                    dx = ixy(1)-nix
                    dy = ixy(2)-niy
                    dxm = 1.0-dx
                    dym = 1.0-dy
! interpolate the intensity
                    EBSDpattern(ll,mm) = EBSDpattern(ll,mm) + accum_e_MC(ll,mm) * ( master_array(nix,niy) * dxm * dym + &
                    master_array(nixp,niy) * dx * dym + master_array(nix,niyp) * dxm * dy + &
                    master_array(nixp,niyp) * dx * dy )

                end do
            end do

            EBSDpattern = prefactor * EBSDpattern

            if (ebsdnl%binning.ne.1) then
                do ll=1,ebsdnl%numsx,ebsdnl%binning
                    do mm=1,ebsdnl%numsy,ebsdnl%binning
                        binned(ll/ebsdnl%binning+1,mm/ebsdnl%binning+1) = &
                        sum(EBSDpattern(ll:ll+ebsdnl%binning-1,mm:mm+ebsdnl%binning-1))
                    end do
                end do
! and divide by binning^2
                binned = binned * bindx
            else
                binned = EBSDpattern
            end if

            do ll = 1,ebsdnl%numsx/ebsdnl%binning
                do mm = 1,ebsdnl%numsy/ebsdnl%binning
                    imagedictflt((mm-1)*ebsdnl%numsx/ebsdnl%binning+ll) = binned(ll,mm)
                end do
            end do

            do ll = 1,imgwd
                imagedictfltflip((ll-1)*imght+1:ll*imght) = imagedictflt((imgwd-ll)*imght+1:(imgwd-ll+1)*imght)
            end do
            imagedictflt = imagedictfltflip
            dict((pp-1)*L+1:pp*L) = imagedictflt(1:L)/NORM2(imagedictflt(1:L))
            projweight = DOT_PRODUCT(meandict,dict((pp-1)*L+1:pp*L))
            dict((pp-1)*L+1:pp*L) = dict((pp-1)*L+1:pp*L) - projweight*meandict
            dict((pp-1)*L+1:pp*L) = dict((pp-1)*L+1:pp*L)/NORM2(dict((pp-1)*L+1:pp*L))

            eulerarray(1:3,(ii-1)*numdictsingle+pp) = 180.0/cPi*ro2eu(FZarray(1:3,(ii-1)*numdictsingle+pp))
        end do

!$OMP END DO
!$OMP END PARALLEL

    end if

    dicttranspose = 0.0

    do ll = 1,L
        do mm = 1,Nd
            dicttranspose((ll-1)*Nd+mm) = dict((mm-1)*L+ll)
        end do
    end do

    call clEnqueueWriteBuffer(command_queue, cl_dict, cl_bool(.true.), 0_8, size_in_bytes_dict, dicttranspose(1), ierr)
    if(ierr /= CL_SUCCESS) stop 'Error: cannot write to buffer.'

    experimentalloop: do jj = 1,ceiling(float(totnumexpt)/float(numexptsingle))

        expt = 0.0

        if (jj .le. floor(float(totnumexpt)/float(numexptsingle))) then
            do pp = 1,numexptsingle
                read(iunitexpt,rec=(jj-1)*numexptsingle+pp) imageexpt
                imageexptflt = float(imageexpt(26:recordsize))
                do mm = 1,L
                    if (imageexptflt(mm) .lt. 0) imageexptflt(mm) = imageexptflt(mm) + 256.0
                end do
                imageexptflt(1:L) = imageexptflt - meanexpt(1:L)
                expt((pp-1)*L+1:pp*L) = imageexptflt(1:L)/NORM2(imageexptflt(1:L))

            end do

        else if (jj .eq. ceiling(float(totnumexpt)/float(numexptsingle))) then
            do pp = 1,MODULO(totnumexpt,numexptsingle)
                read(iunitexpt,rec=(jj-1)*numexptsingle+pp) imageexpt
                imageexptflt = float(imageexpt(26:recordsize))
                do mm = 1,L
                    if (imageexptflt(mm) .lt. 0) imageexptflt(mm) = imageexptflt(mm) + 256.0
                end do
                imageexptflt(1:L) = imageexptflt - meanexpt(1:L)
                expt((pp-1)*L+1:pp*L) = imageexptflt(1:L)/NORM2(imageexptflt(1:L))

            end do

        end if

        call clEnqueueWriteBuffer(command_queue, cl_expt, cl_bool(.true.), 0_8, size_in_bytes_expt, expt(1), ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot write to buffer.'

        call InnerProdGPU(cl_expt,cl_dict,Ne,Nd,L,result,source,source_length,platform,device,context,command_queue)

        if (floor(float(totnumexpt)/float(numexptsingle)/10.0) .gt. 0) then
            if (mod(jj,floor(float(totnumexpt)/float(numexptsingle)/10.0)) .eq. 0) then
                if (ii .le. floor(float(FZcnt)/float(numdictsingle))) then
                    write(6,'(A,I8,A,I8,A)')'Completed dot product of ',jj*numexptsingle,' experimental pattern with ',&
                    ii*numdictsingle,' dictionary patterns.'
                else if (ii .eq. ceiling(float(FZcnt)/float(numdictsingle))) then
                    write(6,'(A,I8,A,I8,A)')'Completed dot product of ',jj*numexptsingle,'experimental pattern with all ',&
                    FZcnt,'dictionary patterns'
                end if
            end if
        else
            if (mod(jj,2) .eq. 0) then
                if (ii .le. floor(float(FZcnt)/float(numdictsingle))) then
                    write(6,'(A,I8,A,I8,A)')'Completed dot product of ',jj*numexptsingle,' experimental pattern with ',&
                    ii*numdictsingle,' dictionary patterns.'
                else if (ii .eq. ceiling(float(FZcnt)/float(numdictsingle))) then
                    write(6,'(A,I8,A,I8,A)')'Completed dot product of ',jj*numexptsingle,'experimental pattern with all ',&
                    FZcnt,'dictionary patterns'
                end if
            end if
        end if
!print*,maxval(result)
!$OMP PARALLEL PRIVATE(TID,qq,resultarray,indexarray) &
!$OMP& SHARED(nthreads,result,resultmain,indexmain,ii,jj,numdictsingle,numexptsingle,indexlist,resulttmp,indextmp)
        TID = OMP_GET_THREAD_NUM()
        nthreads = OMP_GET_NUM_THREADS()
        if (TID .eq. 0 .and. ii .eq. 1 .and. jj .eq. 1) write(6,'(A20,I3)'),'Number of threads =',nthreads
!$OMP BARRIER
!$OMP DO SCHEDULE(DYNAMIC)
        do qq = 1,numexptsingle

            resultarray(1:numdictsingle) = result((qq-1)*numdictsingle+1:qq*numdictsingle)
            indexarray(1:numdictsingle) = indexlist((ii-1)*numdictsingle+1:ii*numdictsingle)

            call SSORT(resultarray,indexarray,numdictsingle,-2)
            resulttmp(nnk+1:2*nnk,(jj-1)*numexptsingle+qq) = resultarray(1:nnk)
            indextmp(nnk+1:2*nnk,(jj-1)*numexptsingle+qq) = indexarray(1:nnk)

            call SSORT(resulttmp(:,(jj-1)*numexptsingle+qq),indextmp(:,(jj-1)*numexptsingle+qq),2*nnk,-2)

            resultmain(1:nnk,(jj-1)*numexptsingle+qq) = resulttmp(1:nnk,(jj-1)*numexptsingle+qq)

            indexmain(1:nnk,(jj-1)*numexptsingle+qq) = indextmp(1:nnk,(jj-1)*numexptsingle+qq)

        end do
!$OMP END DO
!$OMP END PARALLEL

    end do experimentalloop

end do dictionaryloop

open(unit=iunit,file='Results.ctf',action='write')
call WriteHeader(iunit,'ctf')
!open(unit=iunit,file='euler-1.bin',status='unknown',action='write',&
!form='unformatted',access='stream')
!open(unit=13,file='BestMatch-100-meansub.bin',form='unformatted',action='write',status='unknown',access='stream')
do ii = 1,totnumexpt

!do jj = 1,nnk
!    index = indexarray(jj,ii)
!    samples(1:4,jj) = eu2qu((cPi/180.D0)*eulerangles(index,1:3))
!end do

!call DI_EMforDD(samples, dictlist, nnk, seed, muhat, kappahat,'WAT')
!euler = 180.0/cPi*qu2eu(muhat)

    index = indexmain(1,ii)
    euler = eulerarray(1:3,index)
    write(str1,'(F12.3)') float(MODULO(ii-1,512))
    write(str2,'(F12.3)') float(floor(float(ii-1)/512.0))
    write(str3,'(I2)') 10
    write(str4,'(I2)') 0
    write(str5,'(F12.3)') euler(1)
    write(str6,'(F12.3)') euler(2)
    write(str7,'(F12.3)') euler(3)
    write(str8,'(I8)') index
    write(str9,'(I3)') 255
    write(str10,'(I3)') 255
    write(iunit,'(A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A)')'1',TAB,trim(adjustl(str1)),TAB,&
    trim(adjustl(str2)),TAB,trim(adjustl(str3)),TAB,trim(adjustl(str4)),TAB,trim(adjustl(str5)),&
    TAB,trim(adjustl(str6)),TAB,trim(adjustl(str7)),TAB,trim(adjustl(str8)),TAB,trim(adjustl(str9)),&
    TAB,trim(adjustl(str10))
end do

close(iunit)

end subroutine MasterSubroutine

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
!> @param source the opencl kernel as a character array
!> @param length of character array
!> @param platform opencl platform type
!> @param device opencl device type
!> @param context opencl context type
!> @param command_queue opencl command queue
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

size_in_bytes_result = Ne*Nd*sizeof(result(1))
Wexp = L
Wdict = Nd
localsize = (/16,16/)
globalsize = (/Ne,Nd/)

!=====================
! INITIALIZATION
!=====================


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


! finally get the kernel and release the program
kernel = clCreateKernel(prog, 'InnerProd', ierr)
call clReleaseProgram(prog, ierr)

cl_result = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_result, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory for result.'


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


!print*,result(1)
call clReleaseKernel(kernel, ierr)
call clReleaseMemObject(cl_result, ierr)

end subroutine InnerProdGPU

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

