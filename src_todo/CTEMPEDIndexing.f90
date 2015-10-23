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
! CTEMsoft2013:CTEMPEDIndexing.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMPEDIndexing
!
!> @author Saransh Singh/Marc De Graef, Carnegie Mellon University
!
!> @brief Indexing of PED patterns using the dictionary approach
!
!> @date 03/13/15 SS 1.0 original
!--------------------------------------------------------------------------

program PEDIndexing

use local
use typedefs
use NameListTypedefs
use NameListHandlers
use files
use io
use initializers

IMPLICIT NONE

character(fnlen)                            :: nmldeffile, progname, progdesc
type(PEDKINIndxListType)                    :: pednl
type(unitcell),pointer                      :: cell
type(DynType),save                          :: Dyn
type(gnode),save                            :: rlp
type(reflisttype),pointer                   :: reflist, nexts, rltmpa
logical                                     :: verbose

nmldeffile = 'CTEMPEDIndexing.nml'
progname = 'CTEMPEDIndexing.f90'
progdesc = 'Program to index PED patterns using the dynamically calculated dictionary'

! deal with the command line arguments, if any
call Interpret_Program_Arguments(nmldeffile,2,(/ 0, 40 /), progname)

! deal with the namelist stuff
call GetPEDIndxNameList(nmldeffile,pednl)

! print some information
call CTEMsoft(progname, progdesc)

! perform the dictionary indexing computations
call MasterSubroutine(pednl, progname)


end program PEDIndexing

!--------------------------------------------------------------------------
!
! SUBROUTINE:MasterSubroutine
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief Master subroutine to control dictionary generation, inner products computations, sorting values
!> and indexing of points, all in parallel using openMP
!
!> @param pednl ped indexing namelist pointer
!> @param progname name of the program
!
!> @date 01/27/15  SS 1.0 original
!--------------------------------------------------------------------------
subroutine MasterSubroutine(pednl,progname)

use local
use typedefs
use NameListTypedefs
use NameListHandlers
use files
use dictmod
use others, only: SSORT
use dictmod
use crystal
use initializers
use gvectors
use io
use diffraction
use symmetry
use quaternions
use NameListTypedefs
use constants
use rotations
use so3
use math
use cl
use omp_lib

IMPLICIT NONE

type(PEDKINIndxListType),INTENT(IN)                 :: pednl
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
character(len = source_length)                      :: source

integer(kind=irg)                                   :: Ne,Nd,L,totnumexpt,numdictsingle,numexptsingle,imght,imgwd,nnk
integer(kind=8)                                     :: size_in_bytes_dict,size_in_bytes_expt
integer(kind=1),allocatable                         :: imageexpt(:)
real(kind=sgl),allocatable                          :: imageexptflt(:)
real(kind=sgl),allocatable                          :: result(:),expt(:),dict(:),dicttranspose(:),resultarray(:),&
eulerarray(:,:),resultmain(:,:),resulttmp(:,:)
integer(kind=irg)                                   :: i,j,ii,jj,kk,ll,mm,pp,qq
integer(kind=irg)                                   :: FZcnt, pgnum, io_int(3)
type(FZpointd),pointer                              :: FZlist, FZtmp
real(kind=sgl)                                      :: la,dval,dmin,glen,gmax,io_real(3),om(3,3),k(3),sgmax,FN(3),xgmin,Ig
real(kind=sgl)                                      :: Igmax,maxint,w,ku(3),kp(3),rnmpp,dx,dy,eu(3),x,y,euler(3)
integer(kind=irg)                                   :: gp(3),imh,imk,iml,nref,gg(3),ix,iy,iz,ww,nsize,tdp,sx,sy
type(reflisttype),pointer                           :: reflist, nexts, rltmpa
real(kind=sgl),allocatable                          :: pedpattern(:,:),image(:,:),xx(:,:),yy(:,:),line(:),dot(:,:),imagevector(:)
integer(kind=irg)                                   :: TID,nthreads,index
integer(kind=irg),allocatable                       :: indexlist(:),indexarray(:),indexmain(:,:),indextmp(:,:)
character(fnlen)                                    :: str1,str2,str3,str4,str5,str6,str7,str8,str9,str10
character                                           :: TAB
real(kind=4)                                        :: test(144*144),test1(144*144)
TAB = CHAR(9)

!==========================
! READING FROM NAMELIST
!==========================

Ne = pednl%numexptsingle
Nd = pednl%numdictsingle
L = (pednl%imght*pednl%imgwd)
totnumexpt = pednl%totnumexpt
numdictsingle = pednl%numdictsingle
numexptsingle = pednl%numexptsingle
imght = pednl%imght
imgwd = pednl%imgwd
nnk = pednl%nnk

!==============================
! INITIALIZING SOME PARAMETERS
!==============================

sgmax = 0.50
size_in_bytes_dict = Nd*L*sizeof(dict(1))
size_in_bytes_expt = Ne*L*sizeof(expt(1))

!=========================================
! ALLOCATION AND INITIALIZATION OF ARRAYS
!=========================================

allocate(expt(Ne*L),stat=istat)
if (istat .ne. 0) stop 'Could not allocate array for experimental patterns'

allocate(dict(Nd*L),dicttranspose(Nd*L),stat=istat)
if (istat .ne. 0) stop 'Could not allocate array for dictionary patterns'

allocate(result(Ne*Nd),stat=istat)
if (istat .ne. 0) stop 'Could not allocate array for results'

allocate(imageexpt(imght*imgwd),imageexptflt(imght*imgwd),stat=istat)
if (istat .ne. 0) stop 'Could not allocate array for reading experimental image patterns'

!================================
! INITIALIZATION OF CELL POINTER
!================================

nullify(cell)
allocate(cell)

verbose = .TRUE.
call Initialize_Cell(cell,Dyn,rlp,pednl%xtalname, pednl%dmin, pednl%voltage, verbose)

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
write (*,*) 'N = ',pednl%ncubochoric

call sampleRFZ(pednl%ncubochoric, pgnum, FZcnt, FZlist)

io_int(1) = FZcnt
call WriteValue(' Number of incident beam directions       : ', io_int, 1, "(I8)")

! allocate some arrays

allocate(resultarray(numdictsingle),stat=istat)
if (istat .ne. 0) stop 'could not allocate result arrays'

resultarray = 0.0

allocate(indexarray(numdictsingle),stat=istat)
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

!==================================================
! CALCULATING THE LIST OF ALL POSSIBLE REFLECTIONS
!==================================================

! generation of all potential reflections inside a reciprocal space sphere
! computed from the camera length and the detector size ...

! first set the maximum |g| value that can possibly give rise to a diffracted beam on the detector (diagonal)

gmax = sqrt(2.0) * float(pednl%npix) * pednl%rnmpp
io_real(1) = gmax
call WriteValue(' Length of longest g-vector : ', io_real, 1, "(F8.4)")

! this code is taken from the Initialize_ReflectionList routine, but we do not
! need everything from that routine
! get the size of the lookup table
gp = shape(cell%LUT)
imh = (gp(1)-1)/4
imk = (gp(2)-1)/4
iml = (gp(3)-1)/4

write (*,*) 'shape of LUT = ',shape(cell%LUT)

nullify(reflist)
nullify(rltmpa)
nref = 0

! transmitted beam has excitation error zero
gg = (/ 0,0,0 /)
call AddReflection(rltmpa, reflist, cell, nref, gg)   ! this guarantees that 000 is always the first reflection
rltmpa%xg = 0.0
xgmin = 100000.0
Igmax = 0.0

! now compute |U_g|^2 for all allowed reflections;
ixl: do ix=-imh,imh
    iyl:  do iy=-imk,imk
        izl:   do iz=-iml,iml
            if ((abs(ix)+abs(iy)+abs(iz)).ne.0) then  ! avoid double counting the origin
                gg = (/ ix, iy, iz /)
                glen = CalcLength(cell, float(gg), 'r' )

! find all reflections, ignoring double diffraction spots
                if ((IsGAllowed(cell,gg)).and.(glen.le.gmax).and.(glen.gt.0.0)) then ! allowed by the lattice centering, if any
                    call AddReflection(rltmpa, reflist, cell, nref, gg )
! we'll use the sangle field of the rltail structure to store |Ug|^2; we will also need the extinction distance
                    rltmpa%sangle = cdabs(cell%LUT(ix, iy, iz))**2
                    if (rltmpa%sangle.gt.Igmax) Igmax = rltmpa%sangle
                    rltmpa%xg = 1.0/(cdabs(cell%LUT(ix,iy,iz))*cell%mLambda)
                    if (rltmpa%xg.lt.xgmin) xgmin = rltmpa%xg
                end if ! IsGAllowed
            end if
        end do izl
    end do iyl
end do ixl

io_int(1) = nref
call WriteValue(' Length of the master list of reflections : ', io_int, 1, "(I8)")


!=============================================
!=============================================
! create the coordinate arrays for the Gaussian peaks
rnmpp = 1.0/pednl%rnmpp
ww = 4
tdp = 2*ww+1
allocate(xx(-ww:ww,-ww:ww), yy(-ww:ww,-ww:ww), line(-ww:ww), dot(-ww:ww,-ww:ww))
line = (/ (float(i),i=-ww,ww) /) * rnmpp
xx = spread(line,dim=1,ncopies=2*ww+1)
yy = transpose(xx)


!=============================================
!=============================================
! create the output array
nsize = pednl%npix/2 + ww
allocate(pedpattern(-nsize:nsize,-nsize:nsize))
allocate(image(pednl%npix,pednl%npix),imagevector(imght*imgwd))
maxint = Igmax
write (*,*) ' Maximum intensity = ',maxint


!=====================================================
! I/O FOR EXPERIMENTAL DATA SET
!=====================================================

call Message(' -> opened '//trim(pednl%exptfile)//' for I/O', frm = "(A)" )
open(unit=iunitexpt,file=trim(pednl%exptfile),status='old',form='unformatted',access='direct',recl=imght*imgwd,iostat=ierr)

Fztmp => FZlist

dictionaryloop: do ii = 1,ceiling(float(FZcnt)/float(numdictsingle))
    dict = 0.0
    result = 0.0
    if (ii .le. floor(float(FZcnt)/float(numdictsingle))) then
        do pp = 1,numdictsingle
            pedpattern = 0.0
            image = 0.0
            imagevector = 0.0
! convert the rodrigues vector to a passive rotation matrix.
            om = ro2om(FZtmp%rod)
! multiplication with (0,0,1) produces the normalized beam direction in a
! cartesian reference frame; so now we can compute the excitation errors
! for every reflection and keep only the ones that are sufficiently small
            k = (/ 0.0, 0.0, 1.0 /)
            ku = matmul(om,k)
            FN = ku
            k = ku/sngl(cell%mLambda)
! first we go through the entire reflection list and compute the excitation errors
! those points that satisfy the cutoff are linked via the nexts pointers
            rltmpa => reflist%next
            nexts => rltmpa
            do j=1,nref
                gg = rltmpa%hkl
                rltmpa%sg = Calcsg(cell,float(gg),k,FN)
! should we consider this point any further ? If so, add it to the strong reflection linked list
                if (abs(rltmpa%sg).le.sgmax) then
                    nexts%nexts => rltmpa
                    nexts => rltmpa
                end if
                rltmpa => rltmpa%next
            end do
! then, for each point in the nexts list, we compute the components of k' = k+g+s
! and place them in the proper reference frame; we skip the incident beam since it is
! meaningless in the kinematical approximation
            nexts => reflist%next%nexts
            do
! determine the vector k'
                kp = k + float(nexts%hkl) + nexts%sg*ku
                kp = matmul(transpose(om),kp)
! get the intensity for each point
                w = sngl(cPi)*nexts%sg*pednl%thickness
                if (abs(w).lt.1.0e-6) then
                    Ig = nexts%sangle  ! * (sngl(cPi)*pednl%thickness/nexts%xg)**2
                else
                    Ig = nexts%sangle * (sin(w)/w)**2 ! * (sngl(cPi)*pednl%thickness/nexts%xg)**2
                end if

! determine the spot coordinates on the detector
                x = rnmpp * kp(1)
                y = rnmpp * kp(2)

! and plot that spot as a small Gaussian in the pedpattern array, assuming it falls on the detector.
                if ((abs(x).le.nsize-ww).and.(abs(y).le.nsize-ww)) then
                    sx = nint(x)
                    sy = nint(y)
                    dx = x-sx
                    dy = y-sy
                    dot = (Ig/Igmax)**0.2 * exp(-((xx-dx)**2+(yy-dy)**2)*0.003)
                    pedpattern(sx-ww:sx+ww,sy-ww:sy+ww) = pedpattern(sx-ww:sx+ww,sy-ww:sy+ww) + dot(-ww:ww,-ww:ww)
                end if

! and repeat this until the end of the list
                if (.not. associated(nexts%nexts)) EXIT
                nexts => nexts%nexts
            end do

! save the pedpattern to file
            image(1:pednl%npix,1:pednl%npix) = pedpattern(-nsize+ww:nsize-ww,-nsize+ww:nsize-ww)

            do qq = 1,imght
                imagevector((qq-1)*imgwd+1:qq*imgwd) = image(qq,1:imgwd)
            end do

            dict((pp-1)*L+1:pp*L) = imagevector(1:L)/NORM2(imagevector(1:L))
!if (pp .eq. 1 .and. ii .eq. 1) then
!test = imagevector/NORM2(imagevector)
!open(unit=13,file='testdict.txt',action='write')
!do kk = 1,imght
!do ll = 1,imgwd
!write(13, '(F15.6)', advance='no') image(ll,kk)
!end do
!write(13, *) ''  ! this gives you the line break
!end do
!close(13)
!end if
! reset the nexts linked list and start over
            nexts => reflist%next
            rltmpa => nexts%nexts
            do
                nullify(nexts%nexts)
                if (.not. associated(rltmpa%nexts)) EXIT
                    nexts => rltmpa
                    rltmpa => rltmpa%nexts
            end do
            eulerarray(1:3,(ii-1)*numdictsingle+pp) = 180.0/cPi*ro2eu(FZtmp%rod)
            FZtmp => FZtmp%next                  ! point to the next entry
        end do

    else if (ii .eq. ceiling(float(FZcnt)/float(numdictsingle))) then
        do pp = 1,MODULO(FZcnt,numdictsingle)

            pedpattern = 0.0
            image = 0.0
            imagevector = 0.0
! convert the rodrigues vector to a passive rotation matrix.
            om = ro2om(FZtmp%rod)
! multiplication with (0,0,1) produces the normalized beam direction in a
! cartesian reference frame; so now we can compute the excitation errors
! for every reflection and keep only the ones that are sufficiently small
            k = (/ 0.0, 0.0, 1.0 /)
            ku = matmul(om,k)
            FN = ku
            k = ku/sngl(cell%mLambda)
! first we go through the entire reflection list and compute the excitation errors
! those points that satisfy the cutoff are linked via the nexts pointers
            rltmpa => reflist%next
            nexts => rltmpa
            do j=1,nref
                gg = rltmpa%hkl
                rltmpa%sg = Calcsg(cell,float(gg),k,FN)
! should we consider this point any further ? If so, add it to the strong reflection linked list
                if (abs(rltmpa%sg).le.sgmax) then
                    nexts%nexts => rltmpa
                    nexts => rltmpa
                end if
                rltmpa => rltmpa%next
            end do

! then, for each point in the nexts list, we compute the components of k' = k+g+s
! and place them in the proper reference frame; we skip the incident beam since it is
! meaningless in the kinematical approximation
            nexts => reflist%next%nexts
            do
! determine the vector k'
                kp = k + float(nexts%hkl) + nexts%sg*ku
                kp = matmul(transpose(om),kp)

! get the intensity for each point
                w = sngl(cPi)*nexts%sg*pednl%thickness
                if (abs(w).lt.1.0e-6) then
                    Ig = nexts%sangle  ! * (sngl(cPi)*pednl%thickness/nexts%xg)**2
                else
                    Ig = nexts%sangle * (sin(w)/w)**2 ! * (sngl(cPi)*pednl%thickness/nexts%xg)**2
                end if

! determine the spot coordinates on the detector
                x = rnmpp * kp(1)
                y = rnmpp * kp(2)

! and plot that spot as a small Gaussian in the pedpattern array, assuming it falls on the detector.
                if ((abs(x).le.nsize-ww).and.(abs(y).le.nsize-ww)) then
                    sx = nint(x)
                    sy = nint(y)
                    dx = x-sx
                    dy = y-sy
                    dot = (Ig/Igmax)**0.2 * exp(-((xx-dx)**2+(yy-dy)**2)*0.003)
                    pedpattern(sx-ww:sx+ww,sy-ww:sy+ww) = pedpattern(sx-ww:sx+ww,sy-ww:sy+ww) + dot(-ww:ww,-ww:ww)
                end if

! and repeat this until the end of the list
                if (.not. associated(nexts%nexts)) EXIT
                nexts => nexts%nexts
            end do

! save the pedpattern to file
                image(1:pednl%npix,1:pednl%npix) = pedpattern(-nsize+ww:nsize-ww,-nsize+ww:nsize-ww)
                do qq = 1,imght
                    imagevector((qq-1)*imgwd+1:qq*imgwd) = image(qq,1:imgwd)
                end do

                dict((pp-1)*L+1:pp*L) = imagevector(1:L)/NORM2(imagevector(1:L))
! reset the nexts linked list and start over
                nexts => reflist%next
                rltmpa => nexts%nexts
                do
                    nullify(nexts%nexts)
                    if (.not. associated(rltmpa%nexts)) EXIT
                        nexts => rltmpa
                        rltmpa => rltmpa%nexts
                end do
        eulerarray(1:3,(ii-1)*numdictsingle+pp) = 180.0/cPi*ro2eu(FZtmp%rod)
        FZtmp => FZtmp%next                  ! point to the next entry
        end do

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
                imageexptflt = float(imageexpt)
                do mm = 1,L
                    if (imageexptflt(mm) .lt. 0) imageexptflt(mm) = imageexptflt(mm) + 256.0
                end do

                expt((pp-1)*L+1:pp*L) = imageexptflt(1:L)/NORM2(imageexptflt(1:L))
!if (pp .eq. 383 .and. jj .eq. 1) then
!test1 = imageexptflt/NORM2(imageexptflt)
!open(unit=13,file='testexpt.txt',action='write')
!do kk = 1,imght
!do ll = 1,imgwd
!write(13, '(F15.6)', advance='no') expt((kk-1)*imgwd+ll)
!end do
!write(13, *) ''  ! this gives you the line break
!end do
!close(13)
!end if
            end do



        else if (jj .eq. ceiling(float(totnumexpt)/float(numexptsingle))) then

            do pp = 1,MODULO(totnumexpt,numexptsingle)
                read(iunitexpt,rec=(jj-1)*numexptsingle+pp) imageexpt
                imageexptflt = float(imageexpt)
                do mm = 1,L
                    if (imageexptflt(mm) .lt. 0) imageexptflt(mm) = imageexptflt(mm) + 256.0
                end do

                expt((pp-1)*L+1:pp*L) = imageexptflt(1:L)/NORM2(imageexptflt(1:L))

            end do

        end if

        call clEnqueueWriteBuffer(command_queue, cl_expt, cl_bool(.true.), 0_8, size_in_bytes_expt, expt(1), ierr)
        if(ierr /= CL_SUCCESS) stop 'Error: cannot write to buffer.'

        call InnerProdGPU(cl_expt,cl_dict,Ne,Nd,L,result,source,source_length,platform,device,context,command_queue)

        !if (mod(jj,floor(float(totnumexpt)/float(numexptsingle)/10.0)) .eq. 0) then
        !    if (ii .le. floor(float(FZcnt)/float(numdictsingle))) then
        !        write(6,'(A,I8,A,I8,A)')'Completed dot product of ',jj*numexptsingle,' experimental pattern with ',&
        !        ii*numdictsingle,' dictionary patterns.'
        !    else if (ii .eq. ceiling(float(FZcnt)/float(numdictsingle))) then
        !        write(6,'(A,I8,A,I8,A)')'Completed dot product of ',jj*numexptsingle,'experimental pattern with all ',&
        !        FZcnt,'dictionary patterns'
        !    end if
        !end if

        !write(6,'(A)'),' -> Starting sorting of the dot products'



!$OMP PARALLEL PRIVATE(TID,qq,resultarray,indexarray) &
!$OMP& SHARED(nthreads,result,resultmain,indexmain,ii,jj,numdictsingle,numexptsingle,indexlist,resulttmp,indextmp)
        TID = OMP_GET_THREAD_NUM()
        nthreads = OMP_GET_NUM_THREADS()
        if (TID .eq. 0 .and. ii .eq. 1 .and. jj .eq. 1) write(6,'(A20,I3)'),'Number of threads =',nthreads
!$OMP BARRIER
!$OMP DO SCHEDULE(DYNAMIC)
        do qq = 1,numexptsingle

            resultarray(1:numdictsingle) = result((qq-1)*numdictsingle:qq*numdictsingle)

            indexarray(1:numdictsingle) = indexlist((jj-1)*numdictsingle:jj*numdictsingle)

            call SSORT(resultarray,indexarray,numdictsingle,-2)

            resulttmp(nnk+1:2*nnk,(jj-1)*numexptsingle+qq) = resultarray(1:nnk)
            indextmp(nnk+1:2*nnk,(jj-1)*numexptsingle+qq) = indexarray(1:nnk)

            call SSORT(resulttmp(:,(jj-1)*numexptsingle+qq),indextmp(:,(jj-1)*numexptsingle+qq),2*nnk,-2)

            resultmain(1:nnk,(jj-1)*numexptsingle+qq) = resulttmp(1:nnk,(jj-1)*numexptsingle+qq)
            indexmain(1:nnk,(jj-1)*numexptsingle+qq) = indextmp(1:nnk,(jj-1)*numexptsingle+qq)

        end do
!$OMP END DO
!$OMP END PARALLEL
!print*,resultmain(1:5,1),indexmain(1:5,1)

    end do experimentalloop

    if (ii .le. floor(float(FZcnt)/float(numdictsingle))) then

        write(6,'(A,I8,A,I8,A)'),'Completed dot product of all ',totnumexpt,&
        ' experimental patterns with ',ii*numexptsingle,' dictionary patterns'

    else if (ii .eq. ceiling(float(FZcnt)/float(numdictsingle))) then

        write(6,'(A,I8,A,I8,A)'),'Completed dot product of all ',FZcnt,' dictionary patterns with all ',totnumexpt,&
            ' experimental patterns'


    end if

end do dictionaryloop

open(unit=iunit,file='Results.ctf',action='write')
call WriteHeader(iunit,'ctf')

do ii = 1,totnumexpt

    !do jj = 1,nnk
    !    index = indexarray(jj,ii)
    !    samples(1:4,jj) = eu2qu((cPi/180.D0)*eulerangles(index,1:3))
    !end do

    !call DI_EMforDD(samples, dictlist, nnk, seed, muhat, kappahat,'WAT')
    !euler = 180.0/cPi*qu2eu(muhat)

    index = indexmain(1,ii)
    euler = eulerarray(1:3,index)

    write(str1,'(F12.3)') float(MODULO(ii-1,180))
    write(str2,'(F12.3)') float(floor(float(ii-1)/180.0))
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
write(iunit,'(3A)'),'XCells',Tab,'180'
write(iunit,'(3A)'),'YCells',TAB,'210'
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
