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

program EMMCOpenCL

use local
use files
use NameListTypedefs
use NameListHandlers
use io

IMPLICIT NONE

character(fnlen)                            :: nmldeffile, progname, progdesc
type(MCCLMultiLayerNameListType)                      :: mcnl

nmldeffile = 'EMMCCLMultiLayer.nml'
progname = 'EMMC.f90'
progdesc = 'Monte Carlo backscattered electron simulation'

! deal with the command line arguments, if any
call Interpret_Program_Arguments(nmldeffile,1,(/ 20 /), progname)

! deal with the namelist stuff
call GetMCCLMultiLayerNameList(nmldeffile,mcnl)

! print some information
call EMsoft(progname, progdesc)

! perform a Monte Carlo simulation
call DoMCsimulation(mcnl, progname)

end program EMMCOpenCL

!--------------------------------------------------------------------------
!
! SUBROUTINE:DoMCsimulation
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief Perform the MC simulation
!
!> @param nmlfile namelist file name
!
!> @date 11/29/01  MDG 1.0 original
!> @date 04/08/13  MDG 2.0 rewrite
!> @date 05/14/13  MDG 2.1 replaced IO by namelist file
!> @date 07/23/13  MDG 3.0 complete rewrite
!> @date 07/30/13  MDG 3.1 added Patrick's code for double sample tilt (sigma, omega)
!> @date 09/25/13  MDG 3.2 added a few parameters to the output file
!> @date 03/17/14  MDG 3.3 added a few more for the IDL visualization program
!> @date 06/19/14  MDG 4.0 rewrite with name list handling removed
!> @date 07/23/14  SS  4.1 conversion to OpenCL
!--------------------------------------------------------------------------
subroutine DoMCsimulation(mcnl, progname)

use local
use typedefs
use NameListTypedefs
use initializers
use crystal
use constants
use symmetry
use error
use io
use files
use diffraction, only:CalcWaveLength
use Lambert
use cl

IMPLICIT NONE

type(MCCLMultiLayerNameListType),INTENT(IN)         :: mcnl
character(fnlen),INTENT(IN)             :: progname


type(unitcell),pointer  :: cell
type(DynType)           :: Dyn
type(gnode)             :: rlp

integer(kind=irg)       :: numsy        ! number of Lambert map points along y
integer(kind=irg)       :: numEbins     ! number of energy bins
integer(kind=irg)       :: numzbins     ! number of depth bins
integer(kind=irg)       :: nx           ! no. of pixels
integer(kind=irg)       :: i,j,k,l,ip,istat
integer,parameter       :: k12 = selected_int_kind(15)

real(kind=4)            :: Ze(2)           ! average atomic number
real(kind=4)            :: density(2)      ! density in g/cm^3
real(kind=4)            :: at_wt(2)        ! average atomic weight in g/mole
logical                 :: verbose
real(kind=4),allocatable:: thickness(:),lambdaZ(:)
real(kind=4)          :: dens, avA, avZ, io_real(3), dmin, tf ! used with CalcDensity routine
real(kind=8) , parameter         :: dtoR = 0.01745329251D0 !auxiliary variables
real(kind=4)          :: EkeV, sig, omega ! input values to the kernel. Can only be real kind=4 otherwise values are not properly passed
integer(kind=8)       :: totnum_el, val, totnum_el_out   ! total number of electrons to simulate
integer(kind=4)       :: totnum_thickness
integer(kind=4)       :: globalworkgrpsz, num_el, num_max, steps, prime ! input values to the kernel
integer(kind=8)       size_in_bytes,size_in_bytes_seeds,size_in_bytes_input ! size of arrays passed to kernel. Only accepts kind=8 integers by clCreateBuffer etc., so donot change
integer(kind=8)         :: globalsize(2), localsize(2) ! size of global and local work groups. Again only kind=8 is accepted by clEnqueueNDRangeKernel
character(4)            :: mode
! results from kernel stored here
real(kind=4),allocatable :: Lamresx(:), Lamresy(:), depthres(:), energyres(:)

! final results stored here
integer(kind=4),allocatable :: accum_e(:,:,:), accum_z(:,:,:,:), rnseeds(:), init_seeds(:)
integer(kind=4)         :: idxy(2), iE, px, py, iz, nseeds, randint ! auxiliary variables
real(kind=4)            :: cxyz(3), edis, bse, xy(2) ! auxiliary variables
real(kind=8)            :: delta,rand


! OpenCL variables
type(cl_platform_id)    :: platform
type(cl_device_id)      :: device
type(cl_context)        :: context
type(cl_command_queue)  :: command_queue
type(cl_program)        :: prog
type(cl_kernel)         :: kernel
type(cl_mem)            :: LamX, LamY, depth, energy, seeds, cl_Ze, cl_at_wt, cl_density
type(cl_event)          :: event

character(len = 100)    :: info, filename ! info about the GPU
integer, parameter      :: iunit = 10
integer, parameter              :: source_length = 10000000
character(len = source_length)  :: source
integer(kind=4)       :: num, ierr, irec, io_int(1), ii ! auxiliary variables
numsy = mcnl%numsx
nullify(cell)
allocate(cell)

! get the crystal strucutre from the *.xtal file
verbose = .TRUE.
dmin = 0.05
val = 0
call Initialize_Cell(cell,Dyn,rlp,mcnl%xtalname_film, dmin, sngl(1000.D0*mcnl%EkeV), verbose)
! then calculate density, average atomic number and average atomic weight
call CalcDensity(cell, dens, avZ, avA)
density(1) = dble(dens)
Ze(1) = dble(avZ)
at_wt(1) = dble(avA)
io_real(1:3) = (/ dens, avZ, avA /)
call WriteValue('Density, avZ, avA for film = ',io_real,3,"(2f10.5,',',f10.5)")

call Initialize_Cell(cell,Dyn,rlp,mcnl%xtalname_subs, dmin, sngl(1000.D0*mcnl%EkeV), verbose)
! then calculate density, average atomic number and average atomic weight
call CalcDensity(cell, dens, avZ, avA)
density(2) = dble(dens)
Ze(2) = dble(avZ)
at_wt(2) = dble(avA)
io_real(1:3) = (/ dens, avZ, avA /)
call WriteValue('Density, avZ, avA for substrate = ',io_real,3,"(2f10.5,',',f10.5)")
mode = mcnl%mode

totnum_thickness = (mcnl%filmthickness/mcnl%filmstep)
allocate(thickness(totnum_thickness),stat=istat)

do i = 1,totnum_thickness+1
    thickness(i) = (i-1)*mcnl%filmstep
end do

if (mode .eq. 'full') then
    steps = 600
else if (mode .eq. 'bse1') then
    steps = 1
else
    stop 'Unknown mode specified in namelist file'
end if


EkeV = mcnl%EkeV
sig = mcnl%sig*dtoR
omega = mcnl%omega*dtoR
globalworkgrpsz = mcnl%globalworkgrpsz
num_el = mcnl%num_el ! no. of electron simulation by one work item
num_max = globalworkgrpsz*globalworkgrpsz*num_el ! total simulation in one loop
totnum_el = mcnl%totnum_el ! total number of electrons to simulate
globalsize = (/ mcnl%globalworkgrpsz, mcnl%globalworkgrpsz /)
localsize = (/ mcnl%globalworkgrpsz/10, mcnl%globalworkgrpsz/10 /)

numEbins =  int((mcnl%EkeV-mcnl%Ehistmin)/mcnl%Ebinsize)+1
numzbins =  int(mcnl%depthmax/mcnl%depthstep)+1
nx = (mcnl%numsx-1)/2

allocate(Lamresx(num_max), Lamresy(num_max), depthres(num_max), energyres(num_max), stat=istat)
Lamresx = 0.0
Lamresy = 0.0
depthres = 0.0
energyres = 0.0
size_in_bytes = num_max*sizeof(EkeV)
size_in_bytes_input = sizeof(Ze)
size_in_bytes_seeds = 4*globalworkgrpsz*globalworkgrpsz*sizeof(EkeV)

delta = dble(nx)/LPs%sPio2

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
open(unit = iunit, file = trim(openclpathname)//'/EMMC-multilayer.cl', access='direct', status = 'old', &
        action = 'read', iostat = ierr, recl = 1)
if (ierr /= 0) stop 'Cannot open file EMMC.cl'
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
kernel = clCreateKernel(prog, 'MC', ierr)
call clReleaseProgram(prog, ierr)

open(unit = iunit, file = mcnl%primelist, form='unformatted', status='old')
read(iunit) nseeds
allocate(rnseeds(nseeds))
read(iunit) rnseeds
close(unit=iunit,status='keep')

if (globalworkgrpsz**2 .gt. nseeds) call FatalError('EMMCOpenCL:','insufficient prime numbers')



allocate(init_seeds(4*globalworkgrpsz*globalworkgrpsz),stat=istat)
init_seeds = 0
do i = 1,globalworkgrpsz
    do j = 1,globalworkgrpsz
        do k = 1,4
            init_seeds(4*((i-1)*globalworkgrpsz+j)+k) = rnseeds(4*((i-1)*globalworkgrpsz+j)+k)
        end do
    end do
end do
!print*,init_seeds

! create device memory buffers
LamX = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size_in_bytes, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

LamY = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size_in_bytes, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

depth = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size_in_bytes, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

energy = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size_in_bytes, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

seeds = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

cl_Ze = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_input, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

cl_density = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_input, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

cl_at_wt = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes_input, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory.'

!call init_random_seed()
!call clEnqueueWriteBuffer(command_queue, seeds, cl_bool(.true.), 0_8, size_in_bytes_seeds, init_seeds(1), ierr)
call clEnqueueWriteBuffer(command_queue, cl_Ze, cl_bool(.true.), 0_8, size_in_bytes_input, Ze(1), ierr)
call clEnqueueWriteBuffer(command_queue, cl_density, cl_bool(.true.), 0_8, size_in_bytes_input, density(1), ierr)
call clEnqueueWriteBuffer(command_queue, cl_at_wt, cl_bool(.true.), 0_8, size_in_bytes_input, at_wt(1), ierr)
!call clEnqueueReadBuffer(command_queue, cl_at_wt, cl_bool(.true.), 0_8, size_in_bytes_input, Ze(1), ierr)

thicknessloop: do k = 1,totnum_thickness+1
    allocate(accum_e(numEbins,-nx:nx,-nx:nx),accum_z(numEbins,numzbins,-nx/10:nx/10,-nx/10:nx/10),stat=istat)
    accum_e = 0
    accum_z = 0
    allocate(lambdaZ(1:numzbins),stat=istat)
    lambdaZ = 0.0

    randomnumloop: do l = 1,10
        call RANDOM_NUMBER(rand)
        randint = rand*1e8
        init_seeds = init_seeds+randint
        call clEnqueueWriteBuffer(command_queue, seeds, cl_bool(.true.), 0_8, size_in_bytes_seeds, init_seeds(1), ierr)

        mainloop: do i = 1,(totnum_el/num_max+1)

            tf = thickness(k)


! set the kernel arguments
            call clSetKernelArg(kernel, 0, LamX, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument 0.'

            call clSetKernelArg(kernel, 1, LamY, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument 1.'

            call clSetKernelArg(kernel, 2, EkeV, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument 2.'

            call clSetKernelArg(kernel, 3, globalworkgrpsz, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument 3.'

            call clSetKernelArg(kernel, 4, cl_Ze, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument 4.'

            call clSetKernelArg(kernel, 5, cl_density, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument 5.'

            call clSetKernelArg(kernel, 6, cl_at_wt, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument 6.'

            call clSetKernelArg(kernel, 7, num_el, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument 7.'

            call clSetKernelArg(kernel, 8, seeds, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument 8.'

            call clSetKernelArg(kernel, 9, sig, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument 9.'

            call clSetKernelArg(kernel, 10, omega, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument 10.'

            call clSetKernelArg(kernel, 11, depth, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument 11.'

            call clSetKernelArg(kernel, 12, energy, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument 12.'

            call clSetKernelArg(kernel, 13, steps, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument 13.'

            call clSetKernelArg(kernel, 14, tf, ierr)
            if(ierr /= CL_SUCCESS) stop 'Error: cannot set kernel argument 14.'


!print*, EkeV, globalworkgrpsz, Ze, density, at_wt, num_max, prime, sig, omega
! execute the kernel
            call clEnqueueNDRangeKernel(command_queue, kernel, globalsize, localsize, event, ierr)

! wait for the commands to finish
            call clFinish(command_queue, ierr)

! read the resulting vector from device memory
            call clEnqueueReadBuffer(command_queue, LamX, cl_bool(.true.), 0_8, size_in_bytes, Lamresx(1), ierr)
            call clEnqueueReadBuffer(command_queue, LamY, cl_bool(.true.), 0_8, size_in_bytes, Lamresy(1), ierr)
            call clEnqueueReadBuffer(command_queue, depth, cl_bool(.true.), 0_8, size_in_bytes, depthres(1), ierr)
            call clEnqueueReadBuffer(command_queue, energy, cl_bool(.true.), 0_8, size_in_bytes, energyres(1), ierr)
!    call clEnqueueReadBuffer(command_queue, seeds, cl_bool(.true.), 0_8, size_in_bytes_seeds, init_seeds(1), ierr)

            subloop: do j = 1, num_max

                if ((Lamresx(j) .ne. -10.0) .and. (Lamresy(j) .ne. -10.0) &
                .and. (depthres(j) .ne. 10.0) .and. (energyres(j) .ne. 0.0) &
                .and. .not.isnan(Lamresx(j)) .and. .not.isnan(Lamresy(j))) then
        !& .and. (Lamresx(j) .ne. 0.0) .and. (Lamresy(j) .ne. 0.0) .and. (depthres(j) .ne. 0.0)) then
! and get the nearest pixel [ take into account reversal of coordinate frame (x,y) -> (y,-x) ]
                    !if ((nint(delta*Lamresy(j)) .eq. 0.0) .and. (nint(-delta*Lamresx(j)) .eq. 0.0)) then
                    !val1 = val1 + 1
                !print*,Lamresx(j),Lamresy(j)
                    !end if

                    val = val + 1
                    idxy = (/ nint(delta*Lamresy(j)), nint(-delta*Lamresx(j)) /)

                    if (maxval(abs(idxy)).le.nx) then
! If Ec larger than Emin, then we should count this electron
                        if (energyres(j).gt.mcnl%Ehistmin) then
!print *, energyres(j), depthres(j)*1e9, Lamresx(j), Lamresy(j)

                            iE = nint((energyres(j)-mcnl%Ehistmin)/mcnl%Ebinsize)+1
! first add this electron to the correct exit distance vs. energy bin (coarser than the angular plot)
                            edis = abs(depthres(j))  ! distance from last scattering point to surface along trajectory
                            iz = nint(edis/mcnl%depthstep) +1
                            if ( (iz.gt.0).and.(iz.le.numzbins) ) then

                                px = nint(idxy(1)/10.0)
                                py = nint(idxy(2)/10.0)
                                accum_z(iE,iz,px,py) = accum_z(iE,iz,px,py) + 1

                            end if
! then add it to the modified Lambert accumulator array.
                            accum_e(iE,idxy(1),idxy(2)) = accum_e(iE,idxy(1),idxy(2)) + 1
                        end if
                    end if
                end if
            end do subloop


        end do mainloop


    io_int = l
    print*,''
    print*,'Finished loop number',l,'for film thickness of',thickness(k),'nm.'
    print*,''

    end do randomnumloop

    ii = totnum_el/num_max
    totnum_el_out = (ii+1)*num_max*10

    print*,"Total number of electron for the",l-1,"runs for different random seeds for film thickness",&
        thickness(k),"=",totnum_el_out
    print*,''
    io_int(1) = sum(accum_e)
    print*,"Total number of electrons on detector for film thickness of",thickness(k),"nm =",io_int
    print*,''
    bse = real(val)/real(totnum_el)
    print*, "Backscatter yield = ",bse


    do iz=1,numzbins
        lambdaZ(iz) = float(sum(accum_z(:,iz,:,:)))/float(sum(accum_z(:,:,:,:)))
    end do
    if (int(thickness(k)) .lt. 10) then
        write(filename,"(A6,I1)")"lambda",int(thickness(k))
    else if (int(thickness(k)) .lt. 100) then
        write(filename,"(A6,I2)")"lambda",int(thickness(k))
    else
        write(filename,"(A6,I3)")"lambda",int(thickness(k))
    end if
    open(unit=13,file=trim(filename),action='write')
    write(13,"(F15.6)") lambdaZ
    close(13)

    if ( int(thickness(k)) .lt. 10) then
        write(filename,"(A13,I1)")mcnl%dataname,int(thickness(k))
    else if (int(thickness(k)) .lt. 100) then
        write(filename,"(A13,I2)")mcnl%dataname,int(thickness(k))
    else
        write(filename,"(A13,I3)")mcnl%dataname,int(thickness(k))
    end if
    open(dataunit,file=trim(filename),status='unknown',form='unformatted')

! and here we create the output file
    call Message(' ',"(A)")
! write the program identifier
    write (dataunit) progname
! write the version number
    write (dataunit) scversion
! then the name of the crystal data file
    write (dataunit) mcnl%xtalname_film
! energy information etc...
    write (dataunit) numEbins, numzbins, mcnl%numsx, numsy, totnum_el_out
    write (dataunit) mcnl%EkeV, mcnl%Ehistmin, mcnl%Ebinsize, mcnl%depthmax, mcnl%depthstep
    write (dataunit) mcnl%sig, mcnl%omega
    write (dataunit) mcnl%MCmode
! and here are the actual results
    write (dataunit) accum_e
    write (dataunit) accum_z

    close(dataunit,status='keep')

    !print*, "Total number of incident electrons = ",totnum_el
    !print*, "Number of electrons on detector = ",sum(accum_e)
    !print*, "Maximum electron in depth and energy bin resp. = ",maxval(accum_z),maxval(accum_e)
    deallocate(accum_e,accum_z,lambdaZ)

end do thicknessloop




!print*,sum(accum_e(:,0,0)),maxval(sum(accum_e,1))
!=====================
! RELEASE EVERYTHING
!=====================

call clReleaseKernel(kernel, ierr)
call clReleaseCommandQueue(command_queue, ierr)
call clReleaseContext(context, ierr)
call clReleaseMemObject(LamX, ierr)
call clReleaseMemObject(LamY, ierr)
call clReleaseMemObject(depth, ierr)
call clReleaseMemObject(energy, ierr)
call clReleaseMemObject(seeds, ierr)


end subroutine DoMCsimulation
