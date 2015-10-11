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
use JSONsupport
use io

IMPLICIT NONE

character(fnlen)                        :: nmldeffile, progname, progdesc
type(MCCLNameListType)                  :: mcnl
integer(kind=irg)                       :: res, error_cnt

nmldeffile = 'EMMCOpenCL.nml'
progname = 'EMMCOpenCL.f90'
progdesc = 'Monte Carlo backscattered electron simulation'

! print some information
call EMsoft(progname, progdesc)

! deal with the command line arguments, if any
call Interpret_Program_Arguments(nmldeffile,1,(/ 42 /), progname)

! deal with the namelist stuff, either .nml or .json format
res = index(nmldeffile,'.nml',kind=irg)
if (res.eq.0) then
  call JSONreadMCCLNameList(mcnl, nmldeffile, error_cnt)
else
  call GetMCCLNameList(nmldeffile,mcnl)
end if

! perform a Monte Carlo simulation
call DoMCsimulation(mcnl, progname, nmldeffile)

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
!> @date 03/26/15  MDG 5.0 all output now in HDF5 format 
!> @date 05/05/15  MDG 5.1 removed getenv() call; replaced by global path string
!> @date 09/01/15  MDG 5.2 modifications due to Lambert module changes
!> @date 09/09/15  MDG 5.3 added devid selector (GPU device ID) to namelist
!> @date 10/07/15  MDG 5.4 general cleanup of code in preparation of release 3.0
!--------------------------------------------------------------------------
subroutine DoMCsimulation(mcnl, progname, nmldeffile)

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
use HDF5
use NameListHDFwriters
use HDFsupport
use ISO_C_BINDING

IMPLICIT NONE

type(MCCLNameListType),INTENT(INOUT)    :: mcnl
character(fnlen),INTENT(IN)             :: progname
character(fnlen),INTENT(IN)             :: nmldeffile


type(unitcell),pointer  :: cell
type(DynType)           :: Dyn
type(gnode)             :: rlp

integer(kind=irg)       :: numsy        ! number of Lambert map points along y
integer(kind=irg)       :: numEbins     ! number of energy bins
integer(kind=irg)       :: numzbins     ! number of depth bins
integer(kind=irg)       :: nx           ! no. of pixels
integer(kind=irg)       :: i,j,k,l,ip,istat
integer,parameter       :: k12 = selected_int_kind(15)

real(kind=4)            :: Ze           ! average atomic number
real(kind=4)            :: density      ! density in g/cm^3
real(kind=4)            :: at_wt        ! average atomic weight in g/mole
logical                 :: verbose
real(kind=4)            :: dens, avA, avZ, io_real(3), dmin ! used with CalcDensity routine
real(kind=8) , parameter:: dtoR = 0.01745329251D0 !auxiliary variables
real(kind=4)            :: EkeV, sig, omega ! input values to the kernel. Can only be real kind=4 otherwise values are not properly passed
integer(kind=4)         :: totnum_el     ! total number of electrons to simulate
integer(kind=4)         :: globalworkgrpsz, num_el, num_max, steps, prime ! input values to the kernel
integer(kind=8)         :: size_in_bytes,size_in_bytes_seeds ! size of arrays passed to kernel. Only accepts kind=8 integers by clCreateBuffer etc., so donot change
integer(kind=8)         :: globalsize(2), localsize(2) ! size of global and local work groups. Again only kind=8 is accepted by clEnqueueNDRangeKernel
character(4)            :: mode
! results from kernel stored here
real(kind=4),allocatable:: Lamresx(:), Lamresy(:), depthres(:), energyres(:)

! final results stored here
integer(kind=4),allocatable :: accum_e(:,:,:), accum_z(:,:,:,:), rnseeds(:), init_seeds(:)
integer(kind=4)         :: idxy(2), iE, px, py, iz, nseeds, hdferr ! auxiliary variables
real(kind=4)            :: cxyz(3), edis, bse, xy(2), tstart, tstop ! auxiliary variables
real(kind=8)            :: delta,rand
character(11)           :: dstr
character(15)           :: tstrb
character(15)           :: tstre
logical                 :: f_exists

! OpenCL variables
type(cl_platform_id)    :: platform
type(cl_device_id)      :: device(2)
type(cl_context)        :: context
type(cl_command_queue)  :: command_queue
type(cl_program)        :: prog
type(cl_kernel)         :: kernel
type(cl_mem)            :: LamX, LamY, depth, energy,seeds
type(cl_event)          :: event

character(len = 100)    :: info ! info about the GPU
integer, parameter      :: iunit = 10
integer, parameter      :: source_length = 10000000
character(len = source_length)  :: source
integer(kind=4)         :: num, ierr, irec, io_int(2), val,val1 ! auxiliary variables
character(fnlen)        :: groupname, dataset, instring, dataname


type(HDFobjectStackType),pointer  :: HDF_head

nullify(HDF_head)

call timestamp(datestring=dstr, timestring=tstrb)
call CPU_TIME(tstart)

numsy = mcnl%numsx
nullify(cell)
allocate(cell)

! get the crystal strucutre from the *.xtal file
verbose = .TRUE.
dmin = 0.05
val = 0
val1 = 0
call Initialize_Cell(cell,Dyn,rlp,mcnl%xtalname, dmin, sngl(1000.D0*mcnl%EkeV), verbose)

! then calculate density, average atomic number and average atomic weight
call CalcDensity(cell, dens, avZ, avA)
density = dble(dens)
Ze = dble(avZ)
at_wt = dble(avA)
io_real(1:3) = (/ dens, avZ, avA /)
call WriteValue('Density, avZ, avA = ',io_real,3,"(2f10.5,',',f10.5)")
mode = mcnl%mode

if (mode .eq. 'full') then
    steps = 600
else if (mode .eq. 'bse1') then
    steps = 1
else
    stop 'Unknown program mode specified in namelist file'
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
allocate(accum_e(numEbins,-nx:nx,-nx:nx),accum_z(numEbins,numzbins,-nx/10:nx/10,-nx/10:nx/10),stat=istat)
accum_e = 0
accum_z = 0

allocate(Lamresx(num_max), Lamresy(num_max), depthres(num_max), energyres(num_max), stat=istat)
Lamresx = 0.0
Lamresy = 0.0
depthres = 0.0
energyres = 0.0
size_in_bytes = num_max*sizeof(EkeV)
size_in_bytes_seeds = 4*globalworkgrpsz*globalworkgrpsz*sizeof(EkeV)

! changed by MDG [09/01/15] after extensive modifications to Lambert routines
! old code delta = dble(nx)/LPs%sPio2
delta = dble(nx)

!=====================
! INITIALIZATION
!=====================

! get the platform ID
call clGetPlatformIDs(platform, num, ierr)
if(ierr /= CL_SUCCESS) stop "Fatal error: Cannot get CL platform."

! get the device ID
call clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, device, num, ierr)
if(ierr /= CL_SUCCESS) stop "Fatal error: Cannot get CL device."

! get the device name and print it
call clGetDeviceInfo(device(mcnl%devid), CL_DEVICE_NAME, info, ierr)
call WriteValue("CL device: ", info)

! create the context and the command queue
context = clCreateContext(platform, device(mcnl%devid), ierr)
if(ierr /= CL_SUCCESS) stop "Fatal error: Cannot create context"
command_queue = clCreateCommandQueue(context, device(mcnl%devid), CL_QUEUE_PROFILING_ENABLE, ierr)
if(ierr /= CL_SUCCESS) stop "Fatal error: Cannot create command queue"


!=====================
! BUILD THE KERNEL
!=====================

! read the source file
open(unit = iunit, file = trim(openclpathname)//'EMMC.cl', access='direct', status = 'old', &
        action = 'read', iostat = ierr, recl = 1)
if (ierr /= 0) stop 'Fatal error: Cannot open file EMMC.cl in opencl kernel folder'

source = ''
irec = 1
do
read(unit = iunit, rec = irec, iostat = ierr) source(irec:irec)
if (ierr /= 0) exit
if(irec == source_length) stop 'Fatal error: CL source file is too big'
irec = irec + 1
end do
close(unit=iunit)

! create the program
prog = clCreateProgramWithSource(context, source, ierr)
if(ierr /= CL_SUCCESS) stop 'Fatal error: cannot create program from opencl source.'

! build
call clBuildProgram(prog, '-cl-no-signed-zeros', ierr)

! get the compilation log
call clGetProgramBuildInfo(prog, device(mcnl%devid), CL_PROGRAM_BUILD_LOG, source, irec)
if(len(trim(source)) > 0) print*, trim(source)

if(ierr /= CL_SUCCESS) stop 'Fatal error: program build failed.'

! if we get here, then the kernal was successfully compiled
call WriteValue('EMMCOpenCL: ','Kernel Build Successful....')

! finally get the kernel and release the program
kernel = clCreateKernel(prog, 'MC', ierr)
call clReleaseProgram(prog, ierr)
! allocate device memory

open(unit = iunit, file = trim(randomseedfilename), form='unformatted', status='old')
read(iunit) nseeds
allocate(rnseeds(nseeds))
read(iunit) rnseeds
close(unit=iunit,status='keep')

if (globalworkgrpsz**2 .gt. nseeds) call FatalError('EMMCOpenCL:','insufficient prime numbers; reduce globalworkgrpsz parameter')


allocate(init_seeds(4*globalworkgrpsz*globalworkgrpsz),stat=istat)
init_seeds = 0
do i = 1,globalworkgrpsz
    do j = 1,globalworkgrpsz
        do k = 1,4
            init_seeds(4*((i-1)*globalworkgrpsz+j)+k) = rnseeds(4*((i-1)*globalworkgrpsz+j)+k)
        end do
    end do
end do

! create device memory buffers
LamX = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size_in_bytes, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory for LamX.'

LamY = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size_in_bytes, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory for LamY.'

depth = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size_in_bytes, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory for depth.'

energy = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size_in_bytes, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory for energy.'

seeds = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, ierr)
if(ierr /= CL_SUCCESS) stop 'Error: cannot allocate device memory for seeds.'

!call init_random_seed()
call clEnqueueWriteBuffer(command_queue, seeds, cl_bool(.true.), 0_8, size_in_bytes_seeds, init_seeds(1), ierr)

mainloop: do i = 1,(totnum_el/num_max+1)

! set the kernel arguments
    call clSetKernelArg(kernel, 0, LamX, ierr)
    if(ierr /= CL_SUCCESS) stop 'Error: cannot set LamX kernel argument.'

    call clSetKernelArg(kernel, 1, LamY, ierr)
    if(ierr /= CL_SUCCESS) stop 'Error: cannot set LamY kernel argument.'

    call clSetKernelArg(kernel, 2, EkeV, ierr)
    if(ierr /= CL_SUCCESS) stop 'Error: cannot set EkeV kernel argument.'

    call clSetKernelArg(kernel, 3, globalworkgrpsz, ierr)
    if(ierr /= CL_SUCCESS) stop 'Error: cannot set globalworkgrpsz kernel argument.'

    call clSetKernelArg(kernel, 4, Ze, ierr)
    if(ierr /= CL_SUCCESS) stop 'Error: cannot set Ze kernel argument.'

    call clSetKernelArg(kernel, 5, density, ierr)
    if(ierr /= CL_SUCCESS) stop 'Error: cannot set density kernel argument.'

    call clSetKernelArg(kernel, 6, at_wt, ierr)
    if(ierr /= CL_SUCCESS) stop 'Error: cannot set at_wt kernel argument.'

    call clSetKernelArg(kernel, 7, num_el, ierr)
    if(ierr /= CL_SUCCESS) stop 'Error: cannot set num_el kernel argument.'

    call clSetKernelArg(kernel, 8, seeds, ierr)
    if(ierr /= CL_SUCCESS) stop 'Error: cannot set seeds kernel argument.'

    call clSetKernelArg(kernel, 9, sig, ierr)
    if(ierr /= CL_SUCCESS) stop 'Error: cannot set sig kernel argument.'

    call clSetKernelArg(kernel, 10, omega, ierr)
    if(ierr /= CL_SUCCESS) stop 'Error: cannot set omega kernel argument.'

    call clSetKernelArg(kernel, 11, depth, ierr)
    if(ierr /= CL_SUCCESS) stop 'Error: cannot set depth kernel argument.'

    call clSetKernelArg(kernel, 12, energy, ierr)
    if(ierr /= CL_SUCCESS) stop 'Error: cannot set energy kernel argument.'

    call clSetKernelArg(kernel, 13, steps, ierr)
    if(ierr /= CL_SUCCESS) stop 'Error: cannot set steps kernel argument.'

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
! and get the nearest pixel [ take into account reversal of coordinate frame (x,y) -> (y,-x) ]
            if ((nint(delta*Lamresy(j)) .eq. 0.0) .and. (nint(-delta*Lamresx(j)) .eq. 0.0)) then
                val1 = val1 + 1
            end if

            val = val + 1
            idxy = (/ nint(delta*Lamresy(j)), nint(-delta*Lamresx(j)) /)

            if (maxval(abs(idxy)).le.nx) then
! If Ec larger than Emin, then we should count this electron
                if (energyres(j).gt.mcnl%Ehistmin) then
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

    if (mod(i,50).eq.0) then
        io_int(1) = i*num_max
        call WriteValue(' Total number of electrons generated = ',io_int, 1, "(I15)")
        io_int(1) = sum(accum_e)
        call WriteValue(' Number of electrons on detector       = ',io_int, 1, "(I15)")
    end if


end do mainloop

i = totnum_el/num_max

totnum_el = (i+1)*num_max

! output in .h5 format.

! Initialize FORTRAN interface.
!
call h5open_f(hdferr)
call timestamp(timestring=tstre)

! first of all, if the file exists, then delete it and rewrite it on each energyloop
dataname = trim(EMdatapathname)//trim(mcnl%dataname)
inquire(file=trim(dataname), exist=f_exists)

if (f_exists) then
  open(unit=dataunit, file=trim(dataname), status='old',form='unformatted')
  close(unit=dataunit, status='delete')
end if

! Create a new file using the default properties.
hdferr =  HDF_createFile(dataname, HDF_head)

! write the EMheader to the file
call HDF_writeEMheader(HDF_head, dstr, tstrb, tstre, progname)

! update the duration string
  groupname = 'EMheader'
  hdferr = HDF_openGroup(groupname, HDF_head)

  call CPU_TIME(tstop)
  dataset = 'Duration'
  tstop = tstop - tstart
  hdferr = HDF_writeDatasetFloat(dataset, tstop, HDF_head)

  call HDF_pop(HDF_head)

! create a namelist group to write all the namelist files into
groupname = "NMLfiles"
hdferr = HDF_createGroup(groupname, HDF_head)

! read the text file and write the array to the file
dataset = 'MCOpenCLNML'
hdferr = HDF_writeDatasetTextFile(dataset, nmldeffile, HDF_head)

! leave this group
call HDF_pop(HDF_head)

! create a namelist group to write all the namelist files into
groupname = "NMLparameters"
hdferr = HDF_createGroup(groupname, HDF_head)
call HDFwriteMCCLNameList(HDF_head, mcnl)

! leave this group
call HDF_pop(HDF_head)

! then the remainder of the data in a EMData group
groupname = 'EMData'
hdferr = HDF_createGroup(groupname, HDF_head)

dataset = 'numEbins'
hdferr = HDF_writeDatasetInteger(dataset, numEbins, HDF_head)

dataset = 'numzbins'
hdferr = HDF_writeDatasetInteger(dataset, numzbins, HDF_head)

dataset = 'totnum_el'
hdferr = HDF_writeDatasetInteger(dataset, totnum_el, HDF_head)

dataset = 'accum_e'
hdferr = HDF_writeDatasetIntegerArray3D(dataset, accum_e, numEbins, 2*nx+1, 2*nx+1, HDF_head)

dataset = 'accum_z'
hdferr = HDF_writeDatasetIntegerArray4D(dataset, accum_z, numEbins, numzbins, 2*(nx/10)+1, 2*(nx/10)+1, HDF_head)

call HDF_pop(HDF_head,.TRUE.)

! and close the fortran hdf interface
call h5close_f(hdferr)

! and write some infgormation to the console
bse = real(val)/real(totnum_el)
io_real(1) = bse
call WriteValue("Backscatter yield          = ", io_real, 1)
io_int(1) = totnum_el
call WriteValue("Total number of incident electrons = ", io_int, 1)
io_int(1) = sum(accum_e)
call WriteValue("Number of BSE electrons = ", io_int, 1)
io_int(1) = maxval(accum_z)
io_int(2) = maxval(accum_e)
call WriteValue("Maximum electron in depth and energy bin resp. = ", io_int, 2)
!
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
