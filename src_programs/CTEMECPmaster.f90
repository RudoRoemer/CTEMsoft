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
! PROGRAM: CTEMECPmaster
!
!> @author Marc De Graef/ Saransh Singh, Carnegie Mellon University
!
!> @brief Zone axis electron channeling master patterns
!
!> @date 03/18/10 MDG 1.0 f90
!> @date 08/09/10 MDG 2.0 corrected weight factors and g-vector ordering problem
!> @date 11/18/13 MDG 3.0 major rewrite with new libraries
!> @date 06/27/14 MDG 4.0 removal of all globals; separation of namelist handling from computation
!> @date 06/30/14 MDG 4.1 added OpenMP
!> @DATE 07/23/14 SS  4.2 rewrite master pattern
!--------------------------------------------------------------------------

program CTEMECPmaster

use local
use typedefs
use NameListTypedefs
use NameListHandlers
use files
use io

IMPLICIT NONE

character(fnlen)                        :: nmldeffile, progname, progdesc
type(ECPMasterNameListType)                   :: ecpnl

nmldeffile = 'CTEMECPmaster.nml'
progname = 'CTEMECPmaster.f90'
progdesc = 'Master pattern generation for Electron channeling pattern'

! deal with the command line arguments, if any
call Interpret_Program_Arguments(nmldeffile,2,(/ 0, 40 /), progname)

! deal with the namelist stuff
call GetECPMasterNameList(nmldeffile,ecpnl)

! print some information
call CTEMsoft(progname, progdesc)

! perform the zone axis computations
call ECmasterpattern(ecpnl, progname)

end program CTEMECPmaster

!--------------------------------------------------------------------------
!
! SUBROUTINE:ECmasterpattern
!
!> @author Marc De Graef/Saransh Singh, Carnegie Mellon University
!
!> @brief compute a master electron channeling pattern
!
!> @note This is really very similar EBSDmaster computation, except that
!> the final intensity computation is somewhat different.  We could in
!> principle also include the Kossel pattern computation in this program.
!> This program now also includes the Bethe potential approximation, to
!> hopefully speed things up a little bit...
!
!> @param ecpnl name list structure
!> @param progname program name
!
!> @date 11/18/13  MDG 1.0 major rewrite from older ECP program
!> @date 11/22/13  MDG 1.1 output modified for IDL interface
!> @date 03/04/14  MDG 1.2 added scattering matrix mode
!> @date 06/27/14  MDG 2.0 removal of globals, split of namelist and computation; OpenMP
!> @date 06/30/14  MDG 2.1 debug; found some inconsistent array (de)allocations
!> @date 07/23/14  SS  2.2 conversion to master pattern simulation
!--------------------------------------------------------------------------
subroutine ECmasterpattern(ecpnl, progname)

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
use omp_lib
use MBModule

IMPLICIT NONE

type(ECPMasterNameListType),INTENT(IN)        :: ecpnl
character(fnlen),INTENT(IN)                   :: progname

real(kind=dbl)          :: frac
integer(kind=irg)       :: gzero, istat

integer(kind=irg)       :: numEbins, numzbins, nx, ny, totnum_el ! reading from MC file
real(kind=dbl)          :: EkeV, Ehistmin, Ebinsize, depthmax, depthstep, sig, omega  ! reading from MC file
integer(kind=irg), allocatable :: accum_e(:,:,:), accum_z(:,:,:,:) ! reading from MC file

integer(kind=irg)       :: io_int_sgl(1), io_int(6) ! integer output variable
real(kind=dbl)          :: io_real(5) ! real output variable

integer(kind=irg)       :: i, j, isym, pgnum ! variables for point group and Laue group
integer(kind=irg),parameter     :: LaueTest(11) = (/ 149, 151, 153, 156, 158, 160, 161, 164, 165, 166, 167 /)  ! space groups with 2 or mirror at 30 degrees
integer(kind=irg)       :: npyhex, ijmax, numk, skip ! parameters for calckvectors and calcwavelength subroutine

integer(kind=irg)       :: ga(3), gb(3) ! shortest reciprocal lattice vector for zone axis
real(kind=sgl), allocatable :: thick(:), sr(:,:), lambdaZ(:)
real(kind=dbl)          :: intthick
complex(kind=dbl),allocatable   :: Lgh(:,:),Sgh(:,:),Sghtmp(:,:,:)
complex(kind=dbl),allocatable   :: DynMat(:,:)
complex(kind=dbl)       :: czero

integer(kind=irg)       :: nt, nns, nnw, tots, totw ! thickness array and BetheParameters strong and weak beams
real(kind=sgl)          :: FN(3), kk(3), fnat, kn
integer(kind=irg)       :: numset, nref, ipx, ipy, iequiv(2,12), nequiv, ip, jp, izz, IE, iz, one,ierr
integer(kind=irg),allocatable   :: kij(:,:), nat(:)
real(kind=dbl)          :: res(2)

character(fnlen)        :: oldprogname
character(fnlen)        :: xtalname
character(8)            :: MCscversion
character(4)            :: MCmode
character(6)            :: projtype

logical                 :: verbose, usehex, switchmirror

type(unitcell), pointer         :: cell
type(gnode)                     :: rlp
type(DynType)                   :: Dyn
type(kvectorlist), pointer      :: khead, ktmp ! linked list for incident wave vectors for master list
type(kvectorlist), pointer      :: kheadcone,ktmpcone ! linked list for incident wave vectors for individual pattern
real(kind=dbl),allocatable      :: ecpattern(:,:)
type(BetheParameterType)        :: BetheParameters
type(reflisttype),pointer       :: reflist, firstw,rltmp


gzero = 1
frac = 0.05
!dataunit = 10

allocate(cell)

!=============================================================
!read Monte Carlo output file and extract necessary parameters
! first, we need to load the data from the MC program.
!=============================================================

call Message('opening '//trim(ecpnl%energyfile), frm = "(A)" )

open(dataunit,file=trim(ecpnl%energyfile),status='unknown',form='unformatted')

! lines from CTEMMCCL.f90... these are the things we need to read in...
! write (dataunit) progname
!! write the version number
! write (dataunit) MCscversion
!! then the name of the crystal data file
! write (dataunit) xtalname
!! energy information etc...
! write (dataunit) numEbins, numzbins, numsx, numsy, , totnum_el
! write (dataunit) EkeV, Ehistmin, Ebinsize, depthmax, depthstep
! write (dataunit) sig, omega
! write (dataunit) MCmode
!! and here are the actual results
! write (dataunit) accum_e
! write (dataunit) accum_z

read (dataunit) oldprogname
read (dataunit) MCscversion
read (dataunit) xtalname

read(dataunit) numEbins, numzbins, nx, ny, totnum_el
nx = (nx - 1)/2
ny = (ny - 1)/2

read (dataunit) EkeV, Ehistmin, Ebinsize, depthmax, depthstep
io_real(1:5) = (/ EkeV, Ehistmin, Ebinsize, depthmax, depthstep /)
call WriteValue(' EkeV, Ehistmin, Ebinsize, depthmax, depthstep ',io_real,5,"(4F10.5,',',F10.5)")

read (dataunit) sig, omega
read (dataunit) MCmode

allocate(accum_e(numEbins,-nx:nx,-nx:nx),accum_z(numEbins,numzbins,-nx/10:nx/10,-nx/10:nx/10),stat=istat)

read(dataunit) accum_e
! actually, we do not yet need the accum_e array for ECP. This will be removed with an updated version of the MC code
! but we need to skip it in this unformatted file so that we can read the accum_z array ...
!deallocate(accum_e)

read(dataunit) accum_z    ! we only need this array for the depth integrations
close(dataunit,status='keep')
call Message(' -> completed reading '//trim(ecpnl%energyfile), frm = "(A)")
!=============================================
! completed reading monte carlo file
!=============================================

nullify(cell)
allocate(cell)

! load the crystal structure and compute the Fourier coefficient lookup table
verbose = .TRUE.
call Initialize_Cell(cell,Dyn,rlp, xtalname, ecpnl%dmin, sngl(1000.0*EkeV),verbose)

! determine the point group and Laue group number
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

if(usehex)  npyhex = nint(2.0*float(ecpnl%npx)/sqrt(3.0))
ijmax = float(ecpnl%npx)**2   ! truncation value for beam directions


!=============================================
! generating list of incident wave vectors
!=============================================

! determine all independent incident beam directions (use a linked list starting at khead)
! numk is the total number of k-vectors to be included in this computation;
nullify(khead)
if (usehex) then
call Calckvectors(khead,cell, (/ 0.D0, 0.D0, 1.D0 /), (/ 0.D0, 0.D0, 0.D0 /),0.D0,ecpnl%npx,npyhex,numk, &
isym,ijmax,'RoscaLambert',usehex)
else
! Calckvectors(k,ga,ktmax,npx,npy,numk,isym,ijmax,mapmode,usehex)
call Calckvectors(khead,cell, (/ 0.D0, 0.D0, 1.D0 /), (/ 0.D0, 0.D0, 0.D0 /),0.D0,ecpnl%npx,ecpnl%npx,numk, &
isym,ijmax,'RoscaLambert',usehex)
end if
io_int_sgl(1)=numk

call WriteValue('# independent beam directions to be considered = ', io_int_sgl, 1, "(I8)")

ktmp => khead
czero = cmplx(0.D0, 0.D0)
! force dynamical matrix routine to read new Bethe parameters from file
call Set_Bethe_Parameters(BetheParameters,.TRUE.)

nt = nint((depthmax - ecpnl%startthick)/depthstep)
allocate(thick(nt))
thick = ecpnl%startthick + depthstep * (/ (i-1,i=1,nt) /)

nullify(reflist)
nullify(firstw)

nns = 0
nnw = 0
tots = 0
totw = 0

numset = cell % ATOM_ntype  ! number of special positions in the unit cell
izz = numzbins
allocate(lambdaZ(1:izz),stat=istat)
allocate(nat(numset),stat=istat)
allocate(kij(2,numk),stat=istat)

do iz=1,izz
    lambdaZ(iz) = float(sum(accum_z(:,iz,:,:)))/float(sum(accum_z(:,:,:,:)))
end do

kij(1:2,1) = (/ ktmp%i, ktmp%j /)

do i = 2,numk
    ktmp => ktmp%next
    kij(1:2,i) = (/ ktmp%i, ktmp%j /)
end do

open(unit=11,file="test.txt",action="write")

ktmp => khead

nat = 0
fnat = 1.0/float(sum(cell%numat(1:numset)))
intthick = dble(depthmax)


open(unit=dataunit,file=trim(ecpnl%outname),status='unknown',action='write',form = 'unformatted')
! write the program identifier
write (dataunit) progname
! write the version number
write (dataunit) scversion
! then the name of the crystal data file
write (dataunit) xtalname
! then the name of the corresponding Monte Carlo data file
write (dataunit) ecpnl%energyfile
! energy information and array size
write (dataunit) ecpnl%npx,ecpnl%npx,numset
write (dataunit) EkeV
write (dataunit) ecpnl%dmin
! atom type array for asymmetric unit
write (dataunit) cell%ATOM_type(1:numset)
! is this a regular (square) or hexagonal projection ?
!if (usehex) then
 !   projtype = 'hexago'
 !  write (dataunit) projtype
!else
 !   projtype = 'square'
 !   write (dataunit) projtype
!end if

allocate(sr(2*ecpnl%npx+1,2*ecpnl%npx+1),stat=istat)
sr = 0.0

beamloop: do i = 1, numk
    kk = ktmp%k(1:3)
    FN = ecpnl%fn

    call Initialize_ReflectionList(cell, reflist, BetheParameters, FN, kk, ecpnl%dmin, nref)

! determine strong and weak reflections
    call Apply_BethePotentials(cell, reflist, firstw, BetheParameters, nref, nns, nnw)

    allocate(DynMat(nns,nns))

    call GetDynMat(cell, reflist, firstw, rlp, DynMat, nns, nnw)

! then we need to initialize the Sgh and Lgh arrays
    if (allocated(Sgh)) deallocate(Sgh)
    if (allocated(Lgh)) deallocate(Lgh)
    if (allocated(Sghtmp)) deallocate(Sghtmp)

    allocate(Sghtmp(nns,nns,numset),Lgh(nns,nns),Sgh(nns,nns))

    Sgh = czero
    Lgh = czero
    Sghtmp = czero
    nat = 0

    call CalcSgh(cell,reflist,nns,numset,Sghtmp,nat)
! sum Sghtmp over the sites
    Sgh = sum(Sghtmp,3)


! solve the dynamical eigenvalue equation
    kn = ktmp%kn

    call CalcLgh(DynMat,Lgh,intthick,dble(kn),nns,gzero,depthstep,lambdaZ,izz)
    deallocate(DynMat,Sghtmp)

! and store the resulting values

!print*, iequiv(1,:)

!print*, nequiv, isym
!print*, iequiv(1,:)
    ipx = kij(1,i)
    ipy = kij(2,i)
    if (isym.ne.1) then
        call Apply2DLaueSymmetry(ipx,ipy,isym,iequiv,nequiv)
        iequiv(1,1:nequiv) = iequiv(1,1:nequiv) + ecpnl%npx + 1
        iequiv(2,1:nequiv) = iequiv(2,1:nequiv) + ecpnl%npx + 1
        do ip=1,nequiv
            sr(iequiv(1,ip),iequiv(2,ip)) = real(sum(Lgh(1:nns,1:nns)*Sgh(1:nns,1:nns)))/float(sum(nat))
        end do
    else
        sr(ipx+ecpnl%npx+1,ipy+ecpnl%npx+1) = real(sum(Lgh(1:nns,1:nns)*Sgh(1:nns,1:nns)))/float(sum(nat))
    end if
!print*, sr
    totw = totw + nnw
    tots = tots + nns
    deallocate(Lgh, Sgh)

    if (mod(i,2500).eq.0) then
        io_int(1) = i
        call WriteValue('  completed beam direction ',io_int, 1, "(I8)")
! write(*,*) minval(sr),maxval(sr)
    end if

    call Delete_gvectorlist(reflist)
    !write(dataunit) ktmp%k
! and finally the results array
    !write (dataunit) sr

    ktmp => ktmp%next

end do beamloop


do i=1,2*ecpnl%npx+1
   do j=1,2*ecpnl%npx+1
        write(11, '(F9.6)', advance='no') sr(i,j)
    end do
    write(11, *) ''  ! this gives you the line break
end do

write (dataunit) sr
!deallocate(sr)

close(unit=dataunit,status='keep')
close(unit=11,status='keep')

end subroutine ECmasterpattern