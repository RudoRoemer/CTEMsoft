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
! EMsoft:EMECP.f90
!--------------------------------------------------------------------------
!
! PROGRAM: EMECP
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief Electron channeling patterns from master pattern
!
!> @date 08/26/14 SS 1.0 f90
!> @date 13/10/15 SS 2.0 added detector model+new GetVectorCone routine+OpenMP
!--------------------------------------------------------------------------

program EMECP

use local
use typedefs
use NameListTypedefs
use NameListHandlers
use files
use io

IMPLICIT NONE

character(fnlen)                        :: nmldeffile, progname, progdesc
type(ECPpatternNameListType)            :: ecpnl

nmldeffile = 'EMECP.nml'
progname = 'EMECP.f90'
progdesc = 'Electron channeling patterns from master pattern'

! print some information
call EMsoft(progname, progdesc)

! deal with the command line arguments, if any
call Interpret_Program_Arguments(nmldeffile,2,(/ 0, 40 /), progname)

! deal with the namelist stuff
call GetECPpatternNameList(nmldeffile,ecpnl)

! perform the zone axis computations
call ECpattern(ecpnl, progname)

end program EMECP

!-----------------------------------------------------------------------------------
!
! SUBROUTINE: ECpattern
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief Electron channeling patterns from master pattern
!
!> @date 08/27/14 SS 1.0 f90
!> @date 13/10/15 SS 2.0 added detector model+new GetVectorCone routine+OpenMP+hdf5
!-------------------------------------------------------------------------------------
subroutine ECpattern(ecpnl, progname)

use local
use typedefs
use NameListTypedefs
use crystal
use constants
use symmetry
use Lambert
use initializers
use constants
use error
use io
use files
use ECPmod
use rotations

IMPLICIT NONE

type(ECPpatternNameListType),INTENT(IN) :: ecpnl
character(fnlen),INTENT(IN)             :: progname

type(ECPLargeAccumType),pointer         :: acc
type(ECPMasterType),pointer             :: master
type(IncidentListECP),pointer           :: khead, ktmp
real(kind=dbl),allocatable              :: klist(:,:)

integer(kind=irg)                       :: npx,npy,numset,istat,val
integer(kind=irg),allocatable           :: ATOM_type(:)
real(kind=dbl)                          :: EkeV
real(kind=sgl)                          :: dmin, FN(3)
real(kind=sgl),allocatable              :: sr(:,:)

type(unitcell), pointer                 :: cell
type(gnode)                             :: rlp
type(DynType)                           :: Dyn
logical                                 :: verbose

integer(kind=irg)                       :: numk,nix,niy,i,j,ierr,ipx,ipy
real(kind=dbl)                          :: scl,x,dx,dy,dxm,dym
real(kind=dbl)                          :: dc(3),ixy(2)
real(kind=sgl)                          :: rotmat(3,3)
integer(kind=irg),allocatable           :: kij(:,:)
real(kind=sgl),allocatable              :: ecp(:,:)
real(kind=dbl)                          :: time_start,time_end


!=================================================================
! read Monte Carlo output file and extract necessary parameters
! first, we need to load the data from the output of EMMCOpenCL
! used in the bse1 mode
!=================================================================

call Message('opening '//trim(ecpnl%energyfile), frm = "(A)" )
call  ECPreadMCfile(ecpnl, acc, verbose=.TRUE.)
call Message(' -> completed reading '//trim(ecpnl%energyfile), frm = "(A)")

!=================================================================
! read Master pattern output file and extract necessary parameters
! first, we need to load the data from the ECP master program
!=================================================================

call Message('opening '//trim(ecpnl%masterfile), frm = "(A)" )
call ECPreadMasterfile(ecpnl, master, verbose=.TRUE.)
call Message(' -> completed reading '//trim(ecpnl%masterfile), frm = "(A)")

!=================================================================
! completed reading the file; generating list of incident vectors
!=================================================================

numk = 0
rotmat = eu2om(ecpnl%eu(1:3))

call GetVectorsCone(ecpnl, khead, rotmat, numk)
allocate(kij(2,numk),klist(3,numk),stat=istat)

ktmp => khead
! converting to array for OpenMP parallelization
do i = 1,numk
   klist(1:3,i) = ktmp%k(1:3)
   kij(1:2,i) = (/ktmp%i,ktmp%j/)
   ktmp => ktmp%next
end do

! determine the scale factor for the Lambert interpolation; the square has
! an edge length of 2 x sqrt(pi/2)
scl = float(npx)/LPs%sPio2

allocate(ecp(-ecpnl%npix:ecpnl%npix,-ecpnl%npix:ecpnl%npix),stat=istat)
ecp = 0.0

imageloop: do while(associated(ktmp%next))
    dc = ktmp%k(1:3)
! make sure the third one is positive; if not, switch all
    if (dc(3) .lt. 0.0) dc = -dc
! convert these direction cosines to coordinates in the Rosca-Lambert projection
    ixy = scl * LambertSphereToSquare( dc, istat )
    x = ixy(1)
    ixy(1) = ixy(2)
    ixy(2) = -x
! four-point interpolation (bi-quadratic)
    nix = nint(ixy(1))
    niy = nint(ixy(2))
!print*,nix,niy
    dx = ixy(1)-nix
    dy = ixy(2)-niy
    dxm = 1.0-dx
    dym = 1.0-dy
! interpolate the intensity
    ipx = kij(1,i)
    ipy = kij(2,i)
    ecp(ipx,ipy) =  sr(nix,niy) * dxm * dym + &
                    sr(nix+1,niy) * dx * dym + &
                    sr(nix,niy+1) * dxm * dy + &
                    sr(nix+1,niy+1) * dx * dy
    ktmp => ktmp%next
    i = i + 1
end do imageloop
nullify(ktmp)
call Delete_kvectorlist(khead)
!print*,val
open(unit=12,file="test.txt",action="write")

do i= -ecpnl%npix,ecpnl%npix
    do j= -ecpnl%npix,ecpnl%npix
        write(12, '(F15.6)', advance='no') ecp(i,j)
    end do
    write(12, *) ''  ! this gives you the line break
end do

close(unit=12,status='keep')

end subroutine ECpattern
