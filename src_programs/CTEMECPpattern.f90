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
! CTEMsoft2013:CTEMECP.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMECP
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief Electron channeling patterns from master pattern
!
!> @date 08/26/14 SS 1.0 f90
!--------------------------------------------------------------------------

program CTEMECP

use local
use typedefs
use NameListTypedefs
use NameListHandlers
use files
use io

IMPLICIT NONE

character(fnlen)                        :: nmldeffile, progname, progdesc
type(ECPpatternNameListType)                   :: ecpnl

nmldeffile = 'CTEMECP.nml'
progname = 'CTEMECP.f90'
progdesc = 'Electron channeling patterns from master pattern'

! deal with the command line arguments, if any
call Interpret_Program_Arguments(nmldeffile,2,(/ 0, 40 /), progname)

! deal with the namelist stuff
call GetECPpatternNameList(nmldeffile,ecpnl)

! print some information
call CTEMsoft(progname, progdesc)

! perform the zone axis computations
call ECpattern(ecpnl, progname)

end program CTEMECP

!--------------------------------------------------------------------------
!
! SUBROUTINE: ECpattern
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief Electron channeling patterns from master pattern
!
!> @date 08/27/14 SS 1.0 f90
!--------------------------------------------------------------------------
subroutine ECpattern(ecpnl, progname)

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

IMPLICIT NONE

type(ECPpatternNameListType),INTENT(IN)        :: ecpnl
character(fnlen),INTENT(IN)             :: progname

character(fnlen)                        :: oldprogname
character(fnlen)                        :: xtalname
character(5)                            :: oldscversion
character(fnlen)                        :: energyfile

integer(kind=irg)                       :: npx,npy,numset,istat
integer(kind=irg),allocatable           :: ATOM_type(:)
real(kind=dbl)                          :: EkeV
real(kind=sgl)                          :: dmin
real(kind=sgl),allocatable              :: sr(:,:)

type(unitcell), pointer                 :: cell
type(gnode)                             :: rlp
type(DynType)                           :: Dyn
logical                                 :: verbose

type(kvectorlist),pointer               :: khead,ktmp
integer(kind=irg)                       :: numk,nix,niy,i,j,ierr,ipx,ipy
real(kind=dbl)                          :: scl,x,dx,dy,dxm,dym
real(kind=dbl)                          :: dc(3),ixy(2)
integer(kind=irg),allocatable           :: kij(:,:)
real(kind=sgl),allocatable              :: ecp(:,:)


!=================================================================
! read Master pattern output file and extract necessary parameters
! first, we need to load the data from the ECP master program
!=================================================================

call Message('opening '//trim(ecpnl%masterfile), frm = "(A)" )

open(dataunit,file=trim(ecpnl%masterfile),status='unknown',form='unformatted')

! lines from CTEMECPmaster.f90... these are the things we need to read in...
! write (dataunit) progname
!! write the version number
! write (dataunit) scversion
!! then the name of the crystal data file
! write (dataunit) xtalname
!! write the name of corresponding Monte Carlo data file
! write (dataunit) ecpnl%energyfile
!! energy information etc...
! write (dataunit) ecpnl%npx,ecpnl%npx,numset
! write (dataunit) EkeV
! write (dataunit) ecpnl%dmin
!! atom type array for asymmetric unit
! write (dataunit) cell%ATOM_type(1:numset)
!! finally the masterpattern array
! write (dataunit) sr

read (dataunit) oldprogname
read (dataunit) oldscversion
read (dataunit) xtalname
read (dataunit) energyfile

read (dataunit) npx,npy,numset
read (dataunit) EkeV
read (dataunit) dmin

allocate(ATOM_type(1:numset))
read (dataunit) ATOM_type

allocate(sr(-npx:npx,-npy:npy),stat=istat)
read (dataunit) sr

close(dataunit,status='keep')
call Message(' -> completed reading '//trim(ecpnl%masterfile), frm = "(A)")

!=============================================
! completed reading master pattern file
! proceed to crystallography section
!=============================================

nullify(cell)
allocate(cell)

! load the crystal structure and compute the Fourier coefficient lookup table
verbose = .TRUE.
call Initialize_Cell(cell,Dyn,rlp,xtalname,dmin,sngl(1000.0*EkeV),verbose)

nullify(khead)
nullify(ktmp)
allocate(khead)
numk = 0

call CalckvectorsECP(khead,cell,ecpnl%k,ecpnl%thetac,ecpnl%npix,ecpnl%npix,numk)
allocate(kij(2,numk),stat=istat)

ktmp => khead
kij(1:2,1) = (/ ktmp%i, ktmp%j /)

do i = 2,numk
    ktmp => ktmp%next
    kij(1:2,i) = (/ ktmp%i, ktmp%j /)
end do


! determine the scale factor for the Lambert interpolation; the square has
! an edge length of 2 x sqrt(pi/2)
scl = float(npx)/1.25331413732D0
!print*,"success..."


allocate(ecp(-ecpnl%npix:ecpnl%npix,-ecpnl%npix:ecpnl%npix),stat=istat)
ecp = 0.0

ktmp => khead

imageloop: do i = 1,numk
    dc = ktmp%k(1:3)
    dc = dc/sqrt(sum(dc**2))
! make sure the third one is positive; if not, switch all
    if (dc(3) .lt. 0.0) dc = -dc
! convert these direction cosines to coordinates in the Rosca-Lambert projection
    ixy = scl * LambertSphereToSquare( dc, istat )
    x = ixy(1)
    ixy(1) = ixy(2)
    ixy(2) = -x
! four-point interpolation (bi-quadratic)
    nix = int(npx+ixy(1))-npx
    niy = int(npx+ixy(2))-npx
    dx = ixy(1)-nix
    dy = ixy(2)-niy
    dxm = 1.0-dx
    dym = 1.0-dy
! interpolate the intensity
    ipx = kij(1,i)
    ipy = kij(2,i)
!print*,nix,niy
!print*,"success..."

    ecp(ipx,ipy) =  sr(nix,niy) * dxm * dym + &
                    sr(nix+1,niy) * dx * dym + &
                    sr(nix,niy+1) * dxm * dy + &
                    sr(nix+1,niy+1) * dx * dy

    ktmp => ktmp%next
end do imageloop

nullify(ktmp)
call Delete_kvectorlist(khead)
nullify(khead)

open(unit=12,file="test.txt",action="write")

do i= -ecpnl%npix,ecpnl%npix
    do j= -ecpnl%npix,ecpnl%npix
        write(12, '(F15.6)', advance='no') ecp(i,j)
    end do
    write(12, *) ''  ! this gives you the line break
end do

close(unit=12,status='keep')

end subroutine ECpattern