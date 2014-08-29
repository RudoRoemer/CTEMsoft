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

!--------------------------------------------------------------------------
! CTEMsoft2013:NameListHandlers.f90
!--------------------------------------------------------------------------
!
! PROGRAM: NameListHandlers
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief routines for reading and returning name list type structures
!
!> @date 06/13/14 MDG 1.0 original
!--------------------------------------------------------------------------
module NameListHandlers

use local
use NameListTypedefs

contains

!--------------------------------------------------------------------------
!
! SUBROUTINE:GetKosselNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read namelist file and fill knl structure (used by CTEMKossel.f90)
!
!> @param nmlfile namelist file name
!> @param knl Kossel name list structure
!
!> @date 06/13/14  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine GetKosselNameList(nmlfile, knl)

use error

IMPLICIT NONE

character(fnlen),INTENT(IN)             :: nmlfile
type(KosselNameListType),INTENT(INOUT)  :: knl

integer(kind=irg)       :: stdout
integer(kind=irg)       :: numthick
integer(kind=irg)       :: npix
integer(kind=irg)       :: maxHOLZ
integer(kind=irg)       :: nthreads
integer(kind=irg)       :: k(3)
integer(kind=irg)       :: fn(3)
real(kind=sgl)          :: voltage
real(kind=sgl)          :: dmin
real(kind=sgl)          :: convergence
real(kind=sgl)          :: startthick
real(kind=sgl)          :: thickinc
real(kind=sgl)          :: minten
character(fnlen)        :: xtalname
character(fnlen)        :: outname

namelist /Kossellist/ stdout, xtalname, voltage, k, fn, dmin, convergence, minten, nthreads, &
                              startthick, thickinc, numthick, outname, npix, maxHOLZ

! set the input parameters to default values (except for xtalname, which must be present)
stdout = 6                      ! standard output
numthick = 10                   ! number of increments
npix = 256                      ! output arrays will have size npix x npix
maxHOLZ = 3                     ! output arrays will have size npix x npix
nthreads = 4                    ! default number of threads for OpenMP
k = (/ 0, 0, 1 /)               ! beam direction [direction indices]
fn = (/ 0, 0, 1 /)              ! foil normal [direction indices]
voltage = 200000.0              ! acceleration voltage [V]
dmin = 0.025                    ! smallest d-spacing to include in dynamical matrix [nm]
convergence = 25.0              ! beam convergence angle [mrad]
startthick = 10.0               ! starting thickness [nm]
thickinc = 10.0                 ! thickness increment
minten = 1.0E-5                 ! minimum intensity in diffraction disk to make it into the output file
xtalname = 'undefined'          ! initial value to check that the keyword is present in the nml file
outname = 'Kosselout.data'      ! output filename

! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=Kossellist)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(xtalname).eq.'undefined') then
  call FatalError('CTEMKossel:',' structure file name is undefined in '//nmlfile)
 end if

! if we get here, then all appears to be ok, and we need to fill in the knl fields
knl%stdout = stdout
knl%numthick = numthick
knl%npix = npix
knl%maxHOLZ = maxHOLZ
knl%nthreads = nthreads
knl%k = k
knl%fn = fn
knl%voltage = voltage
knl%dmin = dmin
knl%convergence = convergence
knl%startthick = startthick
knl%thickinc = thickinc
knl%minten = minten
knl%xtalname = xtalname
knl%outname = outname

end subroutine GetKosselNameList


!--------------------------------------------------------------------------
!
! SUBROUTINE:GetMCNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read namelist file and fill mcnl structure (used by CTEMMC.f90)
!
!> @param nmlfile namelist file name
!> @param mcnl Monte Carloname list structure
!
!> @date 06/18/14  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine GetMCNameList(nmlfile, mcnl)

use error

IMPLICIT NONE

character(fnlen),INTENT(IN)             :: nmlfile
type(MCNameListType),INTENT(INOUT)      :: mcnl

integer(kind=irg)       :: stdout
integer(kind=irg)       :: numsx
integer(kind=irg)       :: primeseed
integer(kind=irg)       :: num_el
integer(kind=irg)       :: nthreads
real(kind=dbl)          :: sig
real(kind=dbl)          :: omega
real(kind=dbl)          :: EkeV
real(kind=dbl)          :: Ehistmin
real(kind=dbl)          :: Ebinsize
real(kind=dbl)          :: depthmax
real(kind=dbl)          :: depthstep
character(4)            :: MCmode
character(fnlen)        :: xtalname
character(fnlen)        :: dataname

! define the IO namelist to facilitate passing variables to the program.
namelist  / MCdata / stdout, xtalname, sig, numsx, num_el, primeseed, EkeV, &
                dataname, nthreads, Ehistmin, Ebinsize, depthmax, depthstep, omega, MCmode

! set the input parameters to default values (except for xtalname, which must be present)
stdout = 6
numsx = 1501
primeseed = 932117
num_el = 12500000
nthreads = 1
sig = 70.D0
omega = 0.D0
EkeV = 30.D0
Ehistmin = 5.D0
Ebinsize = 0.5D0
depthmax = 100.D0
depthstep = 1.0D0
MCmode = 'CSDA'
xtalname = 'undefined'
dataname = 'MCoutput.data'

! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=MCdata)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(xtalname).eq.'undefined') then
  call FatalError('CTEMMC:',' structure file name is undefined in '//nmlfile)
 end if

! if we get here, then all appears to be ok, and we need to fill in the mcnl fields
mcnl%stdout = stdout
mcnl%numsx = numsx
mcnl%primeseed = primeseed
mcnl%num_el = num_el
mcnl%nthreads = nthreads
mcnl%sig = sig
mcnl%omega = omega
mcnl%EkeV = EkeV
mcnl%Ehistmin = Ehistmin
mcnl%Ebinsize = Ebinsize
mcnl%depthmax = depthmax
mcnl%depthstep = depthstep
mcnl%MCmode = MCmode
mcnl%xtalname = xtalname
mcnl%dataname = dataname
mcnl%stdout = stdout

end subroutine GetMCNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:GetMCCLNameList
!
!> @author Saransh Singh/Marc De Graef, Carnegie Mellon University
!
!> @brief read namelist file and fill mcnl structure (used by CTEMMCCL.f90)
!
!> @param nmlfile namelist file name
!> @param mcnl Monte Carloname list structure
!
!> @date 06/18/14  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine GetMCCLNameList(nmlfile, mcnl)

use error

IMPLICIT NONE

character(fnlen),INTENT(IN)             :: nmlfile
type(MCCLNameListType),INTENT(INOUT)      :: mcnl

integer(kind=irg)       :: stdout
integer(kind=irg)       :: numsx
integer(kind=irg)       :: globalworkgrpsz
integer(kind=irg)       :: num_el
integer(kind=irg)       :: totnum_el
real(kind=dbl)          :: sig
real(kind=dbl)          :: omega
real(kind=dbl)          :: EkeV
real(kind=dbl)          :: Ehistmin
real(kind=dbl)          :: Ebinsize
real(kind=dbl)          :: depthmax
real(kind=dbl)          :: depthstep
character(4)            :: MCmode
character(fnlen)        :: xtalname
character(fnlen)        :: dataname
character(fnlen)        :: primelist
character(fnlen)        :: mode

! define the IO namelist to facilitate passing variables to the program.
namelist  / MCCLdata / stdout, xtalname, sig, numsx, num_el, globalworkgrpsz, EkeV, &
dataname, primelist, totnum_el, Ehistmin, Ebinsize, depthmax, depthstep, omega, MCmode, mode

! set the input parameters to default values (except for xtalname, which must be present)
stdout = 6
numsx = 1501
globalworkgrpsz = 100
num_el = 10
totnum_el = 100000
sig = 70.D0
omega = 0.D0
EkeV = 30.D0
Ehistmin = 5.D0
Ebinsize = 0.5D0
depthmax = 100.D0
depthstep = 1.0D0
MCmode = 'CSDA'
xtalname = 'undefined'
dataname = 'MCoutput.data'
primelist = 'list.txt'
mode = 'full'

! read the namelist file
open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
read(UNIT=dataunit,NML=MCCLdata)
close(UNIT=dataunit,STATUS='keep')

! check for required entries
if (trim(xtalname).eq.'undefined') then
call FatalError('CTEMMC:',' structure file name is undefined in '//nmlfile)
end if

! if we get here, then all appears to be ok, and we need to fill in the mcnl fields
mcnl%stdout = stdout
mcnl%numsx = numsx
mcnl%globalworkgrpsz = globalworkgrpsz
mcnl%num_el = num_el
mcnl%totnum_el = totnum_el
mcnl%sig = sig
mcnl%omega = omega
mcnl%EkeV = EkeV
mcnl%Ehistmin = Ehistmin
mcnl%Ebinsize = Ebinsize
mcnl%depthmax = depthmax
mcnl%depthstep = depthstep
mcnl%MCmode = MCmode
mcnl%xtalname = xtalname
mcnl%dataname = dataname
mcnl%primelist = primelist
mcnl%mode = mode

end subroutine GetMCCLNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:GetEBSDMasterNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read namelist file and fill mcnl structure (used by CTEMEBSDmaster.f90)
!
!> @param nmlfile namelist file name
!> @param emnl EBSD master name list structure
!
!> @date 06/19/14  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine GetEBSDMasterNameList(nmlfile, emnl)

use error

IMPLICIT NONE

character(fnlen),INTENT(IN)                     :: nmlfile
type(EBSDMasterNameListType),INTENT(INOUT)      :: emnl

integer(kind=irg)       :: stdout
integer(kind=irg)       :: npx
integer(kind=irg)       :: Esel
integer(kind=irg)       :: nthreads
real(kind=sgl)          :: dmin
character(fnlen)        :: energyfile
character(fnlen)        :: outname

! define the IO namelist to facilitate passing variables to the program.
namelist /EBSDmastervars/ dmin,npx,nthreads,outname,energyfile,Esel

! set the input parameters to default values (except for xtalname, which must be present)
stdout = 6
npx = 500                       ! Nx pixels (total = 2Nx+1)
nthreads = 1
Esel = -1                       ! selected energy value for single energy run
dmin = 0.025                    ! smallest d-spacing to include in dynamical matrix [nm]
energyfile = 'undefined'        ! default filename for z_0(E_e) data from CTEMMC Monte Carlo simulations
outname = 'EBSDmasterout.data'  ! default filename for final output

! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=EBSDmastervars)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(energyfile).eq.'undefined') then
  call FatalError('CTEMEBSDmaster:',' energy file name is undefined in '//nmlfile)
 end if

! if we get here, then all appears to be ok, and we need to fill in the emnl fields
emnl%stdout = stdout
emnl%npx = npx
emnl%Esel = Esel
emnl%nthreads = nthreads
emnl%dmin = dmin
emnl%energyfile = energyfile
emnl%outname = outname

end subroutine GetEBSDMasterNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:GetECPMasterNameList
!
!> @author Saransh Singh/Marc De Graef, Carnegie Mellon University
!
!> @brief read namelist file and fill mcnl structure (used by CTEMECPmaster.f90)
!
!> @param nmlfile namelist file name
!> @param emnl ECP master name list structure
!
!> @date 06/19/14  SS 1.0 new routine
!--------------------------------------------------------------------------
subroutine GetECPMasterNameList(nmlfile, ecpnl)

use error

IMPLICIT NONE

character(fnlen),INTENT(IN)                     :: nmlfile
type(ECPMasterNameListType),INTENT(INOUT)      :: ecpnl

integer(kind=irg)       :: stdout
integer(kind=irg)       :: npx
integer(kind=irg)       :: Esel
integer(kind=irg)       :: numthick
real(kind=irg)          :: startthick
real(kind=irg)          :: fn(3)
real(kind=sgl)          :: dmin
real(kind=sgl)          :: zintstep
real(kind=sgl)          :: abcdist(3)
real(kind=sgl)          :: albegadist(3)
real(kind=sgl)          :: thickinc
character(fnlen)        :: compmode
character(fnlen)        :: energyfile
character(fnlen)        :: outname
logical                 :: distort

! define the IO namelist to facilitate passing variables to the program.
namelist /ECPmastervars/ stdout, startthick, dmin, fn, abcdist, albegadist, compmode, &
    distort, outname, energyfile, Esel, npx

! set the input parameters to default values (except for xtalname, which must be present)
stdout = 6
startthick = 2.0
fn = (/0.0, 0.0, 1.0/)
Esel = -1                       ! selected energy value for single energy run
dmin = 0.025                    ! smallest d-spacing to include in dynamical matrix [nm]
npx = 256
abcdist = (/0.4, 0.4, 0.4/)
albegadist = (/90.0, 90.0, 90.0/)
compmode = 'Blochwv'
distort = .FALSE.
energyfile = 'undefined'        ! default filename for z_0(E_e) data from CTEMMC Monte Carlo simulations
outname = 'ECPmasterout.data'  ! default filename for final output

! read the namelist file
open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
read(UNIT=dataunit,NML=ECPmastervars)
close(UNIT=dataunit,STATUS='keep')

! check for required entries
if (trim(energyfile).eq.'undefined') then
call FatalError('CTEMECPmaster:',' energy file name is undefined in '//nmlfile)
end if

! if we get here, then all appears to be ok, and we need to fill in the emnl fields
ecpnl%stdout = stdout
ecpnl%Esel = Esel
ecpnl%npx = npx
ecpnl%startthick = startthick
ecpnl%fn = fn
ecpnl%abcdist = abcdist
ecpnl%albegadist = albegadist
ecpnl%dmin = dmin
ecpnl%compmode = compmode
ecpnl%distort = distort
ecpnl%energyfile = energyfile
ecpnl%outname = outname

end subroutine GetECPMasterNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:GetECPNameList
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief read namelist file and fill mcnl structure (used by CTEMECP.f90)
!
!> @param nmlfile namelist file name
!> @param emnl ECP name list structure
!
!> @date 06/19/14  SS 1.0 new routine
!--------------------------------------------------------------------------
subroutine GetECPpatternNameList(nmlfile,ecpnl)

use error

IMPLICIT NONE

character(fnlen),INTENT(IN)                     :: nmlfile
type(ECPpatternNameListType),INTENT(INOUT)             :: ecpnl

integer(kind=irg)       :: stdout
integer(kind=irg)       :: npix
real(kind=sgl)          :: thetac
real(kind=sgl)          :: k(3)
character(fnlen)        :: masterfile
character(fnlen)        :: outname

! define the IO namelist to facilitate passing variables to the program.
namelist /ECPvars/ stdout, npix, masterfile, outname, thetac, k

! set the input parameters to default values (except for xtalname, which must be present)
stdout = 6
npix = 256
thetac = 5.0
k = (/0.0,0.0,1.0/)
masterfile = 'undefined'        ! default filename for master data from CTEMECPmaster
outname = 'ECP.data'  ! default filename for final output

! read the namelist file
open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
read(UNIT=dataunit,NML=ECPvars)
close(UNIT=dataunit,STATUS='keep')

! check for required entries
if (trim(masterfile).eq.'undefined') then
call FatalError('CTEMECP:',' master file name is undefined in '//nmlfile)
end if

! if we get here, then all appears to be ok, and we need to fill in the emnl fields
ecpnl%stdout = stdout
ecpnl%npix = npix
ecpnl%thetac = thetac
ecpnl%k = k
ecpnl%masterfile = masterfile
ecpnl%outname = outname

end subroutine GetECPpatternNameList

end module NameListHandlers
