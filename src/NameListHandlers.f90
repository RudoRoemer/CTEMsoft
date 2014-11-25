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
! SUBROUTINE:GetKosselMasterNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read namelist file and fill knl structure (used by CTEMKosselmaster.f90)
!
!> @param nmlfile namelist file name
!> @param knl Kossel name list structure
!
!> @date 09/09/14  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine GetKosselMasterNameList(nmlfile, knl)

use error

IMPLICIT NONE

character(fnlen),INTENT(IN)             :: nmlfile
type(KosselMasterNameListType),INTENT(INOUT)  :: knl

integer(kind=irg)       :: stdout
integer(kind=irg)       :: numthick
integer(kind=irg)       :: npix
integer(kind=irg)       :: nthreads
real(kind=sgl)          :: voltage
real(kind=sgl)          :: dmin
real(kind=sgl)          :: startthick
real(kind=sgl)          :: thickinc
character(fnlen)        :: xtalname
character(fnlen)        :: outname

namelist /Kosselmasterlist/ stdout, xtalname, voltage, dmin,  nthreads, &
                              startthick, thickinc, numthick, outname, npix

! set the input parameters to default values (except for xtalname, which must be present)
stdout = 6                      ! standard output
numthick = 10                   ! number of increments
npix = 256                      ! output arrays will have size npix x npix
nthreads = 4                    ! default number of threads for OpenMP
voltage = 200000.0              ! acceleration voltage [V]
dmin = 0.025                    ! smallest d-spacing to include in dynamical matrix [nm]
startthick = 10.0               ! starting thickness [nm]
thickinc = 10.0                 ! thickness increment
xtalname = 'undefined'          ! initial value to check that the keyword is present in the nml file
outname = 'Kosselout.data'      ! output filename

! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=Kosselmasterlist)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(xtalname).eq.'undefined') then
  call FatalError('CTEMKosselMaster:',' structure file name is undefined in '//nmlfile)
 end if

! if we get here, then all appears to be ok, and we need to fill in the knl fields
knl%stdout = stdout
knl%numthick = numthick
knl%npix = npix
knl%nthreads = nthreads
knl%voltage = voltage
knl%dmin = dmin
knl%startthick = startthick
knl%thickinc = thickinc
knl%xtalname = xtalname
knl%outname = outname

end subroutine GetKosselMasterNameList


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
!> @date 06/18/14  SS 1.0 new routine
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
! SUBROUTINE:GetMCCLMultiLayerNameList
!
!> @author Saransh Singh/Marc De Graef, Carnegie Mellon University
!
!> @brief read namelist file and fill mcnl structure (used by CTEMMCCL.f90)
!
!> @param nmlfile namelist file name
!> @param mcnl Monte Carloname list structure
!
!> @date 06/18/14  SS 1.0 new routine
!--------------------------------------------------------------------------
subroutine GetMCCLMultiLayerNameList(nmlfile, mcnl)

use error

IMPLICIT NONE

character(fnlen),INTENT(IN)             :: nmlfile
type(MCCLMultiLayerNameListType),INTENT(INOUT)      :: mcnl

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
real(kind=dbl)          :: filmthickness
real(kind=dbl)          :: filmstep
character(4)            :: MCmode
character(fnlen)        :: xtalname_film
character(fnlen)        :: xtalname_subs
character(fnlen)        :: dataname
character(fnlen)        :: primelist
character(fnlen)        :: mode

! define the IO namelist to facilitate passing variables to the program.
namelist  / MCCLdata / stdout, sig, numsx, num_el, globalworkgrpsz, EkeV, &
        dataname, primelist, totnum_el, Ehistmin, Ebinsize, depthmax, &
        depthstep, omega, MCmode, mode, xtalname_film, xtalname_subs, &
        filmthickness, filmstep


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
xtalname_film = 'undefined'
xtalname_subs = 'undefined'
dataname = 'MCoutput.data'
primelist = 'RandomSeeds.data'
mode = 'full'
filmthickness = 20.D0
filmstep = 2.D0

! read the namelist file
open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
read(UNIT=dataunit,NML=MCCLdata)
close(UNIT=dataunit,STATUS='keep')

! check for required entries
if ((trim(xtalname_film).eq.'undefined') .or. (trim(xtalname_subs).eq.'undefined')) then
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
mcnl%xtalname_film = xtalname_film
mcnl%xtalname_subs = xtalname_subs
mcnl%dataname = dataname
mcnl%primelist = primelist
mcnl%mode = mode
mcnl%filmthickness = filmthickness
mcnl%filmstep = filmstep

end subroutine GetMCCLMultiLayerNameList

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
! SUBROUTINE:GetEBSDNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read namelist file and fill enl structure (used by CTEMEBSD.f90)
!
!> @param nmlfile namelist file name
!> @param enl EBSD name list structure
!
!> @date 06/23/14  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine GetEBSDNameList(nmlfile, enl)

use error

IMPLICIT NONE

character(fnlen),INTENT(IN)               :: nmlfile
type(EBSDNameListType),INTENT(INOUT)      :: enl

integer(kind=irg)       :: stdout
integer(kind=irg)       :: numsx
integer(kind=irg)       :: numsy
integer(kind=irg)       :: binning
integer(kind=irg)       :: nthreads
real(kind=sgl)          :: L
real(kind=sgl)          :: thetac
real(kind=sgl)          :: delta
real(kind=sgl)          :: xpc
real(kind=sgl)          :: ypc
real(kind=sgl)          :: energymin
real(kind=sgl)          :: energymax
real(kind=sgl)          :: gammavalue
real(kind=sgl)          :: axisangle(4)
real(kind=dbl)          :: beamcurrent
real(kind=dbl)          :: dwelltime
character(3)            :: scalingmode
character(3)            :: eulerconvention
character(fnlen)        :: anglefile
character(fnlen)        :: masterfile
character(fnlen)        :: energyfile 
character(fnlen)        :: datafile

! define the IO namelist to facilitate passing variables to the program.
namelist  / EBSDdata / stdout, L, thetac, delta, numsx, numsy, xpc, ypc, anglefile, eulerconvention, masterfile, &
                        energyfile, datafile, beamcurrent, dwelltime, energymin, energymax, binning, gammavalue, &
                        scalingmode, axisangle, nthreads

! set the input parameters to default values (except for xtalname, which must be present)
stdout          = 6
numsx           = 640           ! [dimensionless]
numsy           = 480           ! [dimensionless]
binning         = 1             ! binning mode  (1, 2, 4, or 8)
L               = 20000.0       ! [microns]
nthreads        = 1             ! number of OpenMP threads
thetac          = 0.0           ! [degrees]
delta           = 25.0          ! [microns]
xpc             = 0.0           ! [pixels]
ypc             = 0.0           ! [pixels]
energymin       = 15.0          ! minimum energy to consider
energymax       = 30.0          ! maximum energy to consider
gammavalue      = 1.0           ! gamma factor
axisangle       = (/0.0, 0.0, 0.0, 0.0/)        ! no additional axis angle rotation
beamcurrent     = 14.513D0      ! beam current (actually emission current) in nano ampere
dwelltime       = 100.0D0       ! in microseconds
scalingmode     = 'not'         ! intensity selector ('lin', 'gam', or 'not')
eulerconvention = 'tsl'         ! convention for the first Euler angle ['tsl' or 'hkl']
anglefile       = 'undefined'   ! filename
masterfile      = 'undefined'   ! filename
energyfile      = 'undefined'   ! name of file that contains energy histograms for all scintillator pixels (output from MC program)
datafile        = 'undefined'   ! output file name


! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=EBSDdata)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(energyfile).eq.'undefined') then
  call FatalError('CTEMEBSD:',' energy file name is undefined in '//nmlfile)
 end if

 if (trim(anglefile).eq.'undefined') then
  call FatalError('CTEMEBSD:',' angle file name is undefined in '//nmlfile)
 end if

 if (trim(masterfile).eq.'undefined') then
  call FatalError('CTEMEBSD:',' master pattern file name is undefined in '//nmlfile)
 end if

 if (trim(datafile).eq.'undefined') then
  call FatalError('CTEMEBSD:',' output file name is undefined in '//nmlfile)
 end if

! if we get here, then all appears to be ok, and we need to fill in the emnl fields
enl%stdout = stdout
enl%numsx = numsx
enl%numsy = numsy
enl%binning = binning
enl%L = L
enl%nthreads = nthreads
enl%thetac = thetac
enl%delta = delta
enl%xpc = xpc
enl%ypc = ypc
enl%energymin = energymin
enl%energymax = energymax
enl%gammavalue = gammavalue
enl%axisangle = axisangle
enl%beamcurrent = beamcurrent
enl%dwelltime = dwelltime
enl%scalingmode = scalingmode
enl%eulerconvention = eulerconvention
enl%anglefile = anglefile
enl%masterfile = masterfile
enl%energyfile = energyfile
enl%datafile = datafile

end subroutine GetEBSDNameList

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
!> @date 11/25/14  MDG 2.0 added parameters for film on substrate mode
!--------------------------------------------------------------------------
subroutine GetECPNameList(nmlfile, ecpnl)

use error

IMPLICIT NONE

character(fnlen),INTENT(IN)             :: nmlfile
type(ECPNameListType),INTENT(INOUT)     :: ecpnl

integer(kind=irg)       :: stdout
integer(kind=irg)       :: k(3)
integer(kind=irg)       :: fn(3)
integer(kind=irg)       :: numthick
integer(kind=irg)       :: npix
integer(kind=irg)       :: nthreads
integer(kind=irg)       :: gF(3)
integer(kind=irg)       :: gS(3)
integer(kind=irg)       :: tF(3)
integer(kind=irg)       :: tS(3)
real(kind=sgl)          :: voltage
real(kind=sgl)          :: dmin
real(kind=sgl)          :: ktmax
real(kind=sgl)          :: thetac
real(kind=sgl)          :: startthick
real(kind=sgl)          :: thickinc
real(kind=sgl)          :: zintstep
real(kind=sgl)          :: filmthickness
character(7)            :: compmode
character(fnlen)        :: outname
character(fnlen)        :: xtalname
character(fnlen)        :: xtalname2
character(fnlen)        :: energyfile

! namelist /ECPlist/ stdout, xtalname, voltage, k, fn, dmin, distort, abcdist, albegadist, ktmax, &
namelist /ECPlist/ stdout, xtalname, xtalname2, voltage, k, fn, dmin, ktmax, filmthickness, &
                   startthick, thickinc, nthreads, numthick, npix, outname, thetac, compmode, zintstep, &
                   gF, gS, tF, tS, energyfile

! default values
stdout = 6                              ! standard output
k = (/ 0, 0, 1 /)                       ! beam direction [direction indices]
fn = (/ 0, 0, 1 /)                      ! foil normal [direction indices]
gF = (/ 0, 0, 0 /)                      ! plane normal in film
gS = (/ 0, 0, 0 /)                      ! plane normal in substrate
tF = (/ 0, 0, 0 /)                      ! direction in film
tS = (/ 0, 0, 0 /)                      ! direction in substrate
numthick = 10                           ! number of increments
npix = 256                              ! output arrays will have size npix x npix
nthreads = 1                            ! number of OpenMP threads
voltage = 30000.0                       ! acceleration voltage [V]
dmin = 0.025                            ! smallest d-spacing to include in dynamical matrix [nm]
ktmax = 0.0                             ! beam convergence in units of |g_a|
thetac = 0.0                            ! beam convergence in mrad (either ktmax or thetac must be given)
startthick = 2.0                        ! starting thickness [nm]
thickinc = 2.0                          ! thickness increment
zintstep = 1.0                          ! integration step size for ScatMat mode
filmthickness = 0.0                     ! 0.0 if there is no film
compmode = 'Blochwv'                    ! 'Blochwv' or 'ScatMat' solution mode (Bloch is default)
outname = 'ecp.data'                    ! output filename
xtalname = 'undefined'                  ! initial value to check that the keyword is present in the nml file
xtalname2 = 'undefined'                 ! initial value for substrate structure name
energyfile = 'undefined'
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=ECPlist)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(xtalname).eq.'undefined') then
  call FatalError('CTEMEECP:',' crystal file name is undefined in '//nmlfile)
 end if

ecpnl%stdout = stdout
ecpnl%k = k
ecpnl%fn = fn
ecpnl%gF = gF
ecpnl%gS = gS
ecpnl%tF = tF
ecpnl%tS = tS
ecpnl%numthick = numthick
ecpnl%npix = npix
ecpnl%nthreads = nthreads
ecpnl%voltage = voltage
ecpnl%dmin = dmin
ecpnl%ktmax = ktmax
ecpnl%thetac = thetac
ecpnl%startthick = startthick
ecpnl%thickinc = thickinc
ecpnl%zintstep = zintstep
ecpnl%filmthickness = filmthickness
ecpnl%compmode = compmode
ecpnl%outname = outname
ecpnl%xtalname = xtalname
ecpnl%xtalname2 = xtalname2
ecpnl%energyfile = energyfile

end subroutine GetECPNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:GetLACBEDNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read namelist file and fill lacbednl structure (used by CTEMLACBED.f90)
!
!> @param nmlfile namelist file name
!> @param lacbednl LACBED name list structure
!
!> @date 07/01/14  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine GetLACBEDNameList(nmlfile, lacbednl)

use error

IMPLICIT NONE

character(fnlen),INTENT(IN)             :: nmlfile
type(LACBEDNameListType),INTENT(INOUT)  :: lacbednl

integer(kind=irg)       :: stdout
integer(kind=irg)       :: k(3)
integer(kind=irg)       :: fn(3)
integer(kind=irg)       :: maxHOLZ
integer(kind=irg)       :: numthick
integer(kind=irg)       :: npix
integer(kind=irg)       :: nthreads
real(kind=sgl)          :: voltage
real(kind=sgl)          :: dmin
real(kind=sgl)          :: convergence
real(kind=sgl)          :: startthick
real(kind=sgl)          :: thickinc
real(kind=sgl)          :: minten
character(fnlen)        :: xtalname
character(fnlen)        :: outname

namelist /inputlist/ stdout, xtalname, voltage, k, fn, dmin, convergence, minten, &
                              nthreads, startthick, thickinc, numthick, outname, npix, maxHOLZ

stdout = 6                      ! standard output
k = (/ 0, 0, 1 /)               ! beam direction [direction indices]
fn = (/ 0, 0, 1 /)              ! foil normal [direction indices]
maxHOLZ = 2                     ! maximum HOLZ layer index to be used for the output file; note that his number
                                ! does not affect the actual computations; it only determines which reflection 
                                ! families will end up in the output file
numthick = 10                   ! number of increments
npix = 256                      ! output arrays will have size npix x npix
nthreads = 1                    ! number of computational threads
voltage = 200000.0              ! acceleration voltage [V]
dmin = 0.025                    ! smallest d-spacing to include in dynamical matrix [nm]
convergence = 25.0              ! beam convergence angle [mrad]
startthick = 10.0               ! starting thickness [nm]
thickinc = 10.0                 ! thickness increment
minten = 1.0E-5                 ! minimum intensity in diffraction disk to make it into the output file
xtalname = 'undefined'          ! initial value to check that the keyword is present in the nml file
outname = 'lacbedout.data'      ! output filename

! read the namelist file
open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
read(UNIT=dataunit,NML=inputlist)
close(UNIT=dataunit,STATUS='keep')

! check for required entries
if (trim(xtalname).eq.'undefined') then
  call FatalError('CTEMLACBED:',' structure file name is undefined in '//nmlfile)
end if

lacbednl%stdout = stdout
lacbednl%k = k
lacbednl%fn = fn
lacbednl%maxHOLZ = maxHOLZ
lacbednl%numthick = numthick
lacbednl%npix = npix
lacbednl%nthreads = nthreads
lacbednl%voltage = voltage
lacbednl%dmin = dmin
lacbednl%convergence = convergence
lacbednl%startthick = startthick
lacbednl%thickinc = thickinc
lacbednl%minten = minten
lacbednl%xtalname = xtalname
lacbednl%outname = outname

end subroutine GetLACBEDNameList


!--------------------------------------------------------------------------
!
! SUBROUTINE:GetECPpatternNameList
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief read namelist file and fill mcnl structure (used by CTEMECPpattern.f90)
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

! set the input parameters to default values (except for masterfile, which must be present)
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

!--------------------------------------------------------------------------
!
! SUBROUTINE:GetPEDNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read namelist file and fill pednl structure (used by CTEMPED.f90)
!
!> @param nmlfile namelist file name
!> @param pednl PED name list structure
!
!> @date 07/09/14 MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine GetPEDNameList(nmlfile,pednl)

use error

IMPLICIT NONE

character(fnlen),INTENT(IN)                     :: nmlfile
type(PEDNameListType),INTENT(INOUT)             :: pednl

integer(kind=irg)       :: stdout
integer(kind=irg)       :: k(3)
integer(kind=irg)       :: fn(3)
integer(kind=irg)       :: precsample
integer(kind=irg)       :: precazimuthal
integer(kind=irg)       :: npix
real(kind=sgl)          :: voltage
real(kind=sgl)          :: dmin
real(kind=sgl)          :: precangle
real(kind=sgl)          :: prechalfwidth
real(kind=sgl)          :: thickness
real(kind=sgl)          :: camlen
character(5)            :: filemode
character(fnlen)        :: xtalname
character(fnlen)        :: outname

! define the IO namelist to facilitate passing variables to the program.
namelist /inputlist/ stdout, xtalname, voltage, k, fn, dmin, precangle, prechalfwidth, precsample, precazimuthal, &
                              thickness,  outname, npix, camlen, filemode

! set the input parameters to default values (except for xtalname, which must be present)
xtalname = 'undefined'          ! initial value to check that the keyword is present in the nml file
stdout = 6                      ! standard output
voltage = 200000.0              ! acceleration voltage [V]
k = (/ 0, 0, 1 /)               ! beam direction [direction indices]
fn = (/ 0, 0, 1 /)              ! foil normal [direction indices]
dmin = 0.025                    ! smallest d-spacing to include in dynamical matrix [nm]
precangle = 10.472              ! beam precession angle [mrad]; default = 0.6 degrees
prechalfwidth = 0.25            ! beam half width in the tilt direction [mrad]
precsample = 10                 ! number of samples (concentric circles) in beam half width (total = 2*precsample + 1)
precazimuthal = 360             ! number of azimuthal samples for each precession circle
thickness = 10.0                ! sample thickness [nm]
filemode = 'total'              ! 'total' mode or 'eachp'
npix = 256                      ! output arrays will have size npix x npix
outname = 'pedout.data'         ! output filename
camlen = 1000.0                 ! camera length [mm]


! read the namelist file
open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
read(UNIT=dataunit,NML=inputlist)
close(UNIT=dataunit,STATUS='keep')

! check for required entries
if (trim(xtalname).eq.'undefined') then
call FatalError('CTEMPED:',' crystal structure file name is undefined in '//nmlfile)
end if

! if we get here, then all appears to be ok, and we need to fill in the pednl fields
pednl%xtalname = xtalname
pednl%stdout = stdout
pednl%voltage = voltage
pednl%k = k
pednl%fn = fn
pednl%dmin = dmin
pednl%precangle = precangle
pednl%prechalfwidth = prechalfwidth
pednl%precsample = precsample
pednl%precazimuthal = precazimuthal
pednl%thickness = thickness
pednl%filemode = filemode
pednl%npix = npix
pednl%outname = outname
pednl%camlen = camlen

end subroutine GetPEDNameList


!--------------------------------------------------------------------------
!
! SUBROUTINE:GetECCINameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read namelist file and fill eccinl structure (used by CTEMECCI.f90)
!
!> @param nmlfile namelist file name
!> @param eccinl ECCI name list structure
!
!> @date 10/04/14 MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine GetECCINameList(nmlfile,eccinl)

use error

IMPLICIT NONE

character(fnlen),INTENT(IN)                     :: nmlfile
type(ECCINameListType),INTENT(INOUT)            :: eccinl

integer(kind=irg)       :: stdout
integer(kind=irg)       :: nthreads,i
integer(kind=irg)       :: k(3)
integer(kind=irg)       :: nktstep
integer(kind=irg)       :: DF_npix
integer(kind=irg)       :: DF_npiy
integer(kind=irg)       :: numYdisl
integer(kind=irg)       :: numdisl
integer(kind=irg)       :: numsf
real(kind=sgl)          :: voltage
real(kind=sgl)          :: dkt
real(kind=sgl)          :: ktmax
real(kind=sgl)          :: lauec(2)
real(kind=sgl)          :: lauec2(2)
real(kind=sgl)          :: dmin
real(kind=sgl)          :: DF_L
real(kind=sgl)          :: DF_slice
character(4)            :: dispmode
character(4)            :: summode
character(5)            :: progmode
character(fnlen)        :: xtalname
character(fnlen)        :: foilnmlfile
character(fnlen)        :: dispfile
character(fnlen)        :: dataname
character(fnlen)        :: ECPname
character(fnlen)        :: dislYname(3*maxdefects)
character(fnlen)        :: dislname(3*maxdefects)
character(fnlen)        :: sfname(maxdefects)
character(fnlen)        :: sgname
character(fnlen)        :: apbname
character(fnlen)        :: incname
character(fnlen)        :: voidname


! define the IO namelist to facilitate passing variables to the program.
namelist / ECCIlist / DF_L, DF_npix, DF_npiy, DF_slice, dmin, sgname, incname, stdout, &
                      voidname, numdisl, dislname, numYdisl, dislYname, numsf, sfname, &
		      progmode, dispfile, ktmax, dkt, ECPname, summode, lauec, lauec2, &
	   	      dispmode, nthreads, xtalname, voltage, k, nktstep, &
	   	      dataname, foilnmlfile, apbname

! set the input parameters to default values (except for xtalname, which must be present)
stdout = 6
nthreads = 1
k = (/ 0,0,1 /)
nktstep = 20
DF_npix = 256
DF_npiy = 256
numYdisl = 0
numdisl = 0
numsf = 0
voltage = 30000.
dkt = 0.1
ktmax = 5.0
lauec = (/ 0.0, 0.0 /)
lauec2 = (/ 0.0, 0.0 /)
dmin = 0.1
DF_L = 1.0
DF_slice = 1.0
dispmode = 'not'
summode = 'diag'
progmode = 'array'
xtalname = 'undefined'
foilnmlfile = 'FOIL_rundata.nml'
dispfile = 'displacements.data'
dataname = 'ECCIout.data'
ECPname = 'undefined'
dislYname = (/ ('',i=1,3*maxdefects) /)
dislname = (/ ('',i=1,3*maxdefects) /)
sfname = (/ ('',i=1,maxdefects) /)
sgname = 'nofile'
apbname = 'none'
incname = 'none'   
voidname = 'none'

! read the namelist file
open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
read(UNIT=dataunit,NML=ECCIlist)
close(UNIT=dataunit,STATUS='keep')

! check for required entries
if (trim(xtalname).eq.'undefined') then
  call FatalError('CTEMECCI:',' crystal structure file name is undefined in '//nmlfile)
end if

! make sure the ECPname variable has been properly defined
if (trim(ECPname).eq.'undefined') then
  call FatalError('CTEMECCI:',' ECP pattern file name is undefined in '//nmlfile)
end if

! if we get here, then all appears to be ok, and we need to fill in the emnl fields
eccinl%stdout = stdout
eccinl%nthreads = nthreads
eccinl%k = k
eccinl%nktstep = nktstep
eccinl%DF_npix = DF_npix
eccinl%DF_npiy = DF_npiy
eccinl%numYdisl = numYdisl
eccinl%numdisl = numdisl
eccinl%numsf = numsf
eccinl%voltage = voltage
eccinl%dkt = dkt
eccinl%ktmax = ktmax
eccinl%lauec = lauec
eccinl%lauec2 = lauec2
eccinl%dmin = dmin
eccinl%DF_L = DF_L
eccinl%DF_slice = DF_slice
eccinl%dispmode = dispmode
eccinl%summode = summode
eccinl%progmode = progmode
eccinl%xtalname = xtalname
eccinl%foilnmlfile = foilnmlfile
eccinl%dispfile = dispfile
eccinl%dataname = dataname
eccinl%ECPname = ECPname
eccinl%dislYname = dislYname
eccinl%dislname = dislname
eccinl%sfname = sfname
eccinl%sgname = sgname
eccinl%apbname = apbname
eccinl%incname = incname
eccinl%voidname = voidname

end subroutine GetECCINameList



end module NameListHandlers
