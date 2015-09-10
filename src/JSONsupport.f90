! ###################################################################
! Copyright (c) 2013-2015, Marc De Graef/Carnegie Mellon University
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
! EMsoft:JSONsupport.f90
!--------------------------------------------------------------------------
!
! MODULE: JSONsupport
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief routines for conversion between json and nml files and reading of json files
!
!> @date  08/11/15 MDG 1.0 original
!> @date  08/12/15 MDG 1.1 added all routines currently also in NameListHDFwriters.f90
!> @date  08/12/15 MDG 1.2 replaced all the json_failed stuff by short routine JSON_failtest
!--------------------------------------------------------------------------
module JSONsupport

use local
use typedefs
use NameListTypedefs
use json_module
use, intrinsic :: iso_fortran_env , only: error_unit, wp => real64

IMPLICIT NONE

contains


!--------------------------------------------------------------------------
!
! SUBROUTINE:JSON_failtest
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief executes the json_fail routine; mostly to shorten the remaining code a little
!
!> @param error_cnt error counter
!
!> @date 08/12/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSON_failtest(error_cnt)

IMPLICIT NONE

integer(kind=irg),INTENT(INOUT)         :: error_cnt

if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

end subroutine JSON_failtest

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSON_writeNMLintegers
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write a series of integer namelist entries to a json structure 
!
!> @param inp json structure pointer
!> @param io_int list of integers
!> @param intlist list of string descriptors
!> @param n_int number of entries
!> @param error_cnt error counter
!
!> @date 08/11/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSON_writeNMLintegers(inp, io_int, intlist, n_int, error_cnt)

IMPLICIT NONE

type(json_value),INTENT(INOUT),pointer                :: inp
integer(kind=irg),INTENT(IN)                          :: io_int(n_int)
character(20),INTENT(IN)                              :: intlist(n_int)
integer(kind=irg),INTENT(IN)                          :: n_int
integer(kind=irg),INTENT(INOUT)                       :: error_cnt

integer(kind=irg)                                     :: i
character(fnlen)                                      :: dataset

do i=1,n_int
  dataset = intlist(i)
  call json_add(inp, dataset, io_int(i)); call JSON_failtest(error_cnt)
end do

end subroutine JSON_writeNMLintegers

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSON_writeNMLreals
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write a series of real namelist entries to a json structure
!
!> @param inp pointer to json_value 
!> @param io_real list of reals
!> @param reallist list of string descriptors
!> @param n_real number of entries
!> @param error_cnt error counter
!
!> @date 08/11/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSON_writeNMLreals(inp, io_real, reallist, n_real, error_cnt)

IMPLICIT NONE

type(json_value),INTENT(INOUT),pointer                :: inp
real(kind=sgl),INTENT(IN)                             :: io_real(n_real)
character(20),INTENT(IN)                              :: reallist(n_real)
integer(kind=irg),INTENT(IN)                          :: n_real

integer(kind=irg)                                     :: hdferr, i
character(fnlen)                                      :: dataset
integer(kind=irg),INTENT(INOUT)                       :: error_cnt

do i=1,n_real
  dataset = reallist(i)
  call json_add(inp, dataset, dble(io_real(i))); call JSON_failtest(error_cnt)
end do

end subroutine JSON_writeNMLreals

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSON_writeNMLdoubles
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write a series of double namelist entries to a json structure
!
!> @param inp pointer to json_value 
!> @param io_real list ofadoubles 
!> @param reallist list of string descriptors
!> @param n_real number of entries
!> @param error_cnt error counter
!
!> @date 08/11/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSON_writeNMLdoubles(inp, io_real, reallist, n_real, error_cnt)

IMPLICIT NONE

type(json_value),INTENT(INOUT),pointer                :: inp
real(kind=dbl),INTENT(IN)                             :: io_real(n_real)
character(20),INTENT(IN)                              :: reallist(n_real)
integer(kind=irg),INTENT(IN)                          :: n_real

integer(kind=irg)                                     :: hdferr, i
character(fnlen)                                      :: dataset
integer(kind=irg),INTENT(INOUT)                       :: error_cnt

do i=1,n_real
  dataset = reallist(i)
  call json_add(inp, dataset, io_real(i)); call JSON_failtest(error_cnt)
end do

end subroutine JSON_writeNMLdoubles

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSON_initpointers
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief initialize the necessary pointers to write a namelist json file
!
!> @param inp pointer to json_value 
!> @param io_real list of reals
!> @param reallist list of string descriptors
!> @param n_real number of entries
!> @param error_cnt error counter
!
!> @date 08/11/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSON_initpointers(p, inp, jsonname, namelistname, error_cnt)

IMPLICIT NONE

type(json_value),INTENT(INOUT),pointer  :: p, inp
character(fnlen),INTENT(IN)             :: jsonname, namelistname
integer(kind=irg),INTENT(INOUT)         :: error_cnt

! initialize the json state variables
error_cnt = 0
call json_initialize(); call JSON_failtest(error_cnt)

! create the json root pointer
call json_create_object(p,trim(jsonname)); call JSON_failtest(error_cnt)

! we'll use the namelist name to configure the inp structure and add it to p
call json_create_object(inp,trim(namelistname)); call JSON_failtest(error_cnt) 
call json_add(p, inp); call JSON_failtest(error_cnt)

end subroutine JSON_initpointers

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSON_cleanuppointers
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief clean up the pointers and write the json file
!
!> @param p pointer to json_value 
!> @param inp pointer to json_value 
!> @param jsonname json output file name
!> @param error_cnt error counter
!
!> @date 08/11/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSON_cleanuppointers(p, inp, jsonname, error_cnt)

use io 

IMPLICIT NONE

type(json_value),INTENT(INOUT),pointer  :: p, inp
character(fnlen),INTENT(IN)             :: jsonname
integer(kind=irg),INTENT(INOUT)         :: error_cnt

! get rid of inp
nullify(inp)

! write the json file
open(unit=dataunit, file=trim(jsonname), status='REPLACE')
call json_print(p,dataunit); call JSON_failtest(error_cnt)
close(dataunit)

! final cleanup
call json_destroy(p); call JSON_failtest(error_cnt)

end subroutine JSON_cleanuppointers


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! from here on we have the Namelist->json conversion routines for all the 
! namelists defined in the NameListTypedefs.f90 file.
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONwriteKosselNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist file into json file
!
!> @param knl Kossel name list structure
!
!> @date 08/11/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONwriteKosselNameList(knl, jsonname, error_cnt)

use ISO_C_BINDING

IMPLICIT NONE

type(KosselNameListType),INTENT(IN)                   :: knl
character(fnlen),INTENT(IN)                           :: jsonname
integer(kind=irg),INTENT(INOUT)                       :: error_cnt

type(json_value),pointer                              :: p, inp

integer(kind=irg),parameter                           :: n_int = 5, n_real = 6
integer(kind=irg)                                     :: hdferr,  io_int(n_int)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
character(fnlen)                                      :: dataset, namelistname
character(fnlen,kind=c_char)                          :: line2(1)

! initialize the json state variables
namelistname = 'Kossellist'
call JSON_initpointers(p, inp, jsonname, namelistname, error_cnt)

! then we need to add all the necessary fields to the inp structure

! write all the single integers
io_int = (/ knl%stdout, knl%numthick, knl%npix, knl%maxHOLZ, knl%nthreads /)
intlist(1) = 'stdout'
intlist(2) = 'numthick'
intlist(3) = 'npix'
intlist(4) = 'maxHOLZ'
intlist(5) = 'nthreads'
call JSON_writeNMLintegers(inp, io_int, intlist, n_int, error_cnt)

! integer vectors
dataset = 'k'
call json_add(inp, dataset, knl%k); call JSON_failtest(error_cnt)

dataset = 'fn'
call json_add(inp, dataset, knl%fn); call JSON_failtest(error_cnt)

! write all the single reals
io_real = (/ knl%voltage, knl%dmin, knl%convergence, knl%startthick, knl%thickinc, knl%minten /)
reallist(1) = 'voltage'
reallist(2) = 'dmin'
reallist(3) = 'convergence'
reallist(4) = 'startthick'
reallist(5) = 'thickinc'
reallist(6) = 'minten'
call JSON_writeNMLreals(inp, io_real, reallist, n_real, error_cnt)

! write all the strings
dataset = 'xtalname'
call json_add(inp, dataset, knl%xtalname); call JSON_failtest(error_cnt)

dataset = 'outname'
call json_add(inp, dataset, knl%outname); call JSON_failtest(error_cnt)

! and then we write the file and clean up
call JSON_cleanuppointers(p, inp, jsonname, error_cnt)

end subroutine JSONwriteKosselNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONwriteKosselMasterNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist file into json file
!
!> @param knl Kossel name list structure
!> @param jsonname output file name
!> @param error_cnt total number of errors encountered by json routines
!
!> @date 08/11/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONwriteKosselMasterNameList(knl, jsonname, error_cnt)

use ISO_C_BINDING

IMPLICIT NONE

type(KosselMasterNameListType),INTENT(IN)             :: knl
character(fnlen),INTENT(IN)                           :: jsonname
integer(kind=irg),INTENT(INOUT)                       :: error_cnt

type(json_value),pointer                              :: p, inp

integer(kind=irg),parameter                           :: n_int = 4, n_real = 5
integer(kind=irg)                                     :: io_int(n_int)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
character(fnlen)                                      :: dataset, namelistname
character(fnlen,kind=c_char)                          :: line2(1)

! initialize the json state variables
namelistname = 'Kosselmasterlist'
call JSON_initpointers(p, inp, jsonname, namelistname, error_cnt)

! then we need to add all the necessary fields to the inp structure

! write all the single integers
io_int = (/ knl%stdout, knl%numthick, knl%npix, knl%nthreads /)
intlist(1) = 'stdout'
intlist(2) = 'numthick'
intlist(3) = 'npix'
intlist(4) = 'nthreads' 
call JSON_writeNMLintegers(inp, io_int, intlist, n_int, error_cnt)

! write all the single reals
io_real = (/ knl%voltage, knl%dmin, knl%startthick, knl%thickinc, knl%tfraction /)
reallist(1) = 'voltage'
reallist(2) = 'dmin'
reallist(3) = 'startthick'
reallist(4) = 'thickinc'
reallist(5) = 'tfraction' 
call JSON_writeNMLreals(inp, io_real, reallist, n_real, error_cnt)

! write all the strings
dataset = 'Kosselmode'
call json_add(inp, dataset, knl%Kosselmode); call JSON_failtest(error_cnt)

dataset = 'xtalname'
call json_add(inp, dataset, knl%xtalname); call JSON_failtest(error_cnt)

dataset = 'outname'
call json_add(inp, dataset, knl%outname); call JSON_failtest(error_cnt)

! and then we write the file and clean up
call JSON_cleanuppointers(p, inp, jsonname, error_cnt)

end subroutine JSONwriteKosselMasterNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONwriteMCNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist file into JSON file
!
!> @param mcnl Monte Carlo name list structure
!> @param jsonname output file name
!> @param error_cnt total number of errors encountered by json routines
!
!> @date 08/11/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONwriteMCNameList(mcnl, jsonname, error_cnt)

use ISO_C_BINDING

IMPLICIT NONE

type(MCNameListType),INTENT(INOUT)                    :: mcnl
character(fnlen),INTENT(IN)                           :: jsonname
integer(kind=irg),INTENT(INOUT)                       :: error_cnt

type(json_value),pointer                              :: p, inp

integer(kind=irg),parameter                           :: n_int = 5, n_real = 7
integer(kind=irg)                                     :: io_int(n_int)
real(kind=dbl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
character(fnlen)                                      :: dataset, sval(1), namelistname
character(fnlen,kind=c_char)                          :: line2(1)

! initialize the json state variables
namelistname = 'MCdata'
call JSON_initpointers(p, inp, jsonname, namelistname, error_cnt)

! write all the single integers
io_int = (/ mcnl%stdout, mcnl%numsx, mcnl%primeseed, mcnl%num_el, mcnl%nthreads /)
intlist(1) = 'stdout'
intlist(2) = 'numsx'
intlist(3) = 'primeseed'
intlist(4) = 'num_el'
intlist(5) = 'nthreads'
call JSON_writeNMLintegers(inp, io_int, intlist, n_int, error_cnt)

! write all the single doubles
io_real = (/ mcnl%sig, mcnl%omega, mcnl%EkeV, mcnl%Ehistmin, mcnl%Ebinsize, mcnl%depthmax, mcnl%depthstep /)
reallist(1) = 'sig'
reallist(2) = 'omega'
reallist(3) = 'EkeV'
reallist(4) = 'Ehistmin'
reallist(5) = 'Ebinsize'
reallist(6) = 'depthmax'
reallist(7) = 'depthstep'
call JSON_writeNMLdoubles(inp, io_real, reallist, n_real, error_cnt)

! write all the strings
dataset = 'MCmode'
call json_add(inp, dataset, mcnl%MCmode); call JSON_failtest(error_cnt)

dataset = 'xtalname'
call json_add(inp, dataset, mcnl%xtalname); call JSON_failtest(error_cnt)

dataset = 'dataname'
call json_add(inp, dataset, mcnl%dataname); call JSON_failtest(error_cnt)

! and then we write the file and clean up
call JSON_cleanuppointers(p, inp, jsonname, error_cnt)

end subroutine JSONwriteMCNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONwriteMCCLNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist to JSON file
!> @param jsonname output file name
!> @param error_cnt total number of errors encountered by json routines
!
!> @param mcnl Monte Carlo name list structure
!
!> @date 03/21/15 MDG 1.0 new routine
!> @date 09/09/15 MDG 1.1 added devid
!--------------------------------------------------------------------------
subroutine JSONwriteMCCLNameList(mcnl, jsonname, error_cnt)

use ISO_C_BINDING

IMPLICIT NONE

type(MCCLNameListType),INTENT(INOUT)                  :: mcnl
character(fnlen),INTENT(IN)                           :: jsonname
integer(kind=irg),INTENT(INOUT)                       :: error_cnt

type(json_value),pointer                              :: p, inp

integer(kind=irg),parameter                           :: n_int = 6, n_real = 7
integer(kind=irg)                                     :: io_int(n_int)
real(kind=dbl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
character(fnlen)                                      :: dataset, sval(1), namelistname
character(fnlen,kind=c_char)                          :: line2(1)

! initialize the json state variables
namelistname = 'MCCLdata'
call JSON_initpointers(p, inp, jsonname, namelistname, error_cnt)

! write all the single integers
io_int = (/ mcnl%stdout, mcnl%numsx, mcnl%globalworkgrpsz, mcnl%num_el, mcnl%totnum_el, mcnl%devid /)
intlist(1) = 'stdout'
intlist(2) = 'numsx'
intlist(3) = 'globalworkgrpsz'
intlist(4) = 'num_el'
intlist(5) = 'totnum_el'
intlist(6) = 'devid'
call JSON_writeNMLintegers(inp, io_int, intlist, n_int, error_cnt)

! write all the single doubles
io_real = (/ mcnl%sig, mcnl%omega, mcnl%EkeV, mcnl%Ehistmin, mcnl%Ebinsize, mcnl%depthmax, mcnl%depthstep /)
reallist(1) = 'sig'
reallist(2) = 'omega'
reallist(3) = 'EkeV'
reallist(4) = 'Ehistmin'
reallist(5) = 'Ebinsize'
reallist(6) = 'depthmax'
reallist(7) = 'depthstep'
call JSON_writeNMLdoubles(inp, io_real, reallist, n_real, error_cnt)

! write all the strings
dataset = 'MCmode'
call json_add(inp, dataset, mcnl%MCmode); call JSON_failtest(error_cnt)

dataset = 'xtalname'
call json_add(inp, dataset, mcnl%xtalname); call JSON_failtest(error_cnt)

dataset = 'dataname'
call json_add(inp, dataset, mcnl%dataname); call JSON_failtest(error_cnt)

dataset = 'mode'
call json_add(inp, dataset, mcnl%mode); call JSON_failtest(error_cnt)

! and then we write the file and clean up
call JSON_cleanuppointers(p, inp, jsonname, error_cnt)

end subroutine JSONwriteMCCLNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONwriteMCCLMultiLayerNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist to JSON file
!
!> @param mcnl Monte Carlo name list structure
!> @param jsonname output file name
!> @param error_cnt total number of errors encountered by json routines
!
!> @date 08/11/15 MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONwriteMCCLMultiLayerNameList(mcnl, jsonname, error_cnt)

use ISO_C_BINDING

IMPLICIT NONE

type(MCCLMultiLayerNameListType),INTENT(INOUT)        :: mcnl
character(fnlen),INTENT(IN)                           :: jsonname
integer(kind=irg),INTENT(INOUT)                       :: error_cnt

type(json_value),pointer                              :: p, inp

integer(kind=irg),parameter                           :: n_int = 5, n_real = 9
integer(kind=irg)                                     :: io_int(n_int)
real(kind=dbl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
character(fnlen)                                      :: dataset, namelistname
character(fnlen,kind=c_char)                          :: line2(1)

! initialize the json state variables
namelistname = 'MCCLdata'
call JSON_initpointers(p, inp, jsonname, namelistname, error_cnt)

! write all the single integers
io_int = (/ mcnl%stdout, mcnl%numsx, mcnl%globalworkgrpsz, mcnl%num_el, mcnl%totnum_el /)
intlist(1) = 'stdout'
intlist(2) = 'numsx'
intlist(3) = 'globalworkgrpsz'
intlist(4) = 'num_el'
intlist(5) = 'totnum_el'
call JSON_writeNMLintegers(inp, io_int, intlist, n_int, error_cnt)

! write all the single doubles
io_real = (/ mcnl%sig, mcnl%omega, mcnl%EkeV, mcnl%Ehistmin, mcnl%Ebinsize, mcnl%depthmax, mcnl%depthstep, &
             mcnl%filmthickness, mcnl%filmstep /)
reallist(1) = 'sig'
reallist(2) = 'omega'
reallist(3) = 'EkeV'
reallist(4) = 'Ehistmin'
reallist(5) = 'Ebinsize'
reallist(6) = 'depthmax'
reallist(7) = 'depthstep'
reallist(8) = 'filmthickness'
reallist(9) = 'filmstep'
call JSON_writeNMLdoubles(inp, io_real, reallist, n_real, error_cnt)

! write all the strings
dataset = 'MCmode'
call json_add(inp, dataset, mcnl%MCmode); call JSON_failtest(error_cnt)

dataset = 'xtalname_film'
call json_add(inp, dataset, mcnl%xtalname_film); call JSON_failtest(error_cnt)

dataset = 'xtalname_subs'
call json_add(inp, dataset, mcnl%xtalname_subs); call JSON_failtest(error_cnt)

dataset = 'dataname'
call json_add(inp, dataset, mcnl%dataname); call JSON_failtest(error_cnt)

dataset = 'mode'
call json_add(inp, dataset, mcnl%mode); call JSON_failtest(error_cnt)

! and then we write the file and clean up
call JSON_cleanuppointers(p, inp, jsonname, error_cnt)

end subroutine JSONwriteMCCLMultiLayerNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONwriteEBSDMasterNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist to JSON file
!
!> @param emnl EBSD master name list structure
!> @param jsonname output file name
!> @param error_cnt total number of errors encountered by json routines
!
!> @date 08/12/15 MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONwriteEBSDMasterNameList(emnl, jsonname, error_cnt)

use ISO_C_BINDING

IMPLICIT NONE

type(EBSDMasterNameListType),INTENT(INOUT)            :: emnl
character(fnlen),INTENT(IN)                           :: jsonname
integer(kind=irg),INTENT(INOUT)                       :: error_cnt

type(json_value),pointer                              :: p, inp

integer(kind=irg),parameter                           :: n_int = 4, n_real = 1
integer(kind=irg)                                     :: io_int(n_int)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
character(fnlen)                                      :: dataset, namelistname
character(fnlen,kind=c_char)                          :: line2(1)

! initialize the json state variables
namelistname = 'EBSDmastervars'
call JSON_initpointers(p, inp, jsonname, namelistname, error_cnt)

! write all the single integers
io_int = (/ emnl%stdout, emnl%npx, emnl%Esel, emnl%nthreads /)
intlist(1) = 'stdout'
intlist(2) = 'npx'
intlist(3) = 'Esel'
intlist(4) = 'nthreads'
call JSON_writeNMLintegers(inp, io_int, intlist, n_int, error_cnt)

! write a single real
dataset = 'dmin'
call json_add(inp, dataset, dble(emnl%dmin)); call JSON_failtest(error_cnt)

! write all the strings
dataset = 'outname'
call json_add(inp, dataset, emnl%outname); call JSON_failtest(error_cnt)

dataset = 'energyfile'
call json_add(inp, dataset, emnl%energyfile); call JSON_failtest(error_cnt)

! and then we write the file and clean up
call JSON_cleanuppointers(p, inp, jsonname, error_cnt)

end subroutine JSONwriteEBSDMasterNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONwriteECPMasterNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist to JSON file
!
!> @param ecpnl ECP master name list structure
!> @param jsonname output file name
!> @param error_cnt total number of errors encountered by json routines
!
!> @date 08/12/15 MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONwriteECPMasterNameList(ecpnl, jsonname, error_cnt)

use ISO_C_BINDING

IMPLICIT NONE

type(ECPMasterNameListType),INTENT(INOUT)             :: ecpnl
character(fnlen),INTENT(IN)                           :: jsonname
integer(kind=irg),INTENT(INOUT)                       :: error_cnt

type(json_value),pointer                              :: p, inp

integer(kind=irg),parameter                           :: n_int = 4, n_real = 2
integer(kind=irg)                                     :: io_int(n_int)
real(kind=dbl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
character(fnlen)                                      :: dataset, namelistname
character(fnlen,kind=c_char)                          :: line2(1)

! initialize the json state variables
namelistname = 'ECPmastervars'
call JSON_initpointers(p, inp, jsonname, namelistname, error_cnt)

! write all the single integers
io_int = (/ ecpnl%stdout, ecpnl%Esel, ecpnl%npx, ecpnl%nthreads /)
intlist(1) = 'stdout'
intlist(2) = 'Esel'
intlist(3) = 'npx'
intlist(4) = 'nthreads'
call JSON_writeNMLintegers(inp, io_int, intlist, n_int, error_cnt)

dataset = 'distort'
call json_add(inp, dataset, ecpnl%distort); call JSON_failtest(error_cnt)

! integer vectors
dataset = 'fn'
call json_add(inp, dataset, dble(ecpnl%fn)); call JSON_failtest(error_cnt)

! write all the single doubles
io_real = (/ ecpnl%dmin, ecpnl%startthick /)
reallist(1) = 'dmin'
reallist(2) = 'startthick' 
call JSON_writeNMLdoubles(inp, io_real, reallist, n_real, error_cnt)

! 3-vectors (real)
dataset = 'abcdist'
call json_add(inp, dataset, dble(ecpnl%abcdist)); call JSON_failtest(error_cnt)

dataset = 'albegadist'
call json_add(inp, dataset, dble(ecpnl%albegadist)); call JSON_failtest(error_cnt)

! write all the strings
dataset = 'outname'
call json_add(inp, dataset, ecpnl%outname); call JSON_failtest(error_cnt)

dataset = 'energyfile'
call json_add(inp, dataset, ecpnl%energyfile); call JSON_failtest(error_cnt)

dataset = 'compmode'
call json_add(inp, dataset, ecpnl%compmode); call JSON_failtest(error_cnt)

! and then we write the file and clean up
call JSON_cleanuppointers(p, inp, jsonname, error_cnt)

end subroutine JSONwriteECPMasterNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONwriteEBSDNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist to JSON file
!
!> @param enl EBSD name list structure
!> @param jsonname output file name
!> @param error_cnt total number of errors encountered by json routines
!
!> @date 08/12/15 MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONwriteEBSDNameList(enl, jsonname, error_cnt)

use ISO_C_BINDING

IMPLICIT NONE

type(EBSDNameListType),INTENT(INOUT)                  :: enl
character(fnlen),INTENT(IN)                           :: jsonname
integer(kind=irg),INTENT(INOUT)                       :: error_cnt

type(json_value),pointer                              :: p, inp

integer(kind=irg),parameter                           :: n_int = 6, n_real = 8
integer(kind=irg)                                     :: io_int(n_int)
real(kind=sgl)                                        :: io_real(n_real)
real(kind=dbl)                                        :: t(1)
character(20)                                         :: intlist(n_int), reallist(n_real)
character(fnlen)                                      :: dataset, namelistname
character(fnlen,kind=c_char)                          :: line2(1)

! initialize the json state variables
namelistname = 'EBSDdata'
call JSON_initpointers(p, inp, jsonname, namelistname, error_cnt)

! write all the single integers
io_int = (/ enl%stdout, enl%numsx, enl%numsy, enl%binning, enl%nthreads, enl%energyaverage /)
intlist(1) = 'stdout'
intlist(2) = 'numsx'
intlist(3) = 'numsy'
intlist(4) = 'binning'
intlist(5) = 'nthreads'
intlist(6) = 'energyaverage'
call JSON_writeNMLintegers(inp, io_int, intlist, n_int, error_cnt)


! write all the single reals 
io_real = (/ enl%L, enl%thetac, enl%delta, enl%xpc, enl%ypc, enl%energymin, enl%energymax, enl%gammavalue /)
reallist(1) = 'L'
reallist(2) = 'thetac'
reallist(3) = 'delta'
reallist(4) = 'xpc'
reallist(5) = 'ypc'
reallist(6) = 'energymin'
reallist(7) = 'energymax'
reallist(8) = 'gammavalue'
call JSON_writeNMLreals(inp, io_real, reallist, n_real, error_cnt)

! a 4-vector
dataset = 'axisangle'
call json_add(inp, dataset, dble(enl%axisangle)); call JSON_failtest(error_cnt)

! a few doubles
dataset = 'beamcurrent'
call json_add(inp, dataset, enl%beamcurrent); call JSON_failtest(error_cnt)

dataset = 'dwelltime'
call json_add(inp, dataset, enl%dwelltime); call JSON_failtest(error_cnt)

! write all the strings
dataset = 'maskpattern'
call json_add(inp, dataset, enl%maskpattern); call JSON_failtest(error_cnt)

dataset = 'scalingmode'
call json_add(inp, dataset, enl%scalingmode); call JSON_failtest(error_cnt)

dataset = 'eulerconvention'
call json_add(inp, dataset, enl%eulerconvention); call JSON_failtest(error_cnt)

dataset = 'outputformat'
call json_add(inp, dataset, enl%outputformat); call JSON_failtest(error_cnt)

dataset = 'energyfile'
call json_add(inp, dataset, enl%energyfile); call JSON_failtest(error_cnt)

dataset = 'masterfile'
call json_add(inp, dataset, enl%masterfile); call JSON_failtest(error_cnt)

dataset = 'anglefile'
call json_add(inp, dataset, enl%anglefile); call JSON_failtest(error_cnt)

dataset = 'datafile'
call json_add(inp, dataset, enl%datafile); call JSON_failtest(error_cnt)

! and then we write the file and clean up
call JSON_cleanuppointers(p, inp, jsonname, error_cnt)

end subroutine JSONwriteEBSDNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONwriteECPNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist to JSON file
!
!> @param ecpnl ECP namelist structure
!> @param jsonname output file name
!> @param error_cnt total number of errors encountered by json routines
!
!> @date 08/12/15 MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONwriteECPNameList(ecpnl, jsonname, error_cnt)

use ISO_C_BINDING

IMPLICIT NONE

type(ECPNameListType),INTENT(INOUT)                   :: ecpnl
character(fnlen),INTENT(IN)                           :: jsonname
integer(kind=irg),INTENT(INOUT)                       :: error_cnt

type(json_value),pointer                              :: p, inp

integer(kind=irg),parameter                           :: n_int = 4, n_real = 8
integer(kind=irg)                                     :: io_int(n_int)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
character(fnlen)                                      :: dataset, namelistname
character(fnlen,kind=c_char)                          :: line2(1)

! initialize the json state variables
namelistname = 'ECPlist'
call JSON_initpointers(p, inp, jsonname, namelistname, error_cnt)

! write all the single integers
io_int = (/ ecpnl%stdout, ecpnl%numthick, ecpnl%npix, ecpnl%nthreads /)
intlist(1) = 'stdout'
intlist(2) = 'numthick'
intlist(3) = 'npix'
intlist(4) = 'nthreads'
call JSON_writeNMLintegers(inp, io_int, intlist, n_int, error_cnt)

! integer vectors
dataset = 'k'
call json_add(inp, dataset, ecpnl%k); call JSON_failtest(error_cnt)

dataset = 'fn'
call json_add(inp, dataset, ecpnl%fn); call JSON_failtest(error_cnt)

dataset = 'gF'
call json_add(inp, dataset, ecpnl%gF); call JSON_failtest(error_cnt)

dataset = 'gS'
call json_add(inp, dataset, ecpnl%gS); call JSON_failtest(error_cnt)

dataset = 'tF'
call json_add(inp, dataset, ecpnl%tF); call JSON_failtest(error_cnt)

dataset = 'tS'
call json_add(inp, dataset, ecpnl%tS); call JSON_failtest(error_cnt)

! write all the single reals
io_real = (/ ecpnl%voltage, ecpnl%dmin, ecpnl%ktmax, ecpnl%thetac, ecpnl%startthick, ecpnl%thickinc, ecpnl%zintstep, &
             ecpnl%filmthickness /)
reallist(1) = 'voltage'
reallist(2) = 'dmin'
reallist(3) = 'ktmax'
reallist(4) = 'thetac'
reallist(5) = 'startthick'
reallist(6) = 'thickinc'
reallist(7) = 'zintstep'
reallist(8) = 'filmthickness'
call JSON_writeNMLreals(inp, io_real, reallist, n_real, error_cnt)

! write all the strings
dataset = 'compmode'
call json_add(inp, dataset, ecpnl%compmode); call JSON_failtest(error_cnt)

dataset = 'energyfile'
call json_add(inp, dataset, ecpnl%energyfile); call JSON_failtest(error_cnt)

dataset = 'outname'
call json_add(inp, dataset, ecpnl%outname); call JSON_failtest(error_cnt)

dataset = 'xtalname'
call json_add(inp, dataset, ecpnl%xtalname); call JSON_failtest(error_cnt)

dataset = 'xtalname2'
call json_add(inp, dataset, ecpnl%xtalname2); call JSON_failtest(error_cnt)

! and then we write the file and clean up
call JSON_cleanuppointers(p, inp, jsonname, error_cnt)

end subroutine JSONwriteECPNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONwriteLACBEDNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist to JSON file
!
!> @param lacbednl LACBED name list structure
!> @param jsonname output file name
!> @param error_cnt total number of errors encountered by json routines
!
!> @date 08/12/15 MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONwriteLACBEDNameList(lacbednl, jsonname, error_cnt)

use ISO_C_BINDING

IMPLICIT NONE

type(LACBEDNameListType),INTENT(INOUT)                :: lacbednl
character(fnlen),INTENT(IN)                           :: jsonname
integer(kind=irg),INTENT(INOUT)                       :: error_cnt

type(json_value),pointer                              :: p, inp

integer(kind=irg),parameter                           :: n_int = 5, n_real = 6
integer(kind=irg)                                     :: io_int(n_int)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
character(fnlen)                                      :: dataset, namelistname
character(fnlen,kind=c_char)                          :: line2(1)

! initialize the json state variables
namelistname = 'inputlist'
call JSON_initpointers(p, inp, jsonname, namelistname, error_cnt)

! write all the single integers
io_int = (/ lacbednl%stdout, lacbednl%maxHOLZ, lacbednl%numthick, lacbednl%npix, lacbednl%nthreads /)
intlist(1) = 'stdout'
intlist(2) = 'maxHOLZ'
intlist(3) = 'numthick'
intlist(4) = 'npix'
intlist(5) = 'nthreads'
call JSON_writeNMLintegers(inp, io_int, intlist, n_int, error_cnt)

! vectors
dataset = 'k'
call json_add(inp, dataset, lacbednl%k); call JSON_failtest(error_cnt)

dataset = 'fn'
call json_add(inp, dataset, lacbednl%fn); call JSON_failtest(error_cnt)

! write all the single reals
io_real = (/ lacbednl%voltage, lacbednl%dmin, lacbednl%convergence, lacbednl%startthick, lacbednl%thickinc, lacbednl%minten/)
reallist(1) = 'voltage'
reallist(2) = 'dmin'
reallist(3) = 'convergence'
reallist(4) = 'startthick'
reallist(5) = 'thickinc'
reallist(6) = 'minten'
call JSON_writeNMLreals(inp, io_real, reallist, n_real, error_cnt)

! write all the strings
dataset = 'outname'
call json_add(inp, dataset, lacbednl%outname); call JSON_failtest(error_cnt)

dataset = 'xtalname'
call json_add(inp, dataset, lacbednl%xtalname); call JSON_failtest(error_cnt)

! and then we write the file and clean up
call JSON_cleanuppointers(p, inp, jsonname, error_cnt)

end subroutine JSONwriteLACBEDNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONwriteECPpatternNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist into JSON file
!
!> @param ecpnl ECP name list structure
!> @param jsonname output file name
!> @param error_cnt total number of errors encountered by json routines
!
!> @date 08/12/15 MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONwriteECPpatternNameList(ecpnl, jsonname, error_cnt)

use ISO_C_BINDING

IMPLICIT NONE

type(ECPpatternNameListType),INTENT(INOUT)            :: ecpnl
character(fnlen),INTENT(IN)                           :: jsonname
integer(kind=irg),INTENT(INOUT)                       :: error_cnt

type(json_value),pointer                              :: p, inp

integer(kind=irg),parameter                           :: n_int = 2, n_real = 6
integer(kind=irg)                                     :: io_int(n_int)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
character(fnlen)                                      :: dataset, namelistname
character(fnlen,kind=c_char)                          :: line2(1)

! initialize the json state variables
namelistname = 'ECPvars'
call JSON_initpointers(p, inp, jsonname, namelistname, error_cnt)

! write all the single integers
io_int = (/ ecpnl%stdout, ecpnl%npix /)
intlist(1) = 'stdout'
intlist(2) = 'npix'
call JSON_writeNMLintegers(inp, io_int, intlist, n_int, error_cnt)

! single real
dataset = 'thetac'
call json_add(inp, dataset, dble(ecpnl%thetac)); call JSON_failtest(error_cnt)

! real vector
dataset = 'k'
call json_add(inp, dataset, dble(ecpnl%k)); call JSON_failtest(error_cnt)

! write all the strings
dataset = 'outname'
call json_add(inp, dataset, ecpnl%outname); call JSON_failtest(error_cnt)

dataset = 'masterfile'
call json_add(inp, dataset, ecpnl%masterfile); call JSON_failtest(error_cnt)

! and then we write the file and clean up
call JSON_cleanuppointers(p, inp, jsonname, error_cnt)

end subroutine JSONwriteECPpatternNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONwritePEDKINNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist to JSON file
!
!> @param pednl PED name list structure
!> @param jsonname output file name
!> @param error_cnt total number of errors encountered by json routines
!
!> @date 08/12/15 MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONwritePEDKINNameList(pednl, jsonname, error_cnt)

use ISO_C_BINDING

IMPLICIT NONE

type(PEDKINNameListType),INTENT(INOUT)                :: pednl
character(fnlen),INTENT(IN)                           :: jsonname
integer(kind=irg),INTENT(INOUT)                       :: error_cnt

type(json_value),pointer                              :: p, inp

integer(kind=irg),parameter                           :: n_int = 4, n_real = 4
integer(kind=irg)                                     :: io_int(n_int)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
character(fnlen)                                      :: dataset, namelistname
character(fnlen,kind=c_char)                          :: line2(1)

! initialize the json state variables
namelistname = 'inputlist'
call JSON_initpointers(p, inp, jsonname, namelistname, error_cnt)

! write all the single integers
io_int = (/ pednl%stdout, pednl%npix, pednl%ncubochoric, pednl%nthreads /)
intlist(1) = 'stdout'
intlist(2) = 'npix'
intlist(3) = 'ncubochoric'
intlist(4) = 'nthreads'
call JSON_writeNMLintegers(inp, io_int, intlist, n_int, error_cnt)

! write all the single reals
io_real = (/ pednl%voltage, pednl%thickness, pednl%dmin, pednl%rnmpp /)
reallist(1) = 'voltage'
reallist(2) = 'thickness'
reallist(3) = 'dmin'
reallist(4) = 'rnmpp'
call JSON_writeNMLreals(inp, io_real, reallist, n_real, error_cnt)

! write all the strings
dataset = 'outname'
call json_add(inp, dataset, pednl%outname); call JSON_failtest(error_cnt)

dataset = 'xtalname'
call json_add(inp, dataset, pednl%xtalname); call JSON_failtest(error_cnt)

dataset = 'eulername'
call json_add(inp, dataset, pednl%eulername); call JSON_failtest(error_cnt)

! and then we write the file and clean up
call JSON_cleanuppointers(p, inp, jsonname, error_cnt)

end subroutine JSONwritePEDKINNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONwritePEDNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist to JSON file
!
!> @param pednl PED name list structure
!> @param jsonname output file name
!> @param error_cnt total number of errors encountered by json routines
!
!> @date 08/12/15 MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONwritePEDNameList(pednl, jsonname, error_cnt)

use ISO_C_BINDING

IMPLICIT NONE

type(PEDNameListType),INTENT(INOUT)                   :: pednl
character(fnlen),INTENT(IN)                           :: jsonname
integer(kind=irg),INTENT(INOUT)                       :: error_cnt

type(json_value),pointer                              :: p, inp

integer(kind=irg),parameter                           :: n_int = 5, n_real = 6
integer(kind=irg)                                     :: io_int(n_int)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
character(fnlen)                                      :: dataset, namelistname
character(fnlen,kind=c_char)                          :: line2(1)

! initialize the json state variables
namelistname = 'inputlist'
call JSON_initpointers(p, inp, jsonname, namelistname, error_cnt)

! write all the single integers
io_int = (/ pednl%stdout, pednl%precsample, pednl%precazimuthal, pednl%npix, pednl%nthreads /)
intlist(1) = 'stdout'
intlist(2) = 'precsample'
intlist(3) = 'precazimuthal'
intlist(4) = 'npix'
intlist(5) = 'nthreads'
call JSON_writeNMLintegers(inp, io_int, intlist, n_int, error_cnt)

! vectors
dataset = 'k'
call json_add(inp, dataset, pednl%k); call JSON_failtest(error_cnt)

dataset = 'fn'
call json_add(inp, dataset, pednl%fn); call JSON_failtest(error_cnt)

! single reals
io_real = (/ pednl%voltage, pednl%dmin, pednl%precangle, pednl%prechalfwidth, pednl%thickness, pednl%camlen /)
reallist(1) = 'voltage'
reallist(2) = 'dmin'
reallist(3) = 'precangle'
reallist(4) = 'prechalfwidth'
reallist(5) = 'thickness'
reallist(6) = 'camlen'
call JSON_writeNMLreals(inp, io_real, reallist, n_real, error_cnt)


! write all the strings
dataset = 'outname'
call json_add(inp, dataset, pednl%outname); call JSON_failtest(error_cnt)

dataset = 'xtalname'
call json_add(inp, dataset, pednl%xtalname); call JSON_failtest(error_cnt)

dataset = 'filemode'
call json_add(inp, dataset, pednl%filemode); call JSON_failtest(error_cnt)

! and then we write the file and clean up
call JSON_cleanuppointers(p, inp, jsonname, error_cnt)

end subroutine JSONwritePEDNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONwriteECCINameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist to JSON file
!
!> @param eccinl ECCI name list structure
!> @param jsonname output file name
!> @param error_cnt total number of errors encountered by json routines
!
!> @date 08/12/15 MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONwriteECCINameList(eccinl, jsonname, error_cnt)

use ISO_C_BINDING

IMPLICIT NONE

type(ECCINameListType),INTENT(INOUT)                  :: eccinl
character(fnlen),INTENT(IN)                           :: jsonname
integer(kind=irg),INTENT(INOUT)                       :: error_cnt

type(json_value),pointer                              :: p, inp

integer(kind=irg),parameter                           :: n_int = 8, n_real = 6
integer(kind=irg)                                     :: io_int(n_int)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
integer(kind=irg)                                     :: i
character(fnlen)                                      :: dataset, namelistname
character(fnlen,kind=c_char)                          :: line2(1)

! initialize the json state variables
namelistname = 'ECCIlist'
call JSON_initpointers(p, inp, jsonname, namelistname, error_cnt)

! write all the single integers
io_int = (/ eccinl%stdout, eccinl%nthreads, eccinl%nktstep, eccinl%DF_npix, eccinl%DF_npiy, eccinl%numYdisl, eccinl%numdisl, &
           eccinl%numsf /)
intlist(1) = 'stdout'
intlist(2) = 'nthreads'
intlist(3) = 'nktstep'
intlist(4) = 'DF_npix'
intlist(5) = 'DF_npiy'
intlist(6) = 'numYdisl'
intlist(7) = 'numdisl'
intlist(8) = 'numsf'
call JSON_writeNMLintegers(inp, io_int, intlist, n_int, error_cnt)

! vectors
dataset = 'k'
call json_add(inp, dataset, eccinl%k); call JSON_failtest(error_cnt)

! single reals
io_real = (/ eccinl%voltage, eccinl%dkt, eccinl%ktmax, eccinl%dmin, eccinl%DF_L, eccinl%DF_slice /)
reallist(1) = 'voltage'
reallist(2) = 'dkt'
reallist(3) = 'ktmax'
reallist(4) = 'dmin'
reallist(5) = 'DF_L'
reallist(6) = 'DF_slice'
call JSON_writeNMLreals(inp, io_real, reallist, n_real, error_cnt)


! 2-vectors
dataset = 'lauec'
call json_add(inp, dataset, dble(eccinl%lauec)); call JSON_failtest(error_cnt)

dataset = 'lauec2'
call json_add(inp, dataset, dble(eccinl%lauec2)); call JSON_failtest(error_cnt)

! write all the strings
dataset = 'dispmode'
call json_add(inp, dataset, eccinl%dispmode); call JSON_failtest(error_cnt)

dataset = 'summode'
call json_add(inp, dataset, eccinl%summode); call JSON_failtest(error_cnt)

dataset = 'progmode'
call json_add(inp, dataset, eccinl%progmode); call JSON_failtest(error_cnt)

dataset = 'xtalname'
call json_add(inp, dataset, eccinl%xtalname); call JSON_failtest(error_cnt)

dataset = 'foilnmlfile'
call json_add(inp, dataset, eccinl%foilnmlfile); call JSON_failtest(error_cnt)

dataset = 'dispfile'
call json_add(inp, dataset, eccinl%dispfile); call JSON_failtest(error_cnt)

dataset = 'dataname'
call json_add(inp, dataset, eccinl%dataname); call JSON_failtest(error_cnt)

dataset = 'ECPname'
call json_add(inp, dataset, eccinl%ECPname); call JSON_failtest(error_cnt)

dataset = 'sgname'
call json_add(inp, dataset, eccinl%sgname); call JSON_failtest(error_cnt)

dataset = 'apbname'
call json_add(inp, dataset, eccinl%apbname); call JSON_failtest(error_cnt)

dataset = 'incname'
call json_add(inp, dataset, eccinl%dispmode); call JSON_failtest(error_cnt)

dataset = 'voidname'
call json_add(inp, dataset, eccinl%voidname); call JSON_failtest(error_cnt)


! maxdefects string arrays
dataset = 'sfname'
call json_add(inp, dataset, eccinl%sfname(1:eccinl%numsf)); call JSON_failtest(error_cnt)

! 3*maxdefects string arrays
dataset = 'dislYname'
call json_add(inp, dataset, eccinl%dislYname(1:eccinl%numYdisl)); call JSON_failtest(error_cnt)

dataset = 'dislname'
call json_add(inp, dataset, eccinl%dislname(1:eccinl%numdisl)); call JSON_failtest(error_cnt)

! and then we write the file and clean up
call JSON_cleanuppointers(p, inp, jsonname, error_cnt)


! IN PRINCIPLE, HERE WE SHOULD ALSO READ ALL OF THE DEFECT FILES AND INSERT THEM INTO THE JSON FILE

! TO BE WRITTEN

end subroutine JSONwriteECCINameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONwriteRFZNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist to an JSON file
!
!> @param rfznl RFZ name list structure
!> @param jsonname output file name
!> @param error_cnt total number of errors encountered by json routines
!
!> @date 08/12/15 MDG 1.0 new routine
!> @date 08/18/15 MDG 1.1 added other rotation representations
!--------------------------------------------------------------------------
subroutine JSONwriteRFZNameList(rfznl, jsonname, error_cnt)

use ISO_C_BINDING

IMPLICIT NONE

type(RFZNameListType),INTENT(INOUT)                   :: rfznl
character(fnlen),INTENT(IN)                           :: jsonname
integer(kind=irg),INTENT(INOUT)                       :: error_cnt

type(json_value),pointer                              :: p, inp

integer(kind=irg),parameter                           :: n_int = 2, n_real = 1
integer(kind=irg)                                     :: io_int(n_int)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
integer(kind=irg)                                     :: i
character(fnlen)                                      :: dataset, namelistname
character(fnlen,kind=c_char)                          :: line2(1)

! initialize the json state variables
namelistname = 'RFZlist'
call JSON_initpointers(p, inp, jsonname, namelistname, error_cnt)

! write all the single integers
io_int = (/ rfznl%pgnum, rfznl%nsteps /)
intlist(1) = 'pgnum'
intlist(2) = 'nsteps'
call JSON_writeNMLintegers(inp, io_int, intlist, n_int, error_cnt)

! strings
dataset = 'euoutname'
call json_add(inp, dataset, rfznl%euoutname); call JSON_failtest(error_cnt)

dataset = 'cuoutname'
call json_add(inp, dataset, rfznl%cuoutname); call JSON_failtest(error_cnt)

dataset = 'hooutname'
call json_add(inp, dataset, rfznl%hooutname); call JSON_failtest(error_cnt)

dataset = 'quoutname'
call json_add(inp, dataset, rfznl%quoutname); call JSON_failtest(error_cnt)

dataset = 'rooutname'
call json_add(inp, dataset, rfznl%rooutname); call JSON_failtest(error_cnt)

dataset = 'omoutname'
call json_add(inp, dataset, rfznl%omoutname); call JSON_failtest(error_cnt)

dataset = 'axoutname'
call json_add(inp, dataset, rfznl%axoutname); call JSON_failtest(error_cnt)

! and then we write the file and clean up
call JSON_cleanuppointers(p, inp, jsonname, error_cnt)

end subroutine JSONwriteRFZNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONwriteDictIndxOpenCLNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist to JSON file
!
!> @param rfznl RFZ name list structure
!> @param jsonname output file name
!> @param error_cnt total number of errors encountered by json routines
!
!> @date 08/12/15 MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONwriteDictIndxOpenCLNameList(dictindxnl, jsonname, error_cnt)

use ISO_C_BINDING

use local

IMPLICIT NONE

type(DictIndxOpenCLListType),INTENT(INOUT)            :: dictindxnl
character(fnlen),INTENT(IN)                           :: jsonname
integer(kind=irg),INTENT(INOUT)                       :: error_cnt

type(json_value),pointer                              :: p, inp

integer(kind=irg),parameter                           :: n_int = 7, n_real = 1
integer(kind=irg)                                     :: io_int(n_int)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
integer(kind=irg)                                     :: i
character(fnlen)                                      :: dataset, namelistname
character(fnlen,kind=c_char)                          :: line2(1)

! initialize the json state variables
namelistname = 'DictIndxOpenCLvars'
call JSON_initpointers(p, inp, jsonname, namelistname, error_cnt)

! write all the single integers
io_int = (/ dictindxnl%numexptsingle, dictindxnl%numdictsingle, dictindxnl%totnumdict, dictindxnl%totnumexpt, dictindxnl%imght, &
            dictindxnl%imgwd, dictindxnl%nnk /)
intlist(1) = 'numexptsingle'
intlist(2) = 'numdictsingle'
intlist(3) = 'totnumdict'
intlist(4) = 'totnumexpt'
intlist(5) = 'imght'
intlist(6) = 'imgwd'
intlist(7) = 'nnk'
call JSON_writeNMLintegers(inp, io_int, intlist, n_int, error_cnt)

dataset = 'MeanSubtraction'
call json_add(inp, dataset, dictindxnl%MeanSubtraction); call JSON_failtest(error_cnt)

! strings
dataset = 'exptfile'
call json_add(inp, dataset, dictindxnl%exptfile); call JSON_failtest(error_cnt)

dataset = 'dictfile'
call json_add(inp, dataset, dictindxnl%dictfile); call JSON_failtest(error_cnt)

dataset = 'eulerfile'
call json_add(inp, dataset, dictindxnl%eulerfile); call JSON_failtest(error_cnt)

! and then we write the file and clean up
call JSON_cleanuppointers(p, inp, jsonname, error_cnt)

end subroutine JSONwriteDictIndxOpenCLNameList


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! here we have the json->namelist conversion routines for all the 
! namelists defined in the NameListTypedefs.f90 file.
!
! To convert a json file to a namelist, we first initialize the namelist 
! to the default values, to make sure that any omissions in the json file
! are intercepted.  To do so, we call the NameListHandler routines with 
! the "initonly" optional keyword.
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------


!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONreadInteger
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read integer from json file into namelist structure (with auto missing detection)
!
!> @param json structure
!> @param ep entry path string
!> @param ival integer variable
!> @param dval integer variable (default value)
!
!> @date 08/12/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONreadInteger(json, ep, ival, dval)

use ISO_C_BINDING
use io

IMPLICIT NONE

type(json_file),INTENT(INOUT)           :: json
character(fnlen),INTENT(IN)             :: ep
integer(kind=irg),INTENT(INOUT)         :: ival
integer(kind=irg),INTENT(IN)            :: dval

logical                                 :: found

! if we find the field 'ep' in the file, then we read its corresponding value
! if it is not there, then we return the dval default value
call json%get(ep, ival, found)
if (.not. found) then
  write(error_unit,'(A)') 'WARNING: field '//trim(ep)//' not found in json file; using default value from namelist template'
  ival = dval
end if

end subroutine JSONreadInteger

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONreadIntegerVec
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read integer vector from json file into namelist structure (with auto missing detection)
!
!> @param json structure
!> @param ep entry path string
!> @param ival integer vector variable
!> @param dval integer vector variable (default value)
!
!> @date 08/12/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONreadIntegerVec(json, ep, ivec, dvec, n)

use ISO_C_BINDING
use io

IMPLICIT NONE

type(json_file),INTENT(INOUT)           :: json
character(fnlen),INTENT(IN)             :: ep
integer(kind=irg),INTENT(INOUT)         :: ivec(n)
integer(kind=irg),INTENT(IN)            :: dvec(n)
integer(kind=irg),INTENT(IN)            :: n

logical                                 :: found
integer(kind=irg),dimension(:),allocatable :: rv

! if we find the field 'ep' in the file, then we read its corresponding value
! if it is not there, then we return the dval default value
call json%get(ep, rv, found)
if (.not. found) then
  write(error_unit,'(A)') 'WARNING: field '//trim(ep)//' not found in json file; using default value from namelist template'
  ivec = dvec
else
  ivec = rv
end if

end subroutine JSONreadIntegerVec

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONreadReal
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read single precision real from json file into namelist structure (with auto missing detection)
!
!> @param json structure
!> @param ep entry path string
!> @param rval real variable
!> @param dval real variable (default value)
!
!> @date 08/12/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONreadReal(json, ep, rval, dval)

use ISO_C_BINDING
use io

IMPLICIT NONE

type(json_file),INTENT(INOUT)           :: json
character(fnlen),INTENT(IN)             :: ep
real(kind=sgl),INTENT(INOUT)            :: rval
real(kind=sgl),INTENT(IN)               :: dval

logical                                 :: found
real(kind=dbl)                          :: rv

! if we find the field 'ep' in the file, then we read its corresponding value
! if it is not there, then we return the dval default value
call json%get(ep, rv, found)
if (.not. found) then
  write(error_unit,'(A)') 'WARNING: field '//trim(ep)//' not found in json file; using default value from namelist template'
  rval = dval
else
  rval = sngl(rv)
end if

end subroutine JSONreadReal

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONreadRealVec
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read single precision real vector from json file into namelist structure (with auto missing detection)
!
!> @param json structure
!> @param ep entry path string
!> @param rval real vector variable
!> @param dval real vector variable (default value)
!
!> @date 08/12/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONreadRealVec(json, ep, rvec, dvec, n)

use ISO_C_BINDING
use io

IMPLICIT NONE

type(json_file),INTENT(INOUT)           :: json
character(fnlen),INTENT(IN)             :: ep
real(kind=sgl),INTENT(INOUT)            :: rvec(n)
real(kind=sgl),INTENT(IN)               :: dvec(n)
integer(kind=irg),INTENT(IN)            :: n

logical                                 :: found
real(kind=dbl),dimension(:),allocatable :: rv

! if we find the field 'ep' in the file, then we read its corresponding value
! if it is not there, then we return the dvec default value
call json%get(ep, rv, found)
if (.not. found) then
  write(error_unit,'(A)') 'WARNING: field '//trim(ep)//' not found in json file; using default value from namelist template'
  rvec = dvec
else
  rvec = sngl(rv)
end if

end subroutine JSONreadRealVec

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONreadDouble
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read double precision real from json file into namelist structure (with auto missing detection)
!
!> @param json structure
!> @param ep entry path string
!> @param rval real variable
!> @param dval real variable (default value)
!
!> @date 08/12/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONreadDouble(json, ep, rval, dval)

use ISO_C_BINDING
use io

IMPLICIT NONE

type(json_file),INTENT(INOUT)           :: json
character(fnlen),INTENT(IN)             :: ep
real(kind=dbl),INTENT(INOUT)            :: rval
real(kind=dbl),INTENT(IN)               :: dval

logical                                 :: found
real(kind=dbl)                          :: rv

! if we find the field 'ep' in the file, then we read its corresponding value
! if it is not there, then we return the dval default value
call json%get(ep, rv, found)
if (.not. found) then
  write(error_unit,'(A)') 'WARNING: field '//trim(ep)//' not found in json file; using default value from namelist template'
  rval = dval
else
  rval = rv
end if

end subroutine JSONreadDouble

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONreadDoubleVec
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read double precision real vector from json file into namelist structure (with auto missing detection)
!
!> @param json structure
!> @param ep entry path string
!> @param rval real vector variable
!> @param dval real vector variable (default value)
!
!> @date 08/12/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONreadDoubleVec(json, ep, rvec, dvec, n)

use ISO_C_BINDING
use io

IMPLICIT NONE

type(json_file),INTENT(INOUT)           :: json
character(fnlen),INTENT(IN)             :: ep
real(kind=dbl),INTENT(INOUT)            :: rvec(n)
real(kind=dbl),INTENT(IN)               :: dvec(n)
integer(kind=irg),INTENT(IN)            :: n

logical                                 :: found
real(kind=dbl),dimension(:),allocatable :: rv

! if we find the field 'ep' in the file, then we read its corresponding value
! if it is not there, then we return the dvec default value
call json%get(ep, rv, found)
if (.not. found) then
  write(error_unit,'(A)') 'WARNING: field '//trim(ep)//' not found in json file; using default value from namelist template'
  rvec = dvec
else
  rvec = rv
end if

end subroutine JSONreadDoubleVec

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONreadString
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read string from json file into namelist structure (with auto missing detection)
!
!> @param json structure
!> @param ep entry path string
!> @param sval string
!> @param dval string (default value)
!
!> @date 08/12/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONreadString(json, ep, sval, dval)

use ISO_C_BINDING
use io

IMPLICIT NONE

type(json_file),INTENT(INOUT)           :: json
character(fnlen),INTENT(IN)             :: ep
character(fnlen),INTENT(INOUT)          :: sval
character(fnlen),INTENT(IN)             :: dval

logical                                 :: found
character(kind=CK,len=:),allocatable    :: cval

! if we find the field 'ep' in the file, then we read its corresponding value
! if it is not there, then we return the dval default value
call json%get(ep, cval, found)
if (.not. found) then
  write(error_unit,'(A)') 'WARNING: field '//trim(ep)//' not found in json file; using default value from namelist template'
  sval = dval
else
  sval = trim(cval)
end if

end subroutine JSONreadString

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONreadLogical
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read logical from json file into namelist structure (with auto missing detection)
!
!> @param json structure
!> @param ep entry path string
!> @param sval logical 
!> @param dval logical (default value)
!
!> @date 08/20/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONreadLogical(json, ep, sval, dval)

use ISO_C_BINDING
use io

IMPLICIT NONE

type(json_file),INTENT(INOUT)           :: json
character(fnlen),INTENT(IN)             :: ep
logical,INTENT(INOUT)                   :: sval
logical,INTENT(IN)                      :: dval

logical                                 :: found, cval

! if we find the field 'ep' in the file, then we read its corresponding value
! if it is not there, then we return the dval default value
call json%get(ep, cval, found)
if (.not. found) then
  write(error_unit,'(A)') 'WARNING: field '//trim(ep)//' not found in json file; using default value from namelist template'
  sval = dval
else
  sval = cval
end if

end subroutine JSONreadLogical



!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONreadKosselNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read json file into namelist structure
!
!> @param knl Kossel name list structure
!> @param jsonname input file name
!> @param error_cnt total number of errors encountered by json routines
!
!> @date 08/12/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONreadKosselNameList(knl, jsonname, error_cnt)

use ISO_C_BINDING
use NameListHandlers

IMPLICIT NONE

type(KosselNameListType),INTENT(INOUT)                :: knl
character(fnlen),INTENT(IN)                           :: jsonname
integer(kind=irg),INTENT(INOUT)                       :: error_cnt

type(json_file)                                       :: json    !the JSON structure read from the file:

type(KosselNameListType)                              :: defknl
logical                                               :: init = .TRUE.
character(fnlen)                                      :: nmlfile = '', ep
real(kind=wp)                                         :: rval
character(kind=CK,len=:),allocatable                  :: cval
real(wp),dimension(:),allocatable                     :: rvec

! first of all, open the file and return an error message if it does not exist
error_cnt = 0
call json_initialize(); call JSON_failtest(error_cnt)

! populate the json structure
call json%load_file(filename = trim(jsonname))
if (json_failed()) then    !if there was an error reading the file
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
else
! ok, we got here so we need to initialize the namelist first to its default values (set in NameListHandlers)
  call GetKosselNameList(nmlfile, defknl, initonly=init)

! then we start reading the values in the json file  
  ep = 'Kossellist.stdout'
  call JSONreadInteger(json, ep, knl%stdout, defknl%stdout)
  ep = 'Kossellist.numthick'
  call JSONreadInteger(json, ep, knl%numthick, defknl%numthick)
  ep = 'Kossellist.npix'
  call JSONreadInteger(json, ep, knl%npix, defknl%npix)
  ep = 'Kossellist.maxHOLZ'
  call JSONreadInteger(json, ep, knl%maxHOLZ, defknl%maxHOLZ)
  ep = 'Kossellist.nthreads'
  call JSONreadInteger(json, ep, knl%nthreads, defknl%nthreads)

  ep = 'Kossellist.k'
  call JSONreadIntegerVec(json, ep, knl%k, defknl%k, size(knl%k))
  ep = 'Kossellist.fn'
  call JSONreadIntegerVec(json, ep, knl%fn, defknl%fn, size(knl%fn))

  ep = 'Kossellist.voltage'
  call JSONreadReal(json, ep, knl%voltage, defknl%voltage)
  ep = 'Kossellist.dmin'
  call JSONreadReal(json, ep, knl%dmin, defknl%dmin)
  ep = 'Kossellist.convergence'
  call JSONreadReal(json, ep, knl%convergence, defknl%convergence)
  ep = 'Kossellist.startthick'
  call JSONreadReal(json, ep, knl%startthick, defknl%startthick)
  ep = 'Kossellist.thickinc'
  call JSONreadReal(json, ep, knl%thickinc, defknl%thickinc)
  ep = 'Kossellist.minten'
  call JSONreadReal(json, ep, knl%minten, defknl%minten)

  ep = 'Kossellist.xtalname'
  call JSONreadString(json, ep, knl%xtalname, defknl%xtalname)
  ep = 'Kossellist.outname'
  call JSONreadString(json, ep, knl%outname, defknl%outname)
end if

call json%destroy(); call JSON_failtest(error_cnt)

end subroutine JSONreadKosselNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONreadKosselMasterNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read json file and fill knl structure (used by EMKosselmaster.f90)
!
!> @param knl Kossel name list structure
!> @param jsonname input file name
!> @param error_cnt total number of errors encountered by json routines
!
!> @date 08/19/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONreadKosselMasterNameList(knl, jsonname, error_cnt)

use ISO_C_BINDING
use NameListHandlers

IMPLICIT NONE

type(KosselMasterNameListType),INTENT(INOUT)          :: knl
character(fnlen),INTENT(IN)                           :: jsonname
integer(kind=irg),INTENT(INOUT)                       :: error_cnt

type(json_file)                                       :: json    !the JSON structure read from the file:

type(KosselMasterNameListType)                        :: defknl
logical                                               :: init = .TRUE.
character(fnlen)                                      :: nmlfile = '', ep, s, s2
real(kind=wp)                                         :: rval
character(kind=CK,len=:),allocatable                  :: cval
real(wp),dimension(:),allocatable                     :: rvec

! first of all, open the file and return an error message if it does not exist
error_cnt = 0
call json_initialize(); call JSON_failtest(error_cnt)

! populate the json structure
call json%load_file(filename = trim(jsonname))
if (json_failed()) then    !if there was an error reading the file
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
else
! ok, we got here so we need to initialize the namelist first to its default values (set in NameListHandlers)
  call GetKosselMasterNameList(nmlfile, defknl, initonly=init)

! then we start reading the values in the json file  
  ep = 'Kosselmasterlist.stdout'
  call JSONreadInteger(json, ep, knl%stdout, defknl%stdout)

  ep = 'Kosselmasterlist.numthick'
  call JSONreadInteger(json, ep, knl%numthick, defknl%numthick)
  ep = 'Kosselmasterlist.npix'
  call JSONreadInteger(json, ep, knl%npix, defknl%npix)
  ep = 'Kosselmasterlist.nthreads'
  call JSONreadInteger(json, ep, knl%nthreads, defknl%nthreads)

  ep = 'Kosselmasterlist.voltage'
  call JSONreadReal(json, ep, knl%voltage, defknl%voltage)
  ep = 'Kosselmasterlist.dmin'
  call JSONreadReal(json, ep, knl%dmin, defknl%dmin)
  ep = 'Kosselmasterlist.startthick'
  call JSONreadReal(json, ep, knl%startthick, defknl%startthick)
  ep = 'Kosselmasterlist.thickinc'
  call JSONreadReal(json, ep, knl%thickinc, defknl%thickinc)
  ep = 'Kosselmasterlist.tfraction'
  call JSONreadReal(json, ep, knl%tfraction, defknl%tfraction)

  ep = 'Kosselmasterlist.Kosselmode'
  s = knl%Kosselmode
  s2 = defknl%Kosselmode
  call JSONreadString(json, ep, s, s2)
  ep = 'Kosselmasterlist.xtalname'
  call JSONreadString(json, ep, knl%xtalname, defknl%xtalname)
  ep = 'Kosselmasterlist.outname'
  call JSONreadString(json, ep, knl%outname, defknl%outname)
end if

call json%destroy(); call JSON_failtest(error_cnt)

end subroutine JSONreadKosselMasterNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONreadMCNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read json file and fill mcnl structure (used by EMMC.f90)
!
!> @param mcnl Monte Carloname list structure
!> @param jsonname input file name
!> @param error_cnt total number of errors encountered by json routines
!
!> @date 08/19/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONreadMCNameList(mcnl, jsonname, error_cnt)

use ISO_C_BINDING
use NameListHandlers

IMPLICIT NONE

type(MCNameListType),INTENT(INOUT)                    :: mcnl
character(fnlen),INTENT(IN)                           :: jsonname
integer(kind=irg),INTENT(INOUT)                       :: error_cnt

type(json_file)                                       :: json    !the JSON structure read from the file:

type(MCNameListType)                                  :: defmcnl
logical                                               :: init = .TRUE.
character(fnlen)                                      :: nmlfile = '', ep, s, s2
real(kind=wp)                                         :: rval
character(kind=CK,len=:),allocatable                  :: cval
real(wp),dimension(:),allocatable                     :: rvec

! first of all, open the file and return an error message if it does not exist
error_cnt = 0
call json_initialize(); call JSON_failtest(error_cnt)

! populate the json structure
call json%load_file(filename = trim(jsonname))
if (json_failed()) then    !if there was an error reading the file
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
else
! ok, we got here so we need to initialize the namelist first to its default values (set in NameListHandlers)
  call GetMCNameList(nmlfile, defmcnl, initonly=init)

! then we start reading the values in the json file  
  ep = 'MCdata.stdout'
  call JSONreadInteger(json, ep, mcnl%stdout, defmcnl%stdout)
  ep = 'MCdata.numsx'
  call JSONreadInteger(json, ep, mcnl%numsx, defmcnl%numsx)
  ep = 'MCdata.num_el'
  call JSONreadInteger(json, ep, mcnl%num_el, defmcnl%num_el)
  ep = 'MCdata.primeseeds'
  call JSONreadInteger(json, ep, mcnl%primeseed, defmcnl%primeseed)
  ep = 'MCdata.nthreads'
  call JSONreadInteger(json, ep, mcnl%nthreads, defmcnl%nthreads)

  ep = 'MCdata.sig'
  call JSONreadDouble(json, ep, mcnl%sig, defmcnl%sig)
  ep = 'MCdata.omega'
  call JSONreadDouble(json, ep, mcnl%omega, defmcnl%omega)
  ep = 'MCdata.EkeV'
  call JSONreadDouble(json, ep, mcnl%EkeV, defmcnl%EkeV)
  ep = 'MCdata.Ehistmin'
  call JSONreadDouble(json, ep, mcnl%Ehistmin, defmcnl%Ehistmin)
  ep = 'MCdata.Ebinsize'
  call JSONreadDouble(json, ep, mcnl%Ebinsize, defmcnl%Ebinsize)
  ep = 'MCdata.depthmax'
  call JSONreadDouble(json, ep, mcnl%depthmax, defmcnl%depthmax)
  ep = 'MCdata.depthstep'
  call JSONreadDouble(json, ep, mcnl%depthstep, defmcnl%depthstep)

  ep = 'MCdata.MCmode'
  s = mcnl%MCmode
  s2 = defmcnl%MCmode
  call JSONreadString(json, ep, s, s2)
  ep = 'MCdata.xtalname'
  call JSONreadString(json, ep, mcnl%xtalname, defmcnl%xtalname)
  ep = 'MCdata.dataname'
  call JSONreadString(json, ep, mcnl%dataname, defmcnl%dataname)
end if

call json%destroy(); call JSON_failtest(error_cnt)

end subroutine JSONreadMCNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONreadMCCLNameList
!
!> @author Saransh Singh/Marc De Graef, Carnegie Mellon University
!
!> @brief read namelist file and fill mcnl structure (used by EMMCCL.f90)
!
!> @param mcnl Monte Carlo name list structure
!> @param jsonname input file name
!> @param error_cnt total number of errors encountered by json routines
!
!> @date 08/19/15  MDG 1.0 new routine
!> @date 09/09/15  MDG 1.1 added devid
!--------------------------------------------------------------------------
subroutine JSONreadMCCLNameList(mcnl, jsonname, error_cnt)

use ISO_C_BINDING
use NameListHandlers

IMPLICIT NONE

type(MCCLNameListType),INTENT(INOUT)                  :: mcnl
character(fnlen),INTENT(IN)                           :: jsonname
integer(kind=irg),INTENT(INOUT)                       :: error_cnt

type(json_file)                                       :: json    !the JSON structure read from the file:

type(MCCLNameListType)                                :: defmcnl
logical                                               :: init = .TRUE.
character(fnlen)                                      :: nmlfile = '', ep, s, s2
real(kind=wp)                                         :: rval
character(kind=CK,len=:),allocatable                  :: cval
real(wp),dimension(:),allocatable                     :: rvec

! first of all, open the file and return an error message if it does not exist
error_cnt = 0
call json_initialize(); call JSON_failtest(error_cnt)

! populate the json structure
call json%load_file(filename = trim(jsonname))
if (json_failed()) then    !if there was an error reading the file
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
else
! ok, we got here so we need to initialize the namelist first to its default values (set in NameListHandlers)
  call GetMCCLNameList(nmlfile, defmcnl, initonly=init)

! then we start reading the values in the json file  
  ep = 'MCCLdata.stdout'
  call JSONreadInteger(json, ep, mcnl%stdout, defmcnl%stdout)
  ep = 'MCCLdata.numsx'
  call JSONreadInteger(json, ep, mcnl%numsx, defmcnl%numsx)
  ep = 'MCCLdata.globalworkgrpsz'
  call JSONreadInteger(json, ep, mcnl%globalworkgrpsz, defmcnl%globalworkgrpsz)
  ep = 'MCCLdata.num_el'
  call JSONreadInteger(json, ep, mcnl%num_el, defmcnl%num_el)
  ep = 'MCCLdata.totnum_el'
  call JSONreadInteger(json, ep, mcnl%totnum_el, defmcnl%totnum_el)
  ep = 'MCCLdata.devid'
  call JSONreadInteger(json, ep, mcnl%devid, defmcnl%devid)

  ep = 'MCCLdata.sig'
  call JSONreadDouble(json, ep, mcnl%sig, defmcnl%sig)
  ep = 'MCCLdata.omega'
  call JSONreadDouble(json, ep, mcnl%omega, defmcnl%omega)
  ep = 'MCCLdata.EkeV'
  call JSONreadDouble(json, ep, mcnl%EkeV, defmcnl%EkeV)
  ep = 'MCCLdata.Ehistmin'
  call JSONreadDouble(json, ep, mcnl%Ehistmin, defmcnl%Ehistmin)
  ep = 'MCCLdata.Ebinsize'
  call JSONreadDouble(json, ep, mcnl%Ebinsize, defmcnl%Ebinsize)
  ep = 'MCCLdata.depthmax'
  call JSONreadDouble(json, ep, mcnl%depthmax, defmcnl%depthmax)
  ep = 'MCCLdata.depthstep'
  call JSONreadDouble(json, ep, mcnl%depthstep, defmcnl%depthstep)

  ep = 'MCCLdata.MCmode'
  s = mcnl%MCmode
  s2 = defmcnl%MCmode
  call JSONreadString(json, ep, s, s2)
  ep = 'MCCLdata.xtalname'
  call JSONreadString(json, ep, mcnl%xtalname, defmcnl%xtalname)
  ep = 'MCCLdata.dataname'
  call JSONreadString(json, ep, mcnl%dataname, defmcnl%dataname)
  ep = 'MCCLdata.mode'
  call JSONreadString(json, ep, mcnl%mode, defmcnl%mode)
end if

call json%destroy(); call JSON_failtest(error_cnt)

end subroutine JSONreadMCCLNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONreadMCCLMultiLayerNameList
!
!> @author Saransh Singh/Marc De Graef, Carnegie Mellon University
!
!> @brief read namelist file and fill mcnl structure (used by EMMCCL.f90)
!
!> @param mcnl Monte Carloname list structure
!> @param jsonname input file name
!> @param error_cnt total number of errors encountered by json routines
!
!> @date 08/19/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONreadMCCLMultiLayerNameList(mcnl, jsonname, error_cnt)

use ISO_C_BINDING
use NameListHandlers

IMPLICIT NONE

type(MCCLMultiLayerNameListType),INTENT(INOUT)        :: mcnl
character(fnlen),INTENT(IN)                           :: jsonname
integer(kind=irg),INTENT(INOUT)                       :: error_cnt

type(json_file)                                       :: json    !the JSON structure read from the file:

type(MCCLMultiLayerNameListType)                      :: defmcnl
logical                                               :: init = .TRUE.
character(fnlen)                                      :: nmlfile = '', ep, s, s2
real(kind=wp)                                         :: rval
character(kind=CK,len=:),allocatable                  :: cval
real(wp),dimension(:),allocatable                     :: rvec

! first of all, open the file and return an error message if it does not exist
error_cnt = 0
call json_initialize(); call JSON_failtest(error_cnt)

! populate the json structure
call json%load_file(filename = trim(jsonname))
if (json_failed()) then    !if there was an error reading the file
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
else
! ok, we got here so we need to initialize the namelist first to its default values (set in NameListHandlers)
  call GetMCCLMultiLayerNameList(nmlfile, defmcnl, initonly=init)

! then we start reading the values in the json file  
  ep = 'MCCLdata.stdout'
  call JSONreadInteger(json, ep, mcnl%stdout, defmcnl%stdout)
  ep = 'MCCLdata.numsx'
  call JSONreadInteger(json, ep, mcnl%numsx, defmcnl%numsx)
  ep = 'MCCLdata.globalworkgrpsz'
  call JSONreadInteger(json, ep, mcnl%globalworkgrpsz, defmcnl%globalworkgrpsz)
  ep = 'MCCLdata.num_el'
  call JSONreadInteger(json, ep, mcnl%num_el, defmcnl%num_el)
  ep = 'MCCLdata.totnum_el'
  call JSONreadInteger(json, ep, mcnl%totnum_el, defmcnl%totnum_el)

  ep = 'MCCLdata.sig'
  call JSONreadDouble(json, ep, mcnl%sig, defmcnl%sig)
  ep = 'MCCLdata.omega'
  call JSONreadDouble(json, ep, mcnl%omega, defmcnl%omega)
  ep = 'MCCLdata.EkeV'
  call JSONreadDouble(json, ep, mcnl%EkeV, defmcnl%EkeV)
  ep = 'MCCLdata.Ehistmin'
  call JSONreadDouble(json, ep, mcnl%Ehistmin, defmcnl%Ehistmin)
  ep = 'MCCLdata.Ebinsize'
  call JSONreadDouble(json, ep, mcnl%Ebinsize, defmcnl%Ebinsize)
  ep = 'MCCLdata.depthmax'
  call JSONreadDouble(json, ep, mcnl%depthmax, defmcnl%depthmax)
  ep = 'MCCLdata.depthstep'
  call JSONreadDouble(json, ep, mcnl%depthstep, defmcnl%depthstep)
  ep = 'MCCLdata.filmthickness'
  call JSONreadDouble(json, ep, mcnl%filmthickness, defmcnl%filmthickness)
  ep = 'MCCLdata.filmstep'
  call JSONreadDouble(json, ep, mcnl%filmstep, defmcnl%filmstep)

  ep = 'MCCLdata.MCmode'
  s = mcnl%MCmode
  s2 = defmcnl%MCmode
  call JSONreadString(json, ep, s, s2)
  ep = 'MCCLdata.xtalname_film'
  call JSONreadString(json, ep, mcnl%xtalname_film, defmcnl%xtalname_film)
  ep = 'MCCLdata.xtalname_subs'
  call JSONreadString(json, ep, mcnl%xtalname_subs, defmcnl%xtalname_subs)
  ep = 'MCCLdata.dataname'
  call JSONreadString(json, ep, mcnl%dataname, defmcnl%dataname)
  ep = 'MCCLdata.mode'
  call JSONreadString(json, ep, mcnl%mode, defmcnl%mode)
end if

call json%destroy(); call JSON_failtest(error_cnt)

end subroutine JSONreadMCCLMultiLayerNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONreadEBSDMasterNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read json file and fill emnl structure (used by EMEBSDmaster.f90)
!
!> @param emnl EBSD master name list structure
!> @param jsonname input file name
!> @param error_cnt total number of errors encountered by json routines
!
!> @date 08/19/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONreadEBSDMasterNameList(emnl, jsonname, error_cnt)

use ISO_C_BINDING
use NameListHandlers

IMPLICIT NONE

type(EBSDMasterNameListType),INTENT(INOUT)            :: emnl
character(fnlen),INTENT(IN)                           :: jsonname
integer(kind=irg),INTENT(INOUT)                       :: error_cnt

type(json_file)                                       :: json    !the JSON structure read from the file:

type(EBSDMasterNameListType)                          :: defemnl
logical                                               :: init = .TRUE.
character(fnlen)                                      :: nmlfile = '', ep, s, s2
real(kind=wp)                                         :: rval
character(kind=CK,len=:),allocatable                  :: cval
real(wp),dimension(:),allocatable                     :: rvec

! first of all, open the file and return an error message if it does not exist
error_cnt = 0
call json_initialize(); call JSON_failtest(error_cnt)

! populate the json structure
call json%load_file(filename = trim(jsonname))
if (json_failed()) then    !if there was an error reading the file
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
else
! ok, we got here so we need to initialize the namelist first to its default values (set in NameListHandlers)
  call GetEBSDMasterNameList(nmlfile, defemnl, initonly=init)

! then we start reading the values in the json file  
  ep = 'EBSDmastervars.stdout'
  call JSONreadInteger(json, ep, emnl%stdout, defemnl%stdout)
  ep = 'EBSDmastervars.npx'
  call JSONreadInteger(json, ep, emnl%npx, defemnl%npx)
  ep = 'EBSDmastervars.Esel'
  call JSONreadInteger(json, ep, emnl%Esel, defemnl%Esel)
  ep = 'EBSDmastervars.nthreads'
  call JSONreadInteger(json, ep, emnl%nthreads, defemnl%nthreads)

  ep = 'EBSDmastervars.dmin'
  call JSONreadReal(json, ep, emnl%dmin, defemnl%dmin)

  ep = 'EBSDmastervars.energyfile'
  call JSONreadString(json, ep, emnl%energyfile, defemnl%energyfile)
  ep = 'EBSDmastervars.outname'
  call JSONreadString(json, ep, emnl%outname, defemnl%outname)

  ep = 'EBSDmastervars.restart'
  call JSONreadLogical(json, ep, emnl%restart, defemnl%restart)
end if

call json%destroy(); call JSON_failtest(error_cnt)

end subroutine JSONreadEBSDMasterNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONreadECPMasterNameList
!
!> @author Saransh Singh/Marc De Graef, Carnegie Mellon University
!
!> @brief read json file and fill mcnl structure (used by EMECPmaster.f90)
!
!> @param emnl ECP master name list structure
!> @param jsonname input file name
!> @param error_cnt total number of errors encountered by json routines
!
!> @date 08/20/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONreadECPMasterNameList(ecpnl, jsonname, error_cnt)

use ISO_C_BINDING
use NameListHandlers

IMPLICIT NONE

type(ECPMasterNameListType),INTENT(INOUT)             :: ecpnl
character(fnlen),INTENT(IN)                           :: jsonname
integer(kind=irg),INTENT(INOUT)                       :: error_cnt

type(json_file)                                       :: json    !the JSON structure read from the file:

type(ECPMasterNameListType)                           :: defecpnl
logical                                               :: init = .TRUE.
character(fnlen)                                      :: nmlfile = '', ep, s, s2
real(kind=wp)                                         :: rval
character(kind=CK,len=:),allocatable                  :: cval
real(wp),dimension(:),allocatable                     :: rvec

! first of all, open the file and return an error message if it does not exist
error_cnt = 0
call json_initialize(); call JSON_failtest(error_cnt)

! populate the json structure
call json%load_file(filename = trim(jsonname))
if (json_failed()) then    !if there was an error reading the file
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
else
! ok, we got here so we need to initialize the namelist first to its default values (set in NameListHandlers)
  call GetECPMasterNameList(nmlfile, defecpnl, initonly=init)

! then we start reading the values in the json file  
  ep = 'ECPmastervars.stdout'
  call JSONreadInteger(json, ep, ecpnl%stdout, defecpnl%stdout)
  ep = 'ECPmastervars.npx'
  call JSONreadInteger(json, ep, ecpnl%npx, defecpnl%npx)
  ep = 'ECPmastervars.Esel'
  call JSONreadInteger(json, ep, ecpnl%Esel, defecpnl%Esel)
  ep = 'ECPmastervars.nthreads'
  call JSONreadInteger(json, ep, ecpnl%nthreads, defecpnl%nthreads)

  ep = 'ECPmastervars.startthick'
  call JSONreadReal(json, ep, ecpnl%startthick, defecpnl%startthick)
  ep = 'ECPmastervars.dmin'
  call JSONreadReal(json, ep, ecpnl%dmin, defecpnl%dmin)

  ep = 'ECPmastervars.fn'
  call JSONreadRealVec(json, ep, ecpnl%fn, defecpnl%fn, size(ecpnl%fn))
  ep = 'ECPmastervars.abcdist'
  call JSONreadRealVec(json, ep, ecpnl%abcdist, defecpnl%abcdist, size(ecpnl%abcdist))
  ep = 'ECPmastervars.albegadist'
  call JSONreadRealVec(json, ep, ecpnl%albegadist, defecpnl%albegadist, size(ecpnl%albegadist))

  ep = 'ECPmastervars.compmode'
  call JSONreadString(json, ep, ecpnl%compmode, defecpnl%compmode)
  ep = 'ECPmastervars.energyfile'
  call JSONreadString(json, ep, ecpnl%energyfile, defecpnl%energyfile)
  ep = 'ECPmastervars.outname'
  call JSONreadString(json, ep, ecpnl%outname, defecpnl%outname)

  ep = 'ECPmastervars.distort'
  call JSONreadLogical(json, ep, ecpnl%distort, defecpnl%distort)
end if

call json%destroy(); call JSON_failtest(error_cnt)

end subroutine JSONreadECPMasterNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONreadEBSDNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read json file and fill enl structure (used by EMEBSD.f90)
!
!> @param enl EBSD name list structure
!> @param jsonname input file name
!> @param error_cnt total number of errors encountered by json routines
!
!> @date 08/20/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONreadEBSDNameList(enl, jsonname, error_cnt)

use ISO_C_BINDING
use NameListHandlers
use error

IMPLICIT NONE

type(EBSDNameListType),INTENT(INOUT)                  :: enl
character(fnlen),INTENT(IN)                           :: jsonname
integer(kind=irg),INTENT(INOUT)                       :: error_cnt

type(json_file)                                       :: json    !the JSON structure read from the file:

type(EBSDNameListType)                                :: defenl
logical                                               :: init = .TRUE.
character(fnlen)                                      :: nmlfile = '', ep, s, s2
real(kind=wp)                                         :: rval
character(kind=CK,len=:),allocatable                  :: cval
real(wp),dimension(:),allocatable                     :: rvec

! first of all, open the file and return an error message if it does not exist
error_cnt = 0
call json_initialize(); call JSON_failtest(error_cnt)

! populate the json structure
call json%load_file(filename = trim(jsonname))
if (json_failed()) then    !if there was an error reading the file
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
else
! ok, we got here so we need to initialize the namelist first to its default values (set in NameListHandlers)
  call GetEBSDNameList(nmlfile, defenl, initonly=init)

  ep = 'EBSDdata.stdout'
  call JSONreadInteger(json, ep, enl%stdout, defenl%stdout)
  ep = 'EBSDdata.numsx'
  call JSONreadInteger(json, ep, enl%numsx, defenl%numsx)
  ep = 'EBSDdata.numsy'
  call JSONreadInteger(json, ep, enl%numsy, defenl%numsy)
  ep = 'EBSDdata.binning'
  call JSONreadInteger(json, ep, enl%binning, defenl%binning)
  ep = 'EBSDdata.nthreads'
  call JSONreadInteger(json, ep, enl%nthreads, defenl%nthreads)
  ep = 'EBSDdata.energyaverage'
  call JSONreadInteger(json, ep, enl%energyaverage, defenl%energyaverage)

  ep = 'EBSDdata.L'
  call JSONreadReal(json, ep, enl%L, defenl%L)
  ep = 'EBSDdata.thetac'
  call JSONreadReal(json, ep, enl%thetac, defenl%thetac)
  ep = 'EBSDdata.delta'
  call JSONreadReal(json, ep, enl%delta, defenl%delta)
  ep = 'EBSDdata.xpc'
  call JSONreadReal(json, ep, enl%xpc, defenl%xpc)
  ep = 'EBSDdata.ypc'
  call JSONreadReal(json, ep, enl%ypc, defenl%ypc)
  ep = 'EBSDdata.omega'
  call JSONreadReal(json, ep, enl%omega, defenl%omega)
  ep = 'EBSDdata.energymin'
  call JSONreadReal(json, ep, enl%energymin, defenl%energymin)
  ep = 'EBSDdata.energymax'
  call JSONreadReal(json, ep, enl%energymax, defenl%energymax)
  ep = 'EBSDdata.gammavalue'
  call JSONreadReal(json, ep, enl%gammavalue, defenl%gammavalue)

  ep = 'EBSDdata.axisangle'
  call JSONreadRealVec(json, ep, enl%axisangle, defenl%axisangle, size(enl%axisangle))

  ep = 'EBSDdata.beamcurrent'
  call JSONreadDouble(json, ep, enl%beamcurrent, defenl%beamcurrent)
  ep = 'EBSDdata.dwelltime'
  call JSONreadDouble(json, ep, enl%dwelltime, defenl%dwelltime)

  ep = 'EBSDdata.maskpattern'
  s = enl%maskpattern
  s2 = defenl%maskpattern
  call JSONreadString(json, ep, s, s2)
  ep = 'EBSDdata.scalingmode'
  s = enl%scalingmode
  s2 = defenl%scalingmode
  call JSONreadString(json, ep, s, s2)
  ep = 'EBSDdata.eulerconvention'
  s = enl%eulerconvention
  s2 = defenl%eulerconvention
  call JSONreadString(json, ep, s, s2)
  ep = 'EBSDdata.outputformat'
  s = enl%outputformat
  s2 = defenl%outputformat
  call JSONreadString(json, ep, s, s2)

  ep = 'EBSDdata.anglefile'
  call JSONreadString(json, ep, enl%anglefile, defenl%anglefile)
  ep = 'EBSDdata.masterfile'
  call JSONreadString(json, ep, enl%masterfile, defenl%masterfile)
  ep = 'EBSDdata.energyfile'
  call JSONreadString(json, ep, enl%energyfile, defenl%energyfile)
  ep = 'EBSDdata.datafile'
  call JSONreadString(json, ep, enl%datafile, defenl%datafile)
end if

call json%destroy(); call JSON_failtest(error_cnt)

end subroutine JSONreadEBSDNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSONreadEBSDoverlapNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read jsonfile and fill enl structure (used by EMEBSDoverlap.f90)
!
!> @param enl EBSD name list structure
!> @param jsonname input file name
!> @param error_cnt total number of errors encountered by json routines
!
!> @date 08/20/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONreadEBSDoverlapNameList(enl, jsonname, error_cnt)

use ISO_C_BINDING
use NameListHandlers
use error

IMPLICIT NONE

type(EBSDoverlapNameListType),INTENT(INOUT)           :: enl
character(fnlen),INTENT(IN)                           :: jsonname
integer(kind=irg),INTENT(INOUT)                       :: error_cnt

type(json_file)                                       :: json    !the JSON structure read from the file:

type(EBSDoverlapNameListType)                         :: defenl
logical                                               :: init = .TRUE.
character(fnlen)                                      :: nmlfile = '', ep, s, s2
real(kind=wp)                                         :: rval
character(kind=CK,len=:),allocatable                  :: cval
real(wp),dimension(:),allocatable                     :: rvec

! first of all, open the file and return an error message if it does not exist
error_cnt = 0
call json_initialize(); call JSON_failtest(error_cnt)

! populate the json structure
call json%load_file(filename = trim(jsonname))
if (json_failed()) then    !if there was an error reading the file
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
else
! ok, we got here so we need to initialize the namelist first to its default values (set in NameListHandlers)
  call GetEBSDoverlapNameList(nmlfile, defenl, initonly=init)

  ep = 'EBSDdata.stdout'
  call JSONreadInteger(json, ep, enl%stdout, defenl%stdout)

  ep = 'EBSDdata.PatternAxisA'
  call JSONreadIntegerVec(json, ep, enl%PatternAxisA, defenl%PatternAxisA, size(defenl%PatternAxisA))
  ep = 'EBSDdata.HorizontalAxisA'
  call JSONreadIntegerVec(json, ep, enl%HorizontalAxisA, defenl%HorizontalAxisA, size(defenl%HorizontalAxisA))

  ep = 'EBSDdata.tA'
  call JSONreadRealVec(json, ep, enl%tA, defenl%tA, size(enl%tA))
  ep = 'EBSDdata.tB'
  call JSONreadRealVec(json, ep, enl%tB, defenl%tB, size(enl%tB))
  ep = 'EBSDdata.gA'
  call JSONreadRealVec(json, ep, enl%gA, defenl%gA, size(enl%gA))
  ep = 'EBSDdata.gB'
  call JSONreadRealVec(json, ep, enl%gB, defenl%gB, size(enl%gB))

  ep = 'EBSDdata.fracA'
  call JSONreadReal(json, ep, enl%fracA, defenl%fracA)

  ep = 'EBSDdata.masterfileA'
  call JSONreadString(json, ep, enl%masterfileA, defenl%masterfileA)
  ep = 'EBSDdata.masterfileB'
  call JSONreadString(json, ep, enl%masterfileB, defenl%masterfileB)
  ep = 'EBSDdata.datafile'
  call JSONreadString(json, ep, enl%datafile, defenl%datafile)
end if

call json%destroy(); call JSON_failtest(error_cnt)

end subroutine JSONreadEBSDoverlapNameList

! line 870 of NameListHandlers

end module JSONsupport
