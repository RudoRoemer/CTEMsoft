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
!--------------------------------------------------------------------------
module JSONsupport

use local
use typedefs
use NameListTypedefs
use json_module
use, intrinsic :: iso_fortran_env , only: error_unit

IMPLICIT NONE

contains

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
  call json_add(inp, dataset, io_int(i))
  if (json_failed()) then
    call json_print_error_message(error_unit)
    error_cnt = error_cnt + 1
  end if
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
  call json_add(inp, dataset, dble(io_real(i)))
  if (json_failed()) then
    call json_print_error_message(error_unit)
    error_cnt = error_cnt + 1
  end if
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
  call json_add(inp, dataset, io_real(i))
  if (json_failed()) then
    call json_print_error_message(error_unit)
    error_cnt = error_cnt + 1
  end if
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
call json_initialize()
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

! create the json root pointer
call json_create_object(p,trim(jsonname))   
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

! we'll use the namelist name to configure the inp structure and add it to p
call json_create_object(inp,trim(namelistname))  
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if
call json_add(p, inp)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

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
call json_print(p,dataunit)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if
close(dataunit)

! final cleanup
call json_destroy(p)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

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
call json_add(inp, dataset, knl%k)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'fn'
call json_add(inp, dataset, knl%fn)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

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
call json_add(inp, dataset, knl%xtalname)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'outname'
call json_add(inp, dataset, knl%outname)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

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
call json_add(inp, dataset, knl%Kosselmode)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'xtalname'
call json_add(inp, dataset, knl%xtalname)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'outname'
call json_add(inp, dataset, knl%outname)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

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
call json_add(inp, dataset, mcnl%MCmode)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'xtalname'
call json_add(inp, dataset, mcnl%xtalname)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'dataname'
call json_add(inp, dataset, mcnl%dataname)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

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
!--------------------------------------------------------------------------
subroutine JSONwriteMCCLNameList(mcnl, jsonname, error_cnt)

use ISO_C_BINDING

IMPLICIT NONE

type(MCCLNameListType),INTENT(INOUT)                  :: mcnl
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
call json_add(inp, dataset, mcnl%MCmode)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'xtalname'
call json_add(inp, dataset, mcnl%xtalname)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'dataname'
call json_add(inp, dataset, mcnl%dataname)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'mode'
call json_add(inp, dataset, mcnl%mode)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

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
call json_add(inp, dataset, mcnl%MCmode)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'xtalname_film'
call json_add(inp, dataset, mcnl%xtalname_film)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'xtalname_subs'
call json_add(inp, dataset, mcnl%xtalname_subs)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'dataname'
call json_add(inp, dataset, mcnl%dataname)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'mode'
call json_add(inp, dataset, mcnl%mode)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

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
call json_add(inp, dataset, dble(emnl%dmin))
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

! write all the strings
dataset = 'outname'
call json_add(inp, dataset, emnl%outname)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'energyfile'
call json_add(inp, dataset, emnl%energyfile)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

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

integer(kind=irg),parameter                           :: n_int = 3, n_real = 2
integer(kind=irg)                                     :: io_int(n_int)
real(kind=dbl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
character(fnlen)                                      :: dataset, namelistname
character(fnlen,kind=c_char)                          :: line2(1)

! initialize the json state variables
namelistname = 'ECPmastervars'
call JSON_initpointers(p, inp, jsonname, namelistname, error_cnt)

! write all the single integers
io_int = (/ ecpnl%stdout, ecpnl%Esel, ecpnl%npx /)
intlist(1) = 'stdout'
intlist(2) = 'Esel'
intlist(3) = 'npx'
call JSON_writeNMLintegers(inp, io_int, intlist, n_int, error_cnt)

dataset = 'distort'
call json_add(inp, dataset, ecpnl%distort)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

! integer vectors
dataset = 'fn'
call json_add(inp, dataset, dble(ecpnl%fn))
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

! write all the single doubles
io_real = (/ ecpnl%dmin, ecpnl%startthick /)
reallist(1) = 'dmin'
reallist(2) = 'startthick' 
call JSON_writeNMLdoubles(inp, io_real, reallist, n_real, error_cnt)

! 3-vectors (real)
dataset = 'abcdist'
call json_add(inp, dataset, dble(ecpnl%abcdist))
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'albegadist'
call json_add(inp, dataset, dble(ecpnl%albegadist))
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

! write all the strings
dataset = 'outname'
call json_add(inp, dataset, ecpnl%outname)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'energyfile'
call json_add(inp, dataset, ecpnl%energyfile)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'compmode'
call json_add(inp, dataset, ecpnl%compmode)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

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
call json_add(inp, dataset, dble(enl%axisangle))
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

! a few doubles
dataset = 'beamcurrent'
call json_add(inp, dataset, enl%beamcurrent)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'dwelltime'
call json_add(inp, dataset, enl%dwelltime)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

! write all the strings
dataset = 'maskpattern'
call json_add(inp, dataset, enl%maskpattern)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'scalingmode'
call json_add(inp, dataset, enl%scalingmode)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'eulerconvention'
call json_add(inp, dataset, enl%eulerconvention)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'outputformat'
call json_add(inp, dataset, enl%outputformat)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'energyfile'
call json_add(inp, dataset, enl%energyfile)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'masterfile'
call json_add(inp, dataset, enl%masterfile)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'anglefile'
call json_add(inp, dataset, enl%anglefile)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'datafile'
call json_add(inp, dataset, enl%datafile)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

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
call json_add(inp, dataset, ecpnl%k)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'fn'
call json_add(inp, dataset, ecpnl%fn)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'gF'
call json_add(inp, dataset, ecpnl%gF)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'gS'
call json_add(inp, dataset, ecpnl%gS)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'tF'
call json_add(inp, dataset, ecpnl%tF)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'tS'
call json_add(inp, dataset, ecpnl%tS)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

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
call json_add(inp, dataset, ecpnl%compmode)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'energyfile'
call json_add(inp, dataset, ecpnl%energyfile)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'outname'
call json_add(inp, dataset, ecpnl%outname)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'xtalname'
call json_add(inp, dataset, ecpnl%xtalname)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'xtalname2'
call json_add(inp, dataset, ecpnl%xtalname2)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

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
call json_add(inp, dataset, lacbednl%k)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'fn'
call json_add(inp, dataset, lacbednl%fn)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

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
call json_add(inp, dataset, lacbednl%outname)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'xtalname'
call json_add(inp, dataset, lacbednl%xtalname)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

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
call json_add(inp, dataset, dble(ecpnl%thetac))
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

! real vector
dataset = 'k'
call json_add(inp, dataset, dble(ecpnl%k))
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

! write all the strings
dataset = 'outname'
call json_add(inp, dataset, ecpnl%outname)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'masterfile'
call json_add(inp, dataset, ecpnl%masterfile)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

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
call json_add(inp, dataset, pednl%outname)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'xtalname'
call json_add(inp, dataset, pednl%xtalname)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'eulername'
call json_add(inp, dataset, pednl%eulername)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

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
call json_add(inp, dataset, pednl%k)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'fn'
call json_add(inp, dataset, pednl%fn)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

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
call json_add(inp, dataset, pednl%outname)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'xtalname'
call json_add(inp, dataset, pednl%xtalname)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'filemode'
call json_add(inp, dataset, pednl%filemode)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

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
call json_add(inp, dataset, eccinl%k)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

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
call json_add(inp, dataset, dble(eccinl%lauec))
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'lauec2'
call json_add(inp, dataset, dble(eccinl%lauec2))
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

! write all the strings
dataset = 'dispmode'
call json_add(inp, dataset, eccinl%dispmode)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'summode'
call json_add(inp, dataset, eccinl%summode)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'progmode'
call json_add(inp, dataset, eccinl%progmode)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'xtalname'
call json_add(inp, dataset, eccinl%xtalname)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'foilnmlfile'
call json_add(inp, dataset, eccinl%foilnmlfile)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'dispfile'
call json_add(inp, dataset, eccinl%dispfile)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'dataname'
call json_add(inp, dataset, eccinl%dataname)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'ECPname'
call json_add(inp, dataset, eccinl%ECPname)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'sgname'
call json_add(inp, dataset, eccinl%sgname)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'apbname'
call json_add(inp, dataset, eccinl%apbname)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'incname'
call json_add(inp, dataset, eccinl%dispmode)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'voidname'
call json_add(inp, dataset, eccinl%voidname)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if


! maxdefects string arrays
dataset = 'sfname'
call json_add(inp, dataset, eccinl%sfname(1:eccinl%numsf))
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

! 3*maxdefects string arrays
dataset = 'dislYname'
call json_add(inp, dataset, eccinl%dislYname(1:eccinl%numYdisl))
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'dislname'
call json_add(inp, dataset, eccinl%dislname(1:eccinl%numdisl))
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

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
dataset = 'outname'
call json_add(inp, dataset, rfznl%outname)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

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
call json_add(inp, dataset, dictindxnl%MeanSubtraction)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

! strings
dataset = 'exptfile'
call json_add(inp, dataset, dictindxnl%exptfile)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'dictfile'
call json_add(inp, dataset, dictindxnl%dictfile)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

dataset = 'eulerfile'
call json_add(inp, dataset, dictindxnl%eulerfile)
if (json_failed()) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

! and then we write the file and clean up
call JSON_cleanuppointers(p, inp, jsonname, error_cnt)

end subroutine JSONwriteDictIndxOpenCLNameList



end module JSONsupport
