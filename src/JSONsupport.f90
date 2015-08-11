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
! from here on we have the Namlist->json conversion routines for all the 
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
call JSON_writeNMLreals(inp, io_real, reallist, n_real, error_cnt)

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
!> @date 03/21/15 MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine JSONwriteMCCLMultiLayerNameList(mcnl, jsonname, error_cnt)

use ISO_C_BINDING

IMPLICIT NONE

type(MCCLMultiLayerNameListType),INTENT(INOUT)        :: mcnl
character(fnlen),INTENT(IN)                           :: jsonname
integer(kind=irg),INTENT(INOUT)                       :: error_cnt

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


end module JSONsupport
