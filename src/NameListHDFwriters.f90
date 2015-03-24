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
! EMsoft:NameListHDFwriters.f90
!--------------------------------------------------------------------------
!
! PROGRAM: NameListHDFwriters
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief routines for reading and returning name list type structures
!
!> @date 03/20/15 MDG 1.0 original, completed on 3/23/15
!--------------------------------------------------------------------------
module NameListHDFwriters

use local
use typedefs
use NameListTypedefs
use HDF5
use h5lt
use HDFsupport

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! first some auxiliary routines to make things easier later on
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!
! SUBROUTINE:HDF_writeNMLintegers
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write a series of integer namelist entries to an HDF file
!
!> @param HDF_head top of stack pointer
!> @param io_int list of integers
!> @param intlist list of string descriptors
!> @param n_int number of entries
!
!> @date 03/20/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine HDF_writeNMLintegers(HDF_head, io_int, intlist, n_int)



IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head
integer(kind=irg),INTENT(IN)                          :: io_int(n_int)
character(20),INTENT(IN)                              :: intlist(n_int)
integer(kind=irg),INTENT(IN)                          :: n_int

integer(kind=irg)                                     :: rnk,  error, ioval(1), i
integer(HSIZE_T)                                      :: dims(1)

rnk = 1
dims(1) = 1

do i=1,n_int
  ioval(1) = io_int(i)
  call h5ltmake_dataset_f(HDF_head%objectID, trim(intlist(i)), rnk, dims, H5T_NATIVE_INTEGER, ioval, error)
  if (error.ne.0) call HDF_handleError(error,'HDF_writeNMLintegers: unable to create '//trim(intlist(i))//' dataset',.TRUE.)
end do

end subroutine HDF_writeNMLintegers

!--------------------------------------------------------------------------
!
! SUBROUTINE:HDF_writeNMLreals
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write a series of real namelist entries to an HDF file
!
!> @param HDF_head top of stack pointer
!> @param io_real list of reals
!> @param reallist list of string descriptors
!> @param n_real number of entries
!
!> @date 03/20/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine HDF_writeNMLreals(HDF_head, io_real, reallist, n_real)



IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head
real(kind=sgl),INTENT(IN)                             :: io_real(n_real)
character(20),INTENT(IN)                              :: reallist(n_real)
integer(kind=irg),INTENT(IN)                          :: n_real

integer(kind=irg)                                     :: rnk,  error, i
integer(HSIZE_T)                                      :: dims(1)
real(kind=sgl)                                        :: ioval(1)

rnk = 1
dims(1) = 1

do i=1,n_real
  ioval(1) = io_real(i)
  call h5ltmake_dataset_f(HDF_head%objectID, trim(reallist(i)), rnk, dims, H5T_NATIVE_REAL, ioval, error)
  if (error.ne.0) call HDF_handleError(error,'HDF_writeNMLreals: unable to create '//trim(reallist(i))//' dataset',.TRUE.)
end do

end subroutine HDF_writeNMLreals

!--------------------------------------------------------------------------
!
! SUBROUTINE:HDF_writeNMLdbles
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write a series of double precisin namelist entries to an HDF file
!
!> @param HDF_head top of stack pointer
!> @param io_real list of doubles
!> @param reallist list of string descriptors
!> @param n_real number of entries
!
!> @date 03/20/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine HDF_writeNMLdbles(HDF_head, io_real, reallist, n_real)



IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head
real(kind=dbl),INTENT(IN)                             :: io_real(n_real)
character(20),INTENT(IN)                              :: reallist(n_real)
integer(kind=irg),INTENT(IN)                          :: n_real

integer(kind=irg)                                     :: rnk,  error, i
integer(HSIZE_T)                                      :: dims(1)
real(kind=dbl)                                        :: ioval(1)

rnk = 1
dims(1) = 1

do i=1,n_real
  ioval(1) = io_real(i)
  call h5ltmake_dataset_f(HDF_head%objectID, trim(reallist(i)), rnk, dims, H5T_NATIVE_DOUBLE, ioval, error)
  if (error.ne.0) call HDF_handleError(error,'HDF_writeNMLdbles: unable to create '//trim(reallist(i))//' dataset',.TRUE.)
end do

end subroutine HDF_writeNMLdbles


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! then the actual namelist write routines
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!
! SUBROUTINE:HDFwriteKosselNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist file into HDF file
!
!> @param HDF_head top of push stack
!> @param knl Kossel name list structure
!
!> @date 03/20/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine HDFwriteKosselNameList(HDF_head, HDF_tail, knl)



IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(KosselNameListType),INTENT(IN)                   :: knl

integer(HID_T)                                        :: grp_id
integer(kind=irg),parameter                           :: n_int = 5, n_real = 6
integer(kind=irg)                                     :: rnk, error,  io_int(n_int)
integer(HSIZE_T)                                      :: dims(1)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)

! create the group for this namelist
error = HDF_createGroup('KosselNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ knl%stdout, knl%numthick, knl%npix, knl%maxHOLZ, knl%nthreads /)
intlist(1) = 'stdout'
intlist(2) = 'numthick'
intlist(3) = 'npix'
intlist(4) = 'maxHOLZ'
intlist(5) = 'nthreads'
call HDF_writeNMLintegers(HDF_head, io_int, intlist, n_int)

! integer vectors
rnk = 1
dims(1) = 3
call h5ltmake_dataset_f(HDF_head%objectID, 'k', rnk, dims, H5T_NATIVE_INTEGER, knl%k, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteKosselNameList: unable to create k dataset',.TRUE.)

call h5ltmake_dataset_f(HDF_head%objectID, 'fn', rnk, dims, H5T_NATIVE_INTEGER, knl%fn, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteKosselNameList: unable to create fn dataset',.TRUE.)

! write all the single reals
io_real = (/ knl%voltage, knl%dmin, knl%convergence, knl%startthick, knl%thickinc, knl%minten /)
reallist(1) = 'voltage'
reallist(2) = 'dmin'
reallist(3) = 'convergence'
reallist(4) = 'startthick'
reallist(5) = 'thickinc'
reallist(6) = 'minten'
call HDF_writeNMLreals(HDF_head, io_real, reallist, n_real)

! write all the strings
call h5ltmake_dataset_string_f(HDF_head%objectID, 'xtalname', knl%xtalname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteKosselNameList: unable to create xtalname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'outname', knl%outname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteKosselNameList: unable to create outname dataset',.TRUE.)

! and pop this group off the stack
call HDF_pop(HDF_head)

end subroutine HDFwriteKosselNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:HDFwriteKosselMasterNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist file into HDF file
!
!> @param HDF_head top of push stack
!> @param knl Kossel name list structure
!
!> @date 03/21/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine HDFwriteKosselMasterNameList(HDF_head, HDF_tail, knl)



IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(KosselMasterNameListType),INTENT(IN)             :: knl

integer(HID_T)                                        :: grp_id
integer(kind=irg),parameter                           :: n_int = 4, n_real = 5
integer(kind=irg)                                     :: rnk, error,  io_int(n_int)
integer(HSIZE_T)                                      :: dims(1)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)

! create the group for this namelist
error = HDF_createGroup('KosselMasterNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ knl%stdout, knl%numthick, knl%npix, knl%nthreads /)
intlist(1) = 'stdout'
intlist(2) = 'numthick'
intlist(3) = 'npix'
intlist(4) = 'nthreads' 
call HDF_writeNMLintegers(HDF_head, io_int, intlist, n_int)

! write all the single reals
io_real = (/ knl%voltage, knl%dmin, knl%startthick, knl%thickinc, knl%tfraction /)
reallist(1) = 'voltage'
reallist(2) = 'dmin'
reallist(3) = 'startthick'
reallist(4) = 'thickinc'
reallist(5) = 'tfraction' 
call HDF_writeNMLreals(HDF_head, io_real, reallist, n_real)

! write all the strings
call h5ltmake_dataset_string_f(HDF_head%objectID, 'Kosselmode', knl%Kosselmode, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteKosselMasterNameList: unable to create Kosselmode dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'xtalname', knl%xtalname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteKosselMasterNameList: unable to create xtalname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'outname', knl%outname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteKosselMasterNameList: unable to create outname dataset',.TRUE.)

! and pop this group off the stack
call HDF_pop(HDF_head)

end subroutine HDFwriteKosselMasterNameList


!--------------------------------------------------------------------------
!
! SUBROUTINE:HDFwriteMCNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist file into HDF file
!
!> @param HDF_head top of push stack
!> @param mcnl Monte Carlo name list structure
!
!> @date 03/21/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine HDFwriteMCNameList(HDF_head, HDF_tail, mcnl)



IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(MCNameListType),INTENT(INOUT)                    :: mcnl

integer(HID_T)                                        :: grp_id
integer(kind=irg),parameter                           :: n_int = 5, n_real = 7
integer(kind=irg)                                     :: rnk, error,  io_int(n_int)
integer(HSIZE_T)                                      :: dims(1)
real(kind=dbl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)

! create the group for this namelist
error = HDF_createGroup('MCNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ mcnl%stdout, mcnl%numsx, mcnl%primeseed, mcnl%num_el, mcnl%nthreads /)
intlist(1) = 'stdout'
intlist(2) = 'numsx'
intlist(3) = 'primeseed'
intlist(4) = 'num_el'
intlist(5) = 'nthreads'
call HDF_writeNMLintegers(HDF_head, io_int, intlist, n_int)

! write all the single doubles
io_real = (/ mcnl%sig, mcnl%omega, mcnl%EkeV, mcnl%Ehistmin, mcnl%Ebinsize, mcnl%depthmax, mcnl%depthstep /)
reallist(1) = 'sig'
reallist(2) = 'omega'
reallist(3) = 'EkeV'
reallist(4) = 'Ehistmin'
reallist(5) = 'Ebinsize'
reallist(6) = 'depthmax'
reallist(7) = 'depthstep'
call HDF_writeNMLdbles(HDF_head, io_real, reallist, n_real)

! write all the strings
call h5ltmake_dataset_string_f(HDF_head%objectID, 'MCmode', mcnl%MCmode, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteMCNameList: unable to create MCmode dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'xtalname', mcnl%xtalname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteMCNameList: unable to create xtalname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'dataname', mcnl%dataname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteMCNameList: unable to create dataname dataset',.TRUE.)

! and pop this group off the stack
call HDF_pop(HDF_head)

end subroutine HDFwriteMCNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:HDFwriteMCCLNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist to HDF file
!
!> @param HDF_head top of push stack
!> @param mcnl Monte Carlon ame list structure
!
!> @date 03/21/15 MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine HDFwriteMCCLNameList(HDF_head, HDF_tail, mcnl)



IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(MCCLNameListType),INTENT(INOUT)                  :: mcnl

integer(HID_T)                                        :: grp_id
integer(kind=irg),parameter                           :: n_int = 5, n_real = 7
integer(kind=irg)                                     :: rnk, error,  io_int(n_int)
integer(HSIZE_T)                                      :: dims(1)
real(kind=dbl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)

! create the group for this namelist
error = HDF_createGroup('MCCLNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ mcnl%stdout, mcnl%numsx, mcnl%globalworkgrpsz, mcnl%num_el, mcnl%totnum_el /)
intlist(1) = 'stdout'
intlist(2) = 'numsx'
intlist(3) = 'globalworkgrpsz'
intlist(4) = 'num_el'
intlist(5) = 'totnum_el'
call HDF_writeNMLintegers(HDF_head, io_int, intlist, n_int)

! write all the single doubles
io_real = (/ mcnl%sig, mcnl%omega, mcnl%EkeV, mcnl%Ehistmin, mcnl%Ebinsize, mcnl%depthmax, mcnl%depthstep /)
reallist(1) = 'sig'
reallist(2) = 'omega'
reallist(3) = 'EkeV'
reallist(4) = 'Ehistmin'
reallist(5) = 'Ebinsize'
reallist(6) = 'depthmax'
reallist(7) = 'depthstep'
call HDF_writeNMLdbles(HDF_head, io_real, reallist, n_real)

! write all the strings
call h5ltmake_dataset_string_f(HDF_head%objectID, 'MCmode', mcnl%MCmode, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteMCCLNameList: unable to create MCmode dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'xtalname', mcnl%xtalname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteMCCLNameList: unable to create xtalname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'dataname', mcnl%dataname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteMCCLNameList: unable to create dataname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'primelist', mcnl%primelist, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteMCCLNameList: unable to create primelist dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'mode', mcnl%mode, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteMCCLNameList: unable to create mode dataset',.TRUE.)

! and pop this group off the stack
call HDF_pop(HDF_head)

end subroutine HDFwriteMCCLNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:HDFwriteMCCLMultiLayerNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist to HDF file
!
!> @param HDF_head top of push stack
!> @param mcnl Monte Carlo name list structure
!
!> @date 03/21/15 MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine HDFwriteMCCLMultiLayerNameList(HDF_head, HDF_tail, mcnl)



IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(MCCLMultiLayerNameListType),INTENT(INOUT)        :: mcnl

integer(HID_T)                                        :: grp_id
integer(kind=irg),parameter                           :: n_int = 5, n_real = 9
integer(kind=irg)                                     :: rnk, error,  io_int(n_int)
integer(HSIZE_T)                                      :: dims(1)
real(kind=dbl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)

! create the group for this namelist
error = HDF_createGroup('MCCLMultiLayerNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ mcnl%stdout, mcnl%numsx, mcnl%globalworkgrpsz, mcnl%num_el, mcnl%totnum_el /)
intlist(1) = 'stdout'
intlist(2) = 'numsx'
intlist(3) = 'globalworkgrpsz'
intlist(4) = 'num_el'
intlist(5) = 'totnum_el'
call HDF_writeNMLintegers(HDF_head, io_int, intlist, n_int)

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
call HDF_writeNMLdbles(HDF_head, io_real, reallist, n_real)

! write all the strings
call h5ltmake_dataset_string_f(HDF_head%objectID, 'MCmode', mcnl%MCmode, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteMCCLMultiLayerNameList: unable to create MCmode dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'xtalname_film', mcnl%xtalname_film, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteMCCLMultiLayerNameList: unable to create xtalname_film dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'xtalname_subs', mcnl%xtalname_subs, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteMCCLMultiLayerNameList: unable to create xtalname_subs dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'dataname', mcnl%dataname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteMCCLMultiLayerNameList: unable to create dataname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'primelist', mcnl%primelist, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteMCCLMultiLayerNameList: unable to create primelist dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'mode', mcnl%mode, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteMCCLMultiLayerNameList: unable to create mode dataset',.TRUE.)

! and pop this group off the stack
call HDF_pop(HDF_head)

end subroutine HDFwriteMCCLMultiLayerNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:HDFwriteEBSDMasterNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist to HDF file
!
!> @param HDF_head top of push stack
!> @param emnl EBSD master name list structure
!
!> @date 03/21/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine HDFwriteEBSDMasterNameList(HDF_head, HDF_tail, emnl)



IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(EBSDMasterNameListType),INTENT(INOUT)            :: emnl

integer(HID_T)                                        :: grp_id
integer(kind=irg),parameter                           :: n_int = 4, n_real = 1
integer(kind=irg)                                     :: rnk, error,  io_int(n_int)
integer(HSIZE_T)                                      :: dims(1)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)

! create the group for this namelist
error = HDF_createGroup('EBSDMasterNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ emnl%stdout, emnl%npx, emnl%Esel, emnl%nthreads /)
intlist(1) = 'stdout'
intlist(2) = 'npx'
intlist(3) = 'Esel'
intlist(4) = 'nthreads'
call HDF_writeNMLintegers(HDF_head, io_int, intlist, n_int)

! write a single real
io_real(1) = emnl%dmin
rnk = 1
dims(1) = 1
call h5ltmake_dataset_f(HDF_head%objectID, 'dmin', rnk, dims, H5T_NATIVE_REAL, io_real, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteEBSDMasterNameList: unable to create dmin dataset',.TRUE.)

! write all the strings
call h5ltmake_dataset_string_f(HDF_head%objectID, 'outname', emnl%outname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteEBSDMasterNameList: unable to create outname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'energyfile', emnl%energyfile, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteEBSDMasterNameList: unable to create energyfile dataset',.TRUE.)

! and pop this group off the stack
call HDF_pop(HDF_head)

end subroutine HDFwriteEBSDMasterNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:HDFwriteECPMasterNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist to HDF file
!
!> @param HDF_head top of push stack
!> @param ecpnl ECP master name list structure
!
!> @date 03/22/15 MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine HDFwriteECPMasterNameList(HDF_head, HDF_tail, ecpnl)



IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(ECPMasterNameListType),INTENT(INOUT)             :: ecpnl

integer(HID_T)                                        :: grp_id
integer(kind=irg),parameter                           :: n_int = 4, n_real = 2
integer(kind=irg)                                     :: rnk, error,  io_int(n_int), distort
integer(HSIZE_T)                                      :: dims(1)
real(kind=dbl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)

! create the group for this namelist
error = HDF_createGroup('ECPMasterNameList',HDF_head, HDF_tail)

! write all the single integers
! distort is a logical, for which there is no real HDF_NATIVE_anything conversion, so we'll store it as a 1 or 0
if (ecpnl%distort) then 
  distort = 1
else 
  distort = 0
end if
io_int = (/ ecpnl%stdout, ecpnl%Esel, ecpnl%npx, distort /)
intlist(1) = 'stdout'
intlist(2) = 'Esel'
intlist(3) = 'npx'
intlist(4) = 'distort'
call HDF_writeNMLintegers(HDF_head, io_int, intlist, n_int)

! integer vectors
rnk = 1
dims(1) = 3
call h5ltmake_dataset_f(HDF_head%objectID, 'fn', rnk, dims, H5T_NATIVE_INTEGER, ecpnl%fn, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPMasterNameList: unable to create fn dataset',.TRUE.)

! write all the single doubles
io_real = (/ ecpnl%dmin, ecpnl%startthick /)
reallist(1) = 'dmin'
reallist(2) = 'startthick' 
call HDF_writeNMLdbles(HDF_head, io_real, reallist, n_real)

! 3-vectors (real)
rnk = 1
dims(1) = 3
call h5ltmake_dataset_f(HDF_head%objectID, 'abcdist', rnk, dims, H5T_NATIVE_REAL, ecpnl%abcdist, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPMasterNameList: unable to create abcdist dataset',.TRUE.)

call h5ltmake_dataset_f(HDF_head%objectID, 'albegadist', rnk, dims, H5T_NATIVE_REAL, ecpnl%albegadist, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPMasterNameList: unable to create albegadist dataset',.TRUE.)


! write all the strings
call h5ltmake_dataset_string_f(HDF_head%objectID, 'outname', ecpnl%outname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPMasterNameList: unable to create outname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'energyfile', ecpnl%energyfile, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPMasterNameList: unable to create energyfile dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'compmode', ecpnl%compmode, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPMasterNameList: unable to create compmode dataset',.TRUE.)

! and pop this group off the stack
call HDF_pop(HDF_head)

end subroutine HDFwriteECPMasterNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:HDFwriteEBSDNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist to HDF file
!
!> @param HDF_head top of push stack
!> @param enl EBSD name list structure
!
!> @date 03/22/15 MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine HDFwriteEBSDNameList(HDF_head, HDF_tail, enl)



IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(EBSDNameListType),INTENT(INOUT)      :: enl

integer(HID_T)                                        :: grp_id
integer(kind=irg),parameter                           :: n_int = 6, n_real = 8
integer(kind=irg)                                     :: rnk, error,  io_int(n_int), distort
integer(HSIZE_T)                                      :: dims(1)
real(kind=sgl)                                        :: io_real(n_real)
real(kind=dbl)                                        :: t(1)
character(20)                                         :: intlist(n_int), reallist(n_real)

! create the group for this namelist
error = HDF_createGroup('EBSDNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ enl%stdout, enl%numsx, enl%numsy, enl%binning, enl%nthreads, enl%energyaverage /)
intlist(1) = 'stdout'
intlist(2) = 'numsx'
intlist(3) = 'numsy'
intlist(4) = 'binning'
intlist(5) = 'nthreads'
intlist(6) = 'energyaverage'
call HDF_writeNMLintegers(HDF_head, io_int, intlist, n_int)

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
call HDF_writeNMLreals(HDF_head, io_real, reallist, n_real)

! a 4-vector
rnk = 1
dims(1) = 4
call h5ltmake_dataset_f(HDF_head%objectID, 'axisangle', rnk, dims, H5T_NATIVE_REAL, enl%axisangle, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteEBSDNameList: unable to create axisangle dataset',.TRUE.)

! a few doubles
rnk = 1
dims(1) = 1
t(1) = enl%beamcurrent
call h5ltmake_dataset_f(HDF_head%objectID, 'beamcurrent', rnk, dims, H5T_NATIVE_DOUBLE, t, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteEBSDNameList: unable to create beamcurrent dataset',.TRUE.)

t(1) = enl%dwelltime
call h5ltmake_dataset_f(HDF_head%objectID, 'dwelltime', rnk, dims, H5T_NATIVE_DOUBLE, t, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteEBSDNameList: unable to create dwelltime dataset',.TRUE.)

! write all the strings
call h5ltmake_dataset_string_f(HDF_head%objectID, 'maskpattern', enl%maskpattern, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteEBSDNameList: unable to create maskpattern dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'scalingmode', enl%scalingmode, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteEBSDNameList: unable to create scalingmode dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'eulerconvention', enl%eulerconvention, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteEBSDNameList: unable to create eulerconvention dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'outputformat', enl%outputformat, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteEBSDNameList: unable to create outputformat dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'energyfile', enl%energyfile, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteEBSDNameList: unable to create energyfile dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'masterfile', enl%masterfile, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteEBSDNameList: unable to create masterfile dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'anglefile', enl%anglefile, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteEBSDNameList: unable to create anglefile dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'datafile', enl%datafile, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteEBSDNameList: unable to create datafile dataset',.TRUE.)

! and pop this group off the stack
call HDF_pop(HDF_head)

end subroutine HDFwriteEBSDNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:HDFwriteECPNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist to HDF file
!
!> @param HDF_head top of push stack
!> @param ecpnl ECP namelist structure
!
!> @date 03/22/15 MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine HDFwriteECPNameList(HDF_head, HDF_tail, ecpnl)



IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(ECPNameListType),INTENT(INOUT)                   :: ecpnl

integer(HID_T)                                        :: grp_id
integer(kind=irg),parameter                           :: n_int = 4, n_real = 8
integer(kind=irg)                                     :: rnk, error,  io_int(n_int), distort
integer(HSIZE_T)                                      :: dims(1)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)

! create the group for this namelist
error = HDF_createGroup('ECPNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ ecpnl%stdout, ecpnl%numthick, ecpnl%npix, ecpnl%nthreads /)
intlist(1) = 'stdout'
intlist(2) = 'numthick'
intlist(3) = 'npix'
intlist(4) = 'nthreads'
call HDF_writeNMLintegers(HDF_head, io_int, intlist, n_int)

! integer vectors
rnk = 1
dims(1) = 3
call h5ltmake_dataset_f(HDF_head%objectID, 'k', rnk, dims, H5T_NATIVE_INTEGER, ecpnl%k, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPNameList: unable to create k dataset',.TRUE.)

call h5ltmake_dataset_f(HDF_head%objectID, 'fn', rnk, dims, H5T_NATIVE_INTEGER, ecpnl%fn, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPNameList: unable to create fn dataset',.TRUE.)

call h5ltmake_dataset_f(HDF_head%objectID, 'gF', rnk, dims, H5T_NATIVE_INTEGER, ecpnl%gF, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPNameList: unable to create gF dataset',.TRUE.)

call h5ltmake_dataset_f(HDF_head%objectID, 'gS', rnk, dims, H5T_NATIVE_INTEGER, ecpnl%gS, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPNameList: unable to create gS dataset',.TRUE.)

call h5ltmake_dataset_f(HDF_head%objectID, 'tF', rnk, dims, H5T_NATIVE_INTEGER, ecpnl%tF, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPNameList: unable to create tF dataset',.TRUE.)

call h5ltmake_dataset_f(HDF_head%objectID, 'tS', rnk, dims, H5T_NATIVE_INTEGER, ecpnl%tS, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPNameList: unable to create tS dataset',.TRUE.)

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
call HDF_writeNMLreals(HDF_head, io_real, reallist, n_real)

! write all the strings
call h5ltmake_dataset_string_f(HDF_head%objectID, 'compmode', ecpnl%compmode, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPNameList: unable to create compmode dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'energyfile', ecpnl%energyfile, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPNameList: unable to create energyfile dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'outname', ecpnl%outname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPNameList: unable to create outname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'xtalname', ecpnl%xtalname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPNameList: unable to create xtalname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'xtalname2', ecpnl%xtalname2, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPNameList: unable to create xtalname2 dataset',.TRUE.)

! and pop this group off the stack
call HDF_pop(HDF_head)

end subroutine HDFwriteECPNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:HDFwriteLACBEDNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist to HDF file
!
!> @param HDF_head top of push stack
!> @param lacbednl LACBED name list structure
!
!> @date 06/22/15 MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine HDFwriteLACBEDNameList(HDF_head, HDF_tail, lacbednl)



IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(LACBEDNameListType),INTENT(INOUT)                :: lacbednl

integer(HID_T)                                        :: grp_id
integer(kind=irg),parameter                           :: n_int = 5, n_real = 6
integer(kind=irg)                                     :: rnk, error,  io_int(n_int), distort
integer(HSIZE_T)                                      :: dims(1)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)

! create the group for this namelist
error = HDF_createGroup('LACBEDNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ lacbednl%stdout, lacbednl%maxHOLZ, lacbednl%numthick, lacbednl%npix, lacbednl%nthreads /)
intlist(1) = 'stdout'
intlist(2) = 'maxHOLZ'
intlist(3) = 'numthick'
intlist(4) = 'npix'
intlist(5) = 'nthreads'
call HDF_writeNMLintegers(HDF_head, io_int, intlist, n_int)

! vectors
rnk = 1
dims(1) = 3
call h5ltmake_dataset_f(HDF_head%objectID, 'k', rnk, dims, H5T_NATIVE_INTEGER, lacbednl%k, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteLACBEDNameList: unable to create k dataset',.TRUE.)

call h5ltmake_dataset_f(HDF_head%objectID, 'fn', rnk, dims, H5T_NATIVE_INTEGER, lacbednl%fn, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteLACBEDNameList: unable to create fn dataset',.TRUE.)

! write all the single reals
io_real = (/ lacbednl%voltage, lacbednl%dmin, lacbednl%convergence, lacbednl%startthick, lacbednl%thickinc, lacbednl%minten/)
reallist(1) = 'voltage'
reallist(2) = 'dmin'
reallist(3) = 'convergence'
reallist(4) = 'startthick'
reallist(5) = 'thickinc'
reallist(6) = 'minten'
call HDF_writeNMLreals(HDF_head, io_real, reallist, n_real)

! write all the strings
call h5ltmake_dataset_string_f(HDF_head%objectID, 'outname', lacbednl%outname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteLACBEDNameList: unable to create outname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'xtalname', lacbednl%xtalname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteLACBEDNameList: unable to create xtalname dataset',.TRUE.)

! and pop this group off the stack
call HDF_pop(HDF_head)

end subroutine HDFwriteLACBEDNameList


!--------------------------------------------------------------------------
!
! SUBROUTINE:HDFwriteECPpatternNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist into HDF file
!
!> @param HDF_head top of push stack
!> @param ecpnl ECP name list structure
!
!> @date 03/22/15 MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine HDFwriteECPpatternNameList(HDF_head, HDF_tail,ecpnl)



IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(ECPpatternNameListType),INTENT(INOUT)            :: ecpnl

integer(HID_T)                                        :: grp_id
integer(kind=irg),parameter                           :: n_int = 2, n_real = 6
integer(kind=irg)                                     :: rnk, error,  io_int(n_int), distort
integer(HSIZE_T)                                      :: dims(1)
real(kind=sgl)                                        :: io_real(n_real), t(1)
character(20)                                         :: intlist(n_int), reallist(n_real)

! create the group for this namelist
error = HDF_createGroup('ECPpatternNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ ecpnl%stdout, ecpnl%npix /)
intlist(1) = 'stdout'
intlist(2) = 'npix'
call HDF_writeNMLintegers(HDF_head, io_int, intlist, n_int)

! single real
rnk = 1
dims(1) = 1
t(1) = ecpnl%thetac
call h5ltmake_dataset_f(HDF_head%objectID, 'thetac', rnk, dims, H5T_NATIVE_REAL, t, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPpatternNameList: unable to create thetac dataset',.TRUE.)

! real vector
dims(1) = 3
call h5ltmake_dataset_f(HDF_head%objectID, 'k', rnk, dims, H5T_NATIVE_REAL, ecpnl%k, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPpatternNameList: unable to create k dataset',.TRUE.)

! write all the strings
call h5ltmake_dataset_string_f(HDF_head%objectID, 'outname', ecpnl%outname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPpatternNameList: unable to create outname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'masterfile', ecpnl%masterfile, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPpatternNameList: unable to create masterfile dataset',.TRUE.)

! and pop this group off the stack
call HDF_pop(HDF_head)

end subroutine HDFwriteECPpatternNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:HDFwritePEDKINNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist to HDF file
!
!> @param HDF_head top of push stack
!> @param pednl PED name list structure
!
!> @date 03/22/15 MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine HDFwritePEDKINNameList(HDF_head, HDF_tail,pednl)



IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(PEDKINNameListType),INTENT(INOUT)                :: pednl

integer(HID_T)                                        :: grp_id
integer(kind=irg),parameter                           :: n_int = 4, n_real = 4
integer(kind=irg)                                     :: rnk, error,  io_int(n_int), distort
integer(HSIZE_T)                                      :: dims(1)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)

! create the group for this namelist
error = HDF_createGroup('PEDKINNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ pednl%stdout, pednl%npix, pednl%ncubochoric, pednl%nthreads /)
intlist(1) = 'stdout'
intlist(2) = 'npix'
intlist(3) = 'ncubochoric'
intlist(4) = 'nthreads'
call HDF_writeNMLintegers(HDF_head, io_int, intlist, n_int)

! write all the single reals
io_real = (/ pednl%voltage, pednl%thickness, pednl%dmin, pednl%rnmpp /)
reallist(1) = 'voltage'
reallist(2) = 'thickness'
reallist(3) = 'dmin'
reallist(4) = 'rnmpp'
call HDF_writeNMLreals(HDF_head, io_real, reallist, n_real)

! write all the strings
call h5ltmake_dataset_string_f(HDF_head%objectID, 'outname', pednl%outname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwritePEDKINNameList: unable to create outname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'xtalname', pednl%xtalname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwritePEDKINNameList: unable to create xtalname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'eulername', pednl%eulername, error)
if (error.ne.0) call HDF_handleError(error,'HDFwritePEDKINNameList: unable to create eulername dataset',.TRUE.)

! and pop this group off the stack
call HDF_pop(HDF_head)

end subroutine HDFwritePEDKINNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:HDFwritePEDNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist to HDF file
!
!> @param HDF_head top of push stack
!> @param pednl PED name list structure
!
!> @date 03/23/15 MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine HDFwritePEDNameList(HDF_head, HDF_tail,pednl)



IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(PEDNameListType),INTENT(INOUT)                   :: pednl

integer(HID_T)                                        :: grp_id
integer(kind=irg),parameter                           :: n_int = 5, n_real = 6
integer(kind=irg)                                     :: rnk, error,  io_int(n_int), distort
integer(HSIZE_T)                                      :: dims(1)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)

! create the group for this namelist
error = HDF_createGroup('PEDNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ pednl%stdout, pednl%precsample, pednl%precazimuthal, pednl%npix, pednl%nthreads /)
intlist(1) = 'stdout'
intlist(2) = 'precsample'
intlist(3) = 'precazimuthal'
intlist(4) = 'npix'
intlist(5) = 'nthreads'
call HDF_writeNMLintegers(HDF_head, io_int, intlist, n_int)

! vectors
rnk =1
dims(1) = 3
call h5ltmake_dataset_f(HDF_head%objectID, 'k', rnk, dims, H5T_NATIVE_INTEGER, pednl%k, error)
if (error.ne.0) call HDF_handleError(error,'HDFwritePEDNameList: unable to create k dataset',.TRUE.)

call h5ltmake_dataset_f(HDF_head%objectID, 'fn', rnk, dims, H5T_NATIVE_INTEGER, pednl%fn, error)
if (error.ne.0) call HDF_handleError(error,'HDFwritePEDNameList: unable to create fn dataset',.TRUE.)

! single reals
io_real = (/ pednl%voltage, pednl%dmin, pednl%precangle, pednl%prechalfwidth, pednl%thickness, pednl%camlen /)
reallist(1) = 'voltage'
reallist(2) = 'dmin'
reallist(3) = 'precangle'
reallist(4) = 'prechalfwidth'
reallist(5) = 'thickness'
reallist(6) = 'camlen'
call HDF_writeNMLreals(HDF_head, io_real, reallist, n_real)

! write all the strings
call h5ltmake_dataset_string_f(HDF_head%objectID, 'outname', pednl%outname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwritePEDNameList: unable to create outname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'xtalname', pednl%xtalname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwritePEDNameList: unable to create xtalname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'filemode', pednl%filemode, error)
if (error.ne.0) call HDF_handleError(error,'HDFwritePEDNameList: unable to create filemode dataset',.TRUE.)

! and pop this group off the stack
call HDF_pop(HDF_head)

end subroutine HDFwritePEDNameList


!--------------------------------------------------------------------------
!
! SUBROUTINE:HDFwriteECCINameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist to HDF file
!
!> @param HDF_head top of push stack
!> @param eccinl ECCI name list structure
!
!> @date 03/23/15 MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine HDFwriteECCINameList(HDF_head, HDF_tail,eccinl)



IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(ECCINameListType),INTENT(INOUT)                  :: eccinl

integer(HID_T)                                        :: grp_id
integer(HSIZE_T)                                      :: dims(1)
integer(kind=irg),parameter                           :: n_int = 8, n_real = 6
integer(kind=irg)                                     :: rnk, error,  io_int(n_int), distort
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
integer(kind=irg)                                     :: i

! create the group for this namelist
error = HDF_createGroup('ECCINameList',HDF_head, HDF_tail)

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
call HDF_writeNMLintegers(HDF_head, io_int, intlist, n_int)

! vectors
rnk =1
dims(1) = 3
call h5ltmake_dataset_f(HDF_head%objectID, 'k', rnk, dims, H5T_NATIVE_INTEGER, eccinl%k, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECCINameList: unable to create k dataset',.TRUE.)

! single reals
io_real = (/ eccinl%voltage, eccinl%dkt, eccinl%ktmax, eccinl%dmin, eccinl%DF_L, eccinl%DF_slice /)
reallist(1) = 'voltage'
reallist(2) = 'dkt'
reallist(3) = 'ktmax'
reallist(4) = 'dmin'
reallist(5) = 'DF_L'
reallist(6) = 'DF_slice'
call HDF_writeNMLreals(HDF_head, io_real, reallist, n_real)

! 2-vectors
dims(1) = 2
call h5ltmake_dataset_f(HDF_head%objectID, 'lauec', rnk, dims, H5T_NATIVE_REAL, eccinl%lauec, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECCINameList: unable to create lauec dataset',.TRUE.)

call h5ltmake_dataset_f(HDF_head%objectID, 'lauec2', rnk, dims, H5T_NATIVE_REAL, eccinl%lauec2, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECCINameList: unable to create lauec2 dataset',.TRUE.)

! write all the strings
call h5ltmake_dataset_string_f(HDF_head%objectID, 'dispmode', eccinl%dispmode, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECCINameList: unable to create dispmode dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'summode', eccinl%summode, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECCINameList: unable to create summode dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'progmode', eccinl%progmode, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECCINameList: unable to create progmode dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'xtalname', eccinl%xtalname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECCINameList: unable to create xtalname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'foilnmlfile', eccinl%foilnmlfile, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECCINameList: unable to create foilnmlfile dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'dispfile', eccinl%dispfile, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECCINameList: unable to create dispfile dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'dataname', eccinl%dataname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECCINameList: unable to create dataname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'ECPname', eccinl%ECPname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECCINameList: unable to create ECPname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'sgname', eccinl%sgname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECCINameList: unable to create sgname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'apbname', eccinl%apbname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECCINameList: unable to create apbname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'incname', eccinl%incname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECCINameList: unable to create incname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'voidname', eccinl%voidname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECCINameList: unable to create voidname dataset',.TRUE.)

! THESE NEED WORK; PERHAPS WE NEED TO REDEFINE ONE OF THE WRAPPER ROUTINES ?

! maxdefects string arrays
!call h5ltmake_dataset_string_f(HDF_head%objectID, 'sfname', eccinl%sfname, error)
!if (error.ne.0) call HDF_handleError(error,'HDFwriteECCINameList: unable to create sfname dataset',.TRUE.)

! 3*maxdefects string arrays
!call h5ltmake_dataset_string_f(HDF_head%objectID, 'dislYname', eccinl%dislYname, error)
!if (error.ne.0) call HDF_handleError(error,'HDFwriteECCINameList: unable to create dislYname dataset',.TRUE.)

!call h5ltmake_dataset_string_f(HDF_head%objectID, 'dislname', eccinl%dislname, error)
!if (error.ne.0) call HDF_handleError(error,'HDFwriteECCINameList: unable to create dislname dataset',.TRUE.)

! and pop this group off the stack
call HDF_pop(HDF_head)

end subroutine HDFwriteECCINameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:HDFwriteRFZNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist to an HDF file
!
!> @param HDF_head top of push stack
!> @param rfznl RFZ name list structure
!
!> @date 03/23/15 MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine HDFwriteRFZNameList(HDF_head, HDF_tail,rfznl)



IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(RFZNameListType),INTENT(INOUT)                   :: rfznl

integer(HID_T)                                        :: grp_id
integer(kind=irg),parameter                           :: n_int = 2, n_real = 1
integer(kind=irg)                                     :: rnk, error,  io_int(n_int), distort
integer(HSIZE_T)                                      :: dims(1)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
integer(kind=irg)                                     :: i

! create the group for this namelist
error = HDF_createGroup('RFZNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ rfznl%pgnum, rfznl%nsteps /)
intlist(1) = 'pgnum'
intlist(2) = 'nsteps'
call HDF_writeNMLintegers(HDF_head, io_int, intlist, n_int)

! strings
call h5ltmake_dataset_string_f(HDF_head%objectID, 'outname', rfznl%outname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteRFZNameList: unable to create outname dataset',.TRUE.)

! and pop this group off the stack
call HDF_pop(HDF_head)

end subroutine HDFwriteRFZNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:HDFwriteDictIndxOpenCLNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write namelist to HDF file
!
!> @param HDF_head top of push stack
!> @param rfznl RFZ name list structure
!
!> @date 03/23/15 MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine HDFwriteDictIndxOpenCLNameList(HDF_head, HDF_tail,dictindxnl)


use local

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(DictIndxOpenCLListType),INTENT(INOUT)            :: dictindxnl

integer(HID_T)                                        :: grp_id
integer(kind=irg),parameter                           :: n_int = 8, n_real = 1
integer(kind=irg)                                     :: rnk, error,  io_int(n_int), distort
integer(HSIZE_T)                                      :: dims(1)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
integer(kind=irg)                                     :: i

! create the group for this namelist
error = HDF_createGroup('DictIndxOpenCLNameList',HDF_head, HDF_tail)

! logical will be written as an integer 1 or 0
if (dictindxnl%MeanSubtraction) then 
  i = 1
else
  i = 0 
end if
! write all the single integers
io_int = (/ dictindxnl%numexptsingle, dictindxnl%numdictsingle, dictindxnl%totnumdict, dictindxnl%totnumexpt, dictindxnl%imght, &
            dictindxnl%imgwd, dictindxnl%nnk, i /)
intlist(1) = 'numexptsingle'
intlist(2) = 'numdictsingle'
intlist(3) = 'totnumdict'
intlist(4) = 'totnumexpt'
intlist(5) = 'imght'
intlist(6) = 'imgwd'
intlist(7) = 'nnk'
intlist(8) = 'MeanSubtraction'
call HDF_writeNMLintegers(HDF_head, io_int, intlist, n_int)

! strings
call h5ltmake_dataset_string_f(HDF_head%objectID, 'exptfile', dictindxnl%exptfile, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteDictIndxOpenCLNameList: unable to create exptfile dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'dictfile', dictindxnl%dictfile, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteDictIndxOpenCLNameList: unable to create dictfile dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%objectID, 'eulerfile', dictindxnl%eulerfile, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteDictIndxOpenCLNameList: unable to create eulerfile dataset',.TRUE.)

! and pop this group off the stack
call HDF_pop(HDF_head)

end subroutine HDFwriteDictIndxOpenCLNameList

end module NameListHDFwriters
