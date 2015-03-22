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
! CTEMsoft2013:NameListHDFwriters.f90
!--------------------------------------------------------------------------
!
! PROGRAM: NameListHDFwriters
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief routines for reading and returning name list type structures
!
!> @date 06/13/14 MDG 1.0 original
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

use error

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head
integer(kind=irg),INTENT(IN)                          :: io_int(n_int)
character(20),INTENT(IN)                              :: intlist(n_int)
integer(kind=irg),INTENT(IN)                          :: n_int

integer(kind=irg)                                     :: rnk, dims(1), error, ioval(1)

rnk = 1
dims(1) = 1

do i=1,n_int
  ioval(1) = io_int(i)
  call h5ltmake_dataset_f(HDF_head%oID, trim(intlist(i)), rnk, dims, H5T_NATIVE_INTEGER, ioval, error)
  if (error.ne.0) call HDF_handleError(error,'HDF_writeNMLintegers: unable to create'//trim(intlist(i))//' dataset',.TRUE.)
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

use error

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head
real(kind=sgl),INTENT(IN)                             :: io_real(n_real)
character(20),INTENT(IN)                              :: reallist(n_real)
integer(kind=irg),INTENT(IN)                          :: n_real

integer(kind=irg)                                     :: rnk, dims(1), error
real(kind=sgl)                                        :: ioval(1)

rnk = 1
dims(1) = 1

do i=1,n_real
  ioval(1) = io_real(i)
  call h5ltmake_dataset_f(HDF_head%oID, trim(reallist(i)), rnk, dims, H5T_NATIVE_REAL, ioval, error)
  if (error.ne.0) call HDF_handleError(error,'HDF_writeNMLreals: unable to create'//trim(reallist(i))//' dataset',.TRUE.)
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

use error

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head
real(kind=dbl),INTENT(IN)                             :: io_real(n_real)
character(20),INTENT(IN)                              :: reallist(n_real)
integer(kind=irg),INTENT(IN)                          :: n_real

integer(kind=irg)                                     :: rnk, dims(1), error
real(kind=dbl)                                        :: ioval(1)

rnk = 1
dims(1) = 1

do i=1,n_real
  ioval(1) = io_real(i)
  call h5ltmake_dataset_f(HDF_head%oID, trim(reallist(i)), rnk, dims, H5T_NATIVE_DOUBLE, ioval, error)
  if (error.ne.0) call HDF_handleError(error,'HDF_writeNMLdbles: unable to create'//trim(reallist(i))//' dataset',.TRUE.)
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
subroutine HDFwriteKosselNameList(HDF_head, knl)

use error

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head
type(KosselNameListType),INTENT(IN)                   :: knl

integer(HID_T)                                        :: grp_id
integer(kind=irg),parameter                           :: n_int = 5, n_real = 6
integer(kind=irg)                                     :: rnk, error, dims(1), io_int(n_int)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)

! create the group for this namelist
call HDF_createGroup('KosselNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ knl%stdout, knl%numthick, knl%npix, knl%maxHOLZ, knl%nthreads /)
intlist = (/ 'stdout', 'numthick', 'npix', 'maxHOLZ', 'nthreads' /)
call HDF_writeNMLintegers(HDF_head, io_int, intlist, n_int)

! integer vectors
rnk = 1
dims(1) = 3
call h5ltmake_dataset_f(HDF_head%oID, 'k', rnk, dims, H5T_NATIVE_INTEGER, knl%k, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteKosselNameList: unable to create k dataset',.TRUE.)

call h5ltmake_dataset_f(HDF_head%oID, 'fn', rnk, dims, H5T_NATIVE_INTEGER, knl%fn, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteKosselNameList: unable to create fn dataset',.TRUE.)

! write all the single reals
io_real = (/ knl%voltage, knl%dmin, knl%convergence, knl%startthick, knl%thickinc, knl%minten /)
reallist = (/ 'voltage', 'dmin', 'convergence', 'startthick', 'thickinc', 'minten' /)
call HDF_writeNMLreals(HDF_head, io_real, reallist, n_real)

! write all the strings
call h5ltmake_dataset_string_f(HDF_head%oID, 'xtalname', knl%xtalname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteKosselNameList: unable to create xtalname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'outname', knl%outname, error)
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
subroutine HDFwriteKosselMasterNameList(HDF_head, knl)

use error

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head
type(KosselMasterNameListType),INTENT(IN)             :: knl

integer(HID_T)                                        :: grp_id
integer(kind=irg),parameter                           :: n_int = 4, n_real = 5
integer(kind=irg)                                     :: rnk, error, dims(1), io_int(n_int)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)

! create the group for this namelist
call HDF_createGroup('KosselMasterNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ knl%stdout, knl%numthick, knl%npix, knl%nthreads /)
intlist = (/ 'stdout', 'numthick', 'npix', 'nthreads' /) 
call HDF_writeNMLintegers(HDF_head, io_int, intlist, n_int)

! write all the single reals
io_real = (/ knl%voltage, knl%dmin, knl%startthick, knl%thickinc, knl%tfraction /)
reallist = (/ 'voltage', 'dmin', 'startthick', 'thickinc', 'tfraction' /) 
call HDF_writeNMLreals(HDF_head, io_real, reallist, n_real)

! write all the strings
call h5ltmake_dataset_string_f(HDF_head%oID, 'Kosselmode', knl%Kosselmode, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteKosselMasterNameList: unable to create Kosselmode dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'xtalname', knl%xtalname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteKosselMasterNameList: unable to create xtalname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'outname', knl%outname, error)
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
subroutine HDFwriteMCNameList(HDF_head, mcnl)

use error

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head
type(MCNameListType),INTENT(INOUT)                    :: mcnl

integer(HID_T)                                        :: grp_id
integer(kind=irg),parameter                           :: n_int = 5, n_real = 7
integer(kind=irg)                                     :: rnk, error, dims(1), io_int(n_int)
real(kind=dbl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)

! create the group for this namelist
call HDF_createGroup('MCNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ mcnl%stdout, mcnl%numsx, mcnl%primeseed, mcnl%num_el, mcnl%nthreads /)
intlist = (/ 'stdout', 'numsx', 'primeseed', 'num_el', 'nthreads' /)
call HDF_writeNMLintegers(HDF_head, io_int, intlist, n_int)

! write all the single doubles
io_real = (/ mcnl%sig, mcnl%omega, mcnl%EkeV, mcnl%Ehistmin, mcnl%Ebinsize, mcnl%depthmax, mcnl%depthstep /)
reallist = (/ 'sig', 'omega', 'EkeV', 'Ehistmin', 'Ebinsize', 'depthmax', 'depthstep' /)
call HDF_writeNMLdbles(HDF_head, io_real, reallist, n_real)

! write all the strings
call h5ltmake_dataset_string_f(HDF_head%oID, 'MCmode', mcnl%MCmode, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteMCNameList: unable to create MCmode dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'xtalname', mcnl%xtalname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteMCNameList: unable to create xtalname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'dataname', mcnl%dataname, error)
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
subroutine HDFwriteMCCLNameList(nmlfile, mcnl)

use error

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head
type(MCCLNameListType),INTENT(INOUT)                  :: mcnl

integer(HID_T)                                        :: grp_id
integer(kind=irg),parameter                           :: n_int = 5, n_real = 7
integer(kind=irg)                                     :: rnk, error, dims(1), io_int(n_int)
real(kind=dbl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)

! create the group for this namelist
call HDF_createGroup('MCCLNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ mcnl%stdout, mcnl%numsx, mcnl%globalworkgrpsz, mcnl%num_el, mcnl%totnum_el /)
intlist = (/ 'stdout', 'numsx', 'globalworkgrpsz', 'num_el', 'totnum_el' /)
call HDF_writeNMLintegers(HDF_head, io_int, intlist, n_int)

! write all the single doubles
io_real = (/ mcnl%sig, mcnl%omega, mcnl%EkeV, mcnl%Ehistmin, mcnl%Ebinsize, mcnl%depthmax, mcnl%depthstep /)
reallist = (/ 'sig', 'omega', 'EkeV', 'Ehistmin', 'Ebinsize', 'depthmax', 'depthstep' /)
call HDF_writeNMLdbles(HDF_head, io_real, reallist, n_real)

! write all the strings
call h5ltmake_dataset_string_f(HDF_head%oID, 'MCmode', mcnl%MCelmode, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteMCCLNameList: unable to create MCmode dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'xtalname', mcnl%xtalname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteMCCLNameList: unable to create xtalname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'dataname', mcnl%dataname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteMCCLNameList: unable to create dataname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'primelist', mcnl%primelist, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteMCCLNameList: unable to create primelist dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'mode', mcnl%mode, error)
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
subroutine HDFwriteMCCLMultiLayerNameList(HDF_head, mcnl)

use error

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head
type(MCCLMultiLayerNameListType),INTENT(INOUT)        :: mcnl

integer(HID_T)                                        :: grp_id
integer(kind=irg),parameter                           :: n_int = 5, n_real = 9
integer(kind=irg)                                     :: rnk, error, dims(1), io_int(n_int)
real(kind=dbl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)

! create the group for this namelist
call HDF_createGroup('MCCLMultiLayerNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ mcnl%stdout, mcnl%numsx, mcnl%globalworkgrpsz, mcnl%num_el, mcnl%totnum_el /)
intlist = (/ 'stdout', 'numsx', 'globalworkgrpsz', 'num_el', 'totnum_el' /)
call HDF_writeNMLintegers(HDF_head, io_int, intlist, n_int)

! write all the single doubles
io_real = (/ mcnl%sig, mcnl%omega, mcnl%EkeV, mcnl%Ehistmin, mcnl%Ebinsize, mcnl%depthmax, mcnl%depthstep, &
             mcnl%filmthickness, mcnl%filmstep /)
reallist = (/ 'sig', 'omega', 'EkeV', 'Ehistmin', 'Ebinsize', 'depthmax', 'depthstep', 'filmthickness', 'filmstep'
call HDF_writeNMLdbles(HDF_head, io_real, reallist, n_real)

! write all the strings
call h5ltmake_dataset_string_f(HDF_head%oID, 'MCmode', mcnl%MCelmode, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteMCCLMultiLayerNameList: unable to create MCmode dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'xtalname_film', mcnl%xtalname_film, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteMCCLMultiLayerNameList: unable to create xtalname_film dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'xtalname_subs', mcnl%xtalname_subs, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteMCCLMultiLayerNameList: unable to create xtalname_subs dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'dataname', mcnl%dataname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteMCCLMultiLayerNameList: unable to create dataname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'primelist', mcnl%primelist, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteMCCLMultiLayerNameList: unable to create primelist dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'mode', mcnl%mode, error)
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
subroutine HDFwriteEBSDMasterNameList(HDF_head, emnl)

use error

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head
type(EBSDMasterNameListType),INTENT(INOUT)            :: emnl

integer(HID_T)                                        :: grp_id
integer(kind=irg),parameter                           :: n_int = 5, n_real = 9
integer(kind=irg)                                     :: rnk, error, dims(1), io_int(n_int)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)

! create the group for this namelist
call HDF_createGroup('EBSDMasterNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ emnl%stdout, emnl%npx, emnl%Esel, emnl%nthreads /)
intlist = (/ 'stdout', 'npx', 'Esel', 'nthreads' /)
call HDF_writeNMLintegers(HDF_head, io_int, intlist, n_int)

! write a single real
io_real(1) = emnl%dmin
rnk = 1
dims(1) = 1
call h5ltmake_dataset_f(HDF_head%oID, 'dmin', rnk, dims, H5T_NATIVE_REAL, ioval, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteEBSDMasterNameList: unable to create dmin dataset',.TRUE.)

! write all the strings
call h5ltmake_dataset_string_f(HDF_head%oID, 'outname', emnl%outname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteEBSDMasterNameList: unable to create outname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'energyfile', emnl%energyfile, error)
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
subroutine HDFwriteECPMasterNameList(HDF_head, ecpnl)

use error

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head
type(ECPMasterNameListType),INTENT(INOUT)             :: ecpnl

integer(HID_T)                                        :: grp_id
integer(kind=irg),parameter                           :: n_int = 5, n_real = 3
integer(kind=irg)                                     :: rnk, error, dims(1), io_int(n_int), distort
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)

! create the group for this namelist
call HDF_createGroup('ECPMasterNameList',HDF_head, HDF_tail)

! write all the single integers
! distort is a logical, for which there is no real HDF_NATIVE_anything conversion, so we'll store it as a 1 or 0
if (ecpnl%distort) then 
  distort = 1
else 
  distort = 0
end if
io_int = (/ ecpnl%stdout, ecpnl%Esel, ecpnl%npx, ecpnl%startthick, distort /)
intlist = (/ 'stdout', 'Esel', 'npx', 'startthick', 'distort' /)
call HDF_writeNMLintegers(HDF_head, io_int, intlist, n_int)

! integer vectors
rnk = 1
dims(1) = 3
call h5ltmake_dataset_f(HDF_head%oID, 'fn', rnk, dims, H5T_NATIVE_INTEGER, ecpnl%fn, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPMasterNameList: unable to create fn dataset',.TRUE.)

! write all the single doubles
io_real = (/ ecpnl%abcdist, ecpnl%albegadist, ecpnl%dmin /)
reallist = (/ 'abcdist', 'albegadist', 'dmin' /)
call HDF_writeNMLdbles(HDF_head, io_real, reallist, n_real)

! write all the strings
call h5ltmake_dataset_string_f(HDF_head%oID, 'outname', ecpnl%outname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPMasterNameList: unable to create outname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'energyfile', ecpnl%energyfile, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPMasterNameList: unable to create energyfile dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'compmode', ecpnl%compmode, error)
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
subroutine HDFwriteEBSDNameList(HDF_head, enl)

use error

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head
type(EBSDNameListType),INTENT(INOUT)      :: enl

integer(HID_T)                                        :: grp_id
integer(kind=irg),parameter                           :: n_int = 6, n_real = 9
integer(kind=irg)                                     :: rnk, error, dims(1), io_int(n_int), distort
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)

! create the group for this namelist
call HDF_createGroup('EBSDNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ enl%stdout, enl%numsx, enl%numsy, enl%binning, enl%nthreads, enl%energyaverage /)
inlist = (/ 'stdout', 'numsx', 'numsy', 'binning', 'nthreads', 'energyaverage' /)
call HDF_writeNMLintegers(HDF_head, io_int, intlist, n_int)

! write all the single reals 
io_real = (/ enl%L, enl%thetac, enl%delta, enl%xpc, enl%ypc, enl%energymin, enl%energymax, enl%gammavalue, enl%axisangle /)
reallist = (/ 'L', 'thetac', 'delta', 'xpc', 'ypc', 'energymin', 'energymax', 'gammavalue', 'axisangle' /)
call HDF_writeNMLreals(HDF_head, io_real, reallist, n_real)

! a few doubles
rnk = 1
dims(1) = 1
call h5ltmake_dataset_f(HDF_head%oID, 'beamcurrent', rnk, dims, H5T_NATIVE_DOUBLE, enl%beamcurrent, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteEBSDNameList: unable to create beamcurrent dataset',.TRUE.)

call h5ltmake_dataset_f(HDF_head%oID, 'dwelltime', rnk, dims, H5T_NATIVE_DOUBLE, enl%dwelltime, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteEBSDNameList: unable to create dwelltime dataset',.TRUE.)

! write all the strings
call h5ltmake_dataset_string_f(HDF_head%oID, 'maskpattern', enl%maskpattern, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteEBSDNameList: unable to create maskpattern dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'scalingmode', enl%scalingmode, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteEBSDNameList: unable to create scalingmode dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'eulerconvention', enl%eulerconvention, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteEBSDNameList: unable to create eulerconvention dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'outputformat', enl%outputformat, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteEBSDNameList: unable to create outputformat dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'energyfile', ecpnl%energyfile, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteEBSDNameList: unable to create energyfile dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'masterfile', ecpnl%masterfile, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteEBSDNameList: unable to create masterfile dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'anglefile', ecpnl%anglefile, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteEBSDNameList: unable to create anglefile dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'datafile', ecpnl%datafile, error)
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
subroutine HDFwriteECPNameList(HDF_head, ecpnl)

use error

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head
type(ECPNameListType),INTENT(INOUT)                   :: ecpnl

integer(HID_T)                                        :: grp_id
integer(kind=irg),parameter                           :: n_int = 4, n_real = 8
integer(kind=irg)                                     :: rnk, error, dims(1), io_int(n_int), distort
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)

! create the group for this namelist
call HDF_createGroup('ECPNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ ecpnl%stdout, ecpnl%numthick, ecpnl%npix, ecpnl%nthreads /)
intlist = (/ 'stdout', 'numthick', 'npix', 'nthreads' /)
call HDF_writeNMLintegers(HDF_head, io_int, intlist, n_int)

! integer vectors
rnk = 1
dims(1) = 3
call h5ltmake_dataset_f(HDF_head%oID, 'k', rnk, dims, H5T_NATIVE_INTEGER, ecpnl%k, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPNameList: unable to create k dataset',.TRUE.)

call h5ltmake_dataset_f(HDF_head%oID, 'fn', rnk, dims, H5T_NATIVE_INTEGER, ecpnl%fn, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPNameList: unable to create fn dataset',.TRUE.)

call h5ltmake_dataset_f(HDF_head%oID, 'gF', rnk, dims, H5T_NATIVE_INTEGER, ecpnl%gF, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPNameList: unable to create gF dataset',.TRUE.)

call h5ltmake_dataset_f(HDF_head%oID, 'gS', rnk, dims, H5T_NATIVE_INTEGER, ecpnl%gS, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPNameList: unable to create gS dataset',.TRUE.)

call h5ltmake_dataset_f(HDF_head%oID, 'tF', rnk, dims, H5T_NATIVE_INTEGER, ecpnl%tF, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPNameList: unable to create tF dataset',.TRUE.)

call h5ltmake_dataset_f(HDF_head%oID, 'tS', rnk, dims, H5T_NATIVE_INTEGER, ecpnl%tS, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPNameList: unable to create tS dataset',.TRUE.)

! write all the single reals
io_real = (/ ecpnl%voltage, ecpnl%dmin, ecpnl%ktmax, ecpnl%thetac, ecpnl%startthick, ecpnl%thickinc, ecpnl%zintstep, &
             ecpnl%filmthickness /)
reallist = (/ 'voltage', 'dmin', 'ktmax', 'thetac', 'startthick', 'thickinc', 'zintstep', 'filmthickness' /)
call HDF_writeNMLreals(HDF_head, io_real, reallist, n_real)

! write all the strings
call h5ltmake_dataset_string_f(HDF_head%oID, 'compmode', ecpnl%compmode, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPNameList: unable to create compmode dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'energyfile', ecpnl%energyfile, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPNameList: unable to create energyfile dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'outname', ecpnl%outname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPNameList: unable to create outname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'xtalname', ecpnl%xtalname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPNameList: unable to create xtalname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'xtalname2', ecpnl%xtalname2, error)
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
subroutine HDFwriteLACBEDNameList(HDF_head, lacbednl)

use error

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head
type(LACBEDNameListType),INTENT(INOUT)                :: lacbednl

integer(HID_T)                                        :: grp_id
integer(kind=irg),parameter                           :: n_int = 5, n_real = 6
integer(kind=irg)                                     :: rnk, error, dims(1), io_int(n_int), distort
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)

! create the group for this namelist
call HDF_createGroup('LACBEDNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ lacbednl%stdout, lacbednl%maxHOLZ, lacbednl%numthick, lacbednl%npix, lacbednl%nthreads /)
intlist = (/ 'stdout', 'maxHOLZ', 'numthick', 'npix', 'nthreads' /)
call HDF_writeNMLintegers(HDF_head, io_int, intlist, n_int)

! vectors
rnk = 1
dims(1) = 3
call h5ltmake_dataset_f(HDF_head%oID, 'k', rnk, dims, H5T_NATIVE_INTEGER, lacbednl%k, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteLACBEDNameList: unable to create k dataset',.TRUE.)

call h5ltmake_dataset_f(HDF_head%oID, 'fn', rnk, dims, H5T_NATIVE_INTEGER, lacbednl%fn, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteLACBEDNameList: unable to create fn dataset',.TRUE.)

! write all the single reals
io_real = (/ lacbednl%voltage, lacbednl%dmin, lacbednl%convergence, lacbednl%startthick, lacbednl%thickinc, lacbednl%minten/)
reallist = (/ 'voltage', 'dmin', 'convergence', 'startthick', 'thickinc', 'minten' /)
call HDF_writeNMLreals(HDF_head, io_real, reallist, n_real)

! write all the strings
call h5ltmake_dataset_string_f(HDF_head%oID, 'outname', lacbednl%outname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteLACBEDNameList: unable to create outname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'xtalname', lacbednl%xtalname, error)
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
subroutine HDFwriteECPpatternNameList(HDF_head,ecpnl)

use error

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head
type(ECPpatternNameListType),INTENT(INOUT)            :: ecpnl

integer(HID_T)                                        :: grp_id
integer(kind=irg),parameter                           :: n_int = 2, n_real = 6
integer(kind=irg)                                     :: rnk, error, dims(1), io_int(n_int), distort
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)

! create the group for this namelist
call HDF_createGroup('ECPpatternNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ ecpnl%stdout, ecpnl%npix /)
intlist = (/ 'stdout', 'npix' /)
call HDF_writeNMLintegers(HDF_head, io_int, intlist, n_int)

! single real
rnk = 1
dims(1) = 1
call h5ltmake_dataset_f(HDF_head%oID, 'thetac', rnk, dims, H5T_NATIVE_REAL, ecpnl%thetac, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPpatternNameList: unable to create thetac dataset',.TRUE.)

! real vector
dims(1) = 3
call h5ltmake_dataset_f(HDF_head%oID, 'k', rnk, dims, H5T_NATIVE_REAL, ecpnl%k, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPpatternNameList: unable to create k dataset',.TRUE.)

! write all the strings
call h5ltmake_dataset_string_f(HDF_head%oID, 'outname', ecpnl%outname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteECPpatternNameList: unable to create outname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'masterfile', ecpnl%masterfile, error)
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
subroutine HDFwritePEDKINNameList(HDF_head,pednl)

use error

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head
type(PEDKINNameListType),INTENT(INOUT)                :: pednl

integer(HID_T)                                        :: grp_id
integer(kind=irg),parameter                           :: n_int = 4, n_real = 4
integer(kind=irg)                                     :: rnk, error, dims(1), io_int(n_int), distort
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)

! create the group for this namelist
call HDF_createGroup('PEDKINNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ pednl%stdout, pednl%npix, pednl%ncubochoric, pednl%nthreads /)
intlist = (/ 'stdout', 'npix', 'ncubochoric', 'nthreads' /)
call HDF_writeNMLintegers(HDF_head, io_int, intlist, n_int)

! write all the single reals
io_real = (/ pednl%voltage, pednl%thickness, pednl%dmin, pednl%rnmpp /)
reallist = (/ 'voltage', 'thickness', 'dmin', 'rnmpp' /)
call HDF_writeNMLreals(HDF_head, io_real, reallist, n_real)

! write all the strings
call h5ltmake_dataset_string_f(HDF_head%oID, 'outname', pednl%outname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwritePEDKINNameList: unable to create outname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'xtalname', pednl%xtalname, error)
if (error.ne.0) call HDF_handleError(error,'HDFwritePEDKINNameList: unable to create xtalname dataset',.TRUE.)

call h5ltmake_dataset_string_f(HDF_head%oID, 'eulername', pednl%eulername, error)
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
!> @brief read namelist file and fill pednl structure (used by CTEMPED.f90)
!
!> @param nmlfile namelist file name
!> @param pednl PED name list structure
!
!> @date 07/09/14 MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine HDFwritePEDNameList(nmlfile,pednl)

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
integer(kind=irg)       :: nthreads
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
                              thickness,  outname, npix, camlen, filemode, nthreads

! set the input parameters to default values (except for xtalname, which must be present)
xtalname = 'undefined'          ! initial value to check that the keyword is present in the nml file
stdout = 6                      ! standard output
voltage = 200000.0              ! acceleration voltage [V]
k = (/ 0, 0, 1 /)               ! beam direction [direction indices]
fn = (/ 0, 0, 1 /)              ! foil normal [direction indices]
dmin = 0.025                    ! smallest d-spacing to include in dynamical matrix [nm]
precangle = 10.472              ! beam precession angle [mrad]; default = 0.6 degrees
prechalfwidth = 0.25            ! beam half width in the tilt direction [mrad]
nthreads = 1                    ! number of OpenMP threads to start
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
pednl%nthreads = nthreads
pednl%outname = outname
pednl%camlen = camlen

end subroutine HDFwritePEDNameList


!--------------------------------------------------------------------------
!
! SUBROUTINE:HDFwriteECCINameList
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
subroutine HDFwriteECCINameList(nmlfile,eccinl)

use error

IMPLICIT NONE

character(fnlen),INTENT(IN)                     :: nmlfile
type(ECCINameListType),INTENT(INOUT)            :: eccinl

integer(kind=irg)                               :: i

integer(kind=irg)       :: stdout
integer(kind=irg)       :: nthreads
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

end subroutine HDFwriteECCINameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:HDFwriteRFZNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read namelist file and fill rfznl structure (used by CTEMsampleRFZ.f90)
!
!> @param nmlfile namelist file name
!> @param rfznl RFZ name list structure
!
!> @date 12/09/14 MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine HDFwriteRFZNameList(nmlfile,rfznl)

use error

IMPLICIT NONE

character(fnlen),INTENT(IN)                     :: nmlfile
type(RFZNameListType),INTENT(INOUT)             :: rfznl

integer(kind=irg)                               :: pgnum, nsteps
character(fnlen)                                :: outname

! namelist components
namelist / RFZlist / pgnum, nsteps, outname

! initialize to default values
pgnum = 32
nsteps = 50
outname = 'anglefile.txt'

! read the namelist file
open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
read(UNIT=dataunit,NML=RFZlist)
close(UNIT=dataunit,STATUS='keep')

! and copy the variables to the rfznl variable
rfznl%pgnum  = pgnum
rfznl%nsteps = nsteps
rfznl%outname= outname

end subroutine HDFwriteRFZNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:HDFwriteDictIndxOpenCLNameList
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief read namelist file and fill DictIndxOpenCLListType (used by CTEMDictIndxOpenCL.f90)
!
!> @param nmlfile namelist file name
!> @param rfznl RFZ name list structure
!
!> @date 13/01/15 SS 1.0 new routine
!--------------------------------------------------------------------------
subroutine HDFwriteDictIndxOpenCLNameList(nmlfile,dictindxnl)

use error
use local

IMPLICIT NONE

character(fnlen),INTENT(IN)                                 :: nmlfile
type(DictIndxOpenCLListType),INTENT(INOUT)                  :: dictindxnl

integer(kind=irg)                                           :: numexptsingle
integer(kind=irg)                                           :: numdictsingle
integer(kind=irg)                                           :: totnumexpt
integer(kind=irg)                                           :: totnumdict
integer(kind=irg)                                           :: imght
integer(kind=irg)                                           :: imgwd
integer(kind=irg)                                           :: nnk
character(fnlen)                                            :: exptfile
character(fnlen)                                            :: dictfile
character(fnlen)                                            :: eulerfile
logical                                                     :: MeanSubtraction

! define the IO namelist to facilitate passing variables to the program.
namelist /DictIndxOpenCLvars/ numexptsingle, numdictsingle, totnumexpt, totnumdict,&
        imght, imgwd, exptfile, dictfile, eulerfile, nnk, MeanSubtraction

! set some of the input parameters to default values 
numdictsingle = 1024
numexptsingle = 1024
imght = 0
imgwd = 0
nnk = 40
exptfile = 'undefined'
dictfile = 'undefined'
eulerfile = 'undefined'
totnumdict = 0
totnumexpt = 0
MeanSubtraction = .TRUE.

! read the namelist file
open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
read(UNIT=dataunit,NML=DictIndxOpenCLvars)
close(UNIT=dataunit,STATUS='keep')

! check for required entries
if (trim(exptfile).eq.'undefined') then
    call FatalError('CTEMDictIndxOpenCL:',' experimental file name is undefined in '//nmlfile)
end if

if (trim(dictfile).eq.'undefined') then
    call FatalError('CTEMDictIndxOpenCL:',' dictionary file name is undefined in '//nmlfile)
end if

if (trim(eulerfile).eq.'undefined') then
    call FatalError('CTEMDictIndxOpenCL:',' euler angle file name is undefined in '//nmlfile)
end if

if (totnumexpt .eq. 0) then
    call FatalError('CTEMDictIndxOpenCL:',' total number of experimental patterns is undefined in '//nmlfile)
end if

if (totnumdict .eq. 0) then
    call FatalError('CTEMDictIndxOpenCL:',' total number of dictionary patterns is undefined in '//nmlfile)
end if

if (imght .eq. 0) then
    call FatalError('CTEMDictIndxOpenCL:',' height of single pattern is undefined in '//nmlfile)
end if

if (imgwd .eq. 0) then
    call FatalError('CTEMDictIndxOpenCL:',' width of single pattern is undefined in '//nmlfile)
end if

! if we get here, then all appears to be ok, and we need to fill in the emnl fields

dictindxnl%numexptsingle = numexptsingle
dictindxnl%numdictsingle = numdictsingle
dictindxnl%imght = imght
dictindxnl%imgwd = imgwd
dictindxnl%exptfile = exptfile
dictindxnl%dictfile = dictfile
dictindxnl%eulerfile = eulerfile
dictindxnl%totnumdict = totnumdict
dictindxnl%totnumexpt = totnumexpt
dictindxnl%nnk = nnk
dictindxnl%MeanSubtraction = MeanSubtraction

end subroutine HDFwriteDictIndxOpenCLNameList

end module NameListHDFwriters
