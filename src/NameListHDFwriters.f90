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
call h5gcreate_f(HDF_head%oID,'KosselNameList',grp_id,error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteKosselNameList: unable to open KosselNameList group',.TRUE.)

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
call h5gcreate_f(HDF_head%oID,'KosselMasterNameList',grp_id,error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteKosselMasterNameList: unable to open KosselMasterNameList group',.TRUE.)

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
call h5gcreate_f(HDF_head%oID,'MCNameList',grp_id,error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteMCNameList: unable to open KosselNameList group',.TRUE.)

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
call h5gcreate_f(HDF_head%oID,'MCCLNameList',grp_id,error)
if (error.ne.0) call HDF_handleError(error,'HDFwriteMCCLNameList: unable to open MCCLNameList group',.TRUE.)

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
call h5gcreate_f(HDF_head%oID,'MCCLMultiLayerNameList',grp_id,error)
if (error.ne.0) call HDF_handleError(error, &
                     'HDFwriteMCCLMultiLayerNameList: unable to open MCCLMultiLayerNameList group',.TRUE.)

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
call h5gcreate_f(HDF_head%oID,'EBSDMasterNameList',grp_id,error)
if (error.ne.0) call HDF_handleError(error, &
                     'HDFwriteEBSDMasterNameList: unable to open EBSDMasterNameList group',.TRUE.)

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
!> @author Saransh Singh/Marc De Graef, Carnegie Mellon University
!
!> @brief read namelist file and fill mcnl structure (used by CTEMECPmaster.f90)
!
!> @param nmlfile namelist file name
!> @param emnl ECP master name list structure
!
!> @date 06/19/14  SS 1.0 new routine
!--------------------------------------------------------------------------
subroutine HDFwriteECPMasterNameList(nmlfile, ecpnl)

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

end subroutine HDFwriteECPMasterNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:HDFwriteEBSDNameList
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
subroutine HDFwriteEBSDNameList(nmlfile, enl)

use error

IMPLICIT NONE

character(fnlen),INTENT(IN)               :: nmlfile
type(EBSDNameListType),INTENT(INOUT)      :: enl

integer(kind=irg)       :: stdout
integer(kind=irg)       :: numsx
integer(kind=irg)       :: numsy
integer(kind=irg)       :: binning
integer(kind=irg)       :: nthreads
integer(kind=irg)       :: energyaverage
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
character(1)            :: maskpattern
character(3)            :: scalingmode
character(3)            :: eulerconvention
character(3)            :: outputformat
character(fnlen)        :: anglefile
character(fnlen)        :: masterfile
character(fnlen)        :: energyfile 
character(fnlen)        :: datafile

! define the IO namelist to facilitate passing variables to the program.
namelist  / EBSDdata / stdout, L, thetac, delta, numsx, numsy, xpc, ypc, anglefile, eulerconvention, masterfile, &
                        energyfile, datafile, beamcurrent, dwelltime, energymin, energymax, binning, gammavalue, &
                        scalingmode, axisangle, nthreads, outputformat, maskpattern, energyaverage

! set the input parameters to default values (except for xtalname, which must be present)
stdout          = 6
numsx           = 640           ! [dimensionless]
numsy           = 480           ! [dimensionless]
binning         = 1             ! binning mode  (1, 2, 4, or 8)
L               = 20000.0       ! [microns]
nthreads        = 1             ! number of OpenMP threads
energyaverage   = 0             ! apply energy averaging (1) or not (0); useful for dictionary computations
thetac          = 0.0           ! [degrees]
delta           = 25.0          ! [microns]
xpc             = 0.0           ! [pixels]
ypc             = 0.0           ! [pixels]
energymin       = 15.0          ! minimum energy to consider
energymax       = 30.0          ! maximum energy to consider
gammavalue      = 1.0           ! gamma factor
axisangle       = (/0.0, 0.0, 1.0, 0.0/)        ! no additional axis angle rotation
beamcurrent     = 14.513D0      ! beam current (actually emission current) in nano ampere
dwelltime       = 100.0D0       ! in microseconds
maskpattern     = 'n'           ! 'y' or 'n' to include a circular mask
scalingmode     = 'not'         ! intensity selector ('lin', 'gam', or 'not')
eulerconvention = 'tsl'         ! convention for the first Euler angle ['tsl' or 'hkl']
outputformat    = 'gui'         ! output format for 'bin' or 'gui' use
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
enl%energyaverage = energyaverage
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
enl%maskpattern = maskpattern
enl%scalingmode = scalingmode
enl%eulerconvention = eulerconvention
enl%outputformat = outputformat
enl%anglefile = anglefile
enl%masterfile = masterfile
enl%energyfile = energyfile
enl%datafile = datafile

end subroutine HDFwriteEBSDNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:HDFwriteECPNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read namelist file and fill ecpnl structure (used by CTEMECP.f90)
!
!> @param nmlfile namelist file name
!> @param knl Kossel name list structure
!
!> @date 06/13/14  MDG 1.0 new routine
!> @date 11/25/14  MDG 2.0 added parameters for film on substrate mode
!--------------------------------------------------------------------------
subroutine HDFwriteECPNameList(nmlfile, ecpnl)

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

end subroutine HDFwriteECPNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:HDFwriteLACBEDNameList
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
subroutine HDFwriteLACBEDNameList(nmlfile, lacbednl)

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

end subroutine HDFwriteLACBEDNameList


!--------------------------------------------------------------------------
!
! SUBROUTINE:HDFwriteECPpatternNameList
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
subroutine HDFwriteECPpatternNameList(nmlfile,ecpnl)

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

end subroutine HDFwriteECPpatternNameList

!--------------------------------------------------------------------------
!
! SUBROUTINE:HDFwritePEDKINNameList
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read namelist file and fill pednl structure (used by CTEMpedKIN.f90)
!
!> @param nmlfile namelist file name
!> @param pednl PED name list structure
!
!> @date 03/02/15 MDG 1.0 new routine
!--------------------------------------------------------------------------
subroutine HDFwritePEDKINNameList(nmlfile,pednl)

use error

IMPLICIT NONE

character(fnlen),INTENT(IN)                     :: nmlfile
type(PEDKINNameListType),INTENT(INOUT)             :: pednl

integer(kind=irg)       :: stdout
integer(kind=irg)       :: npix
integer(kind=irg)       :: ncubochoric
integer(kind=irg)       :: nthreads
real(kind=sgl)          :: voltage
real(kind=sgl)          :: thickness
real(kind=sgl)          :: rnmpp
real(kind=sgl)          :: dmin
character(fnlen)        :: xtalname
character(fnlen)        :: outname
character(fnlen)        :: eulername

! define the IO namelist to facilitate passing variables to the program.
namelist /inputlist/ stdout, xtalname, voltage, npix, rnmpp, ncubochoric, nthreads, &
                              thickness, outname , dmin, eulername

! set the input parameters to default values (except for xtalname, which must be present)
xtalname = 'undefined'          ! initial value to check that the keyword is present in the nml file
stdout = 6                      ! standard output
voltage = 200000.0              ! acceleration voltage [V]
nthreads = 1                    ! number of OpenMP threads to start
thickness = 10.0                ! sample thickness [nm]
npix = 256                      ! output arrays will have size npix x npix
outname = 'pedout.data'         ! output filename
eulername = 'EulerAngles.txt'   ! output filename
dmin = 0.04                     ! smallest d-spacing [nm]
ncubochoric = 100               ! number of samples along the cubochoric edge length
rnmpp = 0.2                     ! nm^{-1} per pattern pixel

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
pednl%thickness = thickness
pednl%dmin = dmin
pednl%npix = npix
pednl%nthreads = nthreads
pednl%outname = outname
pednl%eulername = eulername
pednl%rnmpp = rnmpp
pednl%ncubochoric = ncubochoric

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
