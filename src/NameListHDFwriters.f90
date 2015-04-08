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
!> @date 03/28/15 MDG 2.0 removing all h5lt calls; replaced with HDFsupport calls
!--------------------------------------------------------------------------
module NameListHDFwriters

use local
use typedefs
use NameListTypedefs
use HDF5
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
subroutine HDF_writeNMLintegers(HDF_head, HDF_tail, io_int, intlist, n_int)

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head
type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_tail
integer(kind=irg),INTENT(IN)                          :: io_int(n_int)
character(20),INTENT(IN)                              :: intlist(n_int)
integer(kind=irg),INTENT(IN)                          :: n_int

integer(kind=irg)                                     :: hdferr, i
character(fnlen)                                      :: dataset

do i=1,n_int
  dataset = intlist(i)
  hdferr = HDF_writeDatasetInteger(dataset, io_int(i), HDF_head, HDF_tail)
  if (hdferr.ne.0) call HDF_handleError(hdferr,'HDF_writeNMLintegers: unable to create '//trim(intlist(i))//' dataset',.TRUE.)
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
subroutine HDF_writeNMLreals(HDF_head, HDF_tail, io_real, reallist, n_real)

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head
type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_tail
real(kind=sgl),INTENT(IN)                             :: io_real(n_real)
character(20),INTENT(IN)                              :: reallist(n_real)
integer(kind=irg),INTENT(IN)                          :: n_real

integer(kind=irg)                                     :: hdferr, i
character(fnlen)                                      :: dataset

do i=1,n_real
  dataset = reallist(i)
  hdferr = HDF_writeDatasetFloat(dataset, io_real(i), HDF_head, HDF_tail)
  if (hdferr.ne.0) call HDF_handleError(hdferr,'HDF_writeNMLreals: unable to create '//trim(reallist(i))//' dataset',.TRUE.)
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
subroutine HDF_writeNMLdbles(HDF_head, HDF_tail, io_real, reallist, n_real)

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head
type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_tail
real(kind=dbl),INTENT(IN)                             :: io_real(n_real)
character(20),INTENT(IN)                              :: reallist(n_real)
integer(kind=irg),INTENT(IN)                          :: n_real

integer(kind=irg)                                     :: hdferr, i
character(fnlen)                                      :: dataset

do i=1,n_real
  dataset = reallist(i)
  hdferr = HDF_writeDatasetDouble(dataset, io_real(i), HDF_head, HDF_tail)
  if (hdferr.ne.0) call HDF_handleError(hdferr,'HDF_writeNMLdbles: unable to create '//trim(reallist(i))//' dataset',.TRUE.)
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

use ISO_C_BINDING

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(KosselNameListType),INTENT(IN)                   :: knl

integer(kind=irg),parameter                           :: n_int = 5, n_real = 6
integer(kind=irg)                                     :: hdferr,  io_int(n_int)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
character(fnlen)                                      :: dataset
character(fnlen,kind=c_char)                          :: line2(1)

! create the group for this namelist
hdferr = HDF_createGroup('KosselNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ knl%stdout, knl%numthick, knl%npix, knl%maxHOLZ, knl%nthreads /)
intlist(1) = 'stdout'
intlist(2) = 'numthick'
intlist(3) = 'npix'
intlist(4) = 'maxHOLZ'
intlist(5) = 'nthreads'
call HDF_writeNMLintegers(HDF_head, HDF_tail, io_int, intlist, n_int)

! integer vectors
dataset = 'k'
hdferr = HDF_writeDatasetIntegerArray1D(dataset, knl%k, 3, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteKosselNameList: unable to create k dataset',.TRUE.)

dataset = 'fn'
hdferr = HDF_writeDatasetIntegerArray1D(dataset, knl%fn, 3, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteKosselNameList: unable to create fn dataset',.TRUE.)

! write all the single reals
io_real = (/ knl%voltage, knl%dmin, knl%convergence, knl%startthick, knl%thickinc, knl%minten /)
reallist(1) = 'voltage'
reallist(2) = 'dmin'
reallist(3) = 'convergence'
reallist(4) = 'startthick'
reallist(5) = 'thickinc'
reallist(6) = 'minten'
call HDF_writeNMLreals(HDF_head, HDF_tail, io_real, reallist, n_real)

! write all the strings
dataset = 'xtalname'
line2(1) = knl%xtalname
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteKosselNameList: unable to create xtalname dataset',.TRUE.)

dataset = 'outname'
line2(1) = knl%outname
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteKosselNameList: unable to create outname dataset',.TRUE.)

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

use ISO_C_BINDING

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(KosselMasterNameListType),INTENT(IN)             :: knl

integer(kind=irg),parameter                           :: n_int = 4, n_real = 5
integer(kind=irg)                                     :: hdferr,  io_int(n_int)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
character(fnlen)                                      :: dataset
character(fnlen,kind=c_char)                          :: line2(1)

! create the group for this namelist
hdferr = HDF_createGroup('KosselMasterNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ knl%stdout, knl%numthick, knl%npix, knl%nthreads /)
intlist(1) = 'stdout'
intlist(2) = 'numthick'
intlist(3) = 'npix'
intlist(4) = 'nthreads' 
call HDF_writeNMLintegers(HDF_head,HDF_tail,  io_int, intlist, n_int)

! write all the single reals
io_real = (/ knl%voltage, knl%dmin, knl%startthick, knl%thickinc, knl%tfraction /)
reallist(1) = 'voltage'
reallist(2) = 'dmin'
reallist(3) = 'startthick'
reallist(4) = 'thickinc'
reallist(5) = 'tfraction' 
call HDF_writeNMLreals(HDF_head, HDF_tail, io_real, reallist, n_real)

! write all the strings
dataset = 'Kosselmode'
line2(1) = knl%Kosselmode
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteKosselMasterNameList: unable to create Kosselmode dataset',.TRUE.)

dataset = 'xtalname'
line2(1) = knl%xtalname
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteKosselMasterNameList: unable to create xtalname dataset',.TRUE.)

dataset = 'outname'
line2(1) = knl%outname
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteKosselMasterNameList: unable to create outname dataset',.TRUE.)

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

use ISO_C_BINDING

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(MCNameListType),INTENT(INOUT)                    :: mcnl

integer(kind=irg),parameter                           :: n_int = 5, n_real = 7
integer(kind=irg)                                     :: hdferr,  io_int(n_int)
real(kind=dbl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
character(fnlen)                                      :: dataset, sval(1)
character(fnlen,kind=c_char)                          :: line2(1)

! create the group for this namelist
hdferr = HDF_createGroup('MCNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ mcnl%stdout, mcnl%numsx, mcnl%primeseed, mcnl%num_el, mcnl%nthreads /)
intlist(1) = 'stdout'
intlist(2) = 'numsx'
intlist(3) = 'primeseed'
intlist(4) = 'num_el'
intlist(5) = 'nthreads'
call HDF_writeNMLintegers(HDF_head, HDF_tail, io_int, intlist, n_int)

! write all the single doubles
io_real = (/ mcnl%sig, mcnl%omega, mcnl%EkeV, mcnl%Ehistmin, mcnl%Ebinsize, mcnl%depthmax, mcnl%depthstep /)
reallist(1) = 'sig'
reallist(2) = 'omega'
reallist(3) = 'EkeV'
reallist(4) = 'Ehistmin'
reallist(5) = 'Ebinsize'
reallist(6) = 'depthmax'
reallist(7) = 'depthstep'
call HDF_writeNMLdbles(HDF_head, HDF_tail, io_real, reallist, n_real)

! write all the strings
dataset = 'MCmode'
sval(1) = mcnl%MCmode
hdferr = HDF_writeDatasetStringArray(dataset, sval, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteMCNameList: unable to create MCmode dataset',.TRUE.)

dataset = 'xtalname'
line2(1) = mcnl%xtalname
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteMCNameList: unable to create xtalname dataset',.TRUE.)

dataset = 'dataname'
line2(1) = mcnl%dataname
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteMCNameList: unable to create dataname dataset',.TRUE.)

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

use ISO_C_BINDING

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(MCCLNameListType),INTENT(INOUT)                  :: mcnl

integer(kind=irg),parameter                           :: n_int = 5, n_real = 7
integer(kind=irg)                                     :: hdferr,  io_int(n_int)
real(kind=dbl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
character(fnlen)                                      :: dataset, sval(1)
character(fnlen,kind=c_char)                          :: line2(1)

! create the group for this namelist
hdferr = HDF_createGroup('MCCLNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ mcnl%stdout, mcnl%numsx, mcnl%globalworkgrpsz, mcnl%num_el, mcnl%totnum_el /)
intlist(1) = 'stdout'
intlist(2) = 'numsx'
intlist(3) = 'globalworkgrpsz'
intlist(4) = 'num_el'
intlist(5) = 'totnum_el'
call HDF_writeNMLintegers(HDF_head, HDF_tail, io_int, intlist, n_int)

! write all the single doubles
io_real = (/ mcnl%sig, mcnl%omega, mcnl%EkeV, mcnl%Ehistmin, mcnl%Ebinsize, mcnl%depthmax, mcnl%depthstep /)
reallist(1) = 'sig'
reallist(2) = 'omega'
reallist(3) = 'EkeV'
reallist(4) = 'Ehistmin'
reallist(5) = 'Ebinsize'
reallist(6) = 'depthmax'
reallist(7) = 'depthstep'
call HDF_writeNMLdbles(HDF_head, HDF_tail, io_real, reallist, n_real)

! write all the strings
dataset = 'MCmode'
sval(1) = mcnl%MCmode
hdferr = HDF_writeDatasetStringArray(dataset, sval, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteMCCLNameList: unable to create MCmode dataset',.TRUE.)

dataset = 'xtalname'
line2(1) = mcnl%xtalname
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteMCCLNameList: unable to create xtalname dataset',.TRUE.)

dataset = 'dataname'
line2(1) = mcnl%dataname
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteMCCLNameList: unable to create dataname dataset',.TRUE.)

dataset = 'primelist'
line2(1) = mcnl%primelist
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteMCCLNameList: unable to create primelist dataset',.TRUE.)

dataset = 'mode'
sval(1) = mcnl%mode
hdferr = HDF_writeDatasetStringArray(dataset, sval, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteMCCLNameList: unable to create mode dataset',.TRUE.)

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

use ISO_C_BINDING

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(MCCLMultiLayerNameListType),INTENT(INOUT)        :: mcnl

integer(kind=irg),parameter                           :: n_int = 5, n_real = 9
integer(kind=irg)                                     :: hdferr,  io_int(n_int)
real(kind=dbl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
character(fnlen)                                      :: dataset
character(fnlen,kind=c_char)                          :: line2(1)

! create the group for this namelist
hdferr = HDF_createGroup('MCCLMultiLayerNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ mcnl%stdout, mcnl%numsx, mcnl%globalworkgrpsz, mcnl%num_el, mcnl%totnum_el /)
intlist(1) = 'stdout'
intlist(2) = 'numsx'
intlist(3) = 'globalworkgrpsz'
intlist(4) = 'num_el'
intlist(5) = 'totnum_el'
call HDF_writeNMLintegers(HDF_head, HDF_tail, io_int, intlist, n_int)

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
call HDF_writeNMLdbles(HDF_head, HDF_tail, io_real, reallist, n_real)

! write all the strings
dataset = 'MCmode'
line2(1) = mcnl%MCmode
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteMCCLMultiLayerNameList: unable to create MCmode dataset',.TRUE.)

dataset = 'xtalname_film'
line2(1) = mcnl%xtalname_film
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteMCCLMultiLayerNameList: unable to create xtalname_film dataset',.TRUE.)

dataset = 'xtalname_subs'
line2(1) = mcnl%xtalname_subs
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteMCCLMultiLayerNameList: unable to create xtalname_subs dataset',.TRUE.)

dataset = 'dataname'
line2(1) = mcnl%dataname
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteMCCLMultiLayerNameList: unable to create dataname dataset',.TRUE.)

dataset = 'primelist'
line2(1) = mcnl%primelist
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteMCCLMultiLayerNameList: unable to create primelist dataset',.TRUE.)

dataset = 'mode'
line2(1) = mcnl%mode
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteMCCLMultiLayerNameList: unable to create mode dataset',.TRUE.)

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

use ISO_C_BINDING

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(EBSDMasterNameListType),INTENT(INOUT)            :: emnl

integer(kind=irg),parameter                           :: n_int = 4, n_real = 1
integer(kind=irg)                                     :: hdferr,  io_int(n_int)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
character(fnlen)                                      :: dataset
character(fnlen,kind=c_char)                          :: line2(1)

! create the group for this namelist
hdferr = HDF_createGroup('EBSDMasterNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ emnl%stdout, emnl%npx, emnl%Esel, emnl%nthreads /)
intlist(1) = 'stdout'
intlist(2) = 'npx'
intlist(3) = 'Esel'
intlist(4) = 'nthreads'
call HDF_writeNMLintegers(HDF_head, HDF_tail, io_int, intlist, n_int)

! write a single real
dataset = 'dmin'
hdferr = HDF_writeDatasetFloat(dataset, emnl%dmin, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteEBSDMasterNameList: unable to create dmin dataset',.TRUE.)

! write all the strings
dataset = 'outname'
line2(1) = emnl%outname
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteEBSDMasterNameList: unable to create outname dataset',.TRUE.)

dataset = 'energyfile'
line2(1) = emnl%energyfile
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteEBSDMasterNameList: unable to create energyfile dataset',.TRUE.)

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

use ISO_C_BINDING

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(ECPMasterNameListType),INTENT(INOUT)             :: ecpnl

integer(kind=irg),parameter                           :: n_int = 4, n_real = 2
integer(kind=irg)                                     :: hdferr, io_int(n_int), distort
real(kind=dbl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
character(fnlen)                                      :: dataset
character(fnlen,kind=c_char)                          :: line2(1)

! create the group for this namelist
hdferr = HDF_createGroup('ECPMasterNameList',HDF_head, HDF_tail)

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
call HDF_writeNMLintegers(HDF_head, HDF_tail, io_int, intlist, n_int)

! integer vectors
dataset = 'fn'
hdferr = HDF_writeDatasetFloatArray1D(dataset, ecpnl%fn, 3, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECPMasterNameList: unable to create fn dataset',.TRUE.)

! write all the single doubles
io_real = (/ ecpnl%dmin, ecpnl%startthick /)
reallist(1) = 'dmin'
reallist(2) = 'startthick' 
call HDF_writeNMLdbles(HDF_head, HDF_tail, io_real, reallist, n_real)

! 3-vectors (real)
dataset = 'abcdist'
hdferr = HDF_writeDatasetFloatArray1D(dataset, ecpnl%abcdist, 3, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECPMasterNameList: unable to create abcdist dataset',.TRUE.)

dataset = 'albegadist'
hdferr = HDF_writeDatasetFloatArray1D(dataset, ecpnl%albegadist, 3, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECPMasterNameList: unable to create albegadist dataset',.TRUE.)

! write all the strings
dataset = 'outname'
line2(1) = ecpnl%outname
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECPMasterNameList: unable to create outname dataset',.TRUE.)

dataset = 'energyfile'
line2(1) = ecpnl%energyfile
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECPMasterNameList: unable to create energyfile dataset',.TRUE.)

dataset = 'compmode'
line2(1) = ecpnl%compmode
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECPMasterNameList: unable to create compmode dataset',.TRUE.)

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

use ISO_C_BINDING

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(EBSDNameListType),INTENT(INOUT)                  :: enl

integer(kind=irg),parameter                           :: n_int = 6, n_real = 8
integer(kind=irg)                                     :: hdferr,  io_int(n_int)
real(kind=sgl)                                        :: io_real(n_real)
real(kind=dbl)                                        :: t(1)
character(20)                                         :: intlist(n_int), reallist(n_real)
character(fnlen)                                      :: dataset
character(fnlen,kind=c_char)                          :: line2(1)


! create the group for this namelist
hdferr = HDF_createGroup('EBSDNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ enl%stdout, enl%numsx, enl%numsy, enl%binning, enl%nthreads, enl%energyaverage /)
intlist(1) = 'stdout'
intlist(2) = 'numsx'
intlist(3) = 'numsy'
intlist(4) = 'binning'
intlist(5) = 'nthreads'
intlist(6) = 'energyaverage'
call HDF_writeNMLintegers(HDF_head, HDF_tail, io_int, intlist, n_int)

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
call HDF_writeNMLreals(HDF_head, HDF_tail, io_real, reallist, n_real)

! a 4-vector
dataset = 'axisangle'
hdferr = HDF_writeDatasetFloatArray1D(dataset, enl%axisangle, 4, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteEBSDNameList: unable to create axisangle dataset',.TRUE.)

! a few doubles
dataset = 'beamcurrent'
hdferr = HDF_writeDatasetDouble(dataset, enl%beamcurrent, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteEBSDNameList: unable to create beamcurrent dataset',.TRUE.)

dataset = 'dwelltime'
hdferr = HDF_writeDatasetDouble(dataset, enl%dwelltime, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteEBSDNameList: unable to create dwelltime dataset',.TRUE.)

! write all the strings
dataset = 'maskpattern'
line2(1) = trim(enl%maskpattern)
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteEBSDNameList: unable to create maskpattern dataset',.TRUE.)

dataset = 'scalingmode'
line2(1) = trim(enl%scalingmode)
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteEBSDNameList: unable to create scalingmode dataset',.TRUE.)

dataset = 'eulerconvention'
line2(1) = trim(enl%eulerconvention)
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteEBSDNameList: unable to create eulerconvention dataset',.TRUE.)

dataset = 'outputformat'
line2(1) = trim(enl%outputformat)
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteEBSDNameList: unable to create outputformat dataset',.TRUE.)

dataset = 'energyfile'
line2(1) = trim(enl%energyfile)
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteEBSDNameList: unable to create energyfile dataset',.TRUE.)

dataset = 'masterfile'
line2(1) = trim(enl%masterfile)
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteEBSDNameList: unable to create masterfile dataset',.TRUE.)

dataset = 'anglefile'
line2(1) = trim(enl%anglefile)
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteEBSDNameList: unable to create anglefile dataset',.TRUE.)

dataset = 'datafile'
line2(1) = trim(enl%datafile)
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteEBSDNameList: unable to create datafile dataset',.TRUE.)

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

use ISO_C_BINDING

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(ECPNameListType),INTENT(INOUT)                   :: ecpnl

integer(kind=irg),parameter                           :: n_int = 4, n_real = 8
integer(kind=irg)                                     :: hdferr,  io_int(n_int)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
character(fnlen)                                      :: dataset
character(fnlen,kind=c_char)                          :: line2(1)

! create the group for this namelist
hdferr = HDF_createGroup('ECPNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ ecpnl%stdout, ecpnl%numthick, ecpnl%npix, ecpnl%nthreads /)
intlist(1) = 'stdout'
intlist(2) = 'numthick'
intlist(3) = 'npix'
intlist(4) = 'nthreads'
call HDF_writeNMLintegers(HDF_head, HDF_tail, io_int, intlist, n_int)

! integer vectors
dataset = 'k'
hdferr = HDF_writeDatasetIntegerArray1D(dataset, ecpnl%k, 3, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECPNameList: unable to create k dataset',.TRUE.)

dataset = 'fn'
hdferr = HDF_writeDatasetIntegerArray1D(dataset, ecpnl%fn, 3, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECPNameList: unable to create fn dataset',.TRUE.)

dataset = 'gF'
hdferr = HDF_writeDatasetIntegerArray1D(dataset, ecpnl%gF, 3, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECPNameList: unable to create gF dataset',.TRUE.)

dataset = 'gS'
hdferr = HDF_writeDatasetIntegerArray1D(dataset, ecpnl%gS, 3, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECPNameList: unable to create gS dataset',.TRUE.)

dataset = 'tF'
hdferr = HDF_writeDatasetIntegerArray1D(dataset, ecpnl%tF, 3, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECPNameList: unable to create tF dataset',.TRUE.)

dataset = 'tS'
hdferr = HDF_writeDatasetIntegerArray1D(dataset, ecpnl%tS, 3, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECPNameList: unable to create tS dataset',.TRUE.)

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
call HDF_writeNMLreals(HDF_head, HDF_tail, io_real, reallist, n_real)

! write all the strings
dataset = 'compmode'
line2(1) = ecpnl%compmode
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECPNameList: unable to create compmode dataset',.TRUE.)

dataset = 'energyfile'
line2(1) = ecpnl%energyfile
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECPNameList: unable to create energyfile dataset',.TRUE.)

dataset = 'outname'
line2(1) = ecpnl%outname
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECPNameList: unable to create outname dataset',.TRUE.)

dataset = 'xtalname'
line2(1) = ecpnl%xtalname
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECPNameList: unable to create xtalname dataset',.TRUE.)

dataset = 'xtalname2'
line2(1) = ecpnl%xtalname2
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECPNameList: unable to create xtalname2 dataset',.TRUE.)

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

use ISO_C_BINDING

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(LACBEDNameListType),INTENT(INOUT)                :: lacbednl

integer(kind=irg),parameter                           :: n_int = 5, n_real = 6
integer(kind=irg)                                     :: hdferr,  io_int(n_int)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
character(fnlen)                                      :: dataset
character(fnlen,kind=c_char)                          :: line2(1)

! create the group for this namelist
hdferr = HDF_createGroup('LACBEDNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ lacbednl%stdout, lacbednl%maxHOLZ, lacbednl%numthick, lacbednl%npix, lacbednl%nthreads /)
intlist(1) = 'stdout'
intlist(2) = 'maxHOLZ'
intlist(3) = 'numthick'
intlist(4) = 'npix'
intlist(5) = 'nthreads'
call HDF_writeNMLintegers(HDF_head, HDF_tail, io_int, intlist, n_int)

! vectors
dataset = 'k'
hdferr = HDF_writeDatasetIntegerArray1D(dataset, lacbednl%k, 3, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteLACBEDNameList: unable to create k dataset',.TRUE.)

dataset = 'fn'
hdferr = HDF_writeDatasetIntegerArray1D(dataset, lacbednl%fn, 3, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteLACBEDNameList: unable to create fn dataset',.TRUE.)

! write all the single reals
io_real = (/ lacbednl%voltage, lacbednl%dmin, lacbednl%convergence, lacbednl%startthick, lacbednl%thickinc, lacbednl%minten/)
reallist(1) = 'voltage'
reallist(2) = 'dmin'
reallist(3) = 'convergence'
reallist(4) = 'startthick'
reallist(5) = 'thickinc'
reallist(6) = 'minten'
call HDF_writeNMLreals(HDF_head, HDF_tail, io_real, reallist, n_real)

! write all the strings
dataset = 'outname'
line2(1) = lacbednl%outname
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteLACBEDNameList: unable to create outname dataset',.TRUE.)

dataset = 'xtalname'
line2(1) = lacbednl%xtalname
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteLACBEDNameList: unable to create xtalname dataset',.TRUE.)

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

use ISO_C_BINDING

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(ECPpatternNameListType),INTENT(INOUT)            :: ecpnl

integer(kind=irg),parameter                           :: n_int = 2, n_real = 6
integer(kind=irg)                                     :: hdferr,  io_int(n_int)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
character(fnlen)                                      :: dataset
character(fnlen,kind=c_char)                          :: line2(1)

! create the group for this namelist
hdferr = HDF_createGroup('ECPpatternNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ ecpnl%stdout, ecpnl%npix /)
intlist(1) = 'stdout'
intlist(2) = 'npix'
call HDF_writeNMLintegers(HDF_head, HDF_tail, io_int, intlist, n_int)

! single real
dataset = 'thetac'
hdferr = HDF_writeDatasetFloat(dataset, ecpnl%thetac, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECPpatternNameList: unable to create thetac dataset',.TRUE.)

! real vector
dataset = 'k'
hdferr = HDF_writeDatasetFloatArray1D(dataset, ecpnl%k, 3, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECPpatternNameList: unable to create k dataset',.TRUE.)

! write all the strings
dataset = 'outname'
line2(1) = ecpnl%outname
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECPpatternNameList: unable to create outname dataset',.TRUE.)

dataset = 'masterfile'
line2(1) = ecpnl%masterfile
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECPpatternNameList: unable to create masterfile dataset',.TRUE.)

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

use ISO_C_BINDING

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(PEDKINNameListType),INTENT(INOUT)                :: pednl

integer(kind=irg),parameter                           :: n_int = 4, n_real = 4
integer(kind=irg)                                     :: hdferr,  io_int(n_int)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
character(fnlen)                                      :: dataset
character(fnlen,kind=c_char)                          :: line2(1)

! create the group for this namelist
hdferr = HDF_createGroup('PEDKINNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ pednl%stdout, pednl%npix, pednl%ncubochoric, pednl%nthreads /)
intlist(1) = 'stdout'
intlist(2) = 'npix'
intlist(3) = 'ncubochoric'
intlist(4) = 'nthreads'
call HDF_writeNMLintegers(HDF_head, HDF_tail, io_int, intlist, n_int)

! write all the single reals
io_real = (/ pednl%voltage, pednl%thickness, pednl%dmin, pednl%rnmpp /)
reallist(1) = 'voltage'
reallist(2) = 'thickness'
reallist(3) = 'dmin'
reallist(4) = 'rnmpp'
call HDF_writeNMLreals(HDF_head, HDF_tail, io_real, reallist, n_real)

! write all the strings
dataset = 'outname'
line2(1) = pednl%outname
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwritePEDKINNameList: unable to create outname dataset',.TRUE.)

dataset = 'xtalname'
line2(1) = pednl%xtalname
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwritePEDKINNameList: unable to create xtalname dataset',.TRUE.)

dataset = 'eulername'
line2(1) = pednl%eulername
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwritePEDKINNameList: unable to create eulername dataset',.TRUE.)

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

use ISO_C_BINDING

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(PEDNameListType),INTENT(INOUT)                   :: pednl

integer(kind=irg),parameter                           :: n_int = 5, n_real = 6
integer(kind=irg)                                     :: hdferr,  io_int(n_int)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
character(fnlen)                                      :: dataset
character(fnlen,kind=c_char)                          :: line2(1)

! create the group for this namelist
hdferr = HDF_createGroup('PEDNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ pednl%stdout, pednl%precsample, pednl%precazimuthal, pednl%npix, pednl%nthreads /)
intlist(1) = 'stdout'
intlist(2) = 'precsample'
intlist(3) = 'precazimuthal'
intlist(4) = 'npix'
intlist(5) = 'nthreads'
call HDF_writeNMLintegers(HDF_head, HDF_tail, io_int, intlist, n_int)

! vectors
dataset = 'k'
hdferr = HDF_writeDatasetIntegerArray1D(dataset, pednl%k, 3, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwritePEDNameList: unable to create k dataset',.TRUE.)

dataset = 'fn'
hdferr = HDF_writeDatasetIntegerArray1D(dataset, pednl%fn, 3, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwritePEDNameList: unable to create fn dataset',.TRUE.)

! single reals
io_real = (/ pednl%voltage, pednl%dmin, pednl%precangle, pednl%prechalfwidth, pednl%thickness, pednl%camlen /)
reallist(1) = 'voltage'
reallist(2) = 'dmin'
reallist(3) = 'precangle'
reallist(4) = 'prechalfwidth'
reallist(5) = 'thickness'
reallist(6) = 'camlen'
call HDF_writeNMLreals(HDF_head, HDF_tail, io_real, reallist, n_real)

! write all the strings
dataset = 'outname'
line2(1) = pednl%outname
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwritePEDNameList: unable to create outname dataset',.TRUE.)

dataset = 'xtalname'
line2(1) = pednl%xtalname
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwritePEDNameList: unable to create xtalname dataset',.TRUE.)

dataset = 'filemode'
line2(1) = pednl%filemode
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwritePEDNameList: unable to create filemode dataset',.TRUE.)

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

use ISO_C_BINDING

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(ECCINameListType),INTENT(INOUT)                  :: eccinl

integer(kind=irg),parameter                           :: n_int = 8, n_real = 6
integer(kind=irg)                                     :: hdferr,  io_int(n_int)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
integer(kind=irg)                                     :: i
character(fnlen)                                      :: dataset
character(fnlen,kind=c_char)                          :: line2(1)

! create the group for this namelist
hdferr = HDF_createGroup('ECCINameList',HDF_head, HDF_tail)

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
call HDF_writeNMLintegers(HDF_head, HDF_tail, io_int, intlist, n_int)

! vectors
dataset = 'k'
hdferr = HDF_writeDatasetIntegerArray1D(dataset, eccinl%k, 3, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECCINameList: unable to create k dataset',.TRUE.)

! single reals
io_real = (/ eccinl%voltage, eccinl%dkt, eccinl%ktmax, eccinl%dmin, eccinl%DF_L, eccinl%DF_slice /)
reallist(1) = 'voltage'
reallist(2) = 'dkt'
reallist(3) = 'ktmax'
reallist(4) = 'dmin'
reallist(5) = 'DF_L'
reallist(6) = 'DF_slice'
call HDF_writeNMLreals(HDF_head, HDF_tail, io_real, reallist, n_real)

! 2-vectors
dataset = 'lauec'
hdferr = HDF_writeDatasetFloatArray1D(dataset, eccinl%lauec, 2, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECCINameList: unable to create lauec dataset',.TRUE.)

dataset = 'lauec2'
hdferr = HDF_writeDatasetFloatArray1D(dataset, eccinl%lauec2, 2, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECCINameList: unable to create lauec2 dataset',.TRUE.)

! write all the strings
dataset = 'dispmode'
line2(1) = eccinl%dispmode
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECCINameList: unable to create dispmode dataset',.TRUE.)

dataset = 'summode'
line2(1) = eccinl%dispmode
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECCINameList: unable to create summode dataset',.TRUE.)

dataset = 'progmode'
line2(1) = eccinl%dispmode
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECCINameList: unable to create progmode dataset',.TRUE.)

dataset = 'xtalname'
line2(1) = eccinl%dispmode
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECCINameList: unable to create xtalname dataset',.TRUE.)

dataset = 'foilnmlfile'
line2(1) = eccinl%dispmode
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECCINameList: unable to create foilnmlfile dataset',.TRUE.)

dataset = 'dispfile'
line2(1) = eccinl%dispmode
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECCINameList: unable to create dispfile dataset',.TRUE.)

dataset = 'dataname'
line2(1) = eccinl%dispmode
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECCINameList: unable to create dataname dataset',.TRUE.)

dataset = 'ECPname'
line2(1) = eccinl%dispmode
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECCINameList: unable to create ECPname dataset',.TRUE.)

dataset = 'sgname'
line2(1) = eccinl%dispmode
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECCINameList: unable to create sgname dataset',.TRUE.)

dataset = 'apbname'
line2(1) = eccinl%dispmode
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECCINameList: unable to create apbname dataset',.TRUE.)

dataset = 'incname'
line2(1) = eccinl%dispmode
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECCINameList: unable to create incname dataset',.TRUE.)

dataset = 'voidname'
line2(1) = eccinl%dispmode
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECCINameList: unable to create voidname dataset',.TRUE.)

! maxdefects string arrays
dataset = 'sfname'
hdferr = HDF_writeDatasetStringArray(dataset, eccinl%sfname(1:eccinl%numsf), eccinl%numsf, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECCINameList: unable to create sfname dataset',.TRUE.)

! 3*maxdefects string arrays
dataset = 'dislYname'
hdferr = HDF_writeDatasetStringArray(dataset, eccinl%dislYname(1:eccinl%numYdisl), eccinl%numYdisl, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECCINameList: unable to create dislYname dataset',.TRUE.)

dataset = 'dislname'
hdferr = HDF_writeDatasetStringArray(dataset, eccinl%dislname(1:eccinl%numdisl), eccinl%numdisl, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteECCINameList: unable to create dislname dataset',.TRUE.)

! and pop this group off the stack
call HDF_pop(HDF_head)

! IN PRINCIPLE, HERE WE SHOULD ALSO READ ALL OF THE DEFECT FILES AND INSERT THEM INTO THE HDF FILE

! TO BE WRITTEN



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

use ISO_C_BINDING

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(RFZNameListType),INTENT(INOUT)                   :: rfznl

integer(kind=irg),parameter                           :: n_int = 2, n_real = 1
integer(kind=irg)                                     :: hdferr,  io_int(n_int)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
integer(kind=irg)                                     :: i
character(fnlen)                                      :: dataset
character(fnlen,kind=c_char)                          :: line2(1)

! create the group for this namelist
hdferr = HDF_createGroup('RFZNameList',HDF_head, HDF_tail)

! write all the single integers
io_int = (/ rfznl%pgnum, rfznl%nsteps /)
intlist(1) = 'pgnum'
intlist(2) = 'nsteps'
call HDF_writeNMLintegers(HDF_head, HDF_tail, io_int, intlist, n_int)

! strings
dataset = 'outname'
line2(1) = rfznl%outname
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteRFZNameList: unable to create outname dataset',.TRUE.)

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

use ISO_C_BINDING

use local

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head, HDF_tail
type(DictIndxOpenCLListType),INTENT(INOUT)            :: dictindxnl

integer(kind=irg),parameter                           :: n_int = 8, n_real = 1
integer(kind=irg)                                     :: hdferr,  io_int(n_int)
real(kind=sgl)                                        :: io_real(n_real)
character(20)                                         :: intlist(n_int), reallist(n_real)
integer(kind=irg)                                     :: i
character(fnlen)                                      :: dataset
character(fnlen,kind=c_char)                          :: line2(1)

! create the group for this namelist
hdferr = HDF_createGroup('DictIndxOpenCLNameList',HDF_head, HDF_tail)

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
call HDF_writeNMLintegers(HDF_head, HDF_tail, io_int, intlist, n_int)

! strings
dataset = 'exptfile'
line2(1) = dictindxnl%exptfile
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteDictIndxOpenCLNameList: unable to create exptfile dataset',.TRUE.)

dataset = 'dictfile'
line2(1) = dictindxnl%dictfile
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteDictIndxOpenCLNameList: unable to create dictfile dataset',.TRUE.)

dataset = 'eulerfile'
line2(1) = dictindxnl%eulerfile
hdferr = HDF_writeDatasetStringArray(dataset, line2, 1, HDF_head, HDF_tail)
if (hdferr.ne.0) call HDF_handleError(hdferr,'HDFwriteDictIndxOpenCLNameList: unable to create eulerfile dataset',.TRUE.)

! and pop this group off the stack
call HDF_pop(HDF_head)

end subroutine HDFwriteDictIndxOpenCLNameList

end module NameListHDFwriters
