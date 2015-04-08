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
! EMsoft:HDFsupport.f90
!--------------------------------------------------------------------------
!
! MODULE: HDFsupport
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief HDF5 helper routines
!
!> @details: the EM HDF format has four Groups : 
!>
!> EMheader, created with HDF_writeEMheader
!> 
!> EMNMLfiles, which contains verbatim all the nml files as string arrays
!>
!> EMNMKparameters contains the parsed nml files, one parameter at a time
!>
!> and finally the EMData section, which has all the program output
!>
!
!> @date  03/17/15 MDG 1.0 original
!> @date  03/27/15 MDG 1.1 added integer and real write routines
!> @date  03/29/15 MDG 1.2 removed all h5lt routines
!> @date  03/31/15 MDG 1.3 added support for arrays of c_chars
!> @date  04/07/15 MDG 1.4 added hyperslab routines for char, integer, float and double in 2, 3, and 4D
!> @date  04/08/15 MDG 1.5 removed HDF_tail pointer as it was no longer needed
!> @data  04/08/15 MDG 1.6 added optional overwrite keyword to all HDF_writeDataset... routines
!--------------------------------------------------------------------------
module HDFsupport

use local
use typedefs
use HDF5

private :: HDF_readfromTextfile, HDF_push !, HDF_stackdump

contains

!--------------------------------------------------------------------------
!
! SUBROUTINE:HDF_writeEMheader
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write the EMsoft header information to the HDF file
!
!> @details The EMheader is a group that contains
!> the following basic dataset strings
!>
!> EMsoft version       : scversion from local.f90
!> execution date       : dstr
!> start time           : tstr1
!> end time             : tstr2
!> program name         : prn
!> user name            : username (local.f90) [these can/should be redefined by the user via nml files] 
!> user location        : userlocn (local.f90)
!> user email           : useremail (local.f90)
!> computer name        : read via system call hostnm()
!
!> @note The original version used regular character arrays and had problems when
!> datasets needed to be overwritten.  Because of this, we decided to start using 
!> the fortran2003 extensions, which have a richer data type set and allow for more
!> flexible C-bindings.
!>
!> @param HDF_head pointer to top of push-pop stack
!> @param dstr date string
!> @param tstrb time start string
!> @param tstre time end string
!> @param prn program name
!>
!> @date 03/20/15 MDG 1.0 original
!> @date 03/26/15 MDG 2.0 modified with fortran2003 resources
!--------------------------------------------------------------------------
subroutine HDF_writeEMheader(HDF_head, dstr, tstrb, tstre, prn)

use local
use io
use ISO_C_BINDING

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head
character(11),INTENT(IN)                              :: dstr
character(15),INTENT(IN)                              :: tstrb
character(15),INTENT(IN)                              :: tstre
character(fnlen),INTENT(IN)                           :: prn

integer                                               :: error, istat  ! error flag
integer                                               :: i,ic,nlen 
character(100)                                        :: c
character(fnlen)                                      :: line 
character(fnlen,kind=c_char)                          :: line2(1)


! create and open the EMheader group
istat = HDF_createGroup('EMheader', HDF_head)

! version number /EMheader/Version 'character'
line = 'Version'
line2(1) = scversion
error = HDF_writeDatasetStringArray(line, line2, 1, HDF_head)

! execution data /EMheader/Date 'character'
line = 'Date'
line2(1) = dstr
error = HDF_writeDatasetStringArray(line, line2, 1, HDF_head)

! start time /EMheader/StartTime 'character'
line = 'StartTime'
line2(1) = tstrb
error = HDF_writeDatasetStringArray(line, line2, 1, HDF_head)

! stop time /EMheader/StopTime 'character'
line = 'StopTime'
line2(1) = tstre
error = HDF_writeDatasetStringArray(line, line2, 1, HDF_head)

! program name /EMheader/ProgramName 'character'
line = 'ProgramName'
line2(1) = prn 
error = HDF_writeDatasetStringArray(line, line2, 1, HDF_head)

! user name /EMheader/UserName 'character'
line = 'UserName'
line2(1) = username
error = HDF_writeDatasetStringArray(line, line2, 1, HDF_head)

! user location /EMheader/UserLocation 'character'
line = 'UserLocation'
line2(1) = userlocn
error = HDF_writeDatasetStringArray(line, line2, 1, HDF_head)

! user email /EMheader/UserEmail 'character'
line = 'UserEmail'
line2(1) = useremail
error = HDF_writeDatasetStringArray(line, line2, 1, HDF_head)

! hostname /EMheader/HostName 'character'
call hostnm(c)
! lowercase it
nlen = len(c) 
do i=1,nlen 
   ic = ichar(c(i:i)) 
   if (ic >= 65 .and. ic < 90) c(i:i) = char(ic+32) 
end do 
line = 'HostName'
line2(1) = c
error = HDF_writeDatasetStringArray(line, line2, 1, HDF_head)
if (error.ne.0) call HDF_handleError(error,'HDF_writeEMheader: unable to write HostName',.TRUE.)

! and close this group
call HDF_pop(HDF_head)

end subroutine HDF_writeEMheader

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! to avoid user usage of file, group, dataset, etc, IDs, we use a push-pop
! stack to keep track of the open items and close them again.  The only 
! price to pay for doing things this way, is that the output must be written
! in a particular order at first.
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------


!--------------------------------------------------------------------------
!
! SUBROUTINE:HDF_push
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief push an HDF object to the stack
!
!> @param HDF_head top of the current stack
!> @param oT object type character
!> @param oID object identifier
!> @param oName name
!> @param verbose (optional) 
!
!> @date 03/17/15  MDG 1.0 original
!--------------------------------------------------------------------------
subroutine HDF_push(HDF_head, oT, oID, oName, verbose)

use local
use io

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head
character(LEN=1),INTENT(IN)                           :: oT
integer(HID_T),INTENT(IN)                             :: oID 
character(fnlen),INTENT(IN)                           :: oName
logical,INTENT(IN),OPTIONAL                           :: verbose

type(HDFobjectStackType),pointer                      :: node
integer(kind=irg)                                     :: istat

! if the stack doesn't exist yet, create it.
if (.not.associated(HDF_head)) then 
   allocate(HDF_head,stat=istat)                        ! allocate new value
   if (istat.ne.0) call HDF_handleError(istat,'HDF_push: unable to allocate HDF_head pointer',.TRUE.)
   nullify(HDF_head%next)                               ! nullify next in tail value
   if (PRESENT(verbose)) then 
     if (verbose) call Message('  -> creating HDF_head linked list', frm = "(A)")
   end if
else
   allocate(node,stat=istat)                        ! allocate new value
   if (istat.ne.0) call HDF_handleError(istat,'HDF_push: unable to allocate node pointer',.TRUE.)
   node%next => HDF_head
   HDF_head => node
end if

! set the values
HDF_head % objectType = oT
HDF_head % objectID = oID
HDF_head % objectname = trim(oName)

if (present(verbose)) call HDF_stackdump(HDF_head)

end subroutine HDF_push

!--------------------------------------------------------------------------
!
! SUBROUTINE:HDF_pop
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief pop an HDF object from the stack and close it
!
!> @param HDF_head top of the current stack
!> @param closeall (optional) close all open objects
!
!> @date 03/17/15  MDG 1.0 original
!--------------------------------------------------------------------------
subroutine HDF_pop(HDF_head, closeall, verbose)

use local
use io

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical,INTENT(IN),optional                             :: closeall
logical,INTENT(IN),optional                             :: verbose

integer                                                 :: error, istat
type(HDFobjectStackType),pointer                        :: tmp

nullify(tmp)

if (PRESENT(closeall)) then  
! this would be called if an error arises that forces a complete shutdown of the program, or at the end of a regular program
  do while (associated(HDF_head % next)) 
! close the current object 
    error = HDF_close_level(HDF_head % objectType, HDF_head % objectID)
    if (error.ne.0) then
      call HDF_handleError(istat,'HDF_pop: unable to close requested level for object type '//HDF_head%objectType,.TRUE.)
    end if
! and re-point the stack head
    tmp => HDF_head
    HDF_head => HDF_head % next  
! delete the old entry
    deallocate(tmp)
  end do
  nullify(HDF_head)
else
! close the current object 
  error = HDF_close_level(HDF_head % objectType, HDF_head % objectID)
  if (error.ne.0) then
    call HDF_handleError(istat,'HDF_pop: unable to close requested level for object type '//HDF_head%objectType,.TRUE.)
  end if
! and re-point the stack head
  tmp => HDF_head
  HDF_head => HDF_head % next  
! delete the old entry
  deallocate(tmp)

  if (present(verbose)) call HDF_stackdump(HDF_head)
end if

contains

  function HDF_close_level(oT, oID) result(error)
  
  use local

  IMPLICIT NONE

  character(LEN=1),INTENT(IN)   :: oT
  integer(HID_T),INTENT(IN)     :: oID 
  integer(kind=irg)             :: error

  select case(oT)
  case ('f') 
    call h5fclose_f(oID, error)  ! close the file

  case ('g') 
    call h5gclose_f(oID, error)  ! close the group

  case ('d') 
    call h5dclose_f(oID, error)  ! close the data set

  case ('a') 
    call h5aclose_f(oID, error)  ! close the attribute

  case ('t') 
    call h5tclose_f(oID, error)  ! close the data type

  case ('s') 
    call h5sclose_f(oID, error)  ! close the data space

  case DEFAULT
  end select

end function HDF_close_level

end subroutine HDF_pop

!--------------------------------------------------------------------------
!
! SUBROUTINE:HDF_stackdump
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief print out the entire stack for debugging purposes
!
!> @param HDF_head top of the current stack
!
!> @date 03/19/15  MDG 1.0 original
!--------------------------------------------------------------------------
subroutine HDF_stackdump(HDF_head)

use local
use io

IMPLICIT NONE

type(HDFobjectStackType),INTENT(IN),pointer          :: HDF_head

type(HDFobjectStackType),pointer                     :: tmp
integer(kind=irg)                                    :: io_int(1)

tmp => HDF_head
if (.not.associated(tmp)) then
  call WriteValue('Stack is empty','')
else
  call WriteValue('HDF stack entries','')
  do
    if (.not.associated(tmp)) EXIT
    call WriteValue('','>'//tmp%objectType//'<  >'//trim(tmp%objectName)//'<', frm = "(A$)") 
    io_int(1) = tmp%objectID
    call WriteValue('',io_int,1,frm="(I12)")
    tmp => tmp%next
  end do
end if

end subroutine HDF_stackdump

!--------------------------------------------------------------------------
!
! SUBROUTINE:HDF_handleError
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief deal with an HDF error
!
!> @param error non-zero integer
!> @param OffendingRoutine name of the routine that caused the error + possible message
!> @param NonFatal (optional logical parameter) Fatal or non-fatal error (logical)
!
!> @date 03/17/15  MDG 1.0 original
!--------------------------------------------------------------------------
subroutine HDF_handleError(error,OffendingRoutine, NonFatal)

use io
use local

IMPLICIT NONE

integer(kind=irg),INTENT(IN) :: error              ! returned error code
character(LEN=*),INTENT(IN)  :: OffendingRoutine   ! name of offending routine + message
logical,OPTIONAL,INTENT(IN)  :: NonFatal           ! if true, then report the error but don't stop

integer(kind=irg)            :: io_int(1)

io_int(1) = error
call WriteValue('Error code : ',io_int,1)
call Message('   returned by routine '//OffendingRoutine,frm="(A)")

if (.not.present(NonFatal)) STOP  ! this is not very graceful, but it'll do the job for now ...

end subroutine HDF_handleError

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! from here on, we have basic HDF support routines
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------


!--------------------------------------------------------------------------
!
! SUBROUTINE:HDF_createFile
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Create a new HDF file (this also opens the file)
!
!> @param HDFname filename string
!> @param HDF_head
!
!> @date 03/17/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_createFile(HDFname, HDF_head) result(success)
 
use local

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: HDFname
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
integer(kind=irg)                                       :: success

integer(HID_T)                                          :: file_id ! file identifier
integer                                                 :: error  ! error flag

success = 0

! Create a new file using default properties.
call h5fcreate_f(trim(HDFname), H5F_ACC_TRUNC_F, file_id, error)
if (error.ne.0) then
  call HDF_handleError(error,'HDF_createFile: error creating file')
  success = -1
else
  call HDF_push(HDF_head, 'f', file_id, HDFname)
end if

end function HDF_createFile


!--------------------------------------------------------------------------
!
! SUBROUTINE:HDF_openFile
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Open an HDF file
!
!> @param HDFname filename string
!> @param HDF_head
!> @param readonly (optional) file open mode
!
!> @date 03/17/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_openFile(HDFname, HDF_head, readonly) result(success)

use local

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: HDFname
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical,INTENT(IN),OPTIONAL                             :: readonly
integer(kind=irg)                                       :: success

integer(HID_T)                                          :: file_id ! file identifier
integer                                                 :: error  ! error flag

success = 0
if (present(readonly)) then 
  call H5Fopen_f(trim(HDFname), H5F_ACC_RDONLY_F, file_id, error)
else
  call H5Fopen_f(trim(HDFname), H5F_ACC_RDWR_F, file_id, error)
end if
if (error.ne.0) then 
  call HDF_handleError(error,'HDF_openFile')
  success = -1
else
  call HDF_push(HDF_head, 'f', file_id, HDFname)
end if

end function HDF_openFile


!--------------------------------------------------------------------------
!
! SUBROUTINE:HDF_createGroup
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief create a new group (this also opens the group)
!
!> @param groupname filename string
!> @param HDF_head
!
!> @date 03/17/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_createGroup(groupname, HDF_head) result(success)

IMPLICIT NONE

character(*),INTENT(IN)                                 :: groupname
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
integer(kind=irg)                                       :: success

integer(HID_T)                                          :: group_id!  identifier
integer                                                 :: error  ! error flag

success = 0

call H5Gcreate_f(HDF_head%objectID, groupname, group_id, error)
if (error.ne.0) then
  call HDF_handleError(error,'HDF_createGroup')
  success = -1
else
! and put the group_id onto the HDF_stack
  call HDF_push(HDF_head, 'g', group_id, groupname)
end if

end function HDF_createGroup

!--------------------------------------------------------------------------
!
! SUBROUTINE:HDF_openGroup
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief open an existing group 
!
!> @param groupname filename string
!> @param HDF_head
!
!> @date 03/17/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_openGroup(groupname, HDF_head) result(success)

IMPLICIT NONE

character(*),INTENT(IN)                                 :: groupname
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
integer(kind=irg)                                       :: success

integer(HID_T)                                          :: group_id !  identifier
integer                                                 :: error  ! error flag

success = 0

call H5Gopen_f(HDF_head%objectID, groupname, group_id, error)
if (error.ne.0) then
  call HDF_handleError(error,'HDF_openGroup')
  success = -1
else
! put the group_id onto the HDF_stack
  call HDF_push(HDF_head, 'g', group_id, groupname)
end if

end function HDF_openGroup

!--------------------------------------------------------------------------
!
! SUBROUTINE:HDF_openDataset
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief open an existing dataset
!
!> @param dataname string
!> @param HDF_head
!
!> @date 03/17/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_openDataset(dataname, HDF_head) result(success)

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
integer(kind=irg)                                       :: success

integer(HID_T)                                          :: data_id !  identifier
integer                                                 :: error  ! error flag

success = 0

call H5dopen_f(HDF_head%objectID, dataname, data_id, error)
if (error.ne.0) then
  call HDF_handleError(error,'HDF_openDataset')
  success = -1
else
! put the data_id onto the HDF_stack
  call HDF_push(HDF_head, 'd', data_id, dataname)
end if

end function HDF_openDataset


!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeDatasetTextFile
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write a text file data set to the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param filename of the text file
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeDatasetTextFile(dataname, filename, HDF_head) result(success)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
character(fnlen),INTENT(IN)                             :: filename
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
integer(kind=irg)                                       :: success

character(len=fnlen, KIND=c_char),allocatable, TARGET   :: stringarray(:) 
integer(kind=irg)                                       :: nlines

integer(HSIZE_T)                                        :: dim0 
integer(SIZE_T)                                         :: sdim 
integer(HID_T)                                          :: filetype, space, dset ! Handles
integer                                                 :: hdferr, i, rnk
integer(HSIZE_T), DIMENSION(1:1)                        :: dims
integer(HSIZE_T), DIMENSION(1:2)                        :: maxdims

TYPE(C_PTR), ALLOCATABLE, TARGET                        :: wdata(:)
TYPE(C_PTR)                                             :: f_ptr

success = 0

stringarray =  HDF_readfromTextfile(filename, nlines) 

! first, convert the stringarray to an array of C-pointers, with each string
! terminated by a C_NULL_CHAR.
dims(1) = nlines
allocate(wdata(1:dims(1)))
do i=1,dims(1)
  wdata(i) = C_LOC(stringarray(i))
end do

! then we write this C_ptr to the HDF file in the proper data set

! first create the memory data type (filetype)
call H5Tcopy_f(H5T_STRING, filetype, hdferr)
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetStringArray',.TRUE.)
  success = -1
end if
!
! Create dataspace.
!
rnk = 1
call h5screate_simple_f(rnk, dims, space, hdferr)
!
! Create the dataset and write the variable-length string data to it.
!
call h5dcreate_f(HDF_head%objectID, trim(dataname), filetype, space, dset, hdferr)
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetStringArray:hd5create_f',.TRUE.)
  success = -1
end if
f_ptr = C_LOC(wdata(1))
call h5dwrite_f(dset, filetype, f_ptr, hdferr )
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetStringArray:hd5write_f',.TRUE.)
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)
deallocate(wdata)

! that's it

end function HDF_writeDatasetTextFile

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readfromTextfile
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read a text file and return it in a C_NULL_CHAR terminated string array
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param filename file name (string)
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readfromTextfile(filename,nlines) result(stringarray)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                                  :: filename
integer(kind=irg),INTENT(OUT)                                :: nlines
character(len=fnlen, KIND=c_char), allocatable               :: stringarray(:) 

integer(kind=irg)                                            :: i, j, dt
character(len=fnlen, KIND=c_char), DIMENSION(1)              :: line 

dt = 55
! read the file first to determine the number of lines
open(unit=dt,file=trim(filename),action='read',form='formatted',status='old')
nlines = 0
do
  read (dt,"(A)",end=10) line(1)
  nlines = nlines + 1
end do
10 close(unit=dt,status='keep')

! then re-read the file and store all the lines in the wdata array
allocate(stringarray(1:nlines))
open(unit=dt,file=trim(filename),action='read',form='formatted',status='old')
do i=1,nlines
! initialize the line to null characters before each read
  do j=1,fnlen
    line(1)(j:j) = char(0)
  end do
! read the line
  read (dt,"(A)") line(1)
! find the string length and put the next character equal to C_NULL_CHAR
  j = len(trim(line(1)))+1
! truncate a line if it has more than fnlen characters
  if (j.gt.fnlen) j = fnlen
  line(1)(j:j) = C_NULL_CHAR
! store the line in the array
  stringarray(i) = line(1)
end do
close(unit=dt,status='keep')

end function HDF_readfromTextfile


!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readDatasetStringArray
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read a string array from a data set 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param nlines number of lines in string array
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readDatasetStringArray(dataname, nlines, HDF_head) result(stringarray)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
integer(kind=irg),INTENT(OUT)                           :: nlines
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
character(len=fnlen, KIND=c_char),allocatable, TARGET   :: stringarray(:) 

integer(HID_T)                                          :: filetype, space ! Handles
integer                                                 :: hdferr, i, length
integer(HSIZE_T), DIMENSION(1:1)                        :: dims
integer(HSIZE_T), DIMENSION(1:2)                        :: maxdims

character(len = fnlen, kind=c_char),  POINTER           :: pfstr ! A pointer to a Fortran string
TYPE(C_PTR), DIMENSION(:), ALLOCATABLE, TARGET          :: rdata ! Read buffer
TYPE(C_PTR)                                             :: f_ptr

! Open dataset.
!
hdferr = HDF_openDataset(dataname, HDF_head)
!
! Get the datatype.
!
call H5Dget_type_f(HDF_head%objectID, filetype, hdferr)
!
! Get dataspace and allocate memory for read buffer.
!
call H5Dget_space_f(HDF_head%objectID, space, hdferr)
call H5Sget_simple_extent_dims_f(space, dims, maxdims, hdferr)

ALLOCATE(rdata(1:dims(1)), stringarray(1:dims(1)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1))
call h5dread_f(HDF_head%objectID, H5T_STRING, f_ptr, hdferr)
!
! convert the data to a string array
!
DO i = 1, dims(1)
  call C_F_POINTER(rdata(i), pfstr)
  length = 0
  DO
     IF(pfstr(length+1:length+1).EQ.C_NULL_CHAR.OR.length.GE.fnlen) EXIT
     length = length + 1
  ENDDO
  stringarray(i) = pfstr(1:length)
END DO

nlines = dims(1)

DEALLOCATE(rdata)

call h5sclose_f(space, hdferr)
call H5Tclose_f(filetype, hdferr)
! close the dataset
call HDF_pop(HDF_head)


end function HDF_readDatasetStringArray

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_extractDatasetTextfile
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read a string array from a data set and stores it as a text file 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param textfile name of output text file
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_extractDatasetTextfile(dataname, textfile, HDF_head) result(success)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
character(fnlen),INTENT(IN)                             :: textfile
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
integer(kind=irg)                                       :: success

integer(HID_T)                                          :: filetype, space ! Handles
integer                                                 :: hdferr, i, length, dt
integer(HSIZE_T), DIMENSION(1:1)                        :: dims
integer(HSIZE_T), DIMENSION(1:2)                        :: maxdims

character(len = fnlen, kind=c_char),  POINTER           :: pfstr ! A pointer to a Fortran string
TYPE(C_PTR), DIMENSION(:), ALLOCATABLE, TARGET          :: rdata ! Read buffer
TYPE(C_PTR)                                             :: f_ptr

! Open dataset.
!
hdferr = HDF_openDataset(dataname, HDF_head)
!
! Get the datatype.
!
call H5Dget_type_f(HDF_head%objectID, filetype, hdferr)
!
! Get dataspace and allocate memory for read buffer.
!
call H5Dget_space_f(HDF_head%objectID, space, hdferr)
call H5Sget_simple_extent_dims_f(space, dims, maxdims, hdferr)

ALLOCATE(rdata(1:dims(1)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1))
call h5dread_f(HDF_head%objectID, H5T_STRING, f_ptr, hdferr)
!
! store the datain a textfile
!
dt = 55
open(unit=dt, file=trim(textfile), status='unknown', form='formatted')
DO i = 1, dims(1)
  call C_F_POINTER(rdata(i), pfstr)
  length = 0
  DO
     IF(pfstr(length+1:length+1).EQ.C_NULL_CHAR.OR.length.GE.fnlen) EXIT
     length = length + 1
  ENDDO
  write(dt,"(A)") pfstr(1:length)
END DO
close(unit=dt,status='keep')

DEALLOCATE(rdata)

call h5sclose_f(space, hdferr)
call H5Tclose_f(filetype, hdferr)
! close the dataset
call HDF_pop(HDF_head)

success = 0

end function HDF_extractDatasetTextfile

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeDatasetStringArray
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write a string array
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param filename of the text file
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeDatasetStringArray(dataname, inputarray, nlines, HDF_head, overwrite) result(success)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
character(len=fnlen),INTENT(IN)                         :: inputarray(nlines) 
integer(kind=irg),INTENT(IN)                            :: nlines
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical,INTENT(IN),OPTIONAL                             :: overwrite
integer(kind=irg)                                       :: success

character(len=fnlen, KIND=c_char),TARGET                :: stringarray(nlines) 
integer(HSIZE_T)                                        :: dim0 
integer(SIZE_T)                                         :: sdim 
integer(HID_T)                                          :: filetype, space, dset ! Handles
integer                                                 :: hdferr, i, rnk, l
integer(HSIZE_T), DIMENSION(1:1)                        :: dims
integer(HSIZE_T), DIMENSION(1:2)                        :: maxdims

TYPE(C_PTR), ALLOCATABLE, TARGET                        :: wdata(:)
TYPE(C_PTR)                                             :: f_ptr

success = 0

stringarray = ''
! first, convert the stringarray to an array of C-pointers, with each string
! terminated by a C_NULL_CHAR.
dims(1) = nlines
allocate(wdata(1:dims(1)))
do i=1,dims(1)
  l = len(trim(inputarray(i)))+1
  stringarray(i) = trim(inputarray(i))
  stringarray(i)(l:l) = C_NULL_CHAR
  wdata(i) = C_LOC(stringarray(i)(1:1))
end do

! then we write this C_ptr to the HDF file in the proper data set

! first create the memory data type (filetype)
call H5Tcopy_f(H5T_STRING, filetype, hdferr)
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetStringArray',.TRUE.)
  success = -1
end if
!
! Create dataspace.
!
rnk = 1
call h5screate_simple_f(rnk, dims, space, hdferr)
!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then 
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), filetype, space, dset, hdferr)
end if
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetStringArray:hd5create_f',.TRUE.)
  success = -1
end if
f_ptr = C_LOC(wdata(1))
call h5dwrite_f(dset, filetype, f_ptr, hdferr )
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetStringArray:hd5write_f',.TRUE.)
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

deallocate(wdata)

! that's it

end function HDF_writeDatasetStringArray

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeDatasetCharArray1D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write a 1D array of c_chars    [variable type H5T_STD_U8LE]
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param filename of the text file
!> @param HDF_head
!
!> @date 03/31/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeDatasetCharArray1D(dataname, chararray, dim0, HDF_head, overwrite) result(success)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
!character(len=1, KIND=c_char),TARGET                    :: chararray(dim0) 
character(len=1),TARGET                                 :: chararray(dim0) 
integer(kind=irg),INTENT(IN)                            :: dim0
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical,INTENT(IN),OPTIONAL                             :: overwrite
integer(kind=irg)                                       :: success

integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, i, rnk, l
integer(HSIZE_T), DIMENSION(1:1)                        :: dims

TYPE(C_PTR), dimension(1:1), TARGET                     :: wdata
TYPE(C_PTR)                                             :: f_ptr

success = 0

wdata(1) = C_LOC(chararray(1))
dims(1) = dim0

! then we write this C_ptr to the HDF file in the proper data set

! Create dataspace.
rnk = 1
call h5screate_simple_f(rnk, dims, space, hdferr)
!
! Create the dataset and write the c_char data to it.
!
if (present(overwrite)) then 
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_STD_U8LE, space, dset, hdferr)
end if
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetCharArray1D:hd5create_f',.TRUE.)
  success = -1
end if

call h5dwrite_f(dset, H5T_STD_U8LE, wdata(1), hdferr )
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetCharArray1D:hd5write_f',.TRUE.)
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_writeDatasetCharArray1D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeDatasetCharArray2D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write a 2D array of c_chars    [variable type H5T_STD_U8LE]
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param filename of the text file
!> @param HDF_head
!
!> @date 03/31/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeDatasetCharArray2D(dataname, chararray, dim0, dim1, HDF_head, overwrite) result(success)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
character(len=1),TARGET                                 :: chararray(dim0, dim1) 
integer(kind=irg),INTENT(IN)                            :: dim0
integer(kind=irg),INTENT(IN)                            :: dim1
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical,INTENT(IN),OPTIONAL                             :: overwrite
integer(kind=irg)                                       :: success

integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, i, rnk, l
integer(HSIZE_T), DIMENSION(1:2)                        :: dims

TYPE(C_PTR), dimension(1:1), TARGET                     :: wdata

success = 0

wdata(1) = C_LOC(chararray(1,1))
dims(1:2) = (/ dim0, dim1 /)

! then we write this C_ptr to the HDF file in the proper data set

! Create dataspace.
rnk = 2
call h5screate_simple_f(rnk, dims, space, hdferr)
!
! Create the dataset and write the c_char data to it.
!
if (present(overwrite)) then 
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_STD_U8LE, space, dset, hdferr)
end if 
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetCharArray1D:hd5create_f',.TRUE.)
  success = -1
end if

call h5dwrite_f(dset, H5T_STD_U8LE, wdata(1), hdferr )
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetCharArray1D:hd5write_f',.TRUE.)
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_writeDatasetCharArray2D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeDatasetCharArray3D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write a 3D array of c_chars    [variable type H5T_STD_U8LE]
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param filename of the text file
!> @param HDF_head
!
!> @date 03/31/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeDatasetCharArray3D(dataname, chararray, dim0, dim1, dim2, HDF_head, overwrite) result(success)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
character(len=1),TARGET                                 :: chararray(dim0, dim1, dim2) 
integer(kind=irg),INTENT(IN)                            :: dim0
integer(kind=irg),INTENT(IN)                            :: dim1
integer(kind=irg),INTENT(IN)                            :: dim2
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical,INTENT(IN),OPTIONAL                             :: overwrite
integer(kind=irg)                                       :: success

integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, i, rnk, l
integer(HSIZE_T), DIMENSION(1:3)                        :: dims

TYPE(C_PTR), dimension(1:3), TARGET                     :: wdata

success = 0

wdata(1) = C_LOC(chararray(1,1,1))
dims(1:3) = (/ dim0, dim1, dim2 /)

! then we write this C_ptr to the HDF file in the proper data set

! Create dataspace.
rnk = 3
call h5screate_simple_f(rnk, dims, space, hdferr)
!
! Create the dataset and write the c_char data to it.
!
if (present(overwrite)) then
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_STD_U8LE, space, dset, hdferr)
end if
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetCharArray1D:hd5create_f',.TRUE.)
  success = -1
end if

call h5dwrite_f(dset, H5T_STD_U8LE, wdata(1), hdferr )
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetCharArray1D:hd5write_f',.TRUE.)
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_writeDatasetCharArray3D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeDatasetCharArray4D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write a 4D array of c_chars    [variable type H5T_STD_U8LE]
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param filename of the text file
!> @param HDF_head
!
!> @date 03/31/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeDatasetCharArray4D(dataname, chararray, dim0, dim1, dim2, dim3, HDF_head, overwrite) result(success)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
character(len=1),TARGET                                 :: chararray(dim0, dim1, dim2, dim3) 
integer(kind=irg),INTENT(IN)                            :: dim0
integer(kind=irg),INTENT(IN)                            :: dim1
integer(kind=irg),INTENT(IN)                            :: dim2
integer(kind=irg),INTENT(IN)                            :: dim3
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical,INTENT(IN),OPTIONAL                             :: overwrite
integer(kind=irg)                                       :: success

integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, i, rnk, l
integer(HSIZE_T), DIMENSION(1:4)                        :: dims

TYPE(C_PTR), dimension(1:4), TARGET                     :: wdata

success = 0

wdata(1) = C_LOC(chararray(1,1,1,1))
dims(1:4) = (/ dim0, dim1, dim2, dim3 /)

! then we write this C_ptr to the HDF file in the proper data set

! Create dataspace.
rnk = 4
call h5screate_simple_f(rnk, dims, space, hdferr)
!
! Create the dataset and write the c_char data to it.
!
if (present(overwrite)) then
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_STD_U8LE, space, dset, hdferr)
end if 
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetCharArray1D:hd5create_f',.TRUE.)
  success = -1
end if

call h5dwrite_f(dset, H5T_STD_U8LE, wdata(1), hdferr )
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetCharArray1D:hd5write_f',.TRUE.)
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_writeDatasetCharArray4D



!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeDatasetInteger
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write an integer data set to the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param intval  integer 
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeDatasetInteger(dataname, intval, HDF_head, overwrite) result(success)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
integer(kind=irg),INTENT(IN)                            :: intval
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical,INTENT(IN),OPTIONAL                             :: overwrite
integer(kind=irg)                                       :: success

integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:1)                        :: dims

integer, dimension(1:1), TARGET                         :: wdata
TYPE(C_PTR)                                             :: f_ptr

success = 0

dims(1) = 1
wdata(1) = intval

! get a C pointer to the integer array
f_ptr = C_LOC(wdata(1))

! Create dataspace.
!
rnk = 1
call h5screate_simple_f(rnk, dims, space, hdferr)
!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_STD_I32LE, space, dset, hdferr)
end if 
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetIntegerArray1D:hd5create_f',.TRUE.)
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdferr )
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetIntegerArray1D:hd5write_f',.TRUE.)
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_writeDatasetInteger



!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeDatasetIntegerArray1D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write a 1D integer array data set to the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param intarr 1D integer array
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeDatasetIntegerArray1D(dataname, intarr, dim0, HDF_head, overwrite) result(success)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
integer(kind=irg),INTENT(IN)                            :: intarr(dim0)
integer(kind=irg),INTENT(IN)                            :: dim0
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical,INTENT(IN),OPTIONAL                             :: overwrite
integer(kind=irg)                                       :: success

integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:1)                        :: dims

integer, dimension(1:dim0), TARGET                      :: wdata
TYPE(C_PTR)                                             :: f_ptr

success = 0

dims(1) = dim0
wdata = intarr

! get a C pointer to the integer array
f_ptr = C_LOC(wdata(1))

! Create dataspace.
!
rnk = 1
call h5screate_simple_f(rnk, dims, space, hdferr)
!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_STD_I32LE, space, dset, hdferr)
end if
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetIntegerArray1D:hd5create_f',.TRUE.)
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdferr )
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetIntegerArray1D:hd5write_f',.TRUE.)
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_writeDatasetIntegerArray1D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeDatasetIntegerArray2D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write a 2D integer array data set to the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param intarr 1D integer array
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeDatasetIntegerArray2D(dataname, intarr, dim0, dim1, HDF_head, overwrite) result(success)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
integer(kind=irg),INTENT(IN)                            :: intarr(dim0, dim1)
integer(kind=irg),INTENT(IN)                            :: dim0
integer(kind=irg),INTENT(IN)                            :: dim1
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical,INTENT(IN),OPTIONAL                             :: overwrite
integer(kind=irg)                                       :: success

integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:2)                        :: dims

integer, dimension(1:dim0,1:dim1), TARGET               :: wdata
TYPE(C_PTR)                                             :: f_ptr

success = 0

dims(1:2) = (/ dim0, dim1 /)
wdata = intarr

! get a C pointer to the integer array
f_ptr = C_LOC(wdata(1,1))

! Create dataspace.
!
rnk = 2
call h5screate_simple_f(rnk, dims, space, hdferr)
!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_STD_I32LE, space, dset, hdferr)
end if
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetIntegerArray2D:hd5create_f',.TRUE.)
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdferr )
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetIntegerArray2D:hd5write_f',.TRUE.)
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_writeDatasetIntegerArray2D


!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeDatasetIntegerArray3D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write a 3D integer array data set to the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param intarr 1D integer array
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeDatasetIntegerArray3D(dataname, intarr, dim0, dim1, dim2, HDF_head, overwrite) result(success)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
integer(kind=irg),INTENT(IN)                            :: intarr(dim0, dim1, dim2)
integer(kind=irg),INTENT(IN)                            :: dim0
integer(kind=irg),INTENT(IN)                            :: dim1
integer(kind=irg),INTENT(IN)                            :: dim2
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical,INTENT(IN),OPTIONAL                             :: overwrite
integer(kind=irg)                                       :: success

integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:3)                        :: dims

integer, dimension(1:dim0,1:dim1,1:dim2), TARGET        :: wdata
TYPE(C_PTR)                                             :: f_ptr

success = 0

dims(1:3) = (/ dim0, dim1, dim2 /)
wdata = intarr

! get a C pointer to the integer array
f_ptr = C_LOC(wdata(1,1,1))

! Create dataspace.
!
rnk = 3
call h5screate_simple_f(rnk, dims, space, hdferr)
!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_STD_I32LE, space, dset, hdferr)
end if
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetIntegerArray3D:hd5create_f',.TRUE.)
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdferr )
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetIntegerArray3D:hd5write_f',.TRUE.)
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_writeDatasetIntegerArray3D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeDatasetIntegerArray4D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write a 4D integer array data set to the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param intarr 1D integer array
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeDatasetIntegerArray4D(dataname, intarr, dim0, dim1, dim2, dim3, HDF_head, overwrite) result(success)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
integer(kind=irg),INTENT(IN)                            :: intarr(dim0, dim1, dim2, dim3)
integer(kind=irg),INTENT(IN)                            :: dim0
integer(kind=irg),INTENT(IN)                            :: dim1
integer(kind=irg),INTENT(IN)                            :: dim2
integer(kind=irg),INTENT(IN)                            :: dim3
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical,INTENT(IN),OPTIONAL                             :: overwrite
integer(kind=irg)                                       :: success

integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:4)                        :: dims

integer, dimension(1:dim0,1:dim1,1:dim2,1:dim3), TARGET :: wdata
TYPE(C_PTR)                                             :: f_ptr

success = 0

dims(1:4) = (/ dim0, dim1, dim2, dim3 /)
wdata = intarr

! get a C pointer to the integer array
f_ptr = C_LOC(wdata(1,1,1,1))

! Create dataspace.
!
rnk = 4
call h5screate_simple_f(rnk, dims, space, hdferr)
!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_STD_I32LE, space, dset, hdferr)
end if 
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetIntegerArray4D:hd5create_f',.TRUE.)
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdferr )
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetIntegerArray4D:hd5write_f',.TRUE.)
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_writeDatasetIntegerArray4D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeDatasetFloat
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write a single precision array data set to the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param fltval real 
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeDatasetFloat(dataname, fltval, HDF_head, overwrite) result(success)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
real(kind=sgl),INTENT(IN)                               :: fltval
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical,INTENT(IN),OPTIONAL                             :: overwrite
integer(kind=irg)                                       :: success

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_4)
integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:1)                        :: dims

real(real_kind), dimension(1:1), TARGET                 :: wdata
TYPE(C_PTR)                                             :: f_ptr

success = 0

dims(1) = 1
wdata(1) = fltval

! get a C pointer to the integer array
f_ptr = C_LOC(wdata(1))

! Create dataspace.
!
rnk = 1
call h5screate_simple_f(rnk, dims, space, hdferr)
!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_IEEE_F32LE, space, dset, hdferr)
end if
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetFloatArray1D:hd5create_f',.TRUE.)
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_REAL, f_ptr, hdferr )
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetFloatArray1D:hd5write_f',.TRUE.)
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_writeDatasetFloat

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeDatasetDouble
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write a single precision array data set to the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param fltval real 
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeDatasetDouble(dataname, dblval, HDF_head, overwrite) result(success)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
real(kind=dbl),INTENT(IN)                               :: dblval
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical,INTENT(IN),OPTIONAL                             :: overwrite
integer(kind=irg)                                       :: success

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_8)
integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:1)                        :: dims

real(real_kind), dimension(1:1), TARGET                 :: wdata
TYPE(C_PTR)                                             :: f_ptr

success = 0

dims(1) = 1
wdata(1) = dblval

! get a C pointer to the integer array
f_ptr = C_LOC(wdata(1))

! Create dataspace.
!
rnk = 1
call h5screate_simple_f(rnk, dims, space, hdferr)
!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_IEEE_F64LE, space, dset, hdferr)
end if
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetDouble:hd5create_f',.TRUE.)
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdferr )
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetDouble:hd5write_f',.TRUE.)
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_writeDatasetDouble



!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeDatasetFloatArray1D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write a 1D single precision float array data set to the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param fltarr 1D real array
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeDatasetFloatArray1D(dataname, fltarr, dim0, HDF_head, overwrite) result(success)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
real(kind=sgl),INTENT(IN)                               :: fltarr(dim0)
integer(kind=irg),INTENT(IN)                            :: dim0
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical,INTENT(IN),OPTIONAL                             :: overwrite
integer(kind=irg)                                       :: success

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_4)
integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:1)                        :: dims

real(real_kind), dimension(1:dim0), TARGET              :: wdata
TYPE(C_PTR)                                             :: f_ptr

success = 0

dims(1) = dim0
wdata = fltarr

! get a C pointer to the integer array
f_ptr = C_LOC(wdata(1))

! Create dataspace.
!
rnk = 1
call h5screate_simple_f(rnk, dims, space, hdferr)
!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_IEEE_F32LE, space, dset, hdferr)
end if
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetFloatArray1D:hd5create_f',.TRUE.)
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_REAL, f_ptr, hdferr )
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetFloatArray1D:hd5write_f',.TRUE.)
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_writeDatasetFloatArray1D


!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeDatasetFloatArray2D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write a 2D single precision float array data set to the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param fltarr 2D real array
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeDatasetFloatArray2D(dataname, fltarr, dim0, dim1, HDF_head, overwrite) result(success)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
real(kind=sgl),INTENT(IN)                               :: fltarr(dim0, dim1)
integer(kind=irg),INTENT(IN)                            :: dim0
integer(kind=irg),INTENT(IN)                            :: dim1
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical,INTENT(IN),OPTIONAL                             :: overwrite
integer(kind=irg)                                       :: success

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_4)
integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:2)                        :: dims

real(real_kind), dimension(1:dim0,1:dim1), TARGET       :: wdata
TYPE(C_PTR)                                             :: f_ptr

success = 0

dims(1:2) = (/ dim0, dim1 /)
wdata = fltarr

! get a C pointer to the integer array
f_ptr = C_LOC(wdata(1,1))

! Create dataspace.
!
rnk = 2
call h5screate_simple_f(rnk, dims, space, hdferr)
!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_IEEE_F32LE, space, dset, hdferr)
end if
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetFloatArray1D:hd5create_f',.TRUE.)
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_REAL, f_ptr, hdferr )
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetFloatArray1D:hd5write_f',.TRUE.)
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_writeDatasetFloatArray2D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeDatasetFloatArray3D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write a 3D single precision float array data set to the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param fltarr 3D real array
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeDatasetFloatArray3D(dataname, fltarr, dim0, dim1, dim2, HDF_head, overwrite) result(success)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
real(kind=sgl),INTENT(IN)                               :: fltarr(dim0, dim1, dim2)
integer(kind=irg),INTENT(IN)                            :: dim0
integer(kind=irg),INTENT(IN)                            :: dim1
integer(kind=irg),INTENT(IN)                            :: dim2
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical,INTENT(IN),OPTIONAL                             :: overwrite
integer(kind=irg)                                       :: success

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_4)
integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:3)                        :: dims

real(real_kind), dimension(1:dim0,1:dim1,1:dim2), TARGET :: wdata
TYPE(C_PTR)                                             :: f_ptr

success = 0

dims(1:3) = (/ dim0, dim1, dim2 /)
wdata = fltarr

! get a C pointer to the integer array
f_ptr = C_LOC(wdata(1,1,1))

! Create dataspace.
!
rnk = 3
call h5screate_simple_f(rnk, dims, space, hdferr)
!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_IEEE_F32LE, space, dset, hdferr)
end if
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetFloatArray1D:hd5create_f',.TRUE.)
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_REAL, f_ptr, hdferr )
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetFloatArray1D:hd5write_f',.TRUE.)
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_writeDatasetFloatArray3D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeDatasetFloatArray4D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write a 4D single precision float array data set to the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param fltarr 4D real array
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeDatasetFloatArray4D(dataname, fltarr, dim0, dim1, dim2, dim3, HDF_head, overwrite) result(success)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
real(kind=sgl),INTENT(IN)                               :: fltarr(dim0, dim1, dim2, dim3)
integer(kind=irg),INTENT(IN)                            :: dim0
integer(kind=irg),INTENT(IN)                            :: dim1
integer(kind=irg),INTENT(IN)                            :: dim2
integer(kind=irg),INTENT(IN)                            :: dim3
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical,INTENT(IN),OPTIONAL                             :: overwrite
integer(kind=irg)                                       :: success

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_4)
integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:4)                        :: dims

real(real_kind), dimension(1:dim0,1:dim1,1:dim2,1:dim3), TARGET :: wdata
TYPE(C_PTR)                                             :: f_ptr

success = 0

dims(1:4) = (/ dim0, dim1, dim2, dim3 /)
wdata = fltarr

! get a C pointer to the integer array
f_ptr = C_LOC(wdata(1,1,1,1))

! Create dataspace.
!
rnk = 4
call h5screate_simple_f(rnk, dims, space, hdferr)
!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_IEEE_F32LE, space, dset, hdferr)
end if
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetFloatArray1D:hd5create_f',.TRUE.)
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_REAL, f_ptr, hdferr )
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetFloatArray1D:hd5write_f',.TRUE.)
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_writeDatasetFloatArray4D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeDatasetDoubleArray1D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write a 1D double precision float array data set to the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param dblarr 1D double array
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeDatasetDoubleArray1D(dataname, dblarr, dim0, HDF_head, overwrite) result(success)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
real(kind=dbl),INTENT(IN)                               :: dblarr(dim0)
integer(kind=irg),INTENT(IN)                            :: dim0
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical,INTENT(IN),OPTIONAL                             :: overwrite
integer(kind=irg)                                       :: success

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_8)
integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:1)                        :: dims

real(real_kind), dimension(1:dim0), TARGET              :: wdata
TYPE(C_PTR)                                             :: f_ptr

success = 0

dims(1) = dim0
wdata = dblarr

! get a C pointer to the integer array
f_ptr = C_LOC(wdata(1))

! Create dataspace.
!
rnk = 1
call h5screate_simple_f(rnk, dims, space, hdferr)
!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_IEEE_F64LE, space, dset, hdferr)
end if 
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetFloatArray1D:hd5create_f',.TRUE.)
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdferr )
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetFloatArray1D:hd5write_f',.TRUE.)
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_writeDatasetDoubleArray1D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeDatasetDoubleArray2D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write a 2D double precision float array data set to the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param dblarr 2D double array
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeDatasetDoubleArray2D(dataname, dblarr, dim0, dim1, HDF_head, overwrite) result(success)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
real(kind=dbl),INTENT(IN)                               :: dblarr(dim0, dim1)
integer(kind=irg),INTENT(IN)                            :: dim0
integer(kind=irg),INTENT(IN)                            :: dim1
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical,INTENT(IN),OPTIONAL                             :: overwrite
integer(kind=irg)                                       :: success

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_8)
integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:2)                        :: dims

real(real_kind), dimension(1:dim0,1:dim1), TARGET       :: wdata
TYPE(C_PTR)                                             :: f_ptr

success = 0

dims(1:2) = (/ dim0, dim1 /)
wdata = dblarr

! get a C pointer to the integer array
f_ptr = C_LOC(wdata(1,1))

! Create dataspace.
!
rnk = 2
call h5screate_simple_f(rnk, dims, space, hdferr)
!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_IEEE_F64LE, space, dset, hdferr)
end if
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetFloatArray1D:hd5create_f',.TRUE.)
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdferr )
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetFloatArray1D:hd5write_f',.TRUE.)
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_writeDatasetDoubleArray2D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeDatasetDoubleArray3D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write a 3D double precision float array data set to the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param dblarr 3D double array
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeDatasetDoubleArray3D(dataname, dblarr, dim0, dim1, dim2, HDF_head, overwrite) result(success)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
real(kind=dbl),INTENT(IN)                               :: dblarr(dim0, dim1, dim2)
integer(kind=irg),INTENT(IN)                            :: dim0
integer(kind=irg),INTENT(IN)                            :: dim1
integer(kind=irg),INTENT(IN)                            :: dim2
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical,INTENT(IN),OPTIONAL                             :: overwrite
integer(kind=irg)                                       :: success

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_8)
integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:3)                        :: dims

real(real_kind), dimension(1:dim0,1:dim1,1:dim2), TARGET  :: wdata
TYPE(C_PTR)                                             :: f_ptr

success = 0

dims(1:3) = (/ dim0, dim1, dim2 /)
wdata = dblarr

! get a C pointer to the integer array
f_ptr = C_LOC(wdata(1,1,1))

! Create dataspace.
!
rnk = 3
call h5screate_simple_f(rnk, dims, space, hdferr)
!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_IEEE_F64LE, space, dset, hdferr)
end if
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetFloatArray1D:hd5create_f',.TRUE.)
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdferr )
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetFloatArray1D:hd5write_f',.TRUE.)
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_writeDatasetDoubleArray3D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeDatasetDoubleArray4D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write a 4D double precision float array data set to the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param dblarr 4D double array
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeDatasetDoubleArray4D(dataname, dblarr, dim0, dim1, dim2, dim3, HDF_head, overwrite) result(success)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
real(kind=dbl),INTENT(IN)                               :: dblarr(dim0, dim1, dim2, dim3)
integer(kind=irg),INTENT(IN)                            :: dim0
integer(kind=irg),INTENT(IN)                            :: dim1
integer(kind=irg),INTENT(IN)                            :: dim2
integer(kind=irg),INTENT(IN)                            :: dim3
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical,INTENT(IN),OPTIONAL                             :: overwrite
integer(kind=irg)                                       :: success

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_8)
integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:4)                        :: dims

real(real_kind), dimension(1:dim0,1:dim1,1:dim2,1:dim3), TARGET  :: wdata
TYPE(C_PTR)                                             :: f_ptr

success = 0

dims(1:4) = (/ dim0, dim1, dim2, dim3 /)
wdata = dblarr

! get a C pointer to the integer array
f_ptr = C_LOC(wdata(1,1,1,1))

! Create dataspace.
!
rnk = 4
call h5screate_simple_f(rnk, dims, space, hdferr)
!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_IEEE_F64LE, space, dset, hdferr)
end if
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetFloatArray1D:hd5create_f',.TRUE.)
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdferr )
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetFloatArray1D:hd5write_f',.TRUE.)
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_writeDatasetDoubleArray4D


!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readDatasetCharArray1D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns a 1D char array data set from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param dims dimensions of the array
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readDatasetCharArray1D(dataname, dims, HDF_head) result(rdata)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
integer(HSIZE_T),INTENT(OUT)                            :: dims(1)
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
character(len=1), dimension(:), allocatable, TARGET     :: rdata

integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:2)                        :: maxdims

TYPE(C_PTR)                                             :: f_ptr

! open the data set
call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)

allocate(rdata(1:dims(1)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1))
call h5dread_f( dset, H5T_STD_U8LE, f_ptr, hdferr)
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_readDatasetCharArray1D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readDatasetCharArray2D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns a 2D char array data set from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param dims dimensions of the array
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readDatasetCharArray2D(dataname, dims, HDF_head) result(rdata)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
integer(HSIZE_T),INTENT(OUT)                            :: dims(2)
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
character(len=1), dimension(:,:), allocatable, TARGET   :: rdata

integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:2)                        :: maxdims

TYPE(C_PTR)                                             :: f_ptr

! open the data set
call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)

allocate(rdata(1:dims(1),1:dims(2)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1,1))
call h5dread_f( dset, H5T_STD_U8LE, f_ptr, hdferr)
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_readDatasetCharArray2D


!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readDatasetCharArray3D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns a 3D char array data set from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param dims dimensions of the array
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readDatasetCharArray3D(dataname, dims, HDF_head) result(rdata)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
integer(HSIZE_T),INTENT(OUT)                            :: dims(3)
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
character(len=1), dimension(:,:,:), allocatable, TARGET :: rdata

integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:2)                        :: maxdims

TYPE(C_PTR)                                             :: f_ptr

! open the data set
call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)

allocate(rdata(1:dims(1),1:dims(2),1:dims(3)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1,1,1))
call h5dread_f( dset, H5T_STD_U8LE, f_ptr, hdferr)
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_readDatasetCharArray3D


!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readDatasetCharArray4D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns a 4D char array data set from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param dims dimensions of the array
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readDatasetCharArray4D(dataname, dims, HDF_head) result(rdata)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
integer(HSIZE_T),INTENT(OUT)                            :: dims(4)
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
character(len=1), dimension(:,:,:,:), allocatable, TARGET :: rdata

integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:2)                        :: maxdims

TYPE(C_PTR)                                             :: f_ptr

! open the data set
call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)

allocate(rdata(1:dims(1),1:dims(2),1:dims(3),1:dims(4)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1,1,1,1))
call h5dread_f( dset, H5T_STD_U8LE, f_ptr, hdferr)
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_readDatasetCharArray4D


!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readDatasetInteger
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns an integer data set from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readDatasetInteger(dataname, HDF_head) result(rdata)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
integer,  TARGET                                        :: rdata

integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr

TYPE(C_PTR)                                             :: f_ptr

! open the data set
call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
!
! Read the data.
!
f_ptr = C_LOC(rdata)
call h5dread_f( dset, H5T_NATIVE_INTEGER, f_ptr, hdferr)
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_readDatasetInteger


!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readDatasetIntegerArray1D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns a 1D integer array data set from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param dims dimensions of the array
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readDatasetIntegerArray1D(dataname, dims, HDF_head) result(rdata)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
integer(HSIZE_T),INTENT(OUT)                            :: dims(1)
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
integer, dimension(:), allocatable, TARGET              :: rdata

integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:2)                        :: maxdims

TYPE(C_PTR)                                             :: f_ptr

! open the data set
call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)

allocate(rdata(1:dims(1)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1))
call h5dread_f( dset, H5T_NATIVE_INTEGER, f_ptr, hdferr)
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_readDatasetIntegerArray1D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readDatasetIntegerArray2D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns a 2D integer array data set from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param dims dimensions of the array
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readDatasetIntegerArray2D(dataname, dims, HDF_head) result(rdata)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
integer(HSIZE_T),INTENT(OUT)                            :: dims(2)
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
integer, dimension(:,:), allocatable, TARGET            :: rdata

integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:2)                        :: maxdims

TYPE(C_PTR)                                             :: f_ptr

! open the data set
call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)

allocate(rdata(1:dims(1),1:dims(2)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1,1))
call h5dread_f( dset, H5T_NATIVE_INTEGER, f_ptr, hdferr)
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_readDatasetIntegerArray2D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readDatasetIntegerArray3D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns a 3D integer array data set from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param dims dimensions of the array
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readDatasetIntegerArray3D(dataname, dims, HDF_head) result(rdata)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
integer(HSIZE_T),INTENT(OUT)                            :: dims(3)
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
integer, dimension(:,:,:), allocatable, TARGET          :: rdata

integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:2)                        :: maxdims

TYPE(C_PTR)                                             :: f_ptr

! open the data set
call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)

allocate(rdata(1:dims(1),1:dims(2),1:dims(3)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1,1,1))
call h5dread_f( dset, H5T_NATIVE_INTEGER, f_ptr, hdferr)
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_readDatasetIntegerArray3D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readDatasetIntegerArray4D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns a 4D integer array data set from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param dims dimensions of the array
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readDatasetIntegerArray4D(dataname, dims, HDF_head) result(rdata)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
integer(HSIZE_T),INTENT(OUT)                            :: dims(4)
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
integer, dimension(:,:,:,:), allocatable, TARGET        :: rdata

integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:2)                        :: maxdims

TYPE(C_PTR)                                             :: f_ptr

! open the data set
call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)

allocate(rdata(1:dims(1),1:dims(2),1:dims(3),1:dims(4)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1,1,1,1))
call h5dread_f( dset, H5T_NATIVE_INTEGER, f_ptr, hdferr)
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_readDatasetIntegerArray4D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readDatasetFloat
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns a float data set from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readDatasetFloat(dataname, HDF_head) result(rdata)

use ISO_C_BINDING

IMPLICIT NONE

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_4)

character(fnlen),INTENT(IN)                             :: dataname
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
real(real_kind), TARGET                                 :: rdata

integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr

TYPE(C_PTR)                                             :: f_ptr

! open the data set
call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)

! Read the data.
!
f_ptr = C_LOC(rdata)
call h5dread_f( dset, H5T_NATIVE_REAL, f_ptr, hdferr)
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_readDatasetFloat

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readDatasetFloatArray1D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns a 1D float array data set from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param dims dimensions of the array
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readDatasetFloatArray1D(dataname, dims, HDF_head) result(rdata)

use ISO_C_BINDING

IMPLICIT NONE

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_4)

character(fnlen),INTENT(IN)                             :: dataname
integer(HSIZE_T),INTENT(OUT)                            :: dims(1)
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
real(real_kind), dimension(:), allocatable, TARGET      :: rdata

integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:2)                        :: maxdims

TYPE(C_PTR)                                             :: f_ptr

! open the data set
call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)

allocate(rdata(1:dims(1)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1))
call h5dread_f( dset, H5T_NATIVE_REAL, f_ptr, hdferr)
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_readDatasetFloatArray1D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readDatasetFloatArray2D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns a 2D float array data set from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param dims dimensions of the array
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readDatasetFloatArray2D(dataname, dims, HDF_head) result(rdata)

use ISO_C_BINDING

IMPLICIT NONE

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_4)

character(fnlen),INTENT(IN)                             :: dataname
integer(HSIZE_T),INTENT(OUT)                            :: dims(2)
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
real(real_kind), dimension(:,:), allocatable, TARGET    :: rdata

integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:2)                        :: maxdims

TYPE(C_PTR)                                             :: f_ptr

! open the data set
call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)

allocate(rdata(1:dims(1),1:dims(2)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1,1))
call h5dread_f( dset, H5T_NATIVE_REAL, f_ptr, hdferr)
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_readDatasetFloatArray2D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readDatasetFloatArray3D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns a 3D float array data set from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param dims dimensions of the array
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readDatasetFloatArray3D(dataname, dims, HDF_head) result(rdata)

use ISO_C_BINDING

IMPLICIT NONE

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_4)

character(fnlen),INTENT(IN)                             :: dataname
integer(HSIZE_T),INTENT(OUT)                            :: dims(3)
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
real(real_kind), dimension(:,:,:), allocatable, TARGET  :: rdata

integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:2)                        :: maxdims

TYPE(C_PTR)                                             :: f_ptr

! open the data set
call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)

allocate(rdata(1:dims(1),1:dims(2),1:dims(3)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1,1,1))
call h5dread_f( dset, H5T_NATIVE_REAL, f_ptr, hdferr)
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_readDatasetFloatArray3D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readDatasetFloatArray4D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns a 4D float array data set from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param dims dimensions of the array
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readDatasetFloatArray4D(dataname, dims, HDF_head) result(rdata)

use ISO_C_BINDING

IMPLICIT NONE

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_4)

character(fnlen),INTENT(IN)                             :: dataname
integer(HSIZE_T),INTENT(OUT)                            :: dims(4)
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
real(real_kind), dimension(:,:,:,:), allocatable, TARGET:: rdata

integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:2)                        :: maxdims

TYPE(C_PTR)                                             :: f_ptr

! open the data set
call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)

allocate(rdata(1:dims(1),1:dims(2),1:dims(3),1:dims(4)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1,1,1,1))
call h5dread_f( dset, H5T_NATIVE_REAL, f_ptr, hdferr)
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_readDatasetFloatArray4D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readDatasetDouble
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns a double data set from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readDatasetDouble(dataname, HDF_head) result(rdata)

use ISO_C_BINDING

IMPLICIT NONE

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_8)

character(fnlen),INTENT(IN)                             :: dataname
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
real(real_kind), TARGET                                 :: rdata

integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr

TYPE(C_PTR)                                             :: f_ptr

! open the data set
call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)

! Read the data.
!
f_ptr = C_LOC(rdata)
call h5dread_f( dset, H5T_NATIVE_DOUBLE, f_ptr, hdferr)
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_readDatasetDouble


!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readDatasetDoubleArray1D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns a 1D double array data set from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param dims dimensions of the array
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readDatasetDoubleArray1D(dataname, dims, HDF_head) result(rdata)

use ISO_C_BINDING

IMPLICIT NONE

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_8)

character(fnlen),INTENT(IN)                             :: dataname
integer(HSIZE_T),INTENT(OUT)                            :: dims(1)
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
real(real_kind), dimension(:), allocatable, TARGET      :: rdata

integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:2)                        :: maxdims

TYPE(C_PTR)                                             :: f_ptr

! open the data set
call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)

allocate(rdata(1:dims(1)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1))
call h5dread_f( dset, H5T_NATIVE_DOUBLE, f_ptr, hdferr)
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_readDatasetDoubleArray1D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readDatasetDoubleArray2D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns a 2D double array data set from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param dims dimensions of the array
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readDatasetDoubleArray2D(dataname, dims, HDF_head) result(rdata)

use ISO_C_BINDING

IMPLICIT NONE

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_8)

character(fnlen),INTENT(IN)                             :: dataname
integer(HSIZE_T),INTENT(OUT)                            :: dims(2)
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
real(real_kind), dimension(:,:), allocatable, TARGET    :: rdata

integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:2)                        :: maxdims

TYPE(C_PTR)                                             :: f_ptr

! open the data set
call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)

allocate(rdata(1:dims(1),1:dims(2)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1,1))
call h5dread_f( dset, H5T_NATIVE_DOUBLE, f_ptr, hdferr)
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_readDatasetDoubleArray2D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readDatasetDoubleArray3D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns a 3D double array data set from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param dims dimensions of the array
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readDatasetDoubleArray3D(dataname, dims, HDF_head) result(rdata)

use ISO_C_BINDING

IMPLICIT NONE

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_8)

character(fnlen),INTENT(IN)                             :: dataname
integer(HSIZE_T),INTENT(OUT)                            :: dims(3)
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
real(real_kind), dimension(:,:,:), allocatable, TARGET  :: rdata

integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:2)                        :: maxdims

TYPE(C_PTR)                                             :: f_ptr

! open the data set
call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)

allocate(rdata(1:dims(1),1:dims(2),1:dims(3)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1,1,1))
call h5dread_f( dset, H5T_NATIVE_DOUBLE, f_ptr, hdferr)

!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_readDatasetDoubleArray3D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readDatasetDoubleArray4D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns a 4D double array data set from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param dims dimensions of the array
!> @param HDF_head
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readDatasetDoubleArray4D(dataname, dims, HDF_head) result(rdata)

use ISO_C_BINDING

IMPLICIT NONE

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_8)

character(fnlen),INTENT(IN)                             :: dataname
integer(HSIZE_T),INTENT(OUT)                            :: dims(4)
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
real(real_kind), dimension(:,:,:,:), allocatable, TARGET:: rdata

integer(HID_T)                                          :: space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:2)                        :: maxdims

TYPE(C_PTR)                                             :: f_ptr

! open the data set
call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)

allocate(rdata(1:dims(1),1:dims(2),1:dims(3),1:dims(4)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1,1,1,1))
call h5dread_f( dset, H5T_NATIVE_DOUBLE, f_ptr, hdferr)
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call h5sclose_f(space, hdferr)

! that's it

end function HDF_readDatasetDoubleArray4D

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! hyperslab read and write routines
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeHyperslabCharArray2D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief writes a 2D hyperslab to the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param wdata data to be written 
!> @param hdims original dimensions of complete dataset
!> @param offset offset of the hyperslab
!> @param dim0 ...  dimensions of the hyperslab
!> @param HDF_head
!> @param insert (optional) if present, add to existing dataset
!
!> @date 04/06/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeHyperslabCharArray2D(dataname, wdata, hdims, offset, &
                                       dim0, dim1, HDF_head, insert) result(success)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
character(kind=c_char),INTENT(IN),TARGET                :: wdata(dim0,dim1)
integer(HSIZE_T),INTENT(IN)                             :: hdims(2)
integer(HSIZE_T),INTENT(IN)                             :: offset(2)
integer(HSIZE_T),INTENT(IN)                             :: dim0
integer(HSIZE_T),INTENT(IN)                             :: dim1
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical, OPTIONAL, INTENT(IN)                           :: insert
integer(kind=irg)                                       :: success

integer(HID_T)                                          :: memspace, space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T)                                        :: dims(2)

TYPE(C_PTR)                                             :: f_ptr

success = 0

dims = (/ dim0, dim1 /)
rnk = 2
call h5screate_simple_f(rnk, hdims, space, hdferr)

if (present(insert)) then
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_STD_U8LE, space, dset, hdferr)
end if

f_ptr = C_LOC(wdata(1,1))

call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr)
call h5screate_simple_f(rnk, dims, memspace, hdferr)
call h5dwrite_f(dset, H5T_STD_U8LE, f_ptr, hdferr, memspace, space)

call h5dclose_f(dset, hdferr)

end function HDF_writeHyperslabCharArray2D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeHyperslabCharArray3D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief writes a 3D hyperslab to the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param wdata data to be written 
!> @param hdims original dimensions of complete dataset
!> @param offset offset of the hyperslab
!> @param dim0 ...  dimensions of the hyperslab
!> @param HDF_head
!> @param insert (optional) if present, add to existing dataset
!
!> @date 04/06/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeHyperslabCharArray3D(dataname, wdata, hdims, offset, &
                                       dim0, dim1, dim2, HDF_head, insert) result(success)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
character(kind=c_char),INTENT(IN)                       :: wdata(dim0,dim1,dim2)
integer(HSIZE_T),INTENT(IN)                             :: hdims(3)
integer(HSIZE_T),INTENT(IN)                             :: offset(3)
integer(HSIZE_T),INTENT(IN)                             :: dim0
integer(HSIZE_T),INTENT(IN)                             :: dim1
integer(HSIZE_T),INTENT(IN)                             :: dim2
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical, OPTIONAL, INTENT(IN)                           :: insert
integer(kind=irg)                                       :: success

integer(HID_T)                                          :: memspace, space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T)                                        :: dims(3)

success = 0

dims = (/ dim0, dim1, dim2 /)
rnk = 3
call h5screate_simple_f(rnk, hdims, space, hdferr)

if (present(insert)) then
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_STD_U8LE, space, dset, hdferr)
end if

call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr)
call h5screate_simple_f(rnk, dims, memspace, hdferr)
call h5dwrite_f(dset, H5T_STD_U8LE, wdata, dims, hdferr, memspace, space)

call h5dclose_f(dset, hdferr)

end function HDF_writeHyperslabCharArray3D


!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeHyperslabCharArray4D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief writes a 4D hyperslab to the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param wdata data to be written 
!> @param hdims original dimensions of complete dataset
!> @param offset offset of the hyperslab
!> @param dim0 ...  dimensions of the hyperslab
!> @param HDF_head
!> @param insert (optional) if present, add to existing dataset
!
!> @date 04/06/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeHyperslabCharArray4D(dataname, wdata, hdims, offset, &
                                       dim0, dim1, dim2, dim3, HDF_head, insert) result(success)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
character(kind=c_char),INTENT(IN)                       :: wdata(dim0,dim1,dim2,dim3)
integer(HSIZE_T),INTENT(IN)                             :: hdims(4)
integer(HSIZE_T),INTENT(IN)                             :: offset(4)
integer(HSIZE_T),INTENT(IN)                             :: dim0
integer(HSIZE_T),INTENT(IN)                             :: dim1
integer(HSIZE_T),INTENT(IN)                             :: dim2
integer(HSIZE_T),INTENT(IN)                             :: dim3
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical, OPTIONAL, INTENT(IN)                           :: insert
integer(kind=irg)                                       :: success

integer(HID_T)                                          :: memspace, space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T)                                        :: dims(4)

success = 0

dims = (/ dim0, dim1, dim2, dim3 /)
rnk = 4
call h5screate_simple_f(rnk, hdims, space, hdferr)

if (present(insert)) then
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_STD_U8LE, space, dset, hdferr)
end if

call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr)
call h5screate_simple_f(rnk, dims, memspace, hdferr)
call h5dwrite_f(dset, H5T_STD_U8LE, wdata, dims, hdferr, memspace, space)

call h5dclose_f(dset, hdferr)

end function HDF_writeHyperslabCharArray4D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeHyperslabIntegerArray2D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief writes a 2D hyperslab to the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param wdata data to be written 
!> @param hdims original dimensions of complete dataset
!> @param offset offset of the hyperslab
!> @param dim0 ...  dimensions of the hyperslab
!> @param HDF_head
!> @param insert (optional) if present, add to existing dataset
!
!> @date 04/06/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeHyperslabIntegerArray2D(dataname, wdata, hdims, offset, &
                                          dim0, dim1, HDF_head, insert) result(success)

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
integer(kind=irg),INTENT(IN)                            :: wdata(dim0,dim1)
integer(HSIZE_T),INTENT(IN)                             :: hdims(2)
integer(HSIZE_T),INTENT(IN)                             :: offset(2)
integer(HSIZE_T),INTENT(IN)                             :: dim0
integer(HSIZE_T),INTENT(IN)                             :: dim1
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical, OPTIONAL, INTENT(IN)                           :: insert
integer(kind=irg)                                       :: success

integer(HID_T)                                          :: memspace, space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T)                                        :: dims(2)

success = 0

dims = (/ dim0, dim1 /)
rnk = 2
call h5screate_simple_f(rnk, hdims, space, hdferr)

if (present(insert)) then
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_STD_I32LE, space, dset, hdferr)
end if

call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr)
call h5screate_simple_f(rnk, dims, memspace, hdferr)
call h5dwrite_f(dset, H5T_NATIVE_INTEGER, wdata, dims, hdferr, memspace, space)

call h5dclose_f(dset, hdferr)

end function HDF_writeHyperslabIntegerArray2D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeHyperslabIntegerArray3D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief writes a 3D hyperslab to the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param wdata data to be written 
!> @param hdims original dimensions of complete dataset
!> @param offset offset of the hyperslab
!> @param dim0 ...  dimensions of the hyperslab
!> @param HDF_head
!> @param insert (optional) if present, add to existing dataset
!
!> @date 04/06/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeHyperslabIntegerArray3D(dataname, wdata, hdims, offset, &
                                          dim0, dim1, dim2, HDF_head, insert) result(success)

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
integer(kind=irg),INTENT(IN)                            :: wdata(dim0,dim1,dim2)
integer(HSIZE_T),INTENT(IN)                             :: hdims(3)
integer(HSIZE_T),INTENT(IN)                             :: offset(3)
integer(HSIZE_T),INTENT(IN)                             :: dim0
integer(HSIZE_T),INTENT(IN)                             :: dim1
integer(HSIZE_T),INTENT(IN)                             :: dim2
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical, OPTIONAL, INTENT(IN)                           :: insert
integer(kind=irg)                                       :: success

integer(HID_T)                                          :: memspace, space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T)                                        :: dims(3)

success = 0

dims = (/ dim0, dim1, dim2 /)
rnk = 3
call h5screate_simple_f(rnk, hdims, space, hdferr)

if (present(insert)) then
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_STD_I32LE, space, dset, hdferr)
end if

call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr)
call h5screate_simple_f(rnk, dims, memspace, hdferr)
call h5dwrite_f(dset, H5T_NATIVE_INTEGER, wdata, dims, hdferr, memspace, space)

call h5dclose_f(dset, hdferr)

end function HDF_writeHyperslabIntegerArray3D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeHyperslabIntegerArray4D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief writes a 4D hyperslab to the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param wdata data to be written 
!> @param hdims original dimensions of complete dataset
!> @param offset offset of the hyperslab
!> @param dim0 ...  dimensions of the hyperslab
!> @param HDF_head
!> @param insert (optional) if present, add to existing dataset
!
!> @date 04/06/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeHyperslabIntegerArray4D(dataname, wdata, hdims, offset, &
                                          dim0, dim1, dim2, dim3, HDF_head, insert) result(success)

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
integer(kind=irg),INTENT(IN)                            :: wdata(dim0,dim1,dim2,dim3)
integer(HSIZE_T),INTENT(IN)                             :: hdims(4)
integer(HSIZE_T),INTENT(IN)                             :: offset(4)
integer(HSIZE_T),INTENT(IN)                             :: dim0
integer(HSIZE_T),INTENT(IN)                             :: dim1
integer(HSIZE_T),INTENT(IN)                             :: dim2
integer(HSIZE_T),INTENT(IN)                             :: dim3
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical, OPTIONAL, INTENT(IN)                           :: insert
integer(kind=irg)                                       :: success

integer(HID_T)                                          :: memspace, space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T)                                        :: dims(4)

success = 0

dims = (/ dim0, dim1, dim2, dim3 /)
rnk = 4
call h5screate_simple_f(rnk, hdims, space, hdferr)

if (present(insert)) then
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_STD_I32LE, space, dset, hdferr)
end if

call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr)
call h5screate_simple_f(rnk, dims, memspace, hdferr)
call h5dwrite_f(dset, H5T_NATIVE_INTEGER, wdata, dims, hdferr, memspace, space)

call h5dclose_f(dset, hdferr)

end function HDF_writeHyperslabIntegerArray4D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeHyperslabFloatArray2D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief writes a 2D hyperslab to the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param wdata data to be written 
!> @param hdims original dimensions of complete dataset
!> @param offset offset of the hyperslab
!> @param dim0 ...  dimensions of the hyperslab
!> @param HDF_head
!> @param insert (optional) if present, add to existing dataset
!
!> @date 04/06/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeHyperslabFloatArray2D(dataname, wdata, hdims, offset, &
                                          dim0, dim1, HDF_head, insert) result(success)

IMPLICIT NONE

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_4)

character(fnlen),INTENT(IN)                             :: dataname
real(real_kind),INTENT(IN)                              :: wdata(dim0,dim1)
integer(HSIZE_T),INTENT(IN)                             :: hdims(2)
integer(HSIZE_T),INTENT(IN)                             :: offset(2)
integer(HSIZE_T),INTENT(IN)                             :: dim0
integer(HSIZE_T),INTENT(IN)                             :: dim1
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical, OPTIONAL, INTENT(IN)                           :: insert
integer(kind=irg)                                       :: success

integer(HID_T)                                          :: memspace, space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T)                                        :: dims(2)

success = 0

dims = (/ dim0, dim1 /)
rnk = 2
call h5screate_simple_f(rnk, hdims, space, hdferr)

if (present(insert)) then
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_IEEE_F32LE, space, dset, hdferr)
end if

call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr)
call h5screate_simple_f(rnk, dims, memspace, hdferr)
call h5dwrite_f(dset, H5T_NATIVE_REAL, wdata, dims, hdferr, memspace, space)

call h5dclose_f(dset, hdferr)

end function HDF_writeHyperslabFloatArray2D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeHyperslabFloatArray3D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief writes a 3D hyperslab to the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param wdata data to be written 
!> @param hdims original dimensions of complete dataset
!> @param offset offset of the hyperslab
!> @param dim0 ...  dimensions of the hyperslab
!> @param HDF_head
!> @param insert (optional) if present, add to existing dataset
!
!> @date 04/06/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeHyperslabFloatArray3D(dataname, wdata, hdims, offset, &
                                          dim0, dim1, dim2, HDF_head, insert) result(success)

IMPLICIT NONE

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_4)

character(fnlen),INTENT(IN)                             :: dataname
real(real_kind),INTENT(IN)                              :: wdata(dim0,dim1,dim2)
integer(HSIZE_T),INTENT(IN)                             :: hdims(3)
integer(HSIZE_T),INTENT(IN)                             :: offset(3)
integer(HSIZE_T),INTENT(IN)                             :: dim0
integer(HSIZE_T),INTENT(IN)                             :: dim1
integer(HSIZE_T),INTENT(IN)                             :: dim2
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical, OPTIONAL, INTENT(IN)                           :: insert
integer(kind=irg)                                       :: success

integer(HID_T)                                          :: memspace, space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T)                                        :: dims(3)

success = 0

dims = (/ dim0, dim1, dim2 /)
rnk = 3
call h5screate_simple_f(rnk, hdims, space, hdferr)

if (present(insert)) then
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_IEEE_F32LE, space, dset, hdferr)
end if

call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr)
call h5screate_simple_f(rnk, dims, memspace, hdferr)
call h5dwrite_f(dset, H5T_NATIVE_REAL, wdata, dims, hdferr, memspace, space)

call h5dclose_f(dset, hdferr)

end function HDF_writeHyperslabFloatArray3D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeHyperslabFloatArray4D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief writes a 4D hyperslab to the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param wdata data to be written 
!> @param hdims original dimensions of complete dataset
!> @param offset offset of the hyperslab
!> @param dim0 ...  dimensions of the hyperslab
!> @param HDF_head
!> @param insert (optional) if present, add to existing dataset
!
!> @date 04/06/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeHyperslabFloatArray4D(dataname, wdata, hdims, offset, &
                                          dim0, dim1, dim2, dim3, HDF_head, insert) result(success)

IMPLICIT NONE

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_4)

character(fnlen),INTENT(IN)                             :: dataname
real(real_kind),INTENT(IN)                              :: wdata(dim0,dim1,dim2,dim3)
integer(HSIZE_T),INTENT(IN)                             :: hdims(4)
integer(HSIZE_T),INTENT(IN)                             :: offset(4)
integer(HSIZE_T),INTENT(IN)                             :: dim0
integer(HSIZE_T),INTENT(IN)                             :: dim1
integer(HSIZE_T),INTENT(IN)                             :: dim2
integer(HSIZE_T),INTENT(IN)                             :: dim3
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical, OPTIONAL, INTENT(IN)                           :: insert
integer(kind=irg)                                       :: success

integer(HID_T)                                          :: memspace, space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T)                                        :: dims(4)

success = 0

dims = (/ dim0, dim1, dim2, dim3 /)
rnk = 4
call h5screate_simple_f(rnk, hdims, space, hdferr)

if (present(insert)) then
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_IEEE_F32LE, space, dset, hdferr)
end if

call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr)
call h5screate_simple_f(rnk, dims, memspace, hdferr)
call h5dwrite_f(dset, H5T_NATIVE_REAL, wdata, dims, hdferr, memspace, space)

call h5dclose_f(dset, hdferr)

end function HDF_writeHyperslabFloatArray4D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeHyperslabDoubleArray2D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief writes a 2D hyperslab to the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param wdata data to be written 
!> @param hdims original dimensions of complete dataset
!> @param offset offset of the hyperslab
!> @param dim0 ...  dimensions of the hyperslab
!> @param HDF_head
!> @param insert (optional) if present, add to existing dataset
!
!> @date 04/06/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeHyperslabDoubleArray2D(dataname, wdata, hdims, offset, &
                                          dim0, dim1, HDF_head, insert) result(success)

IMPLICIT NONE

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_8)

character(fnlen),INTENT(IN)                             :: dataname
real(real_kind),INTENT(IN)                              :: wdata(dim0,dim1)
integer(HSIZE_T),INTENT(IN)                             :: hdims(2)
integer(HSIZE_T),INTENT(IN)                             :: offset(2)
integer(HSIZE_T),INTENT(IN)                             :: dim0
integer(HSIZE_T),INTENT(IN)                             :: dim1
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical, OPTIONAL, INTENT(IN)                           :: insert
integer(kind=irg)                                       :: success

integer(HID_T)                                          :: memspace, space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T)                                        :: dims(2)

success = 0

dims = (/ dim0, dim1 /)
rnk = 2
call h5screate_simple_f(rnk, hdims, space, hdferr)

if (present(insert)) then
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_IEEE_F64LE, space, dset, hdferr)
end if

call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr)
call h5screate_simple_f(rnk, dims, memspace, hdferr)
call h5dwrite_f(dset, H5T_NATIVE_REAL, wdata, dims, hdferr, memspace, space)

call h5dclose_f(dset, hdferr)

end function HDF_writeHyperslabDoubleArray2D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeHyperslabDoubleArray3D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief writes a 3D hyperslab to the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param wdata data to be written 
!> @param hdims original dimensions of complete dataset
!> @param offset offset of the hyperslab
!> @param dim0 ...  dimensions of the hyperslab
!> @param HDF_head
!> @param insert (optional) if present, add to existing dataset
!
!> @date 04/06/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeHyperslabDoubleArray3D(dataname, wdata, hdims, offset, &
                                          dim0, dim1, dim2, HDF_head, insert) result(success)

IMPLICIT NONE

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_8)

character(fnlen),INTENT(IN)                             :: dataname
real(real_kind),INTENT(IN)                              :: wdata(dim0,dim1,dim2)
integer(HSIZE_T),INTENT(IN)                             :: hdims(3)
integer(HSIZE_T),INTENT(IN)                             :: offset(3)
integer(HSIZE_T),INTENT(IN)                             :: dim0
integer(HSIZE_T),INTENT(IN)                             :: dim1
integer(HSIZE_T),INTENT(IN)                             :: dim2
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical, OPTIONAL, INTENT(IN)                           :: insert
integer(kind=irg)                                       :: success

integer(HID_T)                                          :: memspace, space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T)                                        :: dims(3)

success = 0

dims = (/ dim0, dim1, dim2 /)
rnk = 3
call h5screate_simple_f(rnk, hdims, space, hdferr)

if (present(insert)) then
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_IEEE_F64LE, space, dset, hdferr)
end if

call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr)
call h5screate_simple_f(rnk, dims, memspace, hdferr)
call h5dwrite_f(dset, H5T_NATIVE_REAL, wdata, dims, hdferr, memspace, space)

call h5dclose_f(dset, hdferr)

end function HDF_writeHyperslabDoubleArray3D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_writeHyperslabDoubleArray4D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief writes a 4D hyperslab to the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param wdata data to be written 
!> @param hdims original dimensions of complete dataset
!> @param offset offset of the hyperslab
!> @param dim0 ...  dimensions of the hyperslab
!> @param HDF_head
!> @param insert (optional) if present, add to existing dataset
!
!> @date 04/06/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeHyperslabDoubleArray4D(dataname, wdata, hdims, offset, &
                                          dim0, dim1, dim2, dim3, HDF_head, insert) result(success)

IMPLICIT NONE

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_8)

character(fnlen),INTENT(IN)                             :: dataname
real(real_kind),INTENT(IN)                              :: wdata(dim0,dim1,dim2,dim3)
integer(HSIZE_T),INTENT(IN)                             :: hdims(4)
integer(HSIZE_T),INTENT(IN)                             :: offset(4)
integer(HSIZE_T),INTENT(IN)                             :: dim0
integer(HSIZE_T),INTENT(IN)                             :: dim1
integer(HSIZE_T),INTENT(IN)                             :: dim2
integer(HSIZE_T),INTENT(IN)                             :: dim3
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
logical, OPTIONAL, INTENT(IN)                           :: insert
integer(kind=irg)                                       :: success

integer(HID_T)                                          :: memspace, space, dset ! Handles
integer                                                 :: hdferr, rnk
integer(HSIZE_T)                                        :: dims(4)

success = 0

dims = (/ dim0, dim1, dim2, dim3 /)
rnk = 4
call h5screate_simple_f(rnk, hdims, space, hdferr)

if (present(insert)) then
  call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
else
  call h5dcreate_f(HDF_head%objectID, trim(dataname), H5T_IEEE_F64LE, space, dset, hdferr)
end if

call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr)
call h5screate_simple_f(rnk, dims, memspace, hdferr)
call h5dwrite_f(dset, H5T_NATIVE_REAL, wdata, dims, hdferr, memspace, space)

call h5dclose_f(dset, hdferr)

end function HDF_writeHyperslabDoubleArray4D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readHyperslabCharArray2D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns a 2D hyperslab from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param offset offset of the hyperslab
!> @param dims dimensions of the hyperslab
!> @param HDF_head
!
!> @date 04/06/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readHyperslabCharArray2D(dataname, offset, dims, HDF_head) result(rdata)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
integer(HSIZE_T),INTENT(IN)                             :: offset(2)
integer(HSIZE_T),INTENT(IN)                             :: dims(2)
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
character(len=1,kind=c_char), dimension(:,:), allocatable, TARGET   :: rdata

integer(HID_T)                                          :: memspace, space, dset ! Handles
integer(HSIZE_T)                                        :: hdims(2), max_dims(2)
integer                                                 :: hdferr, rnk

allocate(rdata(1:dims(1),1:dims(2)))

call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
call h5dget_space_f(dset, space, hdferr)

call h5sget_simple_extent_dims_f(space, hdims, max_dims, hdferr)

rnk = 2
call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr) 
call h5screate_simple_f(rnk, dims, memspace, hdferr)

call h5dread_f(dset, H5T_STD_U8LE, rdata, hdims, hdferr, memspace, space)

call h5sclose_f(space, hdferr)
call h5dclose_f(dset, hdferr)

end function HDF_readHyperslabCharArray2D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readHyperslabCharArray3D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns a 3D hyperslab from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param offset offset of the hyperslab
!> @param dims dimensions of the hyperslab
!> @param HDF_head
!
!> @date 04/06/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readHyperslabCharArray3D(dataname, offset, dims, HDF_head) result(rdata)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
integer(HSIZE_T),INTENT(IN)                             :: offset(3)
integer(HSIZE_T),INTENT(IN)                             :: dims(3)
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
character(len=1,kind=c_char), dimension(:,:,:), allocatable, TARGET   :: rdata

integer(HID_T)                                          :: memspace, space, dset ! Handles
integer(HSIZE_T)                                        :: hdims(3), max_dims(3)
integer                                                 :: hdferr, rnk

allocate(rdata(1:dims(1),1:dims(2),1:dims(3)))

call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
call h5dget_space_f(dset, space, hdferr)

call h5sget_simple_extent_dims_f(space, hdims, max_dims, hdferr)

rnk = 3
call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr) 
call h5screate_simple_f(rnk, dims, memspace, hdferr)

call h5dread_f(dset, H5T_STD_U8LE, rdata, hdims, hdferr, memspace, space)

call h5sclose_f(space, hdferr)
call h5dclose_f(dset, hdferr)

end function HDF_readHyperslabCharArray3D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readHyperslabCharArray4D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns a 4D hyperslab from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param offset offset of the hyperslab
!> @param dims dimensions of the hyperslab
!> @param HDF_head
!
!> @date 04/06/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readHyperslabCharArray4D(dataname, offset, dims, HDF_head) result(rdata)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
integer(HSIZE_T),INTENT(IN)                             :: offset(4)
integer(HSIZE_T),INTENT(IN)                             :: dims(4)
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
character(len=1,kind=c_char), dimension(:,:,:,:), allocatable, TARGET   :: rdata

integer(HID_T)                                          :: memspace, space, dset ! Handles
integer(HSIZE_T)                                        :: hdims(4), max_dims(4)
integer                                                 :: hdferr, rnk

allocate(rdata(1:dims(1),1:dims(2),1:dims(3),1:dims(4)))

call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
call h5dget_space_f(dset, space, hdferr)

call h5sget_simple_extent_dims_f(space, hdims, max_dims, hdferr)

rnk = 4
call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr) 
call h5screate_simple_f(rnk, dims, memspace, hdferr)

call h5dread_f(dset, H5T_STD_U8LE, rdata, hdims, hdferr, memspace, space)

call h5sclose_f(space, hdferr)
call h5dclose_f(dset, hdferr)

end function HDF_readHyperslabCharArray4D


!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readHyperslabIntegerArray2D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns a 2D hyperslab from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param offset offset of the hyperslab
!> @param dims dimensions of the hyperslab
!> @param HDF_head
!
!> @date 04/06/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readHyperslabIntegerArray2D(dataname, offset, dims, HDF_head) result(rdata)

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
integer(HSIZE_T),INTENT(IN)                             :: offset(2)
integer(HSIZE_T),INTENT(IN)                             :: dims(2)
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
integer, dimension(:,:), allocatable, TARGET            :: rdata

integer(HID_T)                                          :: memspace, space, dset ! Handles
integer(HSIZE_T)                                        :: hdims(2), max_dims(2)
integer                                                 :: hdferr, rnk

call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
call h5dget_space_f(dset, space, hdferr)
call h5sget_simple_extent_dims_f(space, hdims, max_dims, hdferr)

rnk = 2
call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr) 
call h5screate_simple_f(rnk, dims, memspace, hdferr)

allocate(rdata(dims(1),dims(2)))
call h5dread_f(dset, H5T_NATIVE_INTEGER, rdata, hdims, hdferr, memspace, space)

call h5sclose_f(space, hdferr)
call h5dclose_f(dset, hdferr)

end function HDF_readHyperslabIntegerArray2D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readHyperslabIntegerArray3D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns a 3D hyperslab from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param offset offset of the hyperslab
!> @param dims dimensions of the hyperslab
!> @param HDF_head
!
!> @date 04/06/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readHyperslabIntegerArray3D(dataname, offset, dims, HDF_head) result(rdata)

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
integer(HSIZE_T),INTENT(IN)                             :: offset(3)
integer(HSIZE_T),INTENT(IN)                             :: dims(3)
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
integer, dimension(:,:,:), allocatable, TARGET          :: rdata

integer(HID_T)                                          :: memspace, space, dset ! Handles
integer(HSIZE_T)                                        :: hdims(2), max_dims(2)
integer                                                 :: hdferr, rnk

call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
call h5dget_space_f(dset, space, hdferr)
call h5sget_simple_extent_dims_f(space, hdims, max_dims, hdferr)

rnk = 3
call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr) 
call h5screate_simple_f(rnk, dims, memspace, hdferr)

allocate(rdata(dims(1),dims(2),dims(3)))
call h5dread_f(dset, H5T_NATIVE_INTEGER, rdata, hdims, hdferr, memspace, space)

call h5sclose_f(space, hdferr)
call h5dclose_f(dset, hdferr)

end function HDF_readHyperslabIntegerArray3D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readHyperslabIntegerArray4D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns a 4D hyperslab from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param offset offset of the hyperslab
!> @param dims dimensions of the hyperslab
!> @param HDF_head
!
!> @date 04/06/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readHyperslabIntegerArray4D(dataname, offset, dims, HDF_head) result(rdata)

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
integer(HSIZE_T),INTENT(IN)                             :: offset(4)
integer(HSIZE_T),INTENT(IN)                             :: dims(4)
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
integer, dimension(:,:,:,:), allocatable, TARGET        :: rdata

integer(HID_T)                                          :: memspace, space, dset ! Handles
integer(HSIZE_T)                                        :: hdims(2), max_dims(2)
integer                                                 :: hdferr, rnk

call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
call h5dget_space_f(dset, space, hdferr)
call h5sget_simple_extent_dims_f(space, hdims, max_dims, hdferr)

rnk = 4
call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr) 
call h5screate_simple_f(rnk, dims, memspace, hdferr)

allocate(rdata(dims(1),dims(2),dims(3),dims(4)))
call h5dread_f(dset, H5T_NATIVE_INTEGER, rdata, hdims, hdferr, memspace, space)

call h5sclose_f(space, hdferr)
call h5dclose_f(dset, hdferr)

end function HDF_readHyperslabIntegerArray4D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readHyperslabFloatArray2D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns a 2D hyperslab from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param offset offset of the hyperslab
!> @param dims dimensions of the hyperslab
!> @param HDF_head
!
!> @date 04/06/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readHyperslabFloatArray2D(dataname, offset, dims, HDF_head) result(rdata)

IMPLICIT NONE

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_4)

character(fnlen),INTENT(IN)                             :: dataname
integer(HSIZE_T),INTENT(IN)                             :: offset(2)
integer(HSIZE_T),INTENT(IN)                             :: dims(2)
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
real(real_kind), dimension(:,:), allocatable, TARGET    :: rdata

integer(HID_T)                                          :: memspace, space, dset ! Handles
integer(HSIZE_T)                                        :: hdims(2), max_dims(2)
integer                                                 :: hdferr, rnk

call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
call h5dget_space_f(dset, space, hdferr)
call h5sget_simple_extent_dims_f(space, hdims, max_dims, hdferr)

rnk = 2
call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr) 
call h5screate_simple_f(rnk, dims, memspace, hdferr)

allocate(rdata(dims(1),dims(2)))
call h5dread_f(dset, H5T_NATIVE_REAL, rdata, hdims, hdferr, memspace, space)

call h5sclose_f(space, hdferr)
call h5dclose_f(dset, hdferr)

end function HDF_readHyperslabFloatArray2D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readHyperslabFloatArray3D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns a 3D hyperslab from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param offset offset of the hyperslab
!> @param dims dimensions of the hyperslab
!> @param HDF_head
!
!> @date 04/06/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readHyperslabFloatArray3D(dataname, offset, dims, HDF_head) result(rdata)

IMPLICIT NONE

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_4)

character(fnlen),INTENT(IN)                             :: dataname
integer(HSIZE_T),INTENT(IN)                             :: offset(3)
integer(HSIZE_T),INTENT(IN)                             :: dims(3)
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
real(real_kind), dimension(:,:,:), allocatable, TARGET  :: rdata

integer(HID_T)                                          :: memspace, space, dset ! Handles
integer(HSIZE_T)                                        :: hdims(2), max_dims(2)
integer                                                 :: hdferr, rnk

call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
call h5dget_space_f(dset, space, hdferr)
call h5sget_simple_extent_dims_f(space, hdims, max_dims, hdferr)

rnk = 3
call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr) 
call h5screate_simple_f(rnk, dims, memspace, hdferr)

allocate(rdata(dims(1),dims(2),dims(3)))
call h5dread_f(dset, H5T_NATIVE_REAL, rdata, hdims, hdferr, memspace, space)

call h5sclose_f(space, hdferr)
call h5dclose_f(dset, hdferr)

end function HDF_readHyperslabFloatArray3D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readHyperslabFloatArray4D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns a 4D hyperslab from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param offset offset of the hyperslab
!> @param dims dimensions of the hyperslab
!> @param HDF_head
!
!> @date 04/06/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readHyperslabFloatArray4D(dataname, offset, dims, HDF_head) result(rdata)

IMPLICIT NONE

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_4)

character(fnlen),INTENT(IN)                             :: dataname
integer(HSIZE_T),INTENT(IN)                             :: offset(4)
integer(HSIZE_T),INTENT(IN)                             :: dims(4)
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
real(real_kind), dimension(:,:,:,:), allocatable, TARGET:: rdata

integer(HID_T)                                          :: memspace, space, dset ! Handles
integer(HSIZE_T)                                        :: hdims(2), max_dims(2)
integer                                                 :: hdferr, rnk

call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
call h5dget_space_f(dset, space, hdferr)
call h5sget_simple_extent_dims_f(space, hdims, max_dims, hdferr)

rnk = 4
call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr) 
call h5screate_simple_f(rnk, dims, memspace, hdferr)

allocate(rdata(dims(1),dims(2),dims(3),dims(4)))
call h5dread_f(dset, H5T_NATIVE_REAL, rdata, hdims, hdferr, memspace, space)

call h5sclose_f(space, hdferr)
call h5dclose_f(dset, hdferr)

end function HDF_readHyperslabFloatArray4D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readHyperslabDoubleArray2D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns a 2D hyperslab from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param offset offset of the hyperslab
!> @param dims dimensions of the hyperslab
!> @param HDF_head
!
!> @date 04/06/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readHyperslabDoubleArray2D(dataname, offset, dims, HDF_head) result(rdata)

IMPLICIT NONE

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_8)

character(fnlen),INTENT(IN)                             :: dataname
integer(HSIZE_T),INTENT(IN)                             :: offset(2)
integer(HSIZE_T),INTENT(IN)                             :: dims(2)
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
real(real_kind), dimension(:,:), allocatable, TARGET    :: rdata

integer(HID_T)                                          :: memspace, space, dset ! Handles
integer(HSIZE_T)                                        :: hdims(2), max_dims(2)
integer                                                 :: hdferr, rnk

call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
call h5dget_space_f(dset, space, hdferr)
call h5sget_simple_extent_dims_f(space, hdims, max_dims, hdferr)

rnk = 2
call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr) 
call h5screate_simple_f(rnk, dims, memspace, hdferr)

allocate(rdata(dims(1),dims(2)))
call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdata, hdims, hdferr, memspace, space)

call h5sclose_f(space, hdferr)
call h5dclose_f(dset, hdferr)

end function HDF_readHyperslabDoubleArray2D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readHyperslabDoubleArray3D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns a 3D hyperslab from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param offset offset of the hyperslab
!> @param dims dimensions of the hyperslab
!> @param HDF_head
!
!> @date 04/06/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readHyperslabDoubleArray3D(dataname, offset, dims, HDF_head) result(rdata)

IMPLICIT NONE

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_8)

character(fnlen),INTENT(IN)                             :: dataname
integer(HSIZE_T),INTENT(IN)                             :: offset(3)
integer(HSIZE_T),INTENT(IN)                             :: dims(3)
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
real(real_kind), dimension(:,:,:), allocatable, TARGET  :: rdata

integer(HID_T)                                          :: memspace, space, dset ! Handles
integer(HSIZE_T)                                        :: hdims(2), max_dims(2)
integer                                                 :: hdferr, rnk

call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
call h5dget_space_f(dset, space, hdferr)
call h5sget_simple_extent_dims_f(space, hdims, max_dims, hdferr)

rnk = 3
call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr) 
call h5screate_simple_f(rnk, dims, memspace, hdferr)

allocate(rdata(dims(1),dims(2),dims(3)))
call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdata, hdims, hdferr, memspace, space)

call h5sclose_f(space, hdferr)
call h5dclose_f(dset, hdferr)

end function HDF_readHyperslabDoubleArray3D

!--------------------------------------------------------------------------
!
! FUNCTION:HDF_readHyperslabDoubleArray4D
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reads and returns a 4D hyperslab from the current file or group ID 
!
!> @note Note that this routine uses fortran-2003 options
!
!> @param dataname dataset name (string)
!> @param offset offset of the hyperslab
!> @param dims dimensions of the hyperslab
!> @param HDF_head
!
!> @date 04/06/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readHyperslabDoubleArray4D(dataname, offset, dims, HDF_head) result(rdata)

IMPLICIT NONE

integer,parameter                                       :: real_kind = SELECTED_REAL_KIND(Fortran_REAL_8)

character(fnlen),INTENT(IN)                             :: dataname
integer(HSIZE_T),INTENT(IN)                             :: offset(4)
integer(HSIZE_T),INTENT(IN)                             :: dims(4)
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
real(real_kind), dimension(:,:,:,:), allocatable, TARGET:: rdata

integer(HID_T)                                          :: memspace, space, dset ! Handles
integer(HSIZE_T)                                        :: hdims(2), max_dims(2)
integer                                                 :: hdferr, rnk

call h5dopen_f(HDF_head%objectID, trim(dataname), dset, hdferr)
call h5dget_space_f(dset, space, hdferr)
call h5sget_simple_extent_dims_f(space, hdims, max_dims, hdferr)

rnk = 4
call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr) 
call h5screate_simple_f(rnk, dims, memspace, hdferr)

allocate(rdata(dims(1),dims(2),dims(3),dims(4)))
call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdata, hdims, hdferr, memspace, space)

call h5sclose_f(space, hdferr)
call h5dclose_f(dset, hdferr)

end function HDF_readHyperslabDoubleArray4D




end module HDFsupport
