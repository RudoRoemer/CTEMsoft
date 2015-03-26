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
!> @details: the EM HDF format has three parts: 
!>
!> EMheader, created with HDF_writeEMheader
!> 
!> EMnamelistfiles, which contains all the variables from all the namelists
!> and is created as follows in the calling program:
!>
!>     create and open the EMnamelistfiles group
!>     call HDF_createGroup('EMnamelistfiles', HDF_head, HDF_tail)
!>
!>     then for each relevant namelist: (from the NameListHDFwriters module)
!>     HDFwriteKosselNameList(HDF_head, knl)
!>     ...
!>
!>     and then close the group
!>     call HDF_pop(HDF_head)
!> 
!> and finally the EMdata section
!>
!
!> @date  03/17/15  MDG 1.0 original
!--------------------------------------------------------------------------
module HDFsupport

use local
use typedefs
use HDF5
use h5lt

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
!>
!> @param HDF_head pointer to top of push-pop stack
!> @param HDF_tail pointer to bottom of push-pop stack
!> @param dstr date string
!> @param tstrb time start string
!> @param tstre time end string
!> @param prn program name
!>
!> @date 03/20/15  MDG 1.0 original
!--------------------------------------------------------------------------
subroutine HDF_writeEMheader(HDF_head, HDF_tail, dstr, tstrb, tstre, prn)

use local
use io

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head
type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_tail
character(11),INTENT(IN)                              :: dstr
character(15),INTENT(IN)                              :: tstrb
character(15),INTENT(IN)                              :: tstre
character(fnlen),INTENT(IN)                           :: prn

integer                                               :: error, istat  ! error flag
integer                                               :: i,ic,nlen 
character(100)                                        :: c


! create and open the EMheader group
istat = HDF_createGroup('EMheader', HDF_head, HDF_tail)

! version number /EMheader/Version 'character'
call h5ltmake_dataset_string_f(HDF_head%objectID, 'Version', scversion, error)
if (error.ne.0) call HDF_handleError(error,'HDF_writeEMheader: unable to write Version',.TRUE.)

! execution data /EMheader/Date 'character'
call h5ltmake_dataset_string_f(HDF_head%objectID, 'Date', dstr, error)
if (error.ne.0) call HDF_handleError(error,'HDF_writeEMheader: unable to write Date',.TRUE.)

! start time /EMheader/StartTime 'character'
call h5ltmake_dataset_string_f(HDF_head%objectID, 'StartTime', tstrb, error)
if (error.ne.0) call HDF_handleError(error,'HDF_writeEMheader: unable to write StartTime',.TRUE.)

! stop time /EMheader/StopTime 'character'
call h5ltmake_dataset_string_f(HDF_head%objectID, 'StopTime', tstre, error)
if (error.ne.0) call HDF_handleError(error,'HDF_writeEMheader: unable to write StopTime',.TRUE.)

! program name /EMheader/ProgramName 'character'
call h5ltmake_dataset_string_f(HDF_head%objectID, 'ProgramName', trim(prn), error)
if (error.ne.0) call HDF_handleError(error,'HDF_writeEMheader: unable to write ProgramName',.TRUE.)

! user name /EMheader/UserName 'character'
call h5ltmake_dataset_string_f(HDF_head%objectID, 'UserName', trim(username), error)
if (error.ne.0) call HDF_handleError(error,'HDF_writeEMheader: unable to write Username',.TRUE.)

! user location /EMheader/UserLocation 'character'
call h5ltmake_dataset_string_f(HDF_head%objectID, 'UserLocation', trim(userlocn), error)
if (error.ne.0) call HDF_handleError(error,'HDF_writeEMheader: unable to write UserLocation',.TRUE.)

! user email /EMheader/UserEmail 'character'
call h5ltmake_dataset_string_f(HDF_head%objectID, 'UserEmail', trim(useremail), error)
if (error.ne.0) call HDF_handleError(error,'HDF_writeEMheader: unable to write UserEmail',.TRUE.)

! hostname /EMheader/HostName 'character'
call hostnm(c)
! lowercase it
nlen = len(c) 
do i=1,nlen 
   ic = ichar(c(i:i)) 
   if (ic >= 65 .and. ic < 90) c(i:i) = char(ic+32) 
end do 
call h5ltmake_dataset_string_f(HDF_head%objectID, 'HostName', trim(c), error)
if (error.ne.0) call HDF_handleError(error,'HDF_writeEMheader: unable to write HostName',.TRUE.)

! and close this group
call HDF_pop(HDF_head)

end subroutine HDF_writeEMheader



! create and open the EMheader group
! call HDF_createGroup('EMnamelistfiles', HDF_head, HDF_tail)
!
! then for each relevant namelist:
! HDFwriteKosselNameList(HDF_head, knl)
! ...
!
! and then close the group
! call HDF_pop(HDF_head)

! loop over the files and write them to the HDF file as individual groups 
!do i=1,numf(1)
!  call HDF_createGroup(trim(nml_list(i)), HDF_head, HDF_tail)
! get the file information
!  call stat(trim(nml_list(i)), values)
!  numb(1) = values(8)
!  allocate(filebuf(numb(1)))
! and read the file
!  open(unit=dataunit,file=trim(nml_list(i)),access='DIRECT',RECL=numb(1))
!  read (dataunit,REC=1) filebuf
!  close(dataunit)
! then create a dataset with this data
!  data_name = trim(nml_list(i))
!  call h5screate_simple_f(rnk, numb, dataspace, error)
!  call h5dcreate_f(HDF_head%objectID, data_name, H5T_NATIVE_CHARACTER, dataspace, dset_id, error)
!  call h5dwrite_f(dset_id, H5T_NATIVE_CHARACTER, filebuf, numb, error)
!  call h5dclose_f(dset_id,error)
!  if (error.ne.0) write (*,*) ' error writing to file '
! and get rid of the string array
!  deallocate(filebuf)
!end do




!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! from here on, we have basic HDF support routines
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
!> @param HDF_tail bottom of the current stack
!> @param oT object type character
!> @param oID object identifier
!> @param verbose (optional) 
!
!> @date 03/17/15  MDG 1.0 original
!--------------------------------------------------------------------------
subroutine HDF_push(HDF_head, HDF_tail, oT, oID, verbose)

use local
use io

IMPLICIT NONE

type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_head
type(HDFobjectStackType),INTENT(INOUT),pointer        :: HDF_tail
character(LEN=1),INTENT(IN)                           :: oT
integer(HID_T),INTENT(IN)                             :: oID 
logical,INTENT(IN),OPTIONAL                           :: verbose

type(HDFobjectStackType),pointer                      :: node
integer(kind=irg)                                     :: istat

! if the stack doesn't exist yet, create it.
if (.not.associated(HDF_tail)) then 
   allocate(HDF_tail,stat=istat)                        ! allocate new value
   if (istat.ne.0) call HDF_handleError(istat,'HDF_push: unable to allocate HDF_stack_tail pointer',.TRUE.)
   nullify(HDF_tail%next)                               ! nullify next in tail value
   HDF_head => HDF_tail                                 ! head points to new value
   if (PRESENT(verbose)) then 
     if (verbose) call Message('  -> creating HDF_stack_tail linker list', frm = "(A)")
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
    call WriteValue('',tmp%objectType//'  '//tmp%objectName, frm = "(A$)") 
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
!
! SUBROUTINE:HDF_createFile
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Create a new HDF file (this also opens the file)
!
!> @param HDFname filename string
!> @param HDF_head
!> @param HDF_tail
!
!> @date 03/17/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_createFile(HDFname, HDF_head, HDF_tail) result(success)
 
use local

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: HDFname
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_tail
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
  call HDF_push(HDF_head, HDF_tail, 'f',file_id)
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
!> @param HDF_tail
!
!> @date 03/17/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_openFile(HDFname, HDF_head, HDF_tail) result(success)

use local

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: HDFname
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_tail
integer(kind=irg)                                       :: success

integer(HID_T)                                          :: file_id ! file identifier
integer                                                 :: error  ! error flag

success = 0

call H5Fopen_f(HDFname, H5F_ACC_RDWR_F, file_id, error);
if (error.ne.0) then 
  call HDF_handleError(error,'HDF_openFile')
  success = -1
else
  call HDF_push(HDF_head, HDF_tail, 'f',file_id)
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
!> @param HDF_tail
!
!> @date 03/17/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_createGroup(groupname, HDF_head, HDF_tail) result(success)

IMPLICIT NONE

character(*),INTENT(IN)                                 :: groupname
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_tail
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
  call HDF_push(HDF_head, HDF_tail, 'g', group_id)
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
!> @param HDF_tail
!
!> @date 03/17/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_openGroup(groupname, HDF_head, HDF_tail) result(success)

IMPLICIT NONE

character(*),INTENT(IN)                                 :: groupname
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_tail
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
  call HDF_push(HDF_head, HDF_tail, 'g', group_id)
end if

end function HDF_openGroup

!--------------------------------------------------------------------------
!
! SUBROUTINE:HDF_createDataset
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief create a new dataset (this also opens the dataset)
!
!> @param dataname string
!> @param HDF_head
!> @param HDF_tail
!
!> @date 03/17/15  MDG 1.0 original
!--------------------------------------------------------------------------
!function HDF_createDataset(dataname, HDF_head, HDF_tail) result(success)
!
!IMPLICIT NONE
!
!character(fnlen),INTENT(IN)                             :: dataname
!type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
!type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_tail
!integer(kind=irg)                                       :: success
!
!integer(HID_T)                                          :: data_id!  identifier
!integer                                                 :: error  ! error flag
!
!success = 0
!
!call H5dcreate_f(HDF_head%objectID, dataname, data_id, error)
!if (error.ne.0) then
!  call HDF_handleError(error,'HDF_createDataset')
!  success = -1
!else
!! and put the data_id onto the HDF_stack
!  call HDF_push(HDF_head, HDF_tail, 'd', data_id)
!end if
!
!end function HDF_createDataset

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
!> @param HDF_tail
!
!> @date 03/17/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_openDataset(dataname, HDF_head, HDF_tail) result(success)

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_tail
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
  call HDF_push(HDF_head, HDF_tail, 'd', data_id)
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
!> @param HDF_tail
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_writeDatasetTextFile(dataname, filename, HDF_head, HDF_tail) result(success)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
character(fnlen),INTENT(IN)                             :: filename
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_tail
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
CALL H5Tcopy_f(H5T_STRING, filetype, hdferr)
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetStringArray',.TRUE.)
  success = -1
end if
!
! Create dataspace.
!
rnk = 1
CALL h5screate_simple_f(rnk, dims, space, hdferr)
!
! Create the dataset and write the variable-length string data to it.
!
CALL h5dcreate_f(HDF_head%objectID, trim(dataname), filetype, space, dset, hdferr)
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetStringArray:hd5create_f',.TRUE.)
  success = -1
end if
f_ptr = C_LOC(wdata(1))
CALL h5dwrite_f(dset, filetype, f_ptr, hdferr )
if (hdferr.ne.0) then
  call HDF_handleError(hdferr,'HDF_writeDatasetStringArray:hd5write_f',.TRUE.)
  success = -1
end if
!
! Close and release resources.
!
CALL h5dclose_f(dset , hdferr)
CALL h5sclose_f(space, hdferr)
deallocate(wdata)

! that's it

end function HDF_writeDatasetTextFile

!--------------------------------------------------------------------------
!
! F?UNCTION:HDF_readfromTextfile
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

dt = 30
! read the file first to determine the number of lines
open(unit=dt,file=trim(filename),form='formatted',status='old')
nlines = 0
do
  read (dt,"(A)",end=10) line(1)
  nlines = nlines + 1
end do
10 close(unit=dt,status='keep')

! then re-read the file and store all the lines in the wdata array
allocate(stringarray(1:nlines))
open(unit=dt,file=trim(filename),form='formatted',status='old')
do i=1,nlines
! initialize the line to null characters before each read
  do j=1,fnlen
    line(1)(j:j) = char(0)
  end do
! read the line
  read (dt,"(A)") line(1)
! find the string length and put the next character equal to C_NULL_CHAR
  j = len(trim(line(1)))+1
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
!> @param nlines number of lines read from file
!> @param HDF_head
!> @param HDF_tail
!
!> @date 03/26/15  MDG 1.0 original
!--------------------------------------------------------------------------
function HDF_readDatasetStringArray(dataname, nlines, HDF_head, HDF_tail) result(stringarray)

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                             :: dataname
integer(kind=irg),INTENT(OUT)                           :: nlines
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_head
type(HDFobjectStackType),INTENT(INOUT),pointer          :: HDF_tail
character(len=fnlen, KIND=c_char),allocatable, TARGET   :: stringarray(:) 

integer(HID_T)                                          :: filetype, space ! Handles
integer                                                 :: hdferr, i, length
integer(HSIZE_T), DIMENSION(1:1)                        :: dims
integer(HSIZE_T), DIMENSION(1:2)                        :: maxdims

character(len = fnlen, kind=c_char),  POINTER           :: pfstr ! A pointer to a Fortran string
TYPE(C_PTR), DIMENSION(:), ALLOCATABLE, TARGET          :: rdata ! Read buffer
TYPE(C_PTR)                                             :: f_ptr


! we'll assume that the file is opened by the calling program

!
! Open dataset.
!
hdferr = HDF_openDataset(dataname, HDF_head, HDF_tail)

!
! Get the datatype.
!
CALL H5Dget_type_f(HDF_head%objectID, filetype, hdferr)

!
! Get dataspace and allocate memory for read buffer.
!
CALL H5Dget_space_f(HDF_head%objectID, space, hdferr)
CALL H5Sget_simple_extent_dims_f(space, dims, maxdims, hdferr)

ALLOCATE(rdata(1:dims(1)), stringarray(1:dims(1)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1))
CALL h5dread_f(HDF_head%objectID, H5T_STRING, f_ptr, hdferr)
!
! convert the data to a string array
!
DO i = 1, dims(1)
  CALL C_F_POINTER(rdata(i), pfstr)
  length = 0
  DO
     IF(pfstr(length+1:length+1).EQ.C_NULL_CHAR.OR.length.GE.fnlen) EXIT
     length = length + 1
  ENDDO
  stringarray(i) = pfstr(1:length)
!     WRITE(*,'(A,"(",I0,"): ",I4," characters : ",A)') trim(DATASET), i, length, data(1:length)
END DO

nlines = dims(1)

DEALLOCATE(rdata)

CALL h5sclose_f(space, hdferr)
CALL H5Tclose_f(filetype, hdferr)
! close the dataset
call HDF_pop(HDF_head)


end function HDF_readDatasetStringArray



end module HDFsupport
