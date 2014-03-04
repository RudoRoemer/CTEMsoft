! ###################################################################
! Copyright (c) 2013, Marc De Graef/Carnegie Mellon University
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
! CTEMsoft2013:HDF5Utilities.f90
!--------------------------------------------------------------------------
!
! MODULE: HDF5Utilities
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief basic library for HDF5 file handling (based on M. Jackson's H5Utilities.cpp library).
!
!
!> @date 4/16/12   MDG 1.0 original
!> @date 2/26/14   MDG 2.0 inclusion in CTEMsoft2013 project and rewrite
!--------------------------------------------------------------------------

module HDF5Utilities

use HDF5
use H5FORTRAN_TYPES
use ISO_C_BINDING
use m_strings
use HDFConstants
use HDFVars

logical   				:: HDFverbose=.FALSE.   ! this can be set to .TRUE. in the main program to cause the routines below to be more HDFverbose

integer,parameter        		:: HERR_T = SELECTED_INT_KIND(9) 


type objectStackType               ! this is a push-pop stack to keep track of the open HDF5 objects
  character(LEN=1)			:: objectType
  integer(HID_T)			:: objectID
  type(objectStackType),pointer 	:: next
end type objectStackType

type(objectStackType),pointer 		:: HDF_stack_head, HDF_stack_tail


public :: HDF_push
public :: HDF_pop




contains



! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!>  Private routine to extend the dde_ar array by one unit
!! @param num The current number of entries in the array
!! @param success Integer parameter returned by the function
function extend_dde_ar(num) result(success)

IMPLICIT NONE

integer,INTENT(IN) :: num           ! number of entries in array
integer :: success   ! returned variable
integer :: ier  ! error flag for allocate command
type(dde_entry_type),allocatable :: dde_ar_local(:)   ! local copy of dde_ar

! first copy the current array
allocate(dde_ar_local(0:num-1),stat=ier)
if (ier.ne.0) then
  success=-1
  return
end if

dde_ar_local = dde_ar

! destroy the old array
deallocate(dde_ar)
! and allocate dde_ar with one more entry
allocate(dde_ar(0:num),stat=ier)
if (ier.ne.0) then
  success=-1
  return
end if

! copy all the old entries back into the array
dde_ar(0:num-1) = dde_ar_local(0:num-1)

! and return 0 to the calling program
deallocate(dde_ar_local)
success = 0

end function extend_dde_ar

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> @brief private routine to extend the dre_ar array by one unit
!! @param num The current number of entries in the array
!! @param success Integer parameter returned by the function
function extend_dre_ar(num) result(success)

IMPLICIT NONE

integer,INTENT(IN) :: num           ! current number of entries in array
integer :: success   ! returned variable
integer :: ier  ! error flag for allocate command
type(dre_entry_type),allocatable :: dre_ar_local(:)   ! local copy of dre_ar

! first copy the current array
allocate(dre_ar_local(0:num-1),stat=ier)
if (ier.ne.0) then
  success=-1
  return
end if

dre_ar_local = dre_ar

! destroy the old array
deallocate(dre_ar)
! and allocate dde_ar with one more entry
allocate(dre_ar(0:num),stat=ier)
if (ier.ne.0) then
  success=-1
  return
end if

! copy all the old entries back into the array
dre_ar(0:num-1) = dre_ar_local(0:num-1)

! and return 0 to the calling program
deallocate(dre_ar_local)
success = 0

end function extend_dre_ar



! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> @brief private routine to extend the dumd_ar array by one unit
!! @param num The current number of entries in the array
!! @param success Integer parameter returned by the function
function extend_dumd_ar(num) result(success)

IMPLICIT NONE

integer,INTENT(IN) :: num           ! current number of entries in array
integer :: success   ! returned variable
integer :: ier  ! error flag for allocate command
type(dumd_entry_type),allocatable :: dumd_ar_local(:)   ! local copy of dmd_ar

! first copy the current array
allocate(dumd_ar_local(0:num-1),stat=ier)
if (ier.ne.0) then
  success=-1
  return
end if

dumd_ar_local = dumd_ar

! destroy the old array
deallocate(dumd_ar)
! and allocate dde_ar with one more entry
allocate(dumd_ar(0:num),stat=ier)
if (ier.ne.0) then
  success=-1
  return
end if

! copy all the old entries back into the array
dumd_ar(0:num-1) = dumd_ar_local(0:num-1)

! and return 0 to the calling program
deallocate(dumd_ar_local)
success = 0

end function extend_dumd_ar





! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> Push an object onto the HDF_stack
!! @param oT  object type character
!! @param oID object ID
subroutine HDF_push(oT, oID)

IMPLICIT NONE

character(LEN=1),INTENT(IN)	:: oT
integer(HID_T),INTENT(IN)	:: oID 
type(objectStackType),pointer	:: node
integer 			:: istat

! if the stack doesn't exist yet, create it.
if (.not.associated(HDF_stack_tail)) then 
   allocate(HDF_stack_tail,stat=istat)           ! allocate new value
   if (istat.ne.0) call H5U_handleError(istat,'HDF_push: unable to allocate HDF_stack_tail pointer',.TRUE.)
   nullify(HDF_stack_tail%next)                  ! nullify next in tail value
   HDF_stack_head => HDF_stack_tail              ! head points to new value
end if

allocate(node,stat=istat)   ! allocate a new node
if (istat.ne.0) call H5U_handleError(istat,'HDF_push: unable to allocate node pointer',.TRUE.)
! set the values
node % objectType = oT
node % objectID = oID
! insert the node into the stack
node % next => HDF_stack_head
! and re-point the head
HDF_stack_head => node

end subroutine HDF_push


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> Pop an object from the HDF_stack and close the corresponding HDF5 object
!! @param closeall  optional; if present, close all the objects in the stack
subroutine HDF_pop(closeall)

IMPLICIT NONE

logical,optional 	:: closeall
integer            	:: error

if (present(closeall)) then  
! this would be called if an error arises that forces a complete shutdown of the program, or at the end of a regular program
  do while (associated(HDF_stack_head % next)) 
! close the current object 
    error = H5U_close_level(HDF_stack_head % objectType, HDF_stack_head % objectID)
! and re-point the stack head
    HDF_stack_head =>HDF_stack_head % next  
  end do
else
! close the current object 
  error = H5U_close_level(HDF_stack_head % objectType, HDF_stack_head % objectID)
! and re-point the stack head
  HDF_stack_head =>HDF_stack_head % next  
end if

contains

  function H5U_close_level(oT, oID) result(error)
  
  IMPLICIT NONE

  character(LEN=1),INTENT(IN) 	:: oT
  integer(HID_T),INTENT(IN) 	:: oID 
  integer            		:: error

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

end function 

end subroutine HDF_pop


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> Get the information from the current object at the top of the HDF_stack
!! @param oT  object type character
!! @param oID object ID
subroutine HDF_get(oT, oID)

IMPLICIT NONE

character(LEN=1),INTENT(OUT)	:: oT
integer(HID_T),INTENT(OUT)  	:: oID 

! simply return the current oT and oID from the top of the stack
if (associated(HDF_stack_head)) then
  oT = HDF_stack_head % objectType
  oID = HDF_stack_head % objectID
else
  call H5U_handleError(-1,'HDF_get: attempting to access unassociated node pointer',.TRUE.)
end if

end subroutine HDF_get

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> Get the information from the current object at the top of the HDF_stack
!! @param oT  object type character
!! @param oID object ID
integer(HID_T) function HDF_getParent()

IMPLICIT NONE

! simply return the current oT and oID from the top of the stack
if (associated(HDF_stack_head)) then
  HDF_getParent = HDF_stack_head % objectID
else
  call H5U_handleError(-1,'HDF_getParent: attempting to access unassociated node pointer',.TRUE.)
end if

end function HDF_getParent



! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> tests whether this is a group or not
logical function H5U_isGroup(gname)

integer		:: printflag, error
type(string) 	:: gname
integer(HID_T) :: grp_id

! turn error messaging off
printflag = 0
call h5eset_auto_f(printflag, error)

call h5gopen_f(HDF_getParent(),char(gname),grp_id,error)

H5U_isGroup = .FALSE.
if (error.eq.0) then 
  call h5gclose_f(grp_id,error)
  H5U_isGroup = .TRUE.
end if

! turn error messaging back on
printflag = 1
call h5eset_auto_f(printflag, error)

end function H5U_isGroup



! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> Read all the Data Dimensions information
function H5U_readDataDimensions() result(success) ! .TRUE. = readonly; .FALSE. = readwrite

IMPLICIT NONE

integer(HID_T)  		:: data_dim_id, dataset_id, model_g_id
integer            			:: DIMID, storage_type, max_corder, error, i
character(LEN=10)     	:: iloc
type(string)			:: success

success= ''

!   Open the Data Model group.
model_g_id = H5U_openGroup(HDF_DATAMODEL,HDF_getParent())

! open the Data Dimensions group
data_dim_id = H5U_openGroup(HDF_DataDimensions,HDF_getParent())  ! Data Model should be the parent object

! find out how many data sets there are in this group
call h5gget_info_f(data_dim_id, storage_type, DIMID, max_corder, error)
if (error.ne.0) then
  call H5U_handleError(error,'H5U_readDataDimensions: error getting number of Data Dimensions',.TRUE.)		! print an error message
  call HDF_pop(.TRUE.)										! close all the open objects
  call h5close_f(error)										! close the Fortran HDF5 interface
  success = 'Program aborted due to : unable to determine # of data dimensions'	! and abort with friendly message
  return
end if

! allocate the data dimensions array
allocate(dde_ar(0:DIMID-1),stat=error)
if (error.ne.0) then
  call H5U_handleError(error,'H5U_readDataDimensions: error allocating dde_ar',.TRUE.)	! print an error message
  call HDF_pop(.TRUE.)										! close all the open objects
  call h5close_f(error)										! close the Fortran HDF5 interface
  success = 'Program aborted due to : unable to allocate dde_ar array'	! and abort with friendly message
  return
end if

do i=0,DIMID-1 
  write (iloc,"(I)") i+1

! open the dataset
  dataset_id = H5U_openDataSet(trim(iloc),data_dim_id)
  
!  get the Name attribute for this dimension
  dde_ar(i) % nm = H5U_getGroupAttributeString(HDF_NAME_TAG,dataset_id)

!  get the Alt Name attribute
  dde_ar(i) % dnm = H5U_getGroupAttributeString(HDF_ALT_NAME_TAG,dataset_id)

!  get the Count attribute
  dde_ar(i) % cnt = H5U_getGroupAttributeInteger(HDF_COUNT_TAG,dataset_id)

!  get the Start Value attribute
  dde_ar(i) % sv = H5U_getGroupAttributeInteger(HDF_START_VALUE_TAG,dataset_id)

!  get the Increment attribute
  dde_ar(i) % inc = H5U_getGroupAttributeInteger(HDF_INCREMENT_TAG,dataset_id)

!  get the Count attribute
  dde_ar(i) % ev = H5U_getGroupAttributeInteger(HDF_END_VALUE_TAG,dataset_id)

!  get the Count attribute
  dde_ar(i) % uni = H5U_getGroupAttributeInteger(HDF_UNIFORM_TAG,dataset_id)

  call HDF_pop  ! close the data set for this dimension
end do

! close the Data Dimensions group
call HDF_pop  

! close the Data Model group
call HDF_pop  

end function H5U_readDataDimensions


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> Scan the Data Record structure
function H5U_readDataRecords() result(success)

IMPLICIT NONE

integer			:: drlevel, res, drec_count, storage_type, max_corder, error
type(string)		:: UIDstring, success
integer(HID_T)		:: group_id, model_g_id


success = ''

!   Open the Data Model group.
model_g_id = H5U_openGroup(HDF_DATAMODEL,HDF_getParent())

! open the Data Records group
group_id = H5U_openGroup(HDF_DataRecords,HDF_getParent())  ! Data Model should be the parent object

! upon first entry of this routine, we must check the Data Records group to see
! how many members it has, using the h5gget_info_f function.  It must have at 
! least one, which may be a dataset or a group.  If it is a data set, then it will
! not have any members itself
call h5gget_info_f(HDF_getParent(),storage_type,drec_count,max_corder,error)
if (error.ne.0) then
  call H5U_handleError(error,'H5U_readDataRecords: h5gget_info_f returned non-zero status',.TRUE.)	! print an error message
  call HDF_pop(.TRUE.)										! close all the open objects
  call h5close_f(error)										! close the Fortran HDF5 interface
  success = 'Program aborted due to : problem in h5gget_info_f'	! and abort with friendly message
  return
end if

! if drec_count is zero, we issue an error message
if (drec_count.eq.0) then
  call H5U_handleError(error,'H5U_readDataRecords: no Data Records found in Data Records group',.TRUE.)	! print an error message
  call HDF_pop(.TRUE.)
  call h5close_f(error)										! close the Fortran HDF5 interface
  success = 'Program aborted due to : Data Records group empty'	! and abort with friendly message
  return
end if

! call the recursive routine 
drlevel = 0
UIDstring = "/"//HDF_DataRecordsPath
error = read_data_records(drlevel,UIDstring,drec_count)
if (error.ne.0) then
  call H5U_handleError(error,'H5U_readDataRecords: recursive function returned non-zero status',.TRUE.)	! print an error message
  call HDF_pop(.TRUE.)										! close all the open objects
  call h5close_f(error)										! close the Fortran HDF5 interface
  success = 'Program aborted due to : problem in recursive Data Records function'	! and abort with friendly message
  return
end if

! and close the Data Records group
call HDF_pop

!  close the Data Model group
call HDF_pop

contains ! recursive private routine

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> Recursively scan the Data Record structure
!! @param drlevel recursive level
!! @param UIDstring  string used to define the data set paths
!! @param drec_count counter of thenumber of objects at this level
	integer recursive function read_data_records(drlevel,UIDstring,drec_count) 
	!
	! this is a recursive routine that parses the Data Records group and fills in the 
	! dre_ar array with the appropriate information
	!
	! This routine is called by H5U_readDataRecords and by itself.
	!
	IMPLICIT NONE
	
	integer,INTENT(INOUT) 		:: drlevel
	type(string),INTENT(INOUT) 	:: UIDstring
	integer,INTENT(IN)			:: drec_count
	type(string)  				:: objname
	type(string),allocatable   		:: UID(:)
	integer 					:: i, drid, sz, GUID, new_reccount ! GUID = Global Unique IDentifier
	character(LEN=80)			:: iloc
	integer(HID_T)				:: local_id
	
	read_data_records = -1
	
	! string array needed to define the relative paths
	allocate(UID(0:drec_count-1))
	do i=0,drec_count-1
	      write (iloc,"(I)") i
	      UID(i) = UIDstring//"/"//trim(iloc)
	end do
	
	! loop over all the objects at this level
	do drid=0,drec_count-1
	! initialize the dre_ar array or extend it
	   if ((drlevel.eq.0).and.(drid.eq.0)) then   ! create the array
	      sz = 1
	      allocate(dre_ar(0:sz-1),stat=error)
	   else  ! extend it by one unit
	      sz = size(dre_ar)
	      error =  extend_dre_ar(sz)
	   end if
	
	! determine whether this is a group or a dataset and open it accordingly
	   if (H5U_isGroup(UID(drid))) then  ! it is a group
	      local_id = H5U_openGroup(char(UID(drid)),HDF_getParent())
	      GUID = H5U_getGroupAttributeInteger(HDF_GUID_TAG,local_id)
	      dre_ar(GUID) % GUID = GUID
	      dre_ar(GUID) % hierarchy = drlevel
	      dre_ar(GUID) % locatorstring = UID(drid)
	      dre_ar(GUID) % entrytype = 'G'
	      dre_ar(GUID) % LUID = H5U_getGroupAttributeInteger(HDF_LUID_TAG,local_id)
	      dre_ar(GUID) % nm = H5U_getGroupAttributeString(HDF_NAME_TAG,local_id)
	      dre_ar(GUID) % dnm = H5U_getGroupAttributeString(HDF_ALT_NAME_TAG,local_id)
	
	! this is a group so it must have elements; check how many, and recursively call this routine
	      call h5gget_info_f(HDF_getParent(),storage_type,new_reccount,max_corder,error)
	      drlevel = drlevel + 1 		! increment the recursion level
	      res = read_data_records(drlevel,UID(drid),new_reccount)
	      read_data_records = res   
	      drlevel = drlevel - 1 		! decrement the recusion level
	      
	      call HDF_pop    ! close the present group
	  else 		! it is a data set
	      local_id =H5U_openDataSet(char(UID(drid)),HDF_getParent())
	      GUID = H5U_getGroupAttributeInteger(HDF_GUID_TAG,local_id)
	      dre_ar(GUID) % GUID = GUID
	      dre_ar(GUID) % hierarchy = drlevel
	      dre_ar(GUID) % locatorstring = UID(drid)
	      dre_ar(GUID) % entrytype = 'D'
	      dre_ar(GUID) % LUID = H5U_getGroupAttributeInteger(HDF_LUID_TAG,local_id)
	      dre_ar(GUID) % nm = H5U_getGroupAttributeString(HDF_NAME_TAG,local_id)
	      dre_ar(GUID) % dnm = H5U_getGroupAttributeString(HDF_ALT_NAME_TAG,local_id)
	! data sets don't have children objects, so need to recurse here   
	   
	      call HDF_pop  ! close the present data set
	  end if 
	end do
	
	! if we get here, everything is ok
	read_data_records = 0
	
	deallocate(UID)
	
	end function read_data_records


end function H5U_readDataRecords  ! function contains read_data_records



! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> Scan the Data Root data set and return the string
function H5U_readDataRoot() result(success)

IMPLICIT NONE

integer(HID_T)				:: local_id, data_type, model_g_id
character(LEN=1024)		:: content
integer(HSIZE_T),dimension(1) :: cdata_dims  ! data identifier
integer					:: error
type(string)				:: success

success = ''

!   Open the Data Model group.
model_g_id = H5U_openGroup(HDF_DATAMODEL,HDF_getParent())

! open the Data Root data set
local_id =H5U_openDataSet(HDF_DataRoot,HDF_getParent())

! get the attribute type
if (HDFbeverbose) call H5U_printMessage('calling h5dget_type_f')
call h5dget_type_f(local_id, data_type, error)
if (error.ne.0) call H5U_handleError(error,'h5dget_type_f in H5U_readDataRoot')

!   get the attribute value
content = ''
if (HDFbeverbose) call H5U_printMessage('calling h5dread_f')
call h5dread_f(local_id, data_type, content, cdata_dims, error)
if (error.ne.0) call H5U_handleError(error,'h5dread_f')

! and return the value
HDFdataRoot = trim(content)

! and close the data set
call HDF_pop

! and the Data Model group
call HDF_pop

end function H5U_readDataRoot



! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> Scan the Required Meta Data attributes 
function H5U_readRequiredMetaData() result(success)

IMPLICIT NONE

integer(HID_T)  		:: rmd_id, metadata_g_id
type(string)			:: success

success = ''

! open the Meta Data group
metadata_g_id = H5U_openGroup(HDF_MetaData,HDF_getParent())

! first we'll read the Required MetaData, so we need to open the data set
rmd_id =H5U_openDataSet(HDF_RequiredMetaData,HDF_getParent())

! these attributes are truly required and the program will abort if they one is missing
dmd_ar % rmdCreator = H5U_getDataAttributeString(HDF_CREATOR_TAG,rmd_id)
dmd_ar % rmdDate = H5U_getDataAttributeString(HDF_DATE_TAG,rmd_id)
dmd_ar % rmdName = H5U_getDataAttributeString(HDF_Name_TAG,rmd_id)
dmd_ar % rmdDescription = H5U_getDataAttributeString(HDF_DESCRIPTION_TAG,rmd_id)

! and these are optional so we should initalize them as empty strings
dmd_ar % rmdPedigree = ''
dmd_ar % rmdOriginalSourceFile = ''
dmd_ar % rmdDistributionRights = ''
dmd_ar % rmdReleaseLimitation = ''
dmd_ar % rmdReleaseNumber = ''

! and now fill them if they are present
dmd_ar % rmdPedigree = H5U_getDataAttributeString(HDF_PEDIGREE_TAG,rmd_id,.TRUE.)  ! Original or Derived
dmd_ar % rmdOriginalSourceFile = H5U_getDataAttributeString(HDF_DERIVED_SRC_TAG,rmd_id,.TRUE.)
dmd_ar % rmdDistributionRights = H5U_getDataAttributeString(HDF_RIGHTS_TAG,rmd_id,.TRUE.)
dmd_ar % rmdReleaseLimitation = H5U_getDataAttributeString(HDF_RELEASE_LIMITATION_TAG,rmd_id,.TRUE.)
dmd_ar % rmdReleaseNumber = H5U_getDataAttributeString(HDF_RELEASE_NUMBER_TAG,rmd_id,.TRUE.)

! close the data set
call HDF_pop 

! and close the Meta Data group 
call HDF_pop

end function H5U_readRequiredMetaData


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> Scan the Required Meta Data attributes 
function H5U_readUserDefinedMetaData() result(success)

IMPLICIT NONE

integer(HID_T)  		:: rmd_id, metadata_g_id
type(string)			:: success
integer				:: attr_num, error

success = ''

! open the Meta Data group
metadata_g_id = H5U_openGroup(HDF_MetaData,HDF_getParent())

! first we'll read the Required MetaData, so we need to open the data set
rmd_id =H5U_openDataSet(HDF_UserMetaData,HDF_getParent())

! then we'll check if there are any attributes in this data set
call h5aget_num_attrs_f(rmd_id, attr_num, error)    ! This routine will need to be replaced in the future !
if (error.ne.0) call H5U_handleError(error,'h5aget_num_attrs_f')

write (*,*) 'There are ',attr_num,' user defined attributes '

! so far there are none in any of my files, so we'll do this part later, after
! we've managed to write a file that contains them...


! close the data set
call HDF_pop 

! close the Meta Data group
call HDF_pop 

end function H5U_readUserDefinedMetaData




! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> Open the HDF file and verify the type and version number; leave Data Model group open upon return
!! @param fileName name of the input file
!! @param rdonly logical read-ony (if .TRUE.) or read-write
integer function H5U_openVerifyFile(fileName,rdonly) ! .TRUE. = readonly; .FALSE. = readwrite

IMPLICIT NONE

character(LEN=*),INTENT(IN)  	:: fileName
logical,INTENT(IN) 		:: rdonly
integer(HID_T)       		:: file_id, root_g_id, model_g_id
real                       	:: version
character(LEN=3)  		:: cdata_out 

H5U_openVerifyFile = 0   ! if this stays zero, then everything is ok

!   Open existing file HDFfileName.
file_id = H5U_openFile(char(HDFfileName), .TRUE.)  ! .TRUE. = readonly; .FALSE. = readwrite

!   Open the root group. (everything belongs to this top level group)
root_g_id = H5U_openGroup("/",file_id)

!   Open the Data Model group.
model_g_id = H5U_openGroup(HDF_DATAMODEL,root_g_id)

!  get the HDF_ModelType attribute for the DATA_MODEL
cdata_out = char(H5U_getGroupAttributeString(HDF_ModelType,model_g_id))

if (HDF_CurrentFileType.ne.cdata_out) then
  call H5U_printMessage('Incorrect HDF file type ')		! print a message
  call HDF_pop(.TRUE.)							! close all open objects
  H5U_openVerifyFile = -1						! and return an error status
  return  
end if 

!  get the HDF_ModelType attribute for the DATA_MODEL
version = H5U_getGroupAttributeReal(HDF_ModelVersion,model_g_id)

if (HDF_CurrentFileVersion.ne.version) then
  call H5U_printMessage('Incorrect file version')              ! print a message
  call HDF_pop(.TRUE.)							! close all open objects
  H5U_openVerifyFile = -1						! and return an error status
  return  
end if 

! and close the Data Model group
call HDF_pop

end function H5U_openVerifyFile







! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> (Somewhat) gracefully handle errors
!! @param error non-zero integer
!! @param OffendingRoutine name of the routine that caused the error
!! @param NonFatal (optional parameter) Fatal or non-fatal error (logical)
subroutine H5U_handleError(error,OffendingRoutine, NonFatal)

IMPLICIT NONE

integer,INTENT(IN)        	:: error            	! returned error code
character(LEN=*),INTENT(IN) 	:: OffendingRoutine	! name of offending routine
logical,OPTIONAL,INTENT(IN)  	:: NonFatal  		! if true, then report the error but don't stop

write (*,*) '  Error code : ',error
write (*,*) '     returned by routine ',OffendingRoutine

if (.not.present(NonFatal)) STOP  ! this is not very graceful, but it'll do the job ...

end subroutine H5U_handleError


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> repot the error, clean up all open links, and abort
!! @param DDreturn error string
!!  @param OffendingRoutine name of offending routine
subroutine H5U_cleanAbort(DDreturn,OffendingRoutine)

IMPLICIT NONE

character(LEN=*),INTENT(IN) 	:: OffendingRoutine   ! name of offending routine
type(string),INTENT(IN)	:: DDreturn
integer				:: error
type(string)			:: success
 
if (DDreturn.ne.'') then
  call H5U_handleError(-1,OffendingRoutine,.TRUE.)					! print an error message
  call HDF_pop(.TRUE.)											! close all the open objects
  call h5close_f(error)											! close the Fortran HDF5 interface
  success = 'Program aborted due to : Error in '//OffendingRoutine  		! and abort with a friendly message
  return
end if

end subroutine H5U_cleanAbort



! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> Print a message 
!! @param mess message string
subroutine H5U_printMessage(mess)

IMPLICIT NONE

character(LEN=*) ::  mess

write (*,"('Message: ',A)") mess

end subroutine H5U_printMessage



! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> Create an HDF5 file
!! @param filename HDF5 filename string
!! @todo Turn off error handling
!! @return file_id file identifier (integer)
function H5U_createFile() result(success)
 
IMPLICIT NONE

integer(HID_T)                	:: file_id ! file identifier
integer				:: error  ! error flag
type(string)			:: success

success = ''

! Create a new file using default properties.
call h5fcreate_f(char(HDFfileName), H5F_ACC_TRUNC_F, file_id, error)
if (error.ne.0) then
  call H5U_handleError(error,'H5U_createFile: error creating file')
  success = 'H5U_createFile: Error creating file'
end if

call HDF_push('f',file_id)


end function H5U_createFile


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> Open an HDF5 file
!! @param filename HDF5 filename string
!! @param readOnly switch for ReadWrite or ReadOnly access mode
!! @todo Turn off error handling
!! @return file_id file identifier (integer)
function H5U_openFile(filename,rdOnly) result(file_id)

IMPLICIT NONE

character(LEN=*), INTENT(IN)  	:: filename  	! file name (input)
logical,INTENT(IN)          	:: rdOnly 	! open access mode
integer(HID_T)                	:: file_id 	! file identifier
integer                    	:: error  	! error flag

file_id = -1

if (HDFbeverbose) call H5U_printMessage('H5U_openFile')

if ( rdOnly ) then
  call H5Fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error);
else 
  call H5Fopen_f(filename, H5F_ACC_RDWR_F, file_id, error);
end if

if (error.ne.0) call H5U_handleError(error,'H5U_openFile')

! and put the file_id onto the HDF_stack
call HDF_push('f',file_id)

end function H5U_openFile


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> Open a group 
!! @param groupname group name string
!! @param parent_id  HID_T for the parent of this group
!! @return group_id group identifier (integer)
function H5U_openGroup(groupname,parent_id) result(group_id)

IMPLICIT NONE

character(LEN=*), INTENT(IN)  ::  groupname  		! group name (input)
integer(HID_T)                :: group_id, parent_id 	! identifiers
integer                       :: error  		! error flag

group_id = -1

if (HDFbeverbose) call H5U_printMessage('H5U_openGroup')

call H5Gopen_f(parent_id, groupname, group_id, error)

if (error.ne.0) call H5U_handleError(error,'H5U_openGroup')

! and put the group_id onto the HDF_stack
call HDF_push('g',group_id)

end function H5U_openGroup


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> Create a group 
!! @param groupname group name string
!! @param parent_id  HID_T for the parent of this group
!! @return group_id group identifier (integer)
function H5U_createGroup(groupname,parent_id) result(group_id)

IMPLICIT NONE

character(LEN=*), INTENT(IN)  ::  groupname  		! group name (input)
integer(HID_T)                :: group_id, parent_id 	!  identifiers
integer                       :: error  		! error flag

group_id = -1

if (HDFbeverbose) call H5U_printMessage('H5U_createGroup')

call H5Gcreate_f(parent_id, groupname, group_id, error)

if (error.ne.0) call H5U_handleError(error,'H5U_createGroup')

! and put the group_id onto the HDF_stack
call HDF_push('g',group_id)

end function H5U_createGroup





! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> Open a data set 
!! @param datasetname data set name string
!! @param parent_id  HID_T for the parent of this group
!! @return group_id group identifier (integer)
function H5U_openDataSet(datasetname,parent_id) result(dataset_id)

IMPLICIT NONE

character(LEN=*), INTENT(IN)  ::  datasetname  		! group name (input)
integer(HID_T)                :: dataset_id, parent_id 	!  identifiers
integer                       :: error  			! error flag

dataset_id = -1

if (HDFbeverbose) call H5U_printMessage('H5U_openDataSet')

call H5Dopen_f(parent_id, datasetname, dataset_id, error)

if (error.ne.0) call H5U_handleError(error,'H5U_openDataSet')

! and put the group_id onto the HDF_stack
call HDF_push('d',dataset_id)

end function H5U_openDataSet



! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> Open an attribute 
!! @param attrname attribute name string
!! @param parent_id  HID_T for the parent of this group
!! @return attr_id attribute identifier (integer)
function H5U_openAttribute(attrname,parent_id,optionalAttribute) result(attr_id)

IMPLICIT NONE

character(LEN=*), INTENT(IN)  	::  attrname  			! attribute name (input)
integer(HID_T)                	:: attr_id, parent_id 		!  identifiers
integer                       	:: error, error2, printflag  	! error flag
logical,optional		:: optionalAttribute

attr_id = -1

if (HDFbeverbose) call H5U_printMessage('H5U_openAttribute')

if (present(optionalAttribute)) then
! turn error messaging off
  printflag = 0
  call h5eset_auto_f(printflag, error2)
end if

call H5Aopen_f(parent_id, attrname, attr_id, error)

if (error.ne.0) then
   call H5U_handleError(error,'H5U_openAttribute',.TRUE.)
else
! and put the group_id onto the HDF_stack
  call HDF_push('a',attr_id)
end if

if (present(optionalAttribute)) then
! turn error messaging off
  printflag = 1
  call h5eset_auto_f(printflag, error2)
end if

end function H5U_openAttribute


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> Create an attribute 
!! @param attrname attribute name string
!! @param parent_id  HID_T for the parent of this group
!! @return attr_id attribute identifier (integer)
function H5U_createAttributeString(attrname,parent_id,inpString,arank,adims) result(success)

IMPLICIT NONE

type(string),INTENT(IN)	:: inpString
character(LEN=*),INTENT(IN)   	:: attrname  ! attribute name (input)
integer,INTENT(IN)		:: arank
integer(HSIZE_T),INTENT(IN)	:: adims(*)
integer(HSIZE_T)		:: inplen
integer(HID_T)                	:: attr_id, parent_id, aspace_id, atype_id !  identifiers
integer                       	:: error  ! error flag
type(string)			:: success
success = ''

if (HDFbeverbose) call H5U_printMessage('H5U_createAttributeString')

! Create data space for the attribute.
call h5screate_simple_f(arank, adims, aspace_id, error)
if (error.ne.0) call H5U_handleError(error,'H5U_createAttributeString: error calling h5screate_simple_f ')
   
! Create datatype for the attribute.
call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
if (error.ne.0) call H5U_handleError(error,'H5U_createAttributeString: error calling h5tcopy_f ')
inplen = len(inpString)
call h5tset_size_f(atype_id, inplen, error)
if (error.ne.0) call H5U_handleError(error,'H5U_createAttributeString: error calling h5tset_size_f ')

! Create dataset attribute.
call h5acreate_f(parent_id, attrname, atype_id, aspace_id, attr_id, error)
if (error.ne.0) call H5U_handleError(error,'H5U_createAttributeString: error calling h5acreate_f ')

! and, finally, write the attribute value
call h5awrite_f(attr_id, atype_id, char(inpString), adims, error)
if (error.ne.0) call H5U_handleError(error,'H5U_createAttributeString: error calling h5awrite_f ')

! Close the attribute.
call h5aclose_f(attr_id, error)
if (error.ne.0) call H5U_handleError(error,'H5U_createAttributeString: error calling h5aclose_f ')


end function H5U_createAttributeString


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> Create an attribute 
!! @param attrname attribute name string
!! @param parent_id  HID_T for the parent of this group
!! @return attr_id attribute identifier (integer)
function H5U_createAttributeReal(attrname,parent_id,inpReal,arank,adims) result(success)

IMPLICIT NONE

real,INTENT(IN)			:: inpReal
character(LEN=*),INTENT(IN)   	:: attrname  ! attribute name (input)
integer,INTENT(IN)		:: arank
integer(HSIZE_T),INTENT(IN)	:: adims(*)
integer(HSIZE_T)		:: inplen
integer(HID_T)                	:: attr_id, parent_id, aspace_id, atype_id !  identifiers
integer                      	:: error  ! error flag
type(string)			:: success
success = ''

if (HDFbeverbose) call H5U_printMessage('H5U_createAttributeReal')

! Create data space for the attribute.
call h5screate_simple_f(arank, adims, aspace_id, error)
if (error.ne.0) call H5U_handleError(error,'H5U_createAttributeReal: error calling h5screate_simple_f ')
   
! Create datatype for the attribute.
call h5tcopy_f(H5T_NATIVE_REAL, atype_id, error)
if (error.ne.0) call H5U_handleError(error,'H5U_createAttributeReal: error calling h5tcopy_f ')
!inplen = len(inpString)
!call h5tset_size_f(atype_id, inplen, error)
!if (error.ne.0) call H5U_handleError(error,'H5U_createAttributeString: error calling h5tset_size_f ')

! Create dataset attribute.
call h5acreate_f(parent_id, attrname, atype_id, aspace_id, attr_id, error)
if (error.ne.0) call H5U_handleError(error,'H5U_createAttributeReal: error calling h5acreate_f ')

! and, finally, write the attribute value
call h5awrite_f(attr_id, atype_id, inpReal, adims, error)
if (error.ne.0) call H5U_handleError(error,'H5U_createAttributeReal: error calling h5awrite_f ')

! Close the attribute.
call h5aclose_f(attr_id, error)
if (error.ne.0) call H5U_handleError(error,'H5U_createAttributeReal: error calling h5aclose_f ')


end function H5U_createAttributeReal


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> Get an attribute but leave the attribute object closed to the calling routine
!! @param attrname attribute name string
!! @param parent_id  HID_T for the parent of this group
!! @return attr_id attribute identifier (integer)
type(string) function H5U_getGroupAttributeString(attrname,parent_id) 

IMPLICIT NONE

character(LEN=*), INTENT(IN)  	::  attrname  ! attribute name (input)
integer(HID_T)                	:: attr_id, parent_id, attr_type !  identifiers
character(LEN=1024)      	:: attribute
integer                      	:: error  ! error flag
integer(HSIZE_T),dimension(1) 	:: cdata_dims  ! data identifier for attribute

! get the attribute id
attr_id = H5U_openAttribute(attrname,parent_id)

! get the attribute type
if (HDFbeverbose) call H5U_printMessage('calling h5aget_type_f')
call h5aget_type_f(attr_id, attr_type, error)
if (error.ne.0) call H5U_handleError(error,'h5aget_type_f in H5U_getGroupAttributeString')

!   get the attribute value
attribute = ''
if (HDFbeverbose) call H5U_printMessage('calling h5aread_f')
call h5aread_f(attr_id, attr_type, attribute, cdata_dims, error)
if (error.ne.0) call H5U_handleError(error,'h5aread_f')

! and return the value
H5U_getGroupAttributeString = trim(attribute)

! and get rid of the attribute object
call HDF_pop

end function H5U_getGroupAttributeString



! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> Get an attribute but leave the attribute object closed to the calling routine
!! @param attrname attribute name string
!! @param parent_id  HID_T for the parent of this group
!! @return attr_id attribute identifier (integer)
function H5U_getGroupAttributeReal(attrname,parent_id) result(attribute)

IMPLICIT NONE

character(LEN=*), INTENT(IN)  ::  attrname  ! attribute name (input)
integer(HID_T)                :: attr_id, parent_id, attr_type !  identifiers
real                          :: attribute
integer                       :: error  ! error flag
integer(HSIZE_T),dimension(1) :: cdata_dims  ! data identifier for attribute

! get the attribute id
attr_id = H5U_openAttribute(attrname,parent_id)

! get the attribute type
if (HDFbeverbose) call H5U_printMessage('calling h5aget_type_f')
call h5aget_type_f(attr_id, attr_type, error)
if (error.ne.0) call H5U_handleError(error,'h5aget_type_f in H5U_getGroupAttributeReal')

!   get the attribute value
if (HDFbeverbose) call H5U_printMessage('calling h5aread_f')
call h5aread_f(attr_id, attr_type, attribute, cdata_dims, error)
if (error.ne.0) call H5U_handleError(error,'h5aread_f')

! and get rid of the attribute object
call HDF_pop

end function H5U_getGroupAttributeReal



! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> Get an attribute but leave the attribute object closed to the calling routine
!! @param attrname attribute name string
!! @param parent_id  HID_T for the parent of this group
!! @return attr_id attribute identifier (integer)
function H5U_getGroupAttributeInteger(attrname,parent_id) result(attribute)

IMPLICIT NONE

character(LEN=*), INTENT(IN)  ::  attrname  ! attribute name (input)
integer(HID_T)                :: attr_id, parent_id, attr_type !  identifiers
integer                       :: attribute
integer                       :: error  ! error flag
integer(HSIZE_T),dimension(1) :: cdata_dims  ! data identifier for attribute

! get the attribute id
attr_id = H5U_openAttribute(attrname,parent_id)

! get the attribute type
if (HDFbeverbose) call H5U_printMessage('calling h5aget_type_f')
call h5aget_type_f(attr_id, attr_type, error)
if (error.ne.0) call H5U_handleError(error,'h5aget_type_f in H5U_getGroupAttributeReal')

!   get the attribute value
if (HDFbeverbose) call H5U_printMessage('calling h5aread_f')
call h5aread_f(attr_id, attr_type, attribute, cdata_dims, error)
if (error.ne.0) call H5U_handleError(error,'h5aread_f')

! and get rid of the attribute object
call HDF_pop

end function H5U_getGroupAttributeInteger



! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> Get an attribute but leave the attribute object closed to the calling routine
!! @param attrname attribute name string
!! @param parent_id  HID_T for the parent of this group
!! @return optoinalAttribute optional logical indicating that this attribute is optional
type(string) function H5U_getDataAttributeString(attrname,parent_id,optionalAttribute) 

IMPLICIT NONE

character(LEN=*), INTENT(IN)  	:: attrname  ! attribute name (input)
integer(HID_T)                	:: attr_id, parent_id, attr_type !  identifiers
character(LEN=1024)      	:: attribute
integer                      	:: error  ! error flag
integer(HSIZE_T),dimension(1) 	:: cdata_dims  ! data identifier for attribute
logical,optional		:: optionalAttribute

if (present(optionalAttribute)) then 
! get the attribute id
  attr_id = H5U_openAttribute(attrname,parent_id,optionalAttribute)
else
! get the attribute id
  attr_id = H5U_openAttribute(attrname,parent_id)
end if

write (*,*) 'Getting ',attrname

if (attr_id.ne.-1) then
! get the attribute type
	if (HDFbeverbose) call H5U_printMessage('calling h5aget_type_f')
	call h5aget_type_f(attr_id, attr_type, error)
	if (error.ne.0) call H5U_handleError(error,'h5aget_type_f in H5U_getDataAttributeString')
	
	!   get the attribute value
	attribute = ''
	if (HDFbeverbose) call H5U_printMessage('calling h5aread_f')
	call h5aread_f(attr_id, attr_type, attribute, cdata_dims, error)
	if (error.ne.0) call H5U_handleError(error,'h5aread_f')
	
	! and return the value
	H5U_getDataAttributeString = trim(attribute)
	
	! and get rid of the attribute object
	call HDF_pop
	
write (*,*) '     -> Found ',trim(attribute)
end if

end function H5U_getDataAttributeString





! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> Close an HDF5 file
!! @param filename HDF5 filename string
!! @param readOnly switch for ReadWrite or ReadOnly access mode
!! @todo Turn off error handling
!! @return error error status (0) if ok, non-zero otherwise
!integer(HERR_T) function H5U_closeFile(file_id)
!
!IMPLICIT NONE
!
!integer(HERR_T)  :: err_t   ! error variable
!integer(HID_T),INTENT(INOUT)  :: file_id ! file identifier
!integer(SIZE_T)   :: num_open  ! number of open identifiers
!integer(SIZE_T)   :: charsRead  ! number of characters read by H5Iget_name_f
!integer               :: error, &  ! error flag
!                             i  ! loop counter
!integer(HID_T),allocatable :: attr_ids(:)
!character(LEN=1024) :: name
!
!err_t = 1
!if (file_id.lt.0) then
!    H5U_closeFile = 1
!    return
!end if
!
!! Get the number of open identifiers of all types
!! except files
!  call H5Fget_obj_count_f(file_id, ior(ior(H5F_OBJ_DATASET_F,H5F_OBJ_GROUP_F),ior(H5F_OBJ_DATATYPE_F,H5F_OBJ_ATTR_F)), &
!                                       num_open, error)
!  if (error.ne.0) call H5U_handleError(error,'H5U_closeFile, at H5Fget_obj_count_f')
!                      
!  if (num_open.gt.0) then
!    if (HDFverbose) write (*,*) 'WARNING: Some IDs weren''t closed. Closing them.' 
!! get the ids
!    allocate(attr_ids(0:num_open-1))
!    attr_ids = 0
!    call H5Fget_obj_ids_f(file_id, ior(ior(H5F_OBJ_DATASET_F,H5F_OBJ_GROUP_F),ior(H5F_OBJ_DATATYPE_F,H5F_OBJ_ATTR_F)), &
!                                     num_open, attr_ids, error)
!    if (error.ne.0) call H5U_handleError(error,'H5U_closeFile, at H5Fget_obj_ids_f')
!    do i=0,num_open-1
!      name = ''
!      call H5Iget_name_f(attr_ids(i), name, 1024_SIZE_T, charsRead, error);
!      if (error.ne.0) call H5U_handleError(error,'H5U_closeFile, at H5Iget_name_f')
!      if (charsRead.lt.0) then
!        if (HDFverbose) write (*,*) 'Error Trying to get the name of an hdf object that was not closed. This is probably pretty bad.' 
!        H5U_closeFile = -1
!        return
!      end if
!      if (HDFverbose) write (*,*) 'H5 Object left open. Id=',attr_ids(i),'; Name= ', trim(name)
!      err_t = H5U_closeHDF5Object(attr_ids(i))
!      if (err_t.ne.0) call H5U_handleError(err_t,'H5U_closeFile, at H5U_closeHDF5Object.')
!    end do
!    deallocate(attr_ids)
!  end if
!  
!  call H5Fclose_f(file_id,err_t)
!  if (err_t.ne.0) call H5U_handleError(err_t,'H5U_closeFile, at H5Fclose_f; error closing HDF5 file ')
!
!  file_id= -1   ! reset the file_id to some negative number to indicate file is closed.
!  H5U_closeFile = err_t
!
!end function H5U_closeFile
!
!
!! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!!
!! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!!> Close and HDF5 object
!!! @param obj_id  object id
!!! @todo Turn off error handling
!!! @return error error status 
!integer(HERR_T) function H5U_closeHDF5Object(obj_id)
!
!IMPLICIT NONE
!
!integer(HID_T),INTENT(IN)  :: obj_id  ! object identifier
!integer(HERR_T) :: err_t  ! error variable
!integer(H5I_TYPE_T) :: obj_type
!integer(SIZE_T)   :: charsRead  ! number of characters read by H5Iget_name_f
!integer               :: error  ! error flag
!character(LEN=1024) :: name
!logical :: nf  ! fatal or non-fatal error
!
!  if (obj_id.lt.0) then  ! Object was not valid.
!    H5U_closeHDF5Object = 0
!    return
!  end if
!
!!  H5I_type_t obj_type; 
!  name = ''
!  call H5Iget_type_f(obj_id, obj_type, error);
!  if (error.ne.0) call H5U_handleError(error,'H5U_closeHDF5Object, after H5Iget_type_f')
!  
!  call H5Iget_name_f( obj_id, name, 1024, charsRead, error);
!  if (error.ne.0) call H5U_handleError(error,'H5U_closeHDF5Object, after H5Fget_name_f')
!  if (charsRead.lt.0) then
!        if (HDFverbose) write (*,*) 'Error Trying to get the name of an hdf object that was not closed. This is probably pretty bad.' 
!        H5U_closeHDF5Object = -1
!        return
!  end if
!
!  nf = .TRUE.    ! none of these error messages are fatal
!  select case(obj_type)
!  case (H5I_FILE_F)
!    call H5Fclose_f(obj_id, err_t)
!    if (err_t.ne.0) call H5U_handleError(err_t,'H5U_closeHDF5Object, case(obj_type) H5Fclose_f',nf)
!  case (H5I_GROUP_F)
!    call H5Gclose_f(obj_id, err_t)
!    if (err_t.ne.0) call H5U_handleError(err_t,'H5U_closeHDF5Object, case(obj_type) H5Gclose_f, H5 Group Object left open.',nf)
!  case (H5I_DATASET_F)
!    call H5Dclose_f(obj_id, err_t)
!    if (err_t.ne.0) call H5U_handleError(err_t,'H5U_closeHDF5Object, case(obj_type) H5Dclose_f, H5 Dataset Object left open.',nf)
! case (H5I_ATTR_F)
!    call H5Aclose_f(obj_id, err_t)
!    if (err_t.ne.0) call H5U_handleError(err_t,'H5U_closeHDF5Object, case(obj_type) H5Aclose_f, H5 Attribute Object left open.',nf)
! case (H5I_DATATYPE_F)
!    call H5Tclose_f(obj_id, err_t)
!    if (err_t.ne.0) call H5U_handleError(err_t,'H5U_closeHDF5Object, case(obj_type) H5Tclose_f, H5 DataType Object left open. ',nf)
! case (H5I_DATASPACE_F)
!    call H5Sclose_f(obj_id, err_t)
!    if (err_t.ne.0) call H5U_handleError(err_t,'H5U_closeHDF5Object, case(obj_type) H5Sclose_f, H5 Data Space Object left open.',nf)
!  case default
!    err_t = -1
!    call H5U_handleError(err_t,'H5U_closeHDF5Object, case(obj_type) H5Sclose_f, unknown HDF object for closing.',nf)
!  end select 
!
!  H5U_closeHDF5Object = err_t
!
!end function H5U_closeHDF5Object


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!> Returns the full path to the object (with leading / removed
!! @param loc_id  location id
!! @todo Turn off error handling
!! @return objPath a string containing the object path 
!character(LEN=*) function H5U_getObjectPath(loc_id)  ! we've dropped the bool trim variable from the C++ version
!
!IMPLICIT NONE
!
!integer(HID_T), INTENT(IN) :: loc_id
!character(LEN=1024)  :: obj_name    ! we'll assume that the name is never longer than 1024 characters
!character(LEN=*)  :: obj_path
!integer(SIZE_T)     :: name_size
!integer                 :: error, &  ! error flag
!                               sloc  ! character location of forward leaning slash
!
!call H5Iget_name_f(loc_id, obj_name, 1024_SIZE_T, name_size, error)
!if (error.ne.0) call H5U_handleError(error,'H5U_getObjectPath, at H5Iget_name_f')
! 
!! next, trim the string, copy into objPath, and remove a possible leading '/' character 
!objPath = obj_name(1:name_size)
!if (objPath.eq.'/') then 
!  objPath = ''
!else
!  sloc = index(objPath,'/')
!  if (sloc.eq.1) then
!    slen = len_trim(objPath)
!    objPath = trim(objPath(2:slen))
!  end if
!end if
!
!H5U_getObjectPath = objPath
!
!end function H5U_getObjectPath


end module HDF5Utilities
