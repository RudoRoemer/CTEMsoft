
program hdftest

  use local
  use HDF5
  use typedefs
  use HDFsupport
  USE ISO_C_BINDING
  
  IMPLICIT NONE

  CHARACTER(fnlen)             :: filename, groupname, dataname
  CHARACTER(fnlen)             :: dataset
  CHARACTER(fnlen)             :: nmlname = "CTEMKossel.nml"

  INTEGER(HSIZE_T)             :: dim0 
  INTEGER(SIZE_T)              :: sdim 
  INTEGER(HID_T)               :: file, filetype, space, dset ! Handles
  INTEGER                      :: hdferr
  INTEGER(HSIZE_T), DIMENSION(1:1) :: dims
  INTEGER(HSIZE_T), DIMENSION(1:2) :: maxdims
  
  TYPE(C_PTR), ALLOCATABLE, TARGET :: wdata(:)
  CHARACTER(len=fnlen, KIND=c_char), DIMENSION(1), TARGET  :: line 
  CHARACTER(len=fnlen, KIND=c_char), ALLOCATABLE, TARGET  :: lines(:) 

  TYPE(C_PTR), DIMENSION(:), ALLOCATABLE, TARGET :: rdata ! Read buffer
  CHARACTER(len = fnlen, kind=c_char),  POINTER  :: data ! A pointer to a Fortran string
  TYPE(C_PTR) :: f_ptr

  INTEGER                           :: i, j, length, nlines

  type(HDFobjectStackType),pointer  :: HDF_head
  type(HDFobjectStackType),pointer  :: HDF_tail

  nullify(HDF_head)
  nullify(HDF_tail)

!
! Initialize FORTRAN interface.
!
CALL h5open_f(hdferr)

! Create a new file using the default properties.
filename = 'test.h5'
hdferr =  HDF_createFile(filename, HDF_head, HDF_tail)

! create a namelist group to write all the namelist files into
groupname = "NMLfiles"
hdferr = HDF_createGroup(groupname, HDF_head, HDF_tail)

! read the text file and write the array to the file
nmlname = "CTEMKossel.nml"
dataset = 'KosselNML'
hdferr = HDF_writeDatasetTextFile(dataset, nmlname, HDF_head, HDF_tail)

! read the text file and write the array to the file
nmlname = "CTEMsampleRFZ.nml"
dataset = 'sampleRFZNML'
hdferr = HDF_writeDatasetTextFile(dataset, nmlname, HDF_head, HDF_tail)

call HDF_pop(HDF_head,.TRUE.)


!===========================================================
  ! Now we begin the read section of this example.

! Open the file using the default properties.
filename = 'test.h5'
hdferr =  HDF_openFile(filename, HDF_head, HDF_tail)

! open the NMLfiles group
hdferr = HDF_OpenGroup(groupname, HDF_head, HDF_tail)

! read a dataset
dataname = 'KosselNML'
lines = HDF_readDatasetStringArray(dataname, nlines, HDF_head, HDF_tail) 
write (*,*) 'data set name : ',trim(dataname),':  number of lines read = ',nlines

do i=1,nlines
  write (*,*) lines(i)
end do
deallocate(lines)

! read a dataset
dataname = 'sampleRFZNML'
lines = HDF_readDatasetStringArray(dataname, nlines, HDF_head, HDF_tail) 
write (*,*) 'data set name : ',trim(dataname),':  number of lines read = ',nlines

do i=1,nlines
  write (*,*) lines(i)
end do
deallocate(lines)

! and close all resources
call HDF_pop(HDF_head,.TRUE.)


! close the Fortran interface
call h5close_f(hdferr)


end program hdftest

