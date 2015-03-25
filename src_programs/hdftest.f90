
program hdftest

  USE HDF5
  USE ISO_C_BINDING
  
  IMPLICIT NONE
  integer,PARAMETER            :: fnlen=132, dataunit = 20

  CHARACTER(LEN=20), PARAMETER :: filename = "test.h5"
  CHARACTER(LEN=3) , PARAMETER :: dataset  = "NML"
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
  CHARACTER(len = fnlen, kind=c_char),  POINTER :: data ! A pointer to a Fortran string
  TYPE(C_PTR) :: f_ptr

  INTEGER :: i, j, length, nlines

! read the namelist file first to determine the number of lines
  open(unit=dataunit,file=trim(nmlname),form='formatted',status='old')
  nlines = 0
  do
    read (dataunit,"(A)",end=10) line(1)
    nlines = nlines + 1
  end do
  10 close(unit=dataunit,status='keep')

! then re-read the file and store all the lines in the wdata array
  dims(1) = nlines
  allocate(wdata(1:nlines), lines(1:nlines))
  open(unit=dataunit,file=trim(nmlname),form='formatted',status='old')
  do i=1,nlines
! initialize the line to null characters before each read
    do j=1,fnlen
      line(1)(j:j) = char(0)
    end do
! read the line
    read (dataunit,"(A)") line(1)
! print it to the console
    write (*,*) '->'//trim(line(1))//'<-', len(trim(line(1))), i
! find the string length and put the next character equal to C_NULL_CHAR
    j = len(trim(line(1)))+1
    line(1)(j:j) = C_NULL_CHAR
! store the line in the array
    lines(i) = line(1)
! and get the pointer to that entry in the array
    wdata(i) = C_LOC(lines(i))
  end do
  close(unit=dataunit,status='keep')

! that's the prep work

  !
  ! Initialize FORTRAN interface.
  !
  CALL h5open_f(hdferr)
  !
  ! Create a new file using the default properties.
  !
  CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file, hdferr)
  !
  ! Create file and memory datatypes.  For this example we will save
  ! the strings as C variable length strings, H5T_STRING is defined
  ! as a variable length string.
  !
  CALL H5Tcopy_f(H5T_STRING, filetype, hdferr)
  !
  ! Create dataspace.
  !
  CALL h5screate_simple_f(1, dims, space, hdferr)
  !
  ! Create the dataset and write the variable-length string data to
  ! it.
  !
  CALL h5dcreate_f(file, dataset, filetype, space, dset, hdferr)

  f_ptr = C_LOC(wdata(1))
  CALL h5dwrite_f(dset, filetype, f_ptr, hdferr )
  !
  ! Close and release resources.
  !
  CALL h5dclose_f(dset , hdferr)
  CALL h5sclose_f(space, hdferr)
  CALL H5Tclose_f(filetype, hdferr)
  CALL h5fclose_f(file , hdferr)
  !
!===========================================================
  ! Now we begin the read section of this example.
  !
  !
  ! Open file and dataset.
  !
  CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file, hdferr)
  CALL h5dopen_f(file, dataset, dset, hdferr)
  !
  ! Get the datatype.
  !
  CALL H5Dget_type_f(dset, filetype, hdferr)
  !
  ! Get dataspace and allocate memory for read buffer.
  !
  CALL H5Dget_space_f(dset, space, hdferr)
  CALL H5Sget_simple_extent_dims_f(space, dims, maxdims, hdferr)

  ALLOCATE(rdata(1:dims(1)))
  !
  ! Read the data.
  !
  f_ptr = C_LOC(rdata(1))
  CALL h5dread_f(dset, H5T_STRING, f_ptr, hdferr)
  !
  ! Output the data to the screen.
  !
  DO i = 1, dims(1)
     CALL C_F_POINTER(rdata(i), data)
     length = 0
     DO
        IF(data(length+1:length+1).EQ.C_NULL_CHAR.OR.length.GE.fnlen) EXIT
        length = length + 1
     ENDDO
     WRITE(*,'(A,"(",I0,"): ",I4," characters : ",A)') DATASET, i, length, data(1:length)
  END DO

  DEALLOCATE(rdata)
  CALL h5dclose_f(dset , hdferr)
  CALL h5sclose_f(space, hdferr)
  CALL H5Tclose_f(filetype, hdferr)
  CALL h5fclose_f(file , hdferr)

end program hdftest

