
program hdftest

use local
use HDFsupport
use hdf5
use h5lt

IMPLICIT NONE

     CHARACTER(LEN=5), PARAMETER :: filename = "atest"    !File name
     CHARACTER(LEN=80) :: fix_filename
     CHARACTER(LEN=9), PARAMETER :: dsetname = "atestdset"        !Dataset name
     CHARACTER(LEN=11), PARAMETER :: aname = "attr_string"   !String Attribute name
     CHARACTER(LEN=14), PARAMETER :: aname2 = "attr_character"!Character Attribute name
     CHARACTER(LEN=11), PARAMETER :: aname3 = "attr_double"   !DOuble Attribute name
     CHARACTER(LEN=9), PARAMETER :: aname4 = "attr_real"      !Real Attribute name
     CHARACTER(LEN=12), PARAMETER :: aname5 = "attr_integer"  !Integer Attribute name
     CHARACTER(LEN=9), PARAMETER :: aname6 = "attr_null"     !Null Attribute name

     !
     !data space rank and dimensions
     !
     INTEGER, PARAMETER :: RANK = 2
     INTEGER, PARAMETER :: NX = 4
     INTEGER, PARAMETER :: NY = 5



     INTEGER(HID_T) :: file_id       ! File identifier
     INTEGER(HID_T) :: dset_id       ! Dataset identifier
     INTEGER(HID_T) :: dataspace     ! Dataspace identifier for dataset

     INTEGER(HID_T) :: attr_id        !String Attribute identifier
     INTEGER(HID_T) :: attr2_id       !Character Attribute identifier
     INTEGER(HID_T) :: attr3_id       !Double Attribute identifier
     INTEGER(HID_T) :: attr4_id       !Real Attribute identifier
     INTEGER(HID_T) :: attr5_id       !Integer Attribute identifier
     INTEGER(HID_T) :: attr6_id       !Null Attribute identifier
     INTEGER(HID_T) :: aspace_id      !String Attribute Dataspace identifier
     INTEGER(HID_T) :: aspace2_id     !Character Attribute Dataspace identifier
     INTEGER(HID_T) :: aspace6_id     !Null Attribute Dataspace identifier
     INTEGER(HID_T) :: dtype_id       !
     INTEGER(HID_T) :: atype_id       !String Attribute Datatype identifier
     INTEGER(HID_T) :: atype2_id      !Character Attribute Datatype identifier
     INTEGER(HID_T) :: atype3_id      !Double Attribute Datatype identifier
     INTEGER(HID_T) :: atype4_id      !Real Attribute Datatype identifier
     INTEGER(HID_T) :: atype5_id      !Integer Attribute Datatype identifier
     INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/2/) ! Attribute dimension
     INTEGER(HSIZE_T), DIMENSION(1) :: adims2 = (/1/) ! Attribute dimension
     INTEGER     ::   arank = 1                      ! Attribure rank
     INTEGER(SIZE_T) :: attrlen    ! Length of the attribute string

     INTEGER(HID_T) :: attr_space     !Returned String Attribute Space identifier
     INTEGER(HID_T) :: attr2_space    !Returned other Attribute Space identifier
     INTEGER(HID_T) :: attr_type      !Returned Attribute Datatype identifier
     INTEGER(HID_T) :: attr2_type      !Returned CHARACTER Attribute Datatype identifier
     INTEGER(HID_T) :: attr3_type      !Returned DOUBLE Attribute Datatype identifier
     INTEGER(HID_T) :: attr4_type      !Returned REAL Attribute Datatype identifier
     INTEGER(HID_T) :: attr5_type      !Returned INTEGER Attribute Datatype identifier
     INTEGER(HID_T) :: attr6_type      !Returned NULL Attribute Datatype identifier
     INTEGER        :: num_attrs      !number of attributes
     INTEGER(HSIZE_T) :: attr_storage   ! attributes storage requirements .MSB.
     CHARACTER(LEN=256) :: attr_name    !buffer to put attr_name
     INTEGER(SIZE_T)    ::  name_size = 80 !attribute name length

     CHARACTER(LEN=35), DIMENSION(2) ::  attr_data  ! String attribute data
     CHARACTER(LEN=35), DIMENSION(2) ::  aread_data ! Buffer to put read back
                                               ! string attr data
     CHARACTER ::  attr_character_data = 'A'
     DOUBLE PRECISION,  DIMENSION(1) ::  attr_double_data = 3.459
     REAL,         DIMENSION(1) ::  attr_real_data = 4.0
     INTEGER,      DIMENSION(1) ::  attr_integer_data = 5


     CHARACTER :: aread_character_data ! variable to put read back Character attr data
     INTEGER, DIMENSION(1)  :: aread_integer_data ! variable to put read back integer attr data
     INTEGER, DIMENSION(1)  :: aread_null_data = 7 ! variable to put read back null attr data
     DOUBLE PRECISION, DIMENSION(1)   :: aread_double_data ! variable to put read back double attr data
     REAL, DIMENSION(1)  :: aread_real_data ! variable to put read back real attr data

     character(len=fnlen)          :: HDFname, datapath,g1,g2,g3,s, progname
     integer                    :: error, total_error, i ! Error flag
     INTEGER(HSIZE_T),allocatable :: data_dims(:)
     INTEGER(HSIZE_T)           :: npoints
     integer                    :: rnk, type_class
     integer(SIZE_T)            :: type_size
     integer(HID_T)             :: grp1_id, grp2_id, grp3_id
     real(kind=8),allocatable   :: buf_dbl1(:), buf_dbl2(:,:), buf_dbl3(:,:,:)
     real(kind=4),allocatable   :: buf_flt1(:)

     character(11)              :: dstring
     character(15)              :: tstring
     character(100)             :: c
     integer                    :: jj,j,ic,nlen, values(13)
     character(100)             :: nml_list(3)
integer                                     :: numf(1)
integer(HSIZE_T)                            :: numb(1)
!character(1),allocatable                    :: filebytes(:)
character,allocatable                    :: filebytes(:)
character(132)                :: line

type(HDFobjectStackType),pointer        :: HDF_head, HDF_tail

!
! Initialize FORTRAN interface.
!
CALL h5open_f(error)
write (*,*) 'Initialize Fortran interface error   = ',error

call timestamp(datestring=dstring, timestring=tstring)

!if (1.eq.0) then
!  fname = 'Nidata.h5'
!  
!  write (*,*) 'calling h5fopen_f to open file '//fname
!       CALL h5fopen_f(fname, H5F_ACC_RDWR_F, file_id, error)
!  
!  ! 
!  datapath = "/Scan 1/EBSD/Data/CI"
!  
!  call h5ltget_dataset_ndims_f(file_id, datapath, rnk, error)
!  write (*,*) ' data dimensionality = ',rnk
!  allocate(data_dims(rnk))
!  
!  call h5ltget_dataset_info_f(file_id, datapath, data_dims, type_class, type_size, error)
!  write(*,*) ' data set info : ',data_dims(1), type_class, type_size, error
!  allocate(buf_flt1(data_dims(1)))
!  
!  call h5ltread_dataset_f(file_id, datapath, H5T_NATIVE_REAL, buf_flt1, data_dims, error)
!  
!  do i=1,10
!    write (*,*) i, buf_flt1(i)
!  end do
!  
!  write (*,*) 'calling h5fclose_f to close file '
!       CALL h5fclose_f(file_id, error)
!  write (*,*) 'error   = ',error
!end if



! now let's do a test an create a new file with this dataset written to it in a different location
HDFname = 'test.h5'

write (*,*) 'opening file ',HDFname
error =  HDF_createFile(HDFname, HDF_head, HDF_tail)
write (*,*) '   error code = ',error

write (*,*) 'Creating EMheader'
progname = 'testprogram.f90'
call HDF_writeEMheader(HDF_head, HDF_tail, dstring, tstring, tstring, progname)
write (*,*) '    error code = ',error

! test the writing of text files
g1 = 'EMinputdata'
write (*,*) 'Creating group ',trim(g1)
error = HDF_createGroup(g1, HDF_head, HDF_tail) 
write (*,*) '    error code = ',error


!nml_list(1) = 'CTEMKossel.nml'
!nml_list(2) = 'CTEMMCOpenCL.nml'
!nml_list(3) = 'CTEMKosselmaster.nml'

!numf = shape(nml_list)
!rnk = 1
!dataunit = 20

! loop over the files and write them to the HDF file as individual data sets
!do i=1,numf(1)
!  g2 = trim(nml_list(i))
!  call h5gcreate_f(grp1_id,g2,grp2_id,error)
! and read the file
!  open(unit=dataunit,file=trim(nml_list(i)),form='formatted',status='old')
!  do 
!    read (dataunit,"(A)",end=10) line
!    write (*,"(A)") trim(line)
!  end do
!  10 close(dataunit)
! then create a dataset with this data
! g2 = trim(nml_list(i))
! call h5screate_simple_f(rnk, numb, dataspace, error)
! call h5dcreate_f(grp1_id, g2, H5T_NATIVE_CHARACTER, dataspace, dset_id, error)
! call h5dwrite_f(dset_id, H5T_NATIVE_CHARACTER, filebytes, numb, error)
! call h5dclose_f(dset_id,error)
! if (error.ne.0) write (*,*) ' error writing to file '
! and get rid of the string array
!   deallocate(filebytes)
!  call h5gclose_f(grp2_id,error)
!end do

!call h5gclose_f(grp1_id,error)

!g1 = 'level1'
!g2 = 'level2'
!g3 = 'CI'
!s = '/'


!call h5gcreate_f(file_id,g1,grp1_id,error)
!write (*,*) 'open group ',trim(g1),' ',error

!call h5gcreate_f(grp1_id,g2,grp2_id,error)
!write (*,*) 'open group ',trim(g2),' ',error

!datapath = trim(s)//trim(g1)//trim(s)//trim(g2)//trim(s)//trim(g3)
!write (*,*) 'datapath = ',datapath

!call h5ltmake_dataset_f(file_id, datapath, rnk, data_dims, H5T_NATIVE_REAL, buf_flt1, error)

!call h5gclose_f(grp2_id,error)
!call h5gclose_f(grp1_id,error)

!call h5fclose_f(file_id,error)

write (*,*) 'closing file'
call HDF_pop(HDF_head, .TRUE.)
call HDF_stackdump(HDF_head)

!
! Close FORTRAN interface.
!
CALL h5close_f(error)
write (*,*) 'close fortran interface error   = ',error


end program hdftest

