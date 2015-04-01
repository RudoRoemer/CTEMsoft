
program hdftest

  use local
  use HDF5
  use typedefs
  use HDFsupport
  use NameListTypedefs
  use NameListHDFwriters
  use ISO_C_BINDING
  
  IMPLICIT NONE

  CHARACTER(fnlen)             :: filename, groupname, dataname
  CHARACTER(fnlen)             :: dataset, textfile
  CHARACTER(fnlen)             :: nmlname = "CTEMKossel.nml"

  INTEGER(SIZE_T)              :: sdim 
  INTEGER(HID_T)               :: file, filetype, space, dset ! Handles
  INTEGER                      :: hdferr
  INTEGER(HSIZE_T), DIMENSION(1:1) :: dims
  INTEGER(HSIZE_T), DIMENSION(1:2) :: dims2
  INTEGER(HSIZE_T), DIMENSION(1:2) :: maxdims
  
! character(len=1,kind=c_char)          :: chararr(256)
  character(len=1)          :: chararr(256)
  character(len=1)          :: chararr2(256,256)

  character(len=1),allocatable :: rdchararr(:)

  TYPE(C_PTR), ALLOCATABLE, TARGET :: wdata(:)
  CHARACTER(len=fnlen, KIND=c_char), DIMENSION(1), TARGET  :: line 
  CHARACTER(len=fnlen, KIND=c_char), ALLOCATABLE, TARGET  :: lines(:) 

  TYPE(C_PTR), DIMENSION(:), ALLOCATABLE, TARGET :: rdata ! Read buffer
  CHARACTER(len = fnlen, kind=c_char),  POINTER  :: data ! A pointer to a Fortran string
  TYPE(C_PTR) :: f_ptr

  INTEGER                           :: i, j, length, nlines
  integer(kind=irg)                 :: intarr(8), dim0, dim1, intarr2(3,3), printflag
  integer*4                         :: i1, i2
  real(kind=sgl)                    :: fltarr(8)
  real(kind=dbl)                    :: dblarr(8)
  integer(kind=irg),allocatable     :: rdintarr(:), rdintarr2(:,:)
  real(kind=sgl),allocatable        :: rdfltarr(:)
  real(kind=dbl),allocatable        :: rddblarr(:)
  logical                           :: HDFstatus

  type(HDFobjectStackType),pointer  :: HDF_head
  type(HDFobjectStackType),pointer  :: HDF_tail

  character(11)                     :: dstr
  character(15)                     :: tstrb
  character(15)                     :: tstre
  character(fnlen)                  :: prn

  type(KosselNameListType)          :: knl
  type(RFZNameListType)             :: rnl

  nullify(HDF_head)
  nullify(HDF_tail)

  rnl%pgnum = 32
  rnl%nsteps = 100
  rnl%outname = 'RFZtest.data'

  knl%stdout = 6
  knl%numthick = 10
  knl%npix = 401
  knl%maxHOLZ = 1
  knl%nthreads = 8
  knl%k = (/ 1, 1, 2 /)
  knl%fn = (/ 1, 2, 2 /)
  knl%voltage = 30.0
  knl%dmin = 0.05
  knl%convergence = 0.5
  knl%startthick = 5.0
  knl%thickinc = 5.0
  knl%minten = 0.2
  knl%xtalname = 'Cu.xtal'
  knl%outname = 'Kosseltest.data'
!
! Initialize FORTRAN interface.
!
CALL h5open_f(hdferr)

printflag = 0
call h5eset_auto_f(printflag, hdferr)

filename = 'doesnotexist.h5'
call h5fis_hdf5_f(trim(filename), HDFstatus, hdferr) 
write (*,*) filename,' HDF5 ? ',HDFstatus, hdferr

filename = 'EMKossel.nml'
call h5fis_hdf5_f(trim(filename), HDFstatus, hdferr) 
write (*,*) filename,' HDF5 ? ',HDFstatus, hdferr

filename = 'test.h5'
call h5fis_hdf5_f(trim(filename), HDFstatus, hdferr) 
write (*,*) filename,' HDF5 ? ',HDFstatus, hdferr



call timestamp(datestring=dstr, timestring=tstrb)
tstre = tstrb

write (*,*) 'date = ', dstr

! Create a new file using the default properties.
hdferr =  HDF_createFile(filename, HDF_head, HDF_tail)
write (*,*) 'file created ',hdferr

! write the EMheader to the file
prn = 'hdftest.f90'
write (*,*) 'creating EMheader'
call HDF_writeEMheader(HDF_head, HDF_tail, dstr, tstrb, tstre, prn)
write (*,*) 'EMheader created '

! create a namelist group to write all the namelist files into
groupname = "NMLfiles"
hdferr = HDF_createGroup(groupname, HDF_head, HDF_tail)

! read the text file and write the array to the file
nmlname = "EMKossel.nml"
dataset = 'KosselNML'
hdferr = HDF_writeDatasetTextFile(dataset, nmlname, HDF_head, HDF_tail)

! read the text file and write the array to the file
nmlname = "CTEMsampleRFZ.nml"
dataset = 'sampleRFZNML'
hdferr = HDF_writeDatasetTextFile(dataset, nmlname, HDF_head, HDF_tail)

call HDF_pop(HDF_head)

! create a NMLparameters group to write all the namelist entries into
groupname = "NMLparameters"
hdferr = HDF_createGroup(groupname, HDF_head, HDF_tail)

call HDFwriteKosselNameList(HDF_head, HDF_tail, knl)
call HDFwriteRFZNameList(HDF_head, HDF_tail, rnl)

! and leave this group
call HDF_pop(HDF_head)

! then the remainder of the data in a EMData group
groupname = 'EMData'
hdferr = HDF_createGroup(groupname, HDF_head, HDF_tail)

allocate(lines(1))
lines(1) = 'This is line 1'
!lines(2) = 'and this is line 2'
dataset = 'StringTest'
hdferr = HDF_writeDatasetStringArray(dataset, lines, 1, HDF_head, HDF_tail)
deallocate(lines)

intarr = (/ 1, 2, 3, 4, 5, 6, 7, 8 /)
dataset = 'intarr1D'
dims = shape(intarr)
dim0 = dims(1)
hdferr = HDF_writeDatasetIntegerArray1D(dataset, intarr, dim0, HDF_head, HDF_tail)

intarr2 = reshape( (/ 1,2,3,4,5,6,7,8,9 /), (/3,3/))
dataset = 'intarr2D'
dims2 = shape(intarr2)
dim0 = dims2(1)
dim1 = dims2(2)
hdferr = HDF_writeDatasetIntegerArray2D(dataset, intarr2, dim0, dim1, HDF_head, HDF_tail)

fltarr = (/ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 /)
dataset = 'fltarr1D'
dims = shape(fltarr)
dim0 = dims(1)
hdferr = HDF_writeDatasetFloatArray1D(dataset, fltarr, dim0, HDF_head, HDF_tail)

dblarr = (/ 1.D0, 2.D0, 3.D0, 4.D0, 5.D0, 6.D0, 7.D0, 8.D0 /)
dataset = 'dblarr1D'
dims = shape(dblarr)
dim0 = dims(1)
hdferr = HDF_writeDatasetDoubleArray1D(dataset, dblarr, dim0, HDF_head, HDF_tail)

do i1=0,255
  chararr(i1+1) = char(i1)
end do

dataset = 'chararray1D'
dim0 = 256
write (*,*) 'about to write chararray'
hdferr = HDF_writeDatasetCharArray1D(dataset, chararr, dim0, HDF_head, HDF_tail)

do i1=0,255
 do i2=0,255
  chararr2(i1+1,i2+1) = char(mod(i1+i2,256))
 end do
end do

dataset = 'chararray2D'
dim0 = 256
write (*,*) 'about to write chararray2'
hdferr = HDF_writeDatasetCharArray2D(dataset, chararr2, dim0, dim0, HDF_head, HDF_tail)


call HDF_pop(HDF_head,.TRUE.)


!===========================================================
  ! Now we begin the read section of this example.

! Open the file using the default properties.
filename = 'test.h5'
hdferr =  HDF_openFile(filename, HDF_head, HDF_tail)

! open the NMLfiles group
groupname = 'NMLfiles'
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

textfile = './testtext.txt'
hdferr = HDF_extractDatasetTextfile(dataname, textfile, HDF_head, HDF_tail)

call HDF_pop(HDF_head)

! next, read one of the integer string arrays

! open the EMData group
groupname = 'EMData'
hdferr = HDF_OpenGroup(groupname, HDF_head, HDF_tail)

dataset = 'StringTest'
lines = HDF_readDatasetStringArray(dataset, nlines, HDF_head, HDF_tail)
do i=1,nlines
  write (*,*) lines(i)
end do
deallocate(lines)

dataname = 'intarr1D'
rdintarr = HDF_readDatasetIntegerArray1D(dataname, dims, HDF_head, HDF_tail)

write (*,*) 'shape of read intarr = ',shape(rdintarr)
write (*,*) rdintarr

dataname = 'intarr2D'
rdintarr2 = HDF_readDatasetIntegerArray2D(dataname, dims2, HDF_head, HDF_tail)

write (*,*) 'shape of read intarr2 = ',shape(rdintarr2)
write (*,*) rdintarr2

dataname = 'fltarr1D'
rdfltarr = HDF_readDatasetFloatArray1D(dataname, dims, HDF_head, HDF_tail)

write (*,*) 'shape of read fltarr = ',shape(rdfltarr)
write (*,*) rdfltarr

dataname = 'dblarr1D'
rddblarr = HDF_readDatasetDoubleArray1D(dataname, dims, HDF_head, HDF_tail)

write (*,*) 'shape of read dblarr = ',shape(rddblarr)
write (*,*) rddblarr

dataname = 'chararray1D'
rdchararr = HDF_readDatasetCharArray1D(dataname, dims, HDF_head, HDF_tail)

write (*,*) 'shape of chararray = ',shape(rdchararr)
write (*,*) rdchararr



! and close all resources
call HDF_pop(HDF_head,.TRUE.)


! close the Fortran interface
call h5close_f(hdferr)


end program hdftest
