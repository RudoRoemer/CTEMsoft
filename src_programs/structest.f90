
program structest

IMPLICIT NONE 

! Define the "kind" parameters for single and double precision reals, 
!> single precision real kind parameter
  integer,parameter                     :: sgl = SELECTED_REAL_KIND(p=6,r=37)   
!> double precision real kind parameter
  integer,parameter                     :: dbl = SELECTED_REAL_KIND(p=13,r=200) 

! Define the "kind" parameters for short and regular integers,
!> short integer kind parameter 
  integer,parameter                     :: ish = SELECTED_INT_KIND(3) 
!> long integer kind parameter  
  integer,parameter                     :: irg = SELECTED_INT_KIND(9)



type intarr4
        integer(kind=irg)       :: array(4)
end type intarr4


type fltarr4
        real(kind=sgl)          :: array(5)
end type fltarr4


type ll
        integer(kind=irg)       :: cntr
        type(ll),pointer        :: next
end type ll


type olio
! some integer variables
        integer(kind=ish)       :: int2
        integer(kind=irg)       :: int4
        integer(kind=irg)       :: int4arr(5)
        type(intarr4)           :: int4array
        type(intarr4),pointer   :: int4arrayptr
! some real variables
        real(kind=sgl)          :: flt4
        real(kind=dbl)          :: flt8
        real(kind=dbl)          :: flt8arr(4)
        type(fltarr4)           :: flt4array
        type(fltarr4),pointer   :: flt4arrayptr
! some character variables
        character(4)            :: s4
        character(5)            :: s4arr(3)
! allocatable arrays
        real(kind=sgl),allocatable      :: fltalloc1D(:)
        real(kind=sgl),allocatable      :: fltalloc2D(:,:)
!  and a linked list
        type(ll),pointer        :: linkedlist
end type olio


! declare a variable of this type
type(olio)                      :: o

! and some auxiliary variables
integer(kind=irg)               :: i, llnum
type(intarr4),pointer           :: i4p
type(fltarr4),pointer           :: f4p
type(ll), pointer               :: llp


!=================================
! fill in all the fields of o
!=================================

! integers
o%int2 = 5_ish
o%int4 = 5_irg
o%int4arr = (/ (i, i=1_irg,5_irg ) /)
o%int4array%array = (/ (10_irg + i, i=1_irg,4_irg ) /)

allocate(i4p)
i4p%array(1:4) = (/ (20_irg + i, i=1_irg,4_irg ) /)
o%int4arrayptr => i4p

! reals
o%flt4 = 25.0_sgl
o%flt8 = 25.0_dbl
o%flt8arr = (/ (dble(i), i=1,4 ) /)
o%flt4array%array = (/ (10.0_sgl + real(i), i=1_irg,5_irg ) /)

allocate(f4p)
f4p%array(1:5) = (/ (20.0_sgl + real(i), i=1_irg,5_irg ) /)
o%flt4arrayptr => f4p

! character variables
o%s4 = 'abcd'
o%s4arr(1) = 'abcde'
o%s4arr(2) = 'fghij'
o%s4arr(3) = 'klmno'

! allocatable arrays
allocate(o%fltalloc1D(3), o%fltalloc2D(2,2))
o%fltalloc1D(1:3) = (/ 10.0, 20.0, 30.0 /)
o%fltalloc2D = reshape( (/ 100.0, 200.0, 300.0, 400.0 /), (/2,2/) )

! linked list
allocate(o%linkedlist)
nullify(o%linkedlist%next)
o%linkedlist%cntr = 1000_irg

llp => o%linkedlist

llnum = 5
do i=1,llnum
  allocate(llp%next)
  llp => llp%next
  nullify(llp%next)
  llp%cntr = 1000_irg + 500_irg * i
end do



!=================================
! and print them
!=================================

write (*,*) 'Integer components:'
write (*,*) '-------------------'

write (*,*) 'int2         = ',o%int2
write (*,*) 'int4         = ',o%int4
write (*,*) 'int4arr      = ',o%int4arr(1:5)
write (*,*) 'int4array    = ',o%int4array%array(1:4)
write (*,*) 'int4arrayptr = ',o%int4arrayptr%array(1:4)


write (*,*) 'Real components:'
write (*,*) '----------------'

write (*,*) 'flt4         = ',o%flt4
write (*,*) 'flt8         = ',o%flt8
write (*,*) 'flt8arr      = ',o%flt8arr(1:4)
write (*,*) 'flt4array    = ',o%flt4array%array(1:5)
write (*,*) 'flt4arrayptr = ',o%flt4arrayptr%array(1:5)

write (*,*) 'Character components:'
write (*,*) '---------------------'

write (*,*) 's4           = '//o%s4
write (*,*) 's4arr        = '//o%s4arr(1)//' '//o%s4arr(2)//' '//o%s4arr(3)


write (*,*) 'Allocatable array components:'
write (*,*) '-----------------------------'

write (*,*) 'fltalloc1D    = ',o%fltalloc1D(1:3)
write (*,*) 'fltalloc2D    = '
write (*,*) o%fltalloc2D(1,1),o%fltalloc2D(1,2)
write (*,*) o%fltalloc2D(2,1),o%fltalloc2D(2,2)

write (*,*) 'Linked list:'
write (*,*) '------------'

llp => o%linkedlist

i = 1

do
  if (.not.associated(llp%next)) EXIT
  write (*,*) 'Item ',i,' -> ',llp%cntr
  llp => llp%next
  i = i+1
end do


end program structest

