!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CTEMsoft:                                                                     !
! Copyright (c) 1996, 2001, 2002  R.A. Vowels/Marc De Graef/CMU                 !
!                                                                               !
!     This program is free software; you can redistribute it and/or modify      !
!     it under the terms of the GNU General Public License as published by      !
!     the Free Software Foundation; either version 2 of the License, or         !
!     (at your option) any later version.                                       !
!                                                                               !
!     This program is distributed in the hope that it will be useful,           !
!     but WITHOUT ANY WARRANTY; without even the implied warranty of            !
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             !
!     GNU General Public License for more details.                              !
!                                                                               !
!     You should have received a copy of the GNU General Public License         !
!     along with this program; if not, write to the Free Software               !
!     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA !
!                                                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ###################################################################
! 
!  FILE: "tiff.f90"
!                                    created: 8/21/2001 
!
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description: TIFF routines, based on the programs in the book
!               "Introduction to Fortran 90/95, Algorithms, and 
!                Structured Programming", by R.A. Vowels (1996)
! 
!  History 
!
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   ?/??/96 RAV 1.0 original
!   8/28/01 MDG 2.0 adapted for grayscale images
!  11/28/01 MDG 2.1 added kind support
! ###################################################################
!
! Here is an example showing how to use the TIFF routines.
!
!program tiff
!
!use TIFF_global
!use TIFF_f90
!
!IMPLICIT NONE
!
!integer   :: i,j
!
!! declare TIFF variables in TIFF_global
! TIFF_nx = 295
! TIFF_ny = 150
! TIFF_filename = "test.tiff"
!! allocate memory for image
! allocate(TIFF_image(0:TIFF_nx-1,0:TIFF_ny-1))
!! fill the image with whatever data you have (between 0 and 255)
! do i=0,TIFF_nx-1
!  do j=0,TIFF_ny-1
!   TIFF_image(i,j) = mod(i*i+j*j,255) 
!  end do
! end do
!! create the file
! call TIFF_Write_File
!
!end program
!
!**************************************************************

!
! The RecordLength parameter may be platform dependent;  you should
! check whether or not this package works;  if it does not work, then
! this is most likely due to an incorrect RecordLength parameter.
!
module TIFF_private

use local

 character(len=256)              :: Record
 integer(kind=irg),save          :: Rec_No=0, L=0
!
! TIFFRecordLength = 64       Compaq TRU64, DEC Alpha, ...
!                  =          MacOSX
!
!
 integer(kind=irg),parameter     :: TIFFRecordLength = 64
end module TIFF_private



module TIFF_global

use local

 integer(kind=irg)               :: TIFF_nx,TIFF_ny
 integer(kind=irg),allocatable   :: TIFF_Image(:,:)
 character(30)                   :: TIFF_filename
end module TIFF_global

!**************************************************************

module TIFF_f90

contains
!
! ###################################################################
! 
!  subroutine TIFF_Write_Byte_Into_Buffer 
!
!                                    created: ?/??/96 
!
!  Author: R.A. Vowels / Marc De Graef
!  
!  Description: write a single byte into a buffer and dump the
!               buffer to file if full
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   ?/??/96 RAV 1.0 original
!   8/28/01 MDG 2.0 commented and change of variable names
! ###################################################################
subroutine TIFF_Write_Byte_Into_Buffer(Byte)

use TIFF_private
IMPLICIT NONE
character(len=1),intent(IN)   :: Byte

! increment byte counter
 L=L+1
! is record full ?
 if (L>len(Record))then  ! yes it is, so write to file
  Rec_No = Rec_No + 1
  write (9,REC=Rec_No) Record
! reset entire record to zero
  Record(L:L)=char(0)
  L = 1
 end if
! add byte to record
 Record(L:L) = Byte
end subroutine TIFF_Write_Byte_Into_Buffer
!
! ###################################################################
! 
!  subroutine TIFF_Write_Word
!
!                                    created: ?/??/96 
!
!  Author: R.A. Vowels / Marc De Graef
!  
!  Description: write a 4-byte word into the buffer
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   ?/??/96 RAV 1.0 original
!   8/28/01 MDG 2.0 commented and change of variable names
! ###################################################################
subroutine TIFF_Write_Word(Word,Length)

use local

IMPLICIT NONE
integer(kind=irg),intent(IN)    :: Word,Length

integer(kind=irg) :: L_Word
integer(kind=irg) :: j
character(len=1)  :: Ch

 L_Word = Word
 do j=1,Length
  Ch = char(iand(L_Word,255))
  call TIFF_write_byte_into_buffer(Ch)
  L_Word = ishft(L_Word,-8)
 end do
end subroutine TIFF_Write_Word
!
! ###################################################################
! 
!  subroutine TIFF_Make_Tag
!
!                                    created: ?/??/96 
!
!  Author: R.A. Vowels / Marc De Graef
!  
!  Description: create a 12 byte Image File Directory Entry
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   ?/??/96 RAV 1.0 original
!   8/28/01 MDG 2.0 commented and change of variable names
! ###################################################################
subroutine TIFF_Make_Tag(Number,Tag_ID, Data_Type,Count,Offset)

use local

IMPLICIT NONE
integer(kind=irg),intent(IN)   :: Number, Tag_ID, Data_Type, Count, Offset

 call TIFF_Write_Word(Tag_ID,2)
 call TIFF_Write_Word(Data_Type,2)
 call TIFF_Write_Word(Count,4)
 call TIFF_Write_Word(Offset,4)
end subroutine TIFF_Make_Tag
!
! ###################################################################
! 
!  subroutine TIFF_Write_File
!
!                                    created: ?/??/96 
!
!  Author: R.A. Vowels / Marc De Graef
!  
!  Description: write the TIFF file to unit 9
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!   ?/??/96 RAV 1.0 original
!   8/28/01 MDG 2.0 commented and change of variable names + grayscale
! ###################################################################
subroutine TIFF_Write_File

use local
use TIFF_global
use TIFF_private

IMPLICIT NONE

integer(kind=irg)    :: I, Row, Col

 Rec_No = 0
 L = 0
! RECL is measured in units of words, not bytes !!!
! This may depend on the platform, and may need to be changed
 open(9,file=TIFF_filename,access="DIRECT",action="WRITE", FORM="UNFORMATTED", RECL=TIFFRecordLength)

! 8 byte header
 call TIFF_Write_Byte_Into_Buffer('I')      ! little endian header
 call TIFF_Write_Byte_Into_Buffer('I')
 call TIFF_Write_Word(42,2)                 ! version number
 call TIFF_Write_Word(8,4)                  ! address of IFD
 
! number of entries in IFD (2 bytes)
 call TIFF_Write_Byte_Into_Buffer(char(11)) 
 call TIFF_Write_Byte_Into_Buffer(char(0))

! image tags (11 of them, 12 bytes each)
 call TIFF_Make_Tag(2,256,3,1,TIFF_nx)               ! ImageWidth
 call TIFF_Make_Tag(3,257,3,1,TIFF_ny)               ! ImageLength
 call TIFF_Make_Tag(4,258,3,1,8)                     ! BitsPerSample
 call TIFF_Make_Tag(5,259,3,1,1)                     ! Compression
 call TIFF_Make_Tag(6,262,3,1,1)                     ! PhotometricInterpretation
 call TIFF_Make_Tag(7,273,4,TIFF_ny,256)             ! StripOffsets
 call TIFF_Make_Tag(9,278,4,1,1)                     ! RowsPerStrip
 call TIFF_Make_Tag(10,279,4,TIFF_ny,256+TIFF_ny*4)  ! StripByteCounts
 call TIFF_Make_Tag(11,282,5,1,146)                  ! Xresolution
 call TIFF_Make_Tag(12,283,5,1,154)                  ! Yresolution
 call TIFF_Make_Tag(14,296,3,1,2)                    ! ResolutionUnit

! end of IFD (4 bytes)
 call TIFF_Write_Word(0,4)                 

! extra values (X and Y resolution, 8 bytes each)
 call TIFF_Write_Word(300,4)               ! x-resolution
 call TIFF_Write_Word(1,4)
 call TIFF_Write_Word(300,4)               ! y-resolution
 call TIFF_Write_Word(1,4)

! pad with zeroes to fill first 256 bytes of file
 do i=L+1,256
  call TIFF_Write_Byte_Into_Buffer(char(0))
 end do

! write strip offsets
! offset = 256 + number of offset entries + number of count entries + stripnumber
 do Row=0,TIFF_ny-1
  call TIFF_Write_Word(256+TIFF_ny*8+Row*TIFF_nx,4)
 end do

! write stripcounts
 do Row=0,TIFF_ny-1
  call TIFF_Write_Word(TIFF_nx,4)
 end do

! write the actual image, one strip at a time
 do Row=0,TIFF_ny-1
  do Col=0,TIFF_nx-1
   call TIFF_Write_Byte_Into_Buffer(char(TIFF_Image(Col,Row)))
  end do
 end do
 
! make sure the last record is actually written to the file
 L = 10000
 call TIFF_Write_Byte_Into_Buffer(char(0))

! close and save file
 close(9,status="KEEP")

 return
end subroutine TIFF_Write_File

end Module TIFF_f90

