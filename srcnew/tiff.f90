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
! CTEMsoft2013:tiff.f90
!--------------------------------------------------------------------------
!
! MODULE: tiff
!
!> @author R.A. Vowels / Marc De Graef, Carnegie Mellon University
!
!> @brief TIFF output routines
!
!> @details  based on the programs in the book
!> "Introduction to Fortran 90/95, Algorithms, and Structured Programming", by R.A. Vowels (1996, page 489-94)\n
!> 
!> Here is an example showing how to use the TIFF routines.\n
!>
!>program tiff\n
!>
!>use TIFF_global\n
!>use TIFF_f90\n
!>\n
!>IMPLICIT NONE\n
!>\n
!>integer   :: i,j\n
!>\n
!>! declare TIFF variables in TIFF_global\n
!> TIFF_nx = 295\n
!> TIFF_ny = 150\n
!> TIFF_filename = "test.tiff"\n
!>! allocate memory for image\n
!> allocate(TIFF_image(0:TIFF_nx-1,0:TIFF_ny-1))\n
!>! fill the image with whatever data you have (between 0 and 255)\n
!> do i=0,TIFF_nx-1\n
!>  do j=0,TIFF_ny-1\n
!>   TIFF_image(i,j) = mod(i*i+j*j,255) \n
!>  end do\n
!> end do\n
!>! create the file\n
!> call TIFF_Write_File\n
!>\n
!>end program\n
!>\n
!> The RecordLength parameter may be platform dependent;  you should
!> check whether or not this package works;  if it does not work, then
!> this is most likely due to an incorrect RecordLength parameter.  
!>
!> Also, the string 'Rekord' is named with a k instead of a c because 'Record' is 
!> an f90 reserved word.
!
!> @date    ?/??/96 RAV 1.0 original
!> @date    8/28/01 MDG 2.0 adapted for grayscale images
!> @date   11/28/01 MDG 2.1 added kind support
!--------------------------------------------------------------------------
module TIFF_f90

use local

IMPLICIT NONE

PRIVATE  :: Rekord, Rec_No, L, TIFFRecordLength, TIFF_Write_Byte_Into_Buffer, TIFF_Write_Word, TIFF_Make_Tag

 character(len=256)                   :: Rekord
 integer(kind=irg)                      :: Rec_No=0, L=0
 integer(kind=irg),parameter     :: TIFFRecordLength = 256

PUBLIC :: TIFF_nx, TIFF_ny, TIFF_Image, TIFF_filename, TIFF_Write_File

 integer(kind=irg)                     :: TIFF_nx,TIFF_ny
 integer(kind=irg),allocatable   :: TIFF_Image(:,:)
 character(100)                         :: TIFF_filename

contains

!--------------------------------------------------------------------------
!
! SUBROUTINE:TIFF_Write_Byte_Into_Buffer
!
!> @author R.A. Vowels /Marc De Graef, Carnegie Mellon University
!
!> @brief write a single byte into a buffer and dump the
!               buffer to file if full
!
!> @param Bite single byte to be added to the buffer

!> @note Renamed 'Byte' to 'Bite' since Byte is an f90 reserved word  
! 
!> @date    ?/??/96 RAV 1.0 original
!> @date    8/28/01 MDG 2.0 commented and change of variable names
!--------------------------------------------------------------------------
subroutine TIFF_Write_Byte_Into_Buffer(Bite)

IMPLICIT NONE

character(len=1),intent(IN)   :: Bite		!< byte variable

! increment byte counter
 L=L+1
 
! is record full ?
 if (L>len(Rekord)) then  ! yes it is, so write to file
  Rec_No = Rec_No + 1
  write (9,REC=Rec_No) Rekord
! reset entire record to zero
  Rekord(1:len(Rekord))=char(0)
  L = 1
 end if

! add byte to record
 Rekord(L:L) = Bite

end subroutine TIFF_Write_Byte_Into_Buffer

!--------------------------------------------------------------------------
!
! SUBROUTINE:TIFF_Write_Word
!
!> @author R.A. Vowels /Marc De Graef, Carnegie Mellon University
!
!> @brief write a 4-byte word into the buffer
!
!> @param Word 4-byte word
!> @param Length length parameter
!
!> @note Renamed 'Byte' to 'Bite' since Byte is an f90 reserved word  
! 
!> @date    ?/??/96 RAV 1.0 original
!> @date    8/28/01 MDG 2.0 commented and change of variable names
!--------------------------------------------------------------------------
subroutine TIFF_Write_Word(Word,Length)

use local

IMPLICIT NONE

integer(kind=irg),intent(IN)    	:: Word		!< 4-byte word
integer(kind=irg),intent(IN)    	:: Length		!< length parameter

integer(kind=irg) 			:: L_Word
integer(kind=irg) 			:: j
character(len=1)  			:: Ch

 L_Word = Word
 do j=1,Length
  Ch = char(iand(L_Word,255))
  call TIFF_Write_Byte_Into_Buffer(Ch)
  L_Word = ishft(L_Word,-8)
 end do
 
end subroutine TIFF_Write_Word

!--------------------------------------------------------------------------
!
! SUBROUTINE:TIFF_Make_Tag
!
!> @author R.A. Vowels /Marc De Graef, Carnegie Mellon University
!
!> @brief create a 12 byte Image File Directory Entry
!
!> @param Numbre Tag number (not used)
!> @param Tag_ID Tag identifier
!> @param Data_Type data type identifier
!> @param Cnt counter
!> @param Offset offset parameter
! 
!> @date    ?/??/96 RAV 1.0 original
!> @date    8/28/01 MDG 2.0 commented and change of variable names
!--------------------------------------------------------------------------
subroutine TIFF_Make_Tag(Numbre,Tag_ID, Data_Type,Cnt,Offset)

use local

IMPLICIT NONE

integer(kind=irg),intent(IN)   :: Numbre		!< Tag number (only used for clarity)
integer(kind=irg),intent(IN)   :: Tag_ID		!< Tag identifier
integer(kind=irg),intent(IN)   :: Data_Type	!< Tag data type
integer(kind=irg),intent(IN)   :: Cnt			!< counter
integer(kind=irg),intent(IN)   :: Offset		!< Offset parameter

 call TIFF_Write_Word(Tag_ID,2)
 call TIFF_Write_Word(Data_Type,2)
 call TIFF_Write_Word(Cnt,4)
 call TIFF_Write_Word(Offset,4)
 
end subroutine TIFF_Make_Tag

!--------------------------------------------------------------------------
!
! SUBROUTINE:TIFF_Write_File
!
!> @author R.A. Vowels /Marc De Graef, Carnegie Mellon University
!
!> @brief write the TIFF file to unit 9
! 
!> @date    ?/??/96 RAV 1.0 original
!> @date    8/28/01 MDG 2.0 commented and change of variable names
!> @date   01/08/10  MDG 3.0 conversion to Intel byte ordering and addition of all required tags
!
!> @note  <b>WARNING</b>: Do not modify this routine at all unless you REALLY understand the TIFF
!>  format!  The file format specification for TIFF can be found at the following URL\n
!>
!>        http://partners.adobe.com/public/developer/en/tiff/TIFF6.pdf
!
!--------------------------------------------------------------------------
subroutine TIFF_Write_File

use local

IMPLICIT NONE

integer(kind=irg)    :: I, Row, Col

 Rec_No = 0
 L = 0

! RECL is measured in units of words, not bytes !!!
! This may depend on the platform, and may need to be changed
 open(9,file=TIFF_filename,access="DIRECT",action="WRITE", FORM="UNFORMATTED", RECL=TIFFRecordLength)
 
! 10 byte header
 call TIFF_Write_Byte_Into_Buffer('I')      ! little endian header
 call TIFF_Write_Byte_Into_Buffer('I')
 call TIFF_Write_Word(42,2)                  ! version number
 call TIFF_Write_Word(10,4)                  ! address of IFD
 call TIFF_Write_Word(0,2)                    ! two empty bytes

! number of entries in IFD (2 bytes)
 call TIFF_Write_Byte_Into_Buffer(char(14)) 
 call TIFF_Write_Byte_Into_Buffer(char(0))

! image tags (14 of them, 12 bytes each)
 call TIFF_Make_Tag(1,254,4,1,0)                     ! new subfile type
 call TIFF_Make_Tag(2,256,3,1,TIFF_nx)           ! ImageWidth
 call TIFF_Make_Tag(3,257,3,1,TIFF_ny)           ! ImageLength
 call TIFF_Make_Tag(4,258,3,1,8)                     ! BitsPerSample
 call TIFF_Make_Tag(5,259,3,1,1)                     ! Compression
 call TIFF_Make_Tag(6,262,3,1,1)                     ! PhotometricInterpretation
 call TIFF_Make_Tag(7,273,4,TIFF_ny,256)       ! StripOffsets
 call TIFF_Make_Tag(8,277,3,1,1)                     ! Samples per pixel
 call TIFF_Make_Tag(9,278,3,1,1)                     ! RowsPerStrip
 call TIFF_Make_Tag(10,279,4,TIFF_ny,256+TIFF_ny*4)  ! StripByteCounts
 call TIFF_Make_Tag(11,282,5,1,184)               ! Xresolution
 call TIFF_Make_Tag(12,283,5,1,192)               ! Yresolution
 call TIFF_Make_Tag(13,284,3,1,1)                   ! planar configuration
 call TIFF_Make_Tag(14,296,3,1,1)                   ! ResolutionUnit

! end of IFD (4 bytes)
 call TIFF_Write_Word(0,4)                 
 
! extra values (X and Y resolution, 8 bytes each; default for screen resolution)
 call TIFF_Write_Word(72,4)               ! x-resolution
 call TIFF_Write_Word(1,4)
 call TIFF_Write_Word(72,4)               ! y-resolution
 call TIFF_Write_Word(1,4)

! pad with zeroes to fill first 256 bytes of file
 do I=L+1,256
  call TIFF_Write_Byte_Into_Buffer(char(0))
 end do

! write strip offsets
! offset = 256 + number of offset entries + number of count entries + stripnumber
 do Row=0,TIFF_ny-1
  call TIFF_Write_Word(256+TIFF_ny*8+Row*TIFF_nx,4)
 end do

! write stripcounts (number of bytes in each strip)
 do Row=0,TIFF_ny-1
  call TIFF_Write_Word(TIFF_nx,4)
 end do

! write the actual image, one strip at a time (start with the top row)
 do Row=TIFF_ny-1,0,-1
  do Col=0,TIFF_nx-1
   call TIFF_Write_Byte_Into_Buffer(char(TIFF_Image(Col,Row)))
  end do
 end do
 
! make sure the last record is actually written to the file
  L=len(Rekord)
  call TIFF_Write_Byte_Into_Buffer(char(0))
 
! close and save file
 close(9,status="KEEP")

end subroutine TIFF_Write_File

end module TIFF_f90