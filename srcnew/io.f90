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
! CTEMsoft2013:io.f90
!--------------------------------------------------------------------------
!
! MODULE: io
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief IO routines
!
!> @details  Several IO routines; the thinking is that this might be the only file that needs
!> to be rewritten if the package IO is modified.  In particular, in this version of the package,
!> all program input, with the exception of the mkxtal.f90 program and some other simple
!> programs, is done by means of namelist files.  This allows the entire package to be
!> controlled from inside the IDL environment or from the standard command line.
! 
!> @date 1/5/99   MDG 1.0 original
!> @date 5/19/01 MDG 2.0 f90 version
!> @date 03/19/13  MDG  3.0 major changes in IO routines (use of interface)
!> @date 05/16/13 MDG 3.1 added stdout as an option to run from IDL
!--------------------------------------------------------------------------

module io

use local

character(132)              	:: mess		!< 132 character string for message output

public

interface ReadValue
	module procedure ReadValueIntShort
	module procedure ReadValueIntLong
	module procedure ReadValueRealSingle
	module procedure ReadValueRealDouble
	module procedure ReadValueString
	module procedure ReadValueStringArray
end interface ReadValue

interface WriteValue
	module procedure WriteValueIntShort
	module procedure WriteValueIntLong
	module procedure WriteValueRealSingle
	module procedure WriteValueRealDouble
	module procedure WriteValueRealComplex
	module procedure WriteValueString
end interface WriteValue

interface Message
	module procedure Message
end interface

contains

!--------------------------------------------------------------------------
!
! SUBROUTINE: Message
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief dump a message to standard output
!
!> @details Simple routine to print a string on the standard output, with optional formatting
!>  instructions
! 
!> @param frm optional string formatting command
!
!> @date 1/5/99   MDG 1.0 original
!> @date 5/19/01 MDG 2.0 f90 version
!> @date 03/19/13 MDG 3.0 made argument optional and introduced default format
!--------------------------------------------------------------------------
subroutine Message(frm)

use local

character(*),OPTIONAL,INTENT(IN)  :: frm	!< optional formatting string

! default format or not ?
if (PRESENT(frm)) then
 write (stdout,fmt=frm) trim(mess)
else    ! default output format: a string with a linefeed before and after ...
 write (stdout,fmt="(/A/)") trim(mess)
end if 

! and reset the mess array
mess(:) = ''

end subroutine Message

! ###################################################################
! 
!  subroutine ReadValueString   
!
!                                    created: 10/13/98 {9:29:46 AM} 
!                                last update: 1/8/98 {9:40:05 AM} 
!  Author: Marc De Graef
!  
!  Description: read a string from standard input
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  03/19/13 MDG 1.0 new routine
! ###################################################################
subroutine ReadValueString(Qstring, in_string, frm)

use local

character(*),INTENT(IN)			:: Qstring, frm
character(*),INTENT(OUT)			:: in_string

! send Qstring to the output
mess = Qstring
call Message("(' ',A,' ',$)")

read (*,fmt=frm) in_string

end subroutine ReadValueString

! ###################################################################
! 
!  subroutine ReadValueStringArray   
!
!                                    created: 10/13/98 {9:29:46 AM} 
!                                last update: 1/8/98 {9:40:05 AM} 
!  Author: Marc De Graef
!  
!  Description: read a series of strings from standard input
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  03/19/13 MDG 1.0 new routine
! ###################################################################
subroutine ReadValueStringArray(Qstring, in_string, num, frm)

use local

character(*),INTENT(IN)			:: Qstring, frm 
character(1),INTENT(OUT)			:: in_string(num)
integer(kind=irg),INTENT(IN)			:: num

! send Qstring to the output
mess = Qstring
call Message("(' ',A,' ',$)")

read (*,fmt=frm) in_string

end subroutine ReadValueStringArray
! ###################################################################
! 
!  subroutine ReadValueIntShort    
!
!                                    created: 10/13/98 {9:29:46 AM} 
!                                last update: 1/8/98 {9:40:05 AM} 
!  Author: Marc De Graef
!  
!  Description: read one or more short integers from standard input
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  03/19/13 MDG 1.0 new routine
! ###################################################################
subroutine ReadValueIntShort(Qstring, in_int, num)

use local

character(*), INTENT(IN)			:: Qstring
integer(kind=ish),INTENT(OUT)			:: in_int(*)
integer(kind=irg),INTENT(IN),OPTIONAL		:: num

! send Qstring to the output
mess = Qstring
call Message("(' ',A,' ',$)")

! one or more than one values expected ?
if (PRESENT(num)) then
  read (*,*) (in_int(i),i=1,num)
else
  read (*,*) in_int(1)
end if
  
end subroutine ReadValueIntShort

! ###################################################################
! 
!  subroutine ReadValueInt    
!
!                                    created: 10/13/98 {9:29:46 AM} 
!                                last update: 1/8/98 {9:40:05 AM} 
!  Author: Marc De Graef
!  
!  Description: read one or more regular integers from standard input
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  03/19/13 MDG 1.0 new routine
! ###################################################################
subroutine ReadValueIntLong(Qstring, in_int, num)

use local

character(*), INTENT(IN)			:: Qstring
integer(kind=irg),INTENT(OUT)			:: in_int(*)
integer(kind=irg),INTENT(IN),OPTIONAL		:: num

! send Qstring to the output
mess = Qstring
call Message("(' ',A,' ',$)")

! one or more than one values expected ?
if (PRESENT(num)) then
  read (*,*) (in_int(i),i=1,num)
else
  read (*,*) in_int(1)
end if

end subroutine ReadValueIntLong

! ###################################################################
! 
!  subroutine ReadValueRealSingle    
!
!                                    created: 10/13/98 {9:29:46 AM} 
!                                last update: 1/8/98 {9:40:05 AM} 
!  Author: Marc De Graef
!  
!  Description: read one or more singe precision reals from standard input
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  03/19/13 MDG 1.0 new routine
! ###################################################################
subroutine ReadValueRealSingle(Qstring, in_real, num)

use local

character(*), INTENT(IN)			:: Qstring
real(kind=sgl),INTENT(OUT)			:: in_real(*)
integer(kind=irg),INTENT(IN),OPTIONAL		:: num

! send Qstring to the output
mess = Qstring
call Message("(' ',A,' ',$)")

! one or more than one values expected ?
if (PRESENT(num)) then
  read (*,*) (in_real(i),i=1,num)
else
  read (*,*) in_real(1)
end if

! the calling routine must now use io_real to get to the value(s)
end subroutine ReadValueRealSingle


! ###################################################################
! 
!  subroutine ReadValueRealDouble
!
!                                    created: 10/13/98 {9:29:46 AM} 
!                                last update: 1/8/98 {9:40:05 AM} 
!  Author: Marc De Graef
!  
!  Description: read one or more double precision reals from standard input
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  03/19/13 MDG 1.0 new routine
! ###################################################################
subroutine ReadValueRealDouble(Qstring, in_real, num)

use local

character(*), INTENT(IN)			:: Qstring
real(kind=dbl),INTENT(OUT)			:: in_real(*)
integer(kind=irg),INTENT(IN),OPTIONAL		:: num

! send Qstring to the output
mess = Qstring
call Message("(' ',A,' ',$)")

! one or more than one values expected ?
if (PRESENT(num)) then
  read (*,*) (in_real(i),i=1,num)
else
  read (*,*) in_real(1)
end if

end subroutine ReadValueRealDouble


! ###################################################################
! 
!  subroutine WriteValueString    
!
!                                    created: 10/13/98 {9:29:46 AM} 
!                                
!  Author: Marc De Graef
!  
!  Description: write  a string to standard output
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  03/19/13 MDG 1.0 new routine
! ###################################################################
subroutine WriteValueString(Qstring, out_string, frm)

use local

character(*), INTENT(IN)			:: Qstring, out_string
character(*),INTENT(IN),OPTIONAL		:: frm


! send Qstring to the output only if it is non-zero length
if (len(Qstring).ne.0) then
  mess = Qstring
  call Message("(' ',A,' ',$)")
end if

mess = trim(out_string)
if (PRESENT(frm)) then 
  call Message(frm)
else
 call Message("(A)")
end if

end subroutine WriteValueString

! ###################################################################
! 
!  subroutine WriteValueIntShort    
!
!                                    created: 10/13/98 {9:29:46 AM} 
!                                
!  Author: Marc De Graef
!  
!  Description: write one or more short integers to standard input
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  03/19/13 MDG 1.0 new routine
! ###################################################################
subroutine WriteValueIntShort(Qstring, out_int, num, frm)

use local

character(*), INTENT(IN)			:: Qstring
character(*),INTENT(IN),OPTIONAL		:: frm
integer(kind=ish),INTENT(IN)			:: out_int(*)
integer(kind=irg),INTENT(IN),OPTIONAL		:: num

! send Qstring to the output
if (len(Qstring).ne.0) then
  mess = Qstring
  call Message("(' ',A,' ',$)")
end if

! one or more than one values expected ?
if (PRESENT(num)) then
 if (PRESENT(frm)) then
  write (stdout,fmt=frm) (out_int(i),i=1,num)
 else
  write (stdout,*) (out_int(i),i=1,num)
 end if
else
 if (PRESENT(frm)) then
  write (stdout,fmt=frm) out_int(1)
 else
  write (stdout,*) out_int(1)
 end if
end if

end subroutine WriteValueIntShort

! ###################################################################
! 
!  subroutine WriteValueIntLong    
!
!                                    created: 10/13/98 {9:29:46 AM} 
!                                
!  Author: Marc De Graef
!  
!  Description: write one or more long integers to standard input
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  03/19/13 MDG 1.0 new routine
! ###################################################################
subroutine WriteValueIntLong(Qstring, out_int, num, frm)

use local

character(*), INTENT(IN)			:: Qstring
character(*),INTENT(IN),OPTIONAL		:: frm
integer(kind=irg),INTENT(IN)			:: out_int(*)
integer(kind=irg),INTENT(IN),OPTIONAL		:: num

! send Qstring to the output
if (len(Qstring).ne.0) then
  mess = Qstring
  call Message("(' ',A,' ',$)")
end if

! one or more than one values expected ?
if (PRESENT(num)) then
 if (PRESENT(frm)) then
  write (stdout,fmt=frm) (out_int(i),i=1,num)
 else
  write (stdout,*) (out_int(i),i=1,num)
 end if
else
 if (PRESENT(frm)) then
  write (stdout,fmt=frm) out_int(1)
 else
  write (stdout,*) out_int(1)
 end if
end if

end subroutine WriteValueIntLong


! ###################################################################
! 
!  subroutine WriteValueRealSingle    
!
!                                    created: 10/13/98 {9:29:46 AM} 
!                                
!  Author: Marc De Graef
!  
!  Description: write one or more single precision reals  to standard input
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  03/19/13 MDG 1.0 new routine
! ###################################################################
subroutine WriteValueRealSingle(Qstring, out_real, num, frm)

use local

character(*), INTENT(IN)			:: Qstring
real(kind=sgl),INTENT(IN)			:: out_real(*)
integer(kind=irg),INTENT(IN),OPTIONAL		:: num
character(*),INTENT(IN),OPTIONAL		:: frm

! send Qstring to the output
if (len(Qstring).ne.0) then
  mess = Qstring
  call Message("(' ',A,' ',$)")
end if

! one or more than one values expected ?
if (PRESENT(num)) then
 if (PRESENT(frm)) then
  write (stdout,fmt=frm) (out_real(i),i=1,num)
 else
  write (stdout,*) (out_real(i),i=1,num)
 end if
else
 if (PRESENT(frm)) then
  write (stdout,fmt=frm) out_real(1)
 else
  write (stdout,*) out_real(1)
 end if
end if

end subroutine WriteValueRealSingle

! ###################################################################
! 
!  subroutine WriteValueRealDouble   
!
!                                    created: 10/13/98 {9:29:46 AM} 
!                                
!  Author: Marc De Graef
!  
!  Description: write one or more double precision reals  to standard input
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  03/19/13 MDG 1.0 new routine
! ###################################################################
subroutine WriteValueRealDouble(Qstring, out_real, num,frm)

use local

character(*), INTENT(IN)			:: Qstring
character(*),INTENT(IN),OPTIONAL		:: frm
real(kind=dbl),INTENT(IN)			:: out_real(*)
integer(kind=irg),INTENT(IN),OPTIONAL		:: num

! send Qstring to the output
if (len(Qstring).ne.0) then
  mess = Qstring
  call Message("(' ',A,' ',$)")
end if

! one or more than one values expected ?
if (PRESENT(num)) then
 if (PRESENT(frm)) then
  write (stdout,fmt=frm) (out_real(i),i=1,num)
 else
  write (stdout,*) (out_real(i),i=1,num)
 end if
else
 if (PRESENT(frm)) then
  write (stdout,fmt=frm) out_real(1)
 else
  write (stdout,*) out_real(1)
 end if
end if

end subroutine WriteValueRealDouble

! ###################################################################
! 
!  subroutine WriteValueRealComplex  
!
!                                    created: 10/13/98 {9:29:46 AM} 
!                                
!  Author: Marc De Graef
!  
!  Description: write one or more complex reals  to standard input
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  03/19/13 MDG 1.0 new routine
! ###################################################################
subroutine WriteValueRealComplex(Qstring, out_cmplx, num,frm)

use local

character(*), INTENT(IN)			:: Qstring
character(*),INTENT(IN),OPTIONAL		:: frm
complex(kind=sgl),INTENT(IN)			:: out_cmplx(*)
integer(kind=irg),INTENT(IN),OPTIONAL		:: num

! send Qstring to the output
if (len(Qstring).ne.0) then
  mess = Qstring
  call Message("(' ',A,' ',$)")
end if

! one or more than one values expected ?
if (PRESENT(num)) then
 if (PRESENT(frm)) then
  write (stdout,fmt=frm) (out_cmplx(i),i=1,num)
 else
  write (stdout,*) (out_cmplx(i),i=1,num)
 end if
else
 if (PRESENT(frm)) then
  write (stdout,fmt=frm) out_cmplx(1)
 else
  write (stdout,*) out_cmplx(1)
 end if
end if

end subroutine WriteValueRealComplex

end module io
