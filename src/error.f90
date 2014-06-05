! ###################################################################
! Copyright (c) 2013-2014, Marc De Graef/Carnegie Mellon University
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
! CTEMsoft2013:error.f90
!--------------------------------------------------------------------------
!
! MODULE: error
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief error handling routines
!
!> @details  Just a simple routine that reports errors and aborts the program
!
!> @date 01/05/99 MDG 1.0 original
!> @date 05/19/01 MDG 2.0 f90 version
!> @date 11/27/01 MDG 2.1 added kind support
!> @date 03/19/13 MDG 3.0 minor changes
!> @date 06/05/14 MDG 4.0 minor changes
!--------------------------------------------------------------------------

module error

contains

!--------------------------------------------------------------------------
!
! SUBROUTINE: FatalError
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Write error message and abort program
!
!> @param var1 routine name string
!> @param var2 explanation string
!> @param stdout optional output unit identifier
!
!> @date   10/13/98 MDG 1.0 original
!> @date    5/19/01 MDG 2.0 f90 version
!> @date   11/27/01 MDG 2.1 added kind support
!> @date   03/19/13 MDG 3.0 updated for new io routines
!> @date   06/05/14 MDG 4.0 added stdout argument
! ###################################################################
subroutine FatalError(var1,var2,stdout)

use io

IMPLICIT NONE

character(*), INTENT(IN)  			:: var1		!< first part of error message (routine name)
character(*), INTENT(IN)  			:: var2		!< second part of error message (brief explanation)
integer(kind=irg)m INTENT(IN),OPTIONAL		:: stdout 	!< output unit identifier

integer(kind=irg)				:: std 
character(*)					:: mess

 std = 6
 if (PRESENT(stdout)) std=stdout

 mess = ' ----> Fatal error in routine '//var1//': '//var2
 call Message(mess, frm='(//A//)', stdout=std) 
 stop '  Progam ended abnormally'

end subroutine


end module

