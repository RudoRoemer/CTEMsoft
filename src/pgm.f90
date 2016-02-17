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
! EMsoft:pgm.f90
!--------------------------------------------------------------------------
!
! MODULE: pgm
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief routine to write an image to a pgm file
!
! 
!> @date 1/5/99   MDG 1.0 original
!> @date 6/9/14   MDG 2.0 removed global variables
!--------------------------------------------------------------------------
module pgm

use local
use files
use io

IMPLICIT NONE

contains

subroutine PGM_Write_File(fname, nx, ny, image) 

character(fnlen),INTENT(IN) 		:: fname
integer(kind=irg),INTENT(IN)   	:: image(nx,ny)
integer(kind=irg),INTENT(IN)  		:: nx,ny

integer(kind=irg)			:: j, stl

stl = len(trim(fname))
call Message('Creating PGM file : '//fname(1:stl), frm = "(A)")

open(unit=dataunit,file=trim(EMsoft_toNativePath(fname)),status='unknown',action='write',form = 'formatted')
write (unit=dataunit,fmt="(A3)") 'P2 '
write (unit=dataunit,fmt="(A)") '# PGM file generated by EMsoft'
write (unit=dataunit,fmt="(I5,I5)") nx,ny
write (unit=dataunit,fmt="(I5)") maxval(image)
do j=ny,1,-1
  write (unit=dataunit,fmt="(I3)") image(1:nx,j)
end do
close(unit=dataunit,status='keep')

end subroutine PGM_Write_File

end module pgm
