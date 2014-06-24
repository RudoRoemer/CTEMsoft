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
! CTEMsoft2013:CTEMmkxtal.f90
!--------------------------------------------------------------------------
!
! PROGRAM: mkxtal 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief create a crystal structure file (very simple program)
!
!> @date   11/23/13 MDG 1.0 original
!--------------------------------------------------------------------------
program CTEMmkxtal

use local
use typedefs
use io
use crystal
use symmetry
use files

IMPLICIT NONE

type(unitcell), pointer         :: cell
character(fnlen)                :: progname, progdesc, fname

 progname = 'CTEMmkxtal.f90'
 progdesc = 'Create a crystal structure file'
 call CTEMsoft(progname, progdesc)

 allocate(cell)
 
 cell%SYM_SGset=0
 call GetLatParm(cell)
 call GetSpaceGroup(cell)
 call GetAsymPos(cell)
 call ReadValue('Enter output file name (*.xtal) ', fname)
 cell%fname = fname
 call SaveData(cell)



end program CTEMmkxtal
       



