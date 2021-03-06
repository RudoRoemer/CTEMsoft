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
! CTEMsoft2013:local.f90
!--------------------------------------------------------------------------
!
! MODULE: local
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief 
!> local definitions of single and double precision, general constants and variables
!
!> @details  
!> defines the kind-parameters for short and long integers, and single/double 
!> precision reals.
! 
!> @todo remove postscript options from this entire package
!
!> @date 1/8/99   MDG 1.0 original
!> @date 5/6/01   MDG 2.0 f90
!> @date 11/27/01  MDG 2.1 added sgl and dbl kinds
!> @date 12/08/01  MDG 2.2 added CTEMsoft subroutine
!> @date 03/19/13  MDG 3.0 rewrite of entire package
!> @date 05/16/13 MDG 3.1 added stdout
!--------------------------------------------------------------------------

module local

! for OpenMP implementations ... 
use omp_lib

!> @note This module must be "use"d by every program, subroutine, and function!

! The entire CTEMsoft package should be processor independent.  This can
! be accomplished by the use of the "kind" parameters.

! Define the "kind" parameters for single and double precision reals, 
!> single precision real kind parameter
  integer,parameter        		:: sgl = SELECTED_REAL_KIND(p=6,r=37)	
!> double precision real kind parameter
  integer,parameter        		:: dbl = SELECTED_REAL_KIND(p=13,r=200)	

! Define the "kind" parameters for short and regular integers,
!> short integer kind parameter 
  integer,parameter        		:: ish = SELECTED_INT_KIND(3) 
!> long integer kind parameter 
  integer,parameter        		:: irg = SELECTED_INT_KIND(9)

!> source code version number
  character(8), parameter  	:: scversion="2.0/2013"
!> source code author name
  character(13), parameter 	:: username="Marc De Graef"	
!> source code author location
  character(26), parameter 	:: userlocn="Carnegie Mellon University"

!> standard string length for filenames
  integer(kind=irg),parameter	:: fnlen=132
  
!> program name string
  character(fnlen)            		:: progname="default"
!> one line program descriptor
  character(fnlen)            		:: progdesc="default"

!> pathname to the namelist template files
  character(fnlen)			:: templatepathname   !='~/.CTEMsoft2013/templatefolder'
  character(fnlen)			:: resourcepathname   !='~/.CTEMsoft2013/resources'
  integer(kind=irg), parameter		:: maxnumtemplates = 256
  character(17)				:: templatecodefilename = 'templatecodes.txt'
  
! input/output information
! psunit    = Postscript output unit number
! dataunit  = Data unit number (for *.xtal files and other in/output)
!> reserved IO unit identifiers for postscript (20) and data (21-23)
  integer(kind=irg), parameter	:: psunit=20, dataunit=21, dataunit2=22, dataunit3=23

! where should standard output go ?  To the terminal (stdout=6) or elsewhere (stdout>10)
! this is not a parameter, but can be changed by each program 
  integer(kind=irg)			:: stdout = 6
  
! parameters governing the size of varous arrays
!> Maximum number of positions in asymmetric unit
  integer, parameter       		:: maxpasym=250		
!>Maximum number of defects of any given type
  integer, parameter       		:: maxdefects=250

! strucdef is a logical variable to determine whether or not a structure has been loaded;
! hexset determines whether or not 4-index hexagonal indices should be used.
!> logical variable to determine whether or not a crystal structuter has been loaded
  logical                  			:: strucdef
!> logical to determine whether to use 3(FALSE) or 4(TRUE) index notation
  logical                  			:: hexset
 
! to be eliminated !!! 
!> where is the postscript viewer on this system ? 
  character(18),parameter  	:: psviewer="/usr/local/bin/gv "
contains


!--------------------------------------------------------------------------
! CTEMsoft2013:local:CTEMsoft.f90
!--------------------------------------------------------------------------
!
! SUBROUTINE: CTEMsoft
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief prints a copyright statement and the program name
!
!> @details prints a copyright statement as well as where the user can find the license information 
!> This is then followed by the program name, a one-line description, and a time stamp.
! 
!> @date  12/08/01 MDG 1.0 original
!> @date  03/19/13 MDG 2.0 minor modifications
!> @date  05/16/13 MDG 2.1 added timestamp and stdout
!--------------------------------------------------------------------------

subroutine CTEMsoft

 write (stdout,"(//1x,'CTEMsoft version ',A8,', Copyright (C) 2001-2013 Marc De Graef/CMU')") scversion
 write (stdout,"(1x,'CTEMsoft comes with ABSOLUTELY NO WARRANTY.')")
 write (stdout,"(1x,'This is free software, and you are welcome to redistribute it')")
 write (stdout,"(1x,'under certain conditions; see Copyright.txt file for details.'//)")

 write (stdout,"(1x,'Program name : ',A)") trim(progname)
 write (stdout,"(1x,A//)") trim(progdesc)

 call timestamp()
 write (stdout,"(1x,//)")
 
end subroutine CTEMsoft


!--------------------------------------------------------------------------
!
! subroutine: timestamp
!
!> @author John Burkardt
!
!> @brief prints the current YMDHMS date as a time stamp.
!
!> @note  This code is distributed under the GNU LGPL license. 
!
!> @date 05/31/01  JB original
!> @date 05/01/13 MDG changed 'm' to 'mo' for month variable, and some other changes
!--------------------------------------------------------------------------
subroutine timestamp ( )

  implicit none

  character ( len = 8 ) ampm
  integer ( kind = irg ) d
  character ( len = 8 ) date
  integer ( kind = irg ) h
  integer ( kind = irg ) mo
  integer ( kind = irg ) mm
  character ( len = 3 ), parameter, dimension(12) :: month = (/ &
    'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec' /)
  integer ( kind = irg ) n
  integer ( kind = irg ) s
  character ( len = 10 ) time
  integer ( kind = irg ) values(8)
  integer ( kind = irg ) y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

  y = values(1)
  mo = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( stdout, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    month(mo), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end subroutine timestamp




end module local
