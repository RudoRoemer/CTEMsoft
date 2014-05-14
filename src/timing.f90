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
! CTEMsoft2013:timing.f90
!--------------------------------------------------------------------------
!
! MODULE: timing
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Provides routines to compute the displacement vector for various dislocations.
! 
!> @date   11/19/01 MDG 1.0 original
!> @date   06/04/13 MDG 2.0 rewrite
!--------------------------------------------------------------------------
module timing

use local
use io

IMPLICIT NONE

real(kind=sgl)     	:: TIME_t_count,TIME_unit_count,TIME_interval,TIME_fraction
integer(kind=irg)   	:: TIME_newcount,TIME_count_rate,TIME_count_max,TIME_count,TIME_old,TIME_loops


contains

!--------------------------------------------------------------------------
!
! SUBROUTINE: Time_reset
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reset time recording
!
!> @date   06/04/01 MDG 1.0 original
!> @date   06/04/13 MDG 2.0 rewrite
!--------------------------------------------------------------------------
subroutine Time_reset

TIME_t_count = 0.0
TIME_unit_count = 0.0
TIME_count = 0
TIME_newcount = 0
TIME_count_rate = 0
TIME_count_max = HUGE(0)
TIME_old = 0
TIME_loops = 0

end subroutine Time_reset

!--------------------------------------------------------------------------
!
! SUBROUTINE: Time_report
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief report time recording
!
!> @date   06/04/01 MDG 1.0 original
!> @date   06/04/13 MDG 2.0 rewrite
!--------------------------------------------------------------------------
subroutine Time_report(interval)

use local
use io

IMPLICIT NONE

real(kind=sgl),intent(IN)   :: interval

  TIME_interval = interval
  TIME_fraction = TIME_interval
  mess = 'Starting computation'; call Message("(/A)")

end subroutine Time_report

!--------------------------------------------------------------------------
!
! SUBROUTINE: Time_start
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief start time recording
!
!> @date   06/04/01 MDG 1.0 original
!> @date   06/04/13 MDG 2.0 rewrite
!--------------------------------------------------------------------------
subroutine Time_start

! start the timing of the computation
 call Time_reset
 call system_clock(TIME_count,TIME_count_rate,TIME_count_max)

end subroutine Time_start

!--------------------------------------------------------------------------
!
! SUBROUTINE: Time_estimate
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief estimare remaining time
!
!> @param numk number of idividual computations
!
!> @date   06/04/01 MDG 1.0 original
!> @date   06/04/13 MDG 2.0 rewrite
!--------------------------------------------------------------------------
subroutine Time_estimate(numk)

use local
use io

IMPLICIT NONE

integer(kind=irg),intent(IN)     	:: numk

integer(kind=irg)      		:: TIME_nc
real(kind=sgl)				:: io_real(1)

! get the current time
 call system_clock(TIME_nc,TIME_count_rate,TIME_count_max)
 TIME_newcount = TIME_nc
 TIME_t_count = float(TIME_newcount-TIME_count)/float(TIME_count_rate)
 TIME_unit_count = TIME_t_count
 io_real(1) = TIME_unit_count
 call WriteValue(' Time for first computation step [s, typically overestimate] :', io_real, 1, "(F10.3)")
 mess = '  Anticipated total computation time :'; call Message("(A$)")
 call PrintTime(TIME_unit_count*float(numk))
 
end subroutine Time_estimate

!--------------------------------------------------------------------------
!
! SUBROUTINE: Time_estimate
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief estimate remaining time
!
!> @param ik current computation
!> @param numk number of idividual computations
!
!> @date   06/04/01 MDG 1.0 original
!> @date   06/04/13 MDG 2.0 rewrite
!--------------------------------------------------------------------------
subroutine Time_remaining(ik,numk)

use local
use io

IMPLICIT NONE

integer(kind=irg),intent(IN)   	:: ik,numk

integer(kind=irg)    			:: TIME_nc, io_int(1)
real(kind=sgl)				:: io_real(1)

 TIME_fraction = TIME_fraction + TIME_interval

! get the current time
 call system_clock(TIME_nc,TIME_count_rate,TIME_count_max)

! correct for the resetting of TIME_nc when TIME_count_max is reached
 if (TIME_nc.lt.TIME_newcount) then     ! we've looped through the entire cycle
   TIME_loops = TIME_loops+1
   TIME_count = 0
 end if 
 TIME_newcount = TIME_nc

! and print it
 TIME_t_count = (float(TIME_loops)*float(TIME_count_max)+float(TIME_newcount-TIME_count))/float(TIME_count_rate)

! reset the time per unit
 TIME_unit_count = TIME_t_count/float(ik)

! print estimated remaining time
 io_int(1) = nint(100.0*TIME_t_count/(TIME_t_count+TIME_unit_count*(float(numk)-float(ik))))
 call WriteValue (' ',io_int, 1, "(1x,I3,' % completed; '$)") 
 io_real(1) = TIME_t_count
 call WriteValue(' Total computation time [s] ', io_real, 1, "(F10.3$)")
 mess = ';  Estimated remaining time : '; call Message("(A$)")
 call PrintTime(TIME_unit_count*(float(numk)-float(ik)))
!
end subroutine Time_remaining

!--------------------------------------------------------------------------
!
! SUBROUTINE: PrintTime
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief print  time
!
!> @param tm time variable
!
!> @date   06/04/01 MDG 1.0 original
!> @date   06/04/13 MDG 2.0 rewrite
!--------------------------------------------------------------------------
subroutine PrintTime(tm)

use local
use io

IMPLICIT NONE

real(kind=sgl),INTENT(IN)		:: tm

integer(kind=irg)    				:: days,hours,minutes,seconds, io_int(4)
real(kind=sgl)       				:: secs

  secs = tm
  days = 0
  hours = 0
  minutes = 0
  if (secs.gt.86400.0) then
    days = int(secs)/86400
    secs = mod(secs,86400.0)
  end if
  if (secs.gt.3600.0) then
    hours = int(secs)/3600
    secs = mod(secs,3600.0)
  end if
  if (secs.gt.60.0) then
    minutes = int(secs)/60
    secs = mod(secs,60.0)
  end if
  seconds = int(secs)
  io_int(1:4) = (/ days, hours, minutes, seconds /)
  call WriteValue(' ',io_int, 4, "(1x,I3,' d,',I3,' h,',I3,' m,',I3,' s')")
end subroutine PrintTime

!--------------------------------------------------------------------------
!
! SUBROUTINE: Time_stop
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief stop time recording
!
!> @param numk total number of computations
!
!> @date   06/04/01 MDG 1.0 original
!> @date   06/04/13 MDG 2.0 rewrite
!--------------------------------------------------------------------------
subroutine Time_stop(numk)

use local
use io

IMPLICIT NONE

integer(kind=irg),INTENT(IN)  :: numk

real(kind=sgl)				:: io_real(1)

  call system_clock(TIME_newcount,TIME_count_rate,TIME_count_max)
  mess = '  Total computation time [s] '; call Message("(A$)")
  call PrintTime((float(TIME_loops)*float(TIME_count_max)+float(TIME_newcount-TIME_count))/float(TIME_count_rate))
  io_real(1)=(float(TIME_loops)*float(TIME_count_max)+float(TIME_newcount-TIME_count))/float(TIME_count_rate)/float(numk)
  call WriteValue(' Time per step/pixel [s] ', io_real, 1, "(F10.3)")
end subroutine Time_stop


end module timing
