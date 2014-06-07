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
! CTEMsoft2013:EDSabscor.f90
!--------------------------------------------------------------------------
!
! PROGRAM: EDSabscor 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief EDS absorption correction factors for four-quadrant detector
!
!> @date 11/06/13 MDG 1.0 initial implementation based on IDL and Mathematica codes
!--------------------------------------------------------------------------
program EDSabscor

use local
use files
use io

IMPLICIT NONE

character(fnlen)	:: cvsfilename
real(kind=sgl)		:: alp, bet		! primary and secondary tilt angles
real(kind=sgl)		:: io_real(2)


mess = 'Absorption Correction Factor Calculation for 4-quadrant Detector'
call Message("(/A/)")

call ReadValue('Enter alpha and beta tilt angles (degrees) : ',io_real,2)
alp = io_real(1)
bet = io_real(2)

call ReadValue('Enter file name for .cvs output of correction factors : ',cvsfilename, '(A)')

mess = 'Computing correction factors now ...'
call Message("(/A/)")

call GetAbsorptionCorrections(alp, bet, cvsfilename)

mess = 'Correction factors stored in file '//trim(cvsfilename)
call Message("(/A/)")

end program EDSabscor

!--------------------------------------------------------------------------
!
! SUBROUTINE:GetAbsorptionCorrections
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief compute the absorption corrections as a function of mu t for 
!> individual quadrants and summed detector mode
!
!> @note The detector parameters are read from a file.
!
!> @param a alpha tilt angle
!> @param b beta tilt angle
!> @param cvsfile output file name (CVS format)
!
!> @date 11/06/13  MDG 1.0 original
!--------------------------------------------------------------------------
subroutine GetAbsorptionCorrections(a,b,cvsfile)

use local
use io
use constants
use error
use files

IMPLICIT NONE

real(kind=sgl),INTENT(IN)	:: a,b
character(fnlen),INTENT(IN)	:: cvsfile

! detector parameters
real(kind=sgl)			:: Rbar, alphad, deltad

! other parameters
integer(kind=irg),parameter	:: dimx = 399, numut = 1000
real(kind=sgl)			:: trunc, eps, cad, sad, cod, cdd, sdd, cc(4), ss(4), apar, bpar, cpar
real(kind=sgl)			:: dtor, ca, sa, cb, sb, rho, Carray(5, numut), dx, dx2, mut(numut), alpha, beta
integer(kind=irg)		:: i, iut

! coordinate arrays
real(kind=sgl)			:: line(0:dimx), xx(0:dimx,0:dimx), yy(0:dimx,0:dimx), sq(0:dimx,0:dimx), &
				mask(0:dimx,0:dimx), num(0:dimx,0:dimx), denom(0:dimx,0:dimx), Et(0:dimx,0:dimx)


dtor = sngl(cPi/180.D0)

! read detector parameters
!call getDetectorParameters(Rbar, alphad, deltad)

Rbar = 0.282
alphad = 22.5*dtor
deltad = 8.0*dtor

! initialize all other parameters
alpha = a*dtor
beta = b*dtor
cad = cos(alphad)
sad = sin(alphad)
cod = 1.0/tan(alphad)
cdd = cos(deltad)
sdd = sin(deltad)

ca = cos(alpha)
sa = sin(alpha)
cb = cos(beta)
sb = sin(beta)

cc = (/ 1.0,1.0,-1.0,-1.0 /) /sqrt(2.0)
ss = (/ 1.0,-1.0,-1.0,1.0 /) /sqrt(2.0)

rho = 1.0/sngl(cPi)/Rbar**2

mut = (/ (i,i=1,numut) /) * 0.005
trunc = 100.0
eps = 0.0001

! generate the coordinate arrays and the mask
line = (/ (i, i=0,dimx) /)
line = (2.0*line/float(dimx)-1.0) * Rbar
do i=0,dimx 
  xx(0:dimx,i) = line(0:dimx)
  yy(i,0:dimx) = line(0:dimx)
end do
sq = sqrt(xx*xx+yy*yy)
where (sq.le.maxval(abs(line)))
  mask = 1.0
elsewhere
  mask = 0.0
end where
sq = sqrt(xx*xx+yy*yy+1.0)
dx = line(3)-line(2)
dx2 = dx*dx


! ok, everything is set to go
do i=1,4
  apar = cc(i)*cb*sa-ss(i)*sb
  bpar = cdd*ca*cb+sdd*apar
  cpar = sad*ca*cb-cad*apar
  apar = ss(i)*cb*sa+cc(i)*sb

  denom = apar*xx+bpar*yy+cpar
  where (denom.le.0.0) 
    num = 0.0
  elsewhere
    num = abs(sq/denom)*mask
  end where 
  where ((num.lt.eps).and.(num.gt.0.0)) num = eps
  where (num.gt.trunc) num = trunc
  
! loop for all values of mu t
  do iut=1,numut 
    Et = 0.0
    where (num.gt.eps)
      Et = (1.0-exp(- mut(iut) * num))/(mut(iut)*num)
    end where
    Et = Et * mask
    where (denom.le.0.0) Et = 0.0
    Carray(i,iut) = sum(rho*Et) * dx2 / (ca*cb)
  end do
end do

do iut=1,numut
  Carray(5,iut) = sum(Carray(1:4,iut))*0.25
end do

! that's it, so now we write everything to a cvs file.
open(unit=dataunit,file=trim(cvsfile),status='unknown',form='formatted')
write (dataunit,"(I6)") numut
do iut=1,numut
  write(dataunit,"(5(f10.6,','),f10.6)") mut(iut),Carray(1:5,iut)
end do
close(unit=dataunit,status='keep')

end subroutine GetAbsorptionCorrections




