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
! CTEMsoft2013:CTEMvrfft.f90
!--------------------------------------------------------------------------
!
! PROGRAM:CTEMvrfft.f90 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief 2D electrostatic lattice potential 
!
!> @details  This program computes the electrostatic lattice
!> potential for a planar cut through the lattice.
!> This version uses the in-place 3-D FFTW functions.  
! 
!> @date 11/26/00   MDG 1.0 original
!> @date  5/22/01 MDG 2.0 f90
!> @date 4/9/13 MDG 3.0 rewrite
!--------------------------------------------------------------------------
program vrfft

use local
use io
use crystalvars
use crystal
use symmetryvars
use symmetry
use files

IMPLICIT NONE

integer(kind=irg)  	:: np(3), io_int(3), i
real(kind=sgl)		:: io_real(1)

 progname = 'CTEMvrfft.f90'
 progdesc = '3D electrostatic lattice potential via fft'
 call CTEMsoft
 

 SG % SYM_reduce=.TRUE.

! read crystal information
 call CrystalData

! generate all atom positions in the fundamental unit cell
 call CalcPositions('v')

! ask the user for the number of pixels along the a-axis
! (need not be a power of 2)
 call ReadValue(' Number of pixels along a, b, c for 3D FFT (even numbers)', io_int, 3)
 do i=1,3
  np(i) = io_int(i)
 end do
 
 mess = 'FFT real space step sizes (nm) '; call Message("(A)")
 io_real(1) = cell % a/float(np(1))
 call WriteValue(' along a : ', io_real, 1, "(f10.5)")
 io_real(1) = cell % b/float(np(2))
 call WriteValue(' along b : ', io_real, 1, "(f10.5)")
 io_real(1) = cell % c/float(np(3))
 call WriteValue(' along c : ', io_real, 1, "(f10.5)")

 call ComputeVreal(np)

end program
!
!
!
subroutine ComputeVreal(np)

use local
use io
use crystalvars
use crystal
use diffraction
use dynamical
use constants
use symmetryvars
use symmetry
use graphics
use postscript
use files

complex(kind=dbl),allocatable	:: Vreal(:,:,:)

! parameters for fftw package
integer,parameter         ::  FFTW_FORWARD=-1,FFTW_BACKWARD=1,  &
                              		FFTW_REAL_TO_COMPLEX=-1,FFTW_COMPLEX_TO_REAL=1, &
                              		FFTW_ESTIMATE=0,FFTW_MEASURE=1, &
                              		FFTW_OUT_OF_PLACE=0,FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16, &
                              		FFTW_THREADSAFE=128

logical,allocatable       	:: z(:,:,:)
complex                   	:: zVg
real(kind=sgl)             	:: cc(3),r(3),Vmod,Vphase,Vpmod,Vpphase
integer(kind=irg)         	:: ind(3),M(3),ih,ik,il,np(3),fcnt,facnt,hra,kra,lra,h,k,l, io_int(1)
logical                   		:: first 
real(kind=sgl)			:: io_real(1)

! initialize the boolean array z
 hra = np(1)/2
 kra = np(2)/2
 lra = np(3)/2
 allocate(z(-hra+1:hra,-kra+1:kra,-lra+1:lra))
 z(-hra+1:hra,-kra+1:kra,-lra+1:lra) = .FALSE.

! initialize the potential array using the fftw alloc routine
 allocate(Vreal(np(1),np(2),np(3)))

! next, fill the reciprocal potential array with Fourier coefficients
 first = .TRUE.
 fcnt=0
 facnt=0
 mess = 'starting computation of Vg''s'; call Message("(A)")

 do h=-hra+1,hra
  ind(1)=h
  do k=-kra+1,kra
   ind(2)=k
   do l=-lra+1,lra
    ind(3)=l

! make sure we have not already done this one in another family
    if (.not.z(h,k,l)) then

! if it is a new one, then determine the entire family, assuming |g|<gmax
     call CalcUcg(ind)
     Vmod = rlp%Vmod
     Vphase = rlp%Vphase
     Vpmod = rlp%Vpmod
     Vpphase = rlp%Vpphase
     zVg=cmplx(Vmod*cos(Vphase),Vmod*sin(Vphase))

! tag all the family members
     call CalcFamily(ind,num,'r')
     do i=1,num
      ih = itmp(i,1)
      ik = itmp(i,2)
      il = itmp(i,3)
      if ((((ih.gt.-hra).and.(ih.le.hra)).and.((ik.gt.-kra).and.(ik.le.kra))).and.((il.gt.-lra).and.(il.le.lra))) z(ih,ik,il)=.TRUE.

! put the values in the correct location in the FFT array Vreal
      if (ih.lt.0) ih=np(1)+ih
      if (ik.lt.0) ik=np(2)+ik
      if (il.lt.0) il=np(3)+il
      ih = ih + 1
      ik = ik + 1
      il = il + 1

! make sure that this point belongs in the array
      if ((((ih.ge.1).and.(ih.le.np(1))).and.((ik.ge.1).and.(ik.le.np(2)))).and.((il.ge.1).and.(il.le.np(3)))) then
        Vreal(ik,ih,il)=zVg
      end if
     end do
     fcnt = fcnt+num
     facnt = facnt+1
    end if
   end do
  end do
  if (mod(h,10).eq.0) write (*,*) 'completed plane ',h
 end do
 io_int(1)=fcnt
 call WriteValue(' Total number of reflections = ', io_int, 1, "(I8)")
 io_int(1)=facnt
 call WriteValue(' Total number of distinct families = ', io_int, 1, "(I8)")
!
! do the 3D inverse discrete fast Fourier transform (fftw)
!
 do i=1,3
  M(i) = np(i)
 end do
 mess = 'initializing FFT'; call Message("(A)")

 call fftwnd_f77_create_plan(plan,3,M,FFTW_BACKWARD,FFTW_ESTIMATE+FFTW_IN_PLACE)
 call fftwnd_f77_one(plan,Vreal,Vreal)
 call fftwnd_f77_destroy_plan(plan)

 mess = 'FFT computed'; call Message("(A)")
!
! and save the entire array in a file, containing coordinates
! (cartesian) and potential value.
!
 open (unit=8,name='vrfft.out',status='unknown',form='unformatted')
 write (8) Vreal
 close (unit=8,status='keep') 
 mess = 'potential saved in file vrfft.out'; call Message("(A)")

end subroutine
