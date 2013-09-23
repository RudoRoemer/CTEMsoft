!----------------------------------------------------------------
! CTEMsoft                                                      !
! Copyright (c) 2001, Marc De Graef/Carnegie Mellon University  !
! All rights reserved.                                          !
!                                                               !
! The source code in this file is protected by a BSD license;   !
! for the complete text of the license see the file BSD.license.!
!                                                               !
!----------------------------------------------------------------
! ###################################################################
!  Introduction to Conventional Transmission Electron Microscopy
!
!  by Marc De Graef
!
!  Cambridge University Press
!  SBN 0521629950 (Paperback)
!  SBN 0521620066 (Hardback)
! 
!  FILE: "vrfft.f90"
!                                    created: 11/26/00 {9:29:46 AM} 
!                                last update: 11/26/00 {10:37:32 AM} 
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description: This program computes the electrostatic lattice
!               potential for a planar cut through the lattice.
!               This version uses the in-place 3-D FFTW functions.
!               The size of the array is limited by dynamic memory.
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  11/26/00 MDG 1.0 original
!   5/22/01 MDG 2.0 f90
! ###################################################################
program vrfft

use local
use io
use crystalvars
use crystal
use symmetryvars
use symmetry
use files

integer  :: np(3)


 progname = 'vrfft.f90'
 progdesc = '3D electrostatic lattice potential via fft'
 call CTEMsoft
 

 SG % SYM_reduce=.TRUE.
! read crystal information
 call CrystalData
! generate all atom positions in the fundamental unit cell
 call CalcPositions('v')
! ask the user for the number of pixels along the a-axis
! (need not be a power of 2)
 mess = 'Number of pixels along a, b, c for 3D FFT (even numbers)'; call GetInt(3)
 do i=1,3
  np(i) = io_int(i)
 end do
 mess = 'FFT real space step sizes (nm) '; call Message("(A)")
 mess = 'along a : '; oi_real(1) = cell % a/float(np(1)); call WriteReal(1,"(f10.5)")
 mess = 'along b : '; oi_real(1) = cell % b/float(np(2)); call WriteReal(1,"(f10.5)")
 mess = 'along c : '; oi_real(1) = cell % c/float(np(3)); call WriteReal(1,"(f10.5)")
 call ComputeVreal(np)
!
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

! parameters for fftw package
integer,parameter         ::  FFTW_FORWARD=-1,FFTW_BACKWARD=1,  &
                              FFTW_REAL_TO_COMPLEX=-1,FFTW_COMPLEX_TO_REAL=1, &
                              FFTW_ESTIMATE=0,FFTW_MEASURE=1, &
                              FFTW_OUT_OF_PLACE=0,FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16, &
                              FFTW_THREADSAFE=128
! other variables
complex,allocatable       :: Vreal(:,:,:)
logical,allocatable       :: z(:,:,:)
complex                   :: zVg
real                      :: cc(3),r(3),Vmod,Vphase,Vpmod,Vpphase
integer(kind=irg)         :: ind(3),M(3),ih,ik,il,np(3),fcnt,facnt,hra,kra,lra,h,k,l
logical                   :: first 

! initialize the boolean array z
 hra = np(1)/2
 kra = np(2)/2
 lra = np(3)/2
 allocate(z(-hra+1:hra,-kra+1:kra,-lra+1:lra))
 z(-hra+1:hra,-kra+1:kra,-lra+1:lra) = .FALSE.
! initialize the potential array
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
 mess = 'Total number of reflections = '; oi_int(1)=fcnt; call WriteInt(1,"(I8)")
 mess = 'Total number of distinct families = '; oi_int(1)=facnt; call WriteInt(1,"(I8)")
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
