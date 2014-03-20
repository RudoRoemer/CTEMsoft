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
! CTEMsoft2013:CTEMTBBFDF.f90
!--------------------------------------------------------------------------
!
! PROGRAM:CTEMTBBFDF 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Two beam bright field dark field image pairs
!
!> @details  there are a few different options for the foil geometry
! 
!> @date   12/11/98 MDG 1.0 original
!> @date    5/27/01 MDG 2.0 f90
!> @date 4/16/13 MDG 3.0 rewrite
!--------------------------------------------------------------------------
program CTEMTBBFDF

use local
use io
use constants
use files
use crystal
use symmetry
use diffraction
use postscript
use dynamical

IMPLICIT NONE

real(kind=sgl),parameter      	:: thr = 1.0E-6
real(kind=sgl),allocatable    	:: BF(:,:), DF(:,:), exe(:,:), z(:,:)
real(kind=sgl)                	:: dz, zmin, zmax
integer(kind=irg)             	:: jdim,zero(3),ind(3), io_int(3), iflag, j,i, ihole, jhole, irad, jrad, nw, &
				   npx, npy, k
real(kind=sgl)                	:: r(200),p(200),xig,xigp,betag,xizero, io_real(3), ss, sm, f1, dd, T0, T1, &
				   zav, zdev, ff, fi, fj, ztot, q, zmi, zma, tps, x0, y0, scl
integer(kind=irg)		:: values(1:8), kk
integer, dimension(:), allocatable :: seed
real(kind=dbl) 			:: zz

 progname = 'CTEMTBBFDF.f90'
 progdesc = 'Two-beam bright field-dark field images, using BFDF.routines'
 call CTEMsoft
                                                  
! compute the relevant parameters for a given crystal
 imanum = 0
 call CrystalData
 call GetVoltage
 call CalcPositions('v')

! initialize all arrays
 call ReadValue(' Enter dimension of image array ', io_int, 1)
 jdim = io_int(1)
 allocate(BF(jdim,jdim))
 allocate(DF(jdim,jdim))
 allocate(exe(jdim,jdim))
 allocate(z(jdim,jdim))

! the following array was declared in module postscript.f90
 allocate(imaint(jdim,jdim))

! extinction distance and absorption length for given plane
 mess = 'Indices for reciprocal lattice vector g'; call Message("(A)")
 call GetIndex(ind,'r')
 call CalcUcg(ind)
 xigp = rlp%xgp
 xig  = rlp%xg
 betag = rlp%Vpphase-rlp%Vphase
 io_real(1) = xig
 call WriteValue(' Extinction distance [nm] = ', io_real, 1, "(f10.4)")
 io_real(1) = xigp
 call WriteValue(' Absorption length   [nm] = ', io_real, 1, "(f10.4)")

! normal absorption length
 zero(1)=0
 zero(2)=0
 zero(3)=0
 call CalcUcg(zero)
 xizero = rlp%xgp
 io_real(1) = xizero
 call WriteValue(' Normal Absorption   [nm] = ', io_real, 1, "(f10.4)")
 
! the following used to be in a separate file (BFDF.routines); in the current version,
! the user is prompted for a selection.

mess = ' Which foil geometry do you wish to use ? '; call Message("(/A)")
mess = '  [1] Straight wedge'; call Message("(A)")
mess = '  [2] Bent wedge'; call Message("(A)")
mess = '  [3] Bent wedge with a hole'; call Message("(A)")
mess = '  [4] Random superposition of cosine waves'; call Message("(A)")
call ReadValue('   -> your choice : ', io_int, 1)
mess = ' '; call Message("(A)")

select case (io_int(1))
  case (1)
! Case 1 straight wedge
! ------
! create the thickness array 
	call ReadValue(' Enter minimum,maximum thickness [nm]  = ', io_real, 2)
 	dz = (io_real(2)-io_real(1))/float(jdim)
 	do j=1,jdim
   	  z(1:jdim,j) = float(j)*dz
 	end do
! create the excitation error array 
 	call ReadValue(' Enter w-parameter (s * xi)  = ', io_real, 1)
 	ss = io_real(1)/xig
  	exe = ss
 case (2)
! Case 2 bent wedge
! ------
! create the thickness array 
	call ReadValue(' Enter minimum,maximum thickness [nm]  = ', io_real, 2)
 	dz = (io_real(2)-io_real(1))/float(jdim)
 	do j=1,jdim
   	  z(1:jdim,j) = float(j)*dz
 	end do
	! create the excitation error array 
 	call ReadValue(' Enter w-parameter (s * xi)  = ', io_real, 1)
	 sm = io_real(1)
	 do i=1,jdim
	   ss = (-sm+2.0*sm*float(i-1)/float(jdim))/xig
	   exe(i,1:jdim) = ss
	 end do
 case(3)
! Case 3 bent foil with circular hole 
! ------
	 call ReadValue(' Enter minimum,maximum thickness [nm]  = ', io_real, 2)
	 zmin = io_real(1)
	 zmax = io_real(2)
	 dz = (zmax-zmin)/float(jdim)
	 ihole = 4*jdim/6
	 jhole = jdim/2
	! outer radius
	 irad = 150
	! inner radius
	 jrad = 40
	 f1 = 1.0/float(irad-jrad)
	 do i=1,jdim
	  do j=1,jdim
	   dd = sqrt(float(i-ihole)**2+float(j-jhole)**2)
	! inside or outside the hole ?
	   z(i,j) = float(j)*dz
	   if (dd.lt.jrad) then
	    z(i,j) = 0.D0
	   else
	    if (dd.lt.irad) then
	     z(i,j)=z(i,j)*f1*(dd-float(jrad))
	    end if
	   end if
	  end do
	 end do
	! create the excitation error array 
 	call ReadValue(' Enter w-parameter (s * xi)  = ', io_real, 1)
	 sm = io_real(1)
	 do i=1,jdim
	   ss = (-sm+2.0*sm*float(i-1)/float(jdim))/xig
	   exe(i,1:jdim) = ss
	 end do
case(4)
! Case 4 Random superposition of cosine waves
! ------
	 mess = ' Average thickness, minimum thickness [nm]  = '; call Message("(A)")
	 call ReadValue(' (minimum thickness may be <0, to create holes)', io_real, 2)
	 zav = io_real(1)
	 zdev = io_real(2)
	 call ReadValue(' How many random waves ? (<100) = ', io_int, 1)
	 nw = io_int(1)
	 nw = min(nw,100)
	 ff = cPi/float(jdim)

	! superimpose a bunch of waves with random periodicities and phases
	  call date_and_time(values=values)
	  kk = 1
	  call random_seed(size=kk)
	  allocate(seed(1:k))
	  seed(:) = values(8)
	  call random_seed(put=seed)
	  do k=1,2*nw
	   call random_number(zz)
	   r(k) = zz*10.0
	   call random_number(zz)
	   p(k) = zz*3.0
	  end do

	 do i=1,jdim
	  fi = float(i)*ff
	  do j=1,jdim
	   fj=float(j)*ff
	   z(i,j) = 0.0 
	   do k=1,nw
	    z(i,j) = z(i,j) + cos(r(2*k)*fi+p(2*k))*cos(r(2*k-1)*(fi+fj)+p(2*k-1)) + &
	                      sin(r(2*k-1)*fj+p(2*k))*sin(r(2*k)*(fi-fj)+p(2*k-1))
	   end do
	  end do
	 end do

	 zmin = minval(z)
	 zmax = maxval(z)
	 ztot = sum(z)/float(jdim)**2
	 io_real(1) = zmin
	 io_real(2) = zmax
	 call WriteValue(' min max after generation : ', io_real, 2, "(1x,2(f8.3,2x))")

	! create the excitation error array 
	 q = (zav-zdev)/(ztot-zmin)
	 do i=1,jdim
	  ss = (-0.5+1.0*float(i-1)/float(jdim))/xig
	  do j=1,jdim
	   exe(i,j) = ss
	   z(i,j) = zdev+q*(z(i,j)-zmin)
	   if (z(i,j).lt.0.0) z(i,j)=0.0
	  end do
	 end do

	 zmi = minval(z)
	 zma = maxval(z)
	 ztot = sum(z)/float(jdim)**2
	 io_real(1) = zmi
	 io_real(2) = zma
	 call WriteValue(' min max after scaling    : ', io_real, 2, "(1x,2(f8.3,2x))")
	 io_real(1)=ztot
	 call WriteValue(' Actual average thickness = ', io_real, 1, "(f10.4)")
	 zmin=zmi 
	 zmax=zma 
case default
	 mess = 'Option not implemented'; call Message("(A)")
end select	


! compute intensities
 T0= SECNDS(0.0)
 mess = 'computation start '; call Message("(A)")
 iflag=0
 do i=1,jdim
  do j=1,jdim
   call TBCalcInten(BF(i,j),DF(i,j),exe(i,j),z(i,j),xig,xigp,xizero,betag)
   if (iflag.eq.0) then 
    if ((BF(i,j).gt.1.D0).or.(DF(i,j).gt.1.D0)) iflag=1
   endif
  end do
 end do
 T1= SECNDS(T0)

 io_real(1)=T1
 call WriteValue(' computation end   -> ', io_real, 1, "(F,' seconds')")
 tps = T1/float(jdim)**2
 io_real(1) = tps 
 call WriteValue(' time per pixel    -> ', io_real, 1, "(F,' seconds')")
 if (iflag.eq.1) then 
  mess = 'WARNING: one or more intensities larger than 1 !'; call Message("(A)")
  mess = 'Normal absorption length is too large'; call Message("(A)")
 endif


! open PostScript file
 call PS_openfile
 PS % pspage = 0
 call PS_newpage(.TRUE.,'Two Beam Bright Field - Dark Field')
 call PS_setfont(PSfonts(2),0.10)
 call PS_text(0.5,1.8,'Input file')
 call PS_text(2.6,1.8,PS % psname)
 call PS_text(0.5,1.6,'Active reflection')
 call PS_textint(2.5,1.6,' ',ind(1))
 call PS_textint(2.6,1.6,' ',ind(2))
 call PS_textint(2.7,1.6,' ',ind(3))
 call PS_text(0.5,1.4,'Extinction distance  [nm]')
 call PS_textvar(2.5,1.4,' ',xig)
 call PS_text(0.5,1.2,'Anomalous absorption length [nm]')
 call PS_textvar(2.5,1.2,' ',xigp)
 call PS_text(0.5,1.0,'Normal absorption length [nm]')
 call PS_textvar(2.5,1.0,' ',xizero)
 call PS_text(0.5,0.8,'Accelerating Voltage [V]')
 call PS_textvar(2.5,0.8,' ',sngl(mAccvol))
 call PS_text(0.5,0.6,'Minimum thickness [nm]')
 call PS_textvar(2.5,0.6,' ',zmin)
 call PS_text(0.5,0.4,'Maximum thickness [nm]')
 call PS_textvar(2.5,0.4,' ',zmax)


! Bright Field image
 imanum = 0
 imaint = int(255.0*BF)
 x0=0.2
 y0=4.0
 npx=jdim
 npy=jdim
 scl=3.0 * float(jdim)/512.0
 call PS_DumpImageDistort(x0,y0,npx,npy,scl,scl)

! Dark Field image
 imaint=int(255.0*DF)
 x0=3.4
 y0=4.0
 npx=jdim
 npy=jdim
 scl=3.0 * float(jdim)/512.0
 call PS_DumpImageDistort(x0,y0,npx,npy,scl,scl)

! close Postscript file
 call PS_closefile

end program CTEMTBBFDF
