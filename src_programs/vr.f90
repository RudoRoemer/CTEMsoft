!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CTEMsoft:  vr.f90                                                             !
! Copyright (c) 2001, 2002  Marc De Graef/Carnegie Mellon University (CMU)      !
!                                                                               !
!     This program is free software; you can redistribute it and/or modify      !
!     it under the terms of the GNU General Public License as published by      !
!     the Free Software Foundation; either version 2 of the License, or         !
!     (at your option) any later version.                                       !
!                                                                               !
!     This program is distributed in the hope that it will be useful,           !
!     but WITHOUT ANY WARRANTY; without even the implied warranty of            !
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             !
!     GNU General Public License for more details.                              !
!                                                                               !
!     You should have received a copy of the GNU General Public License         !
!     along with this program; if not, write to the Free Software               !
!     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA !
!                                                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ###################################################################
!  Introduction to Conventional Transmission Electron Microscopy
!
!  by Marc De Graef
!
!  Cambridge University Press
!  SBN 0521629950 (Paperback)
!  SBN 0521620066 (Hardback)
! 
!  FILE: "vr.f90"
!                                    created: 12/06/98 {9:29:46 AM} 
!                                 
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description: This program computes the electrostatic lattice
!               potential for a planar cut through the lattice.
!               This is the naive version, which does not use
!               the FFT algorithm.  It does include the imaginary
!               part of the potential, using the Weickenmeier-Kohl
!               parametrization of the scattering and form factors.
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  12/06/98 MDG 1.0 original
!  11/26/00 MDG 2.0 modified range of reciprocal space (g_max)
!   5/22/01 MDG 3.0 f90
! ###################################################################
program vr

use local
use io
use crystalvars
use crystal
use diffraction
use constants
use symmetryvars
use symmetry
use graphics
use postscript
use files
use dynamical

character(1)               :: sgn,ans
character(20)              :: vrname,vprname
character(80)              :: line
logical                    :: first,again,nn,more,topbot
logical,allocatable        :: z(:,:,:)
integer(kind=irg),allocatable        :: family(:,:,:),numfam(:)
integer(kind=irg)                    :: h,k,l,totfam,hkl(3),ind(3),hra,kra,lra,ntot,fcnt,i,ii
real(kind=sgl),allocatable           :: Vpg(:,:,:),Vgg(:,:,:),V(:,:),Vp(:,:)
real(kind=sgl)                       :: dx(3),dy(3),twopi,arg,lg,find(3),vmax,vpmax
real(kind=sgl)                       :: rr(4),gg(4),g(3),r(3),cc,Vmod,Vphase,gmax,Vpmod,Vpphase

 progname = 'vr.f90'
 progdesc = '2D section of electrostatic lattice potential'
 call CTEMsoft

 SG % SYM_reduce=.TRUE.
 thr = 1.D-6
 twopi = 2.D0*cPi
! read crystal information
 call CrystalData
 mess = 'Select option 3 to include absorption potential'; call Message("(/A/)")
 call GetVoltage
! generate all atom positions in the fundamental unit cell
 call CalcPositions('v')
! ask the user for the maximum g-vector length to contribute
! to the summation
 mess = 'Enter maximum length of g (nm^-1) : '; call GetReal(1)
 gmax = io_real(1)
! determine the range of reflections inside the sphere with radius
! gmax, along the three main reciprocal directions
  hra = int(gmax/sqrt(cell % rmt(1,1)))
  kra = int(gmax/sqrt(cell % rmt(2,2)))
  lra = int(gmax/sqrt(cell % rmt(3,3)))
! allocate arrays
  allocate(z(-hra:hra,-kra:kra,-lra:lra))
  z(-hra:hra,-kra:kra,-lra:lra) = .FALSE.
  ntot = (2*hra+1)*(2*kra+1)*(2*lra+1)
!  ntot = ntot/32
  write (*,*) 'ntot = ',ntot
  allocate(family(ntot,SG % SYM_numpt,3))
  allocate(numfam(ntot)) 
  allocate(Vpg(ntot,48,2))
  allocate(Vgg(ntot,48,2))
! next, perform the summation over all g-vectors for which |g|<gmax        
 first = .TRUE.
 fcnt = 1
 totfam=0
 do h=-hra,hra
  ind(1)=h
  find(1)=float(h)
  do k=-kra,kra
   ind(2)=k
   find(2)=float(k)
   do l=-lra,lra
    ind(3)=l
    find(3)=float(l)
! make sure we have not already done this one in another family
    if (.not.z(h,k,l)) then
! if it is a new one, then determine the entire family, assuming |g|<gmax
     lg = CalcLength(find,'r')
! if g is too long, then eliminate all family members of g
! from the remainder of the computation.
     if (lg.gt.gmax) then
      call CalcFamily(ind,num,'r')
      do i=1,num
       do j=1,3
        family(fcnt,i,j)=itmp(i,j)
       end do
       z(itmp(i,1),itmp(i,2),itmp(i,3))=.TRUE.
      end do
     else
! if |g|<gmax, then compute the Fourier coefficients
      call CalcFamily(ind,num,'r')
      do i=1,num
        do j=1,3
         family(fcnt,i,j)=itmp(i,j)
        end do
! tag all the family members
        do ii=1,num
         z(itmp(ii,1),itmp(ii,2),itmp(ii,3))=.TRUE.
        end do
        ind(1:3)=itmp(i,1:3)
	call CalcUcg(ind)
	if (rlp%Vmod.ge.thr) then
         Vgg(fcnt,i,1) = rlp%Vmod*cos(rlp%Vphase)
         Vgg(fcnt,i,2) = rlp%Vmod*sin(rlp%Vphase)
         Vpg(fcnt,i,1) = rlp%Vpmod*cos(rlp%Vpphase)
         Vpg(fcnt,i,2) = rlp%Vpmod*sin(rlp%Vpphase)	
! increment family counter
         numfam(fcnt)=num
         totfam=totfam+num-1
         fcnt=fcnt+1
         if (fcnt.gt.ntot) then
          mess = 'Number of families larger than maximum allocated.'; call Message("(A)")
          mess = 'Increase value of ntot parameter in source code.'; call Message("(A)")
          stop
         end if
        end if
      end do
     end if
    end if
   end do
  end do
 end do
 mess = 'Total number of families        = '; oi_int(1)=fcnt; call WriteInt(1,"(I6)")
 mess = 'Total number of family members  = '; oi_int(1)=totfam; call WriteInt(1,"(I6)")
!
! now create the planar section for the potential
!
 mess = 'Enter fractional coordinates for '; call Message("(A)")
 mess = 'lower left corner  : '; call GetReal(3)
 x1 = io_real(1)
 y1 = io_real(2)
 z1 = io_real(3)
 mess = 'lower right corner : '; call GetReal(3)
 x2 = io_real(1)
 y2 = io_real(2)
 z2 = io_real(3)
 mess = 'upper right corner : '; call GetReal(3)
 x3 = io_real(1)
 y3 = io_real(2)
 z3 = io_real(3)
 mess = 'Number of pixels nx,ny : '; call GetInt(2)
 nx = io_int(1)
 ny = io_int(2)
 dx(1)=(x2-x1)/float(nx)
 dx(2)=(y2-y1)/float(nx)
 dx(3)=(z2-z1)/float(nx)
 dy(1)=(x3-x2)/float(ny)
 dy(2)=(y3-y2)/float(ny)
 dy(3)=(z3-z2)/float(ny)
! now loop over all these points and compute V(r) and Vprime(r)
 allocate(V(nx,ny))
 allocate(Vp(nx,ny))
 vmax = 0.0
 vpmax = 0.0
 do i=1,nx
  do j=1,ny
! compute the position r
   r(1)=x1+float(i-1)*dx(1)+float(j-1)*dy(1)
   r(2)=y1+float(i-1)*dx(2)+float(j-1)*dy(2)
   r(3)=z1+float(i-1)*dx(3)+float(j-1)*dy(3)
! sum over all families
   V(i,j)=0.0
   Vp(i,j)=0.0
   do k=1,fcnt-1
    sreal = 0.D0
    simag = 0.D0
    do l=1,numfam(k)
     arg = 0.D0
     do h=1,3
      arg = arg + family(k,l,h)*r(h)
     end do
     cc = cos(twopi*arg)
     ss = sin(twopi*arg)
     V(i,j) = V(i,j)   + Vgg(k,l,1)*cc-Vgg(k,l,2)*ss
     Vp(i,j) = Vp(i,j) + Vpg(k,l,1)*cc-Vpg(k,l,2)*ss
    end do
   end do
   vmax = max(V(i,j),vmax)
   vpmax = max(Vp(i,j),vpmax)
  end do
  if (mod(i,10).eq.0) write (*,*) 'finishing column ',i
 end do
 write (*,*) 'maxima : ',vmax,vpmax
 deallocate(family)
 deallocate(numfam) 
 deallocate(Vpg)
 deallocate(Vgg)
!
! prepare for PostScript output of rendered surfaces
!
! PostScript dimension parameters
 AX % axw = 5.0
 AX % xll = 4.25
 AX % yll = 7.00
 write (*,"(1x,' V(r) temporary PostScript file name : ',$)")
 read (*,"(A20)") vrname
 call axonometry(V,nx,ny,1.0,vrname)

 AX % axw = 5.0
 AX % xll = 4.25
 AX % yll = 2.00
 write (*,"(1x,' V''(r) temporary PostScript file name : ',$)")
 read (*,"(A20)") vprname
 call axonometry(Vp,nx,ny,1.0,vprname)

 mess = 'Combining the two rendered surfaces into a single file'; call Message("(A)")
 call PS_openfile
! first file
 open(UNIT=dataunit,FILE=vrname,STATUS='OLD',FORM='FORMATTED')
 i=0
 ios=0
 do while (ios.ne.-1)
  read (dataunit,fmt="(A)",iostat=ios) line
  write (psunit,*) line
  i = i+1
 end do
 close(dataunit,status='keep')        
 write (*,*) i,' lines transferred'
 write (psunit,*) 'grestore'
! second file
 open(UNIT=dataunit,FILE=vprname,STATUS='OLD',FORM='FORMATTED')
 i=0
 ios=0
 do while (ios.ne.-1)
  read (dataunit,fmt="(A)",iostat=ios) line
  write (psunit,*) line
  i = i+1
 end do
 close(dataunit,status='keep')        
 write (*,*) i,' lines transferred'
 write (psunit,*) 'grestore'
 call PS_closefile
end program
