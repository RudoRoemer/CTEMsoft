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
! EMsoft:EMvr.f90
!--------------------------------------------------------------------------
!
! PROGRAM:EMvr 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief computes a slice of  V(r) using cosines and sines (slow version)
!
!> @details  This program computes the electrostatic lattice potential for a planar cut through the lattice.
!> This is the naive version, which does not use the FFT algorithm.  It does include the imaginary
!> part of the potential, using the Weickenmeier-Kohl parametrization of the scattering and form factors.
! 
!> @date  12/06/98 MDG 1.0 original
!> @date  11/26/00 MDG 2.0 modified range of reciprocal space (g_max)
!> @date   5/22/01 MDG 3.0 f90
!> @date   4/9/13 MDG 3.1 rewrite
!--------------------------------------------------------------------------
program EMvr

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

character(1)               			:: sgn,ans
character(20)              			:: vrname,vprname
character(80)              			:: line
logical                    				:: first,again,nn,more,topbot
logical,allocatable        			:: z(:,:,:)
integer(kind=irg),allocatable        	:: family(:,:,:),numfam(:)
integer(kind=irg)                    		:: h,k,l,totfam,hkl(3),ind(3),hra,kra,lra,ntot,fcnt,i,ii, io_int(2)
real(kind=sgl),allocatable           	:: Vpg(:,:,:),Vgg(:,:,:),V(:,:),Vp(:,:)
real(kind=sgl)                       		:: dx(3),dy(3),twopi,arg,lg,find(3),vmax,vpmax, io_real(3)
real(kind=sgl)                       		:: rr(4),gg(4),g(3),r(3),cc,Vmod,Vphase,gmax,Vpmod,Vpphase

 progname = 'EMvr.f90'
 progdesc = '2D section of electrostatic lattice potential'
 call EMsoft

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
 call ReadValue(' Enter maximum length of g (nm^-1) : ', io_real, 1)
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
!  write (*,*) "ntot = ",ntot  ! PGC debug

  allocate(family(ntot,SG % SYM_NUMpt,3))
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
 io_int(1) = fcnt
 call WriteValue(' Total number of families        = ', io_int, 1, "(I6)")
 io_int(1) = totfam
 call WriteValue(' Total number of family members   = ', io_int, 1, "(I6)")
!
! now create the planar section for the potential
!
 mess = 'Enter fractional coordinates for '; call Message("(A)")
 call ReadValue(' lower left corner  : ', io_real, 3)
 x1 = io_real(1)
 y1 = io_real(2)
 z1 = io_real(3)
 call ReadValue(' lower right corner  : ', io_real, 3)
 x2 = io_real(1)
 y2 = io_real(2)
 z2 = io_real(3)
 call ReadValue(' upper right corner  : ', io_real, 3)
 x3 = io_real(1)
 y3 = io_real(2)
 z3 = io_real(3)
 call ReadValue(' Number of pixels nx,ny : ', io_int, 2)
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
 call ReadValue(' V(r) temporary PostScript file name : ',vrname, "(A)")
 call axonometry(V,nx,ny,1.0,trim(vrname))

 AX % axw = 5.0
 AX % xll = 4.25
 AX % yll = 2.00
 call ReadValue(' V''(r) temporary PostScript file name : ', vprname, "(A)")
 call axonometry(Vp,nx,ny,1.0,trim(vprname))

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
