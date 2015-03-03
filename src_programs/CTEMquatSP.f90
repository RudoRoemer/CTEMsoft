! ###################################################################
! Copyright (c) 2015, Marc De Graef/Carnegie Mellon University
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
! CTEMsoft:CTEMquatSP.f90
!--------------------------------------------------------------------------
!
! Program: CTEMquatSP
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief creates a 3D volume containing a 3D stereographic projection histogram
!
!> @details [initial version written at UCSB, week of 2/16/15]
!> This program takes a binary input file with a list of Euler angle
!> triplets in radians.  The first entry in the file is the number of 
!> triplets as a 4-byte integer, followed by a (3,n) array of 4-byte reals.
!>
!> The program takes namelist file Stereogram.nml that contains the input 
!> parameters as well as the name of the output file;  this file then needs 
!> to be converted to the .mrc format and can be displayed using Chimera.
!> The conversion is controlled form an IDL script that can be found in the 
!> the folder CTEMsoft/IDL/various; in principle, the user can take an 
!> .ang file, dump all the Euler triplets in a data file and have this 
!> program spit out a 3D cube histogram for conversion to .mrc format.
!> There are likely other ways to do all this more efficiently, and completely
!> inside this fortran program, but that will have to wait until the next version.
! 
!> @date    02/17/15 MDG 1.0 original
!> @date    02/19/15 MDG 1.1 tested the use of 3D Gaussian blobs but so far this doesn't work well
!> @date    02/20/15 MDG 1.2 added option to draw Euler space instead of 3D-SP
!> @date    02/25/15 MDG 1.3 added option to generate random distribution
!--------------------------------------------------------------------------
program CTEMquatSP

use local
use typedefs
use dictmod
use quaternions
use constants
use rotations
use so3
use math

IMPLICIT NONE 

integer(kind=irg)               :: nums, seed, FZtype, FZorder, i, j, k, nt, numt, ii, jj, kk, &
                                   dimx, dimy, pgnum, np, fx, fy, fz, npix, nump, bad, nquats, one
real(kind=sgl),allocatable      :: eulers(:,:),cube(:,:,:) 
integer(kind=sgl),allocatable   :: cubeSP(:,:,:), cubeHO(:,:,:), cubeCU(:,:,:), cubeEU(:,:,:)
real(kind=dbl)                  :: x1y1(2), x2y2(2), s1, s2

type(dicttype),pointer          :: dict
real(kind=dbl)                  :: qu(4), c, s, rod(4), n, delta, qq(4)
real(kind=sgl)                  :: x, y, z, xfrac, yfrac, zfrac, tpi, eu(3), ho(3), cu(3), m, ma, sh, sc !, sig
                                
character(fnlen)                :: eulerdatafile, cubefile
character(6)                    :: FZmode               ! 'Random', 'FullSP', 'DrawFZ', 'EulerS' or 'FZonly', 
character(1)                    :: circles,verbose      ! 'y' or 'n'

namelist /Stereogram/ eulerdatafile, cubefile, pgnum, FZmode, verbose, np, npix

eulerdatafile = 'empty'
cubefile = 'empty'
pgnum = 32
FZmode = 'FullSP'
np = 2
npix = 128 
verbose = 'n'

open(UNIT=dataunit,FILE='Stereogram.nml',DELIM='apostrophe',STATUS='old')
read(UNIT=dataunit,NML=Stereogram)
close(UNIT=dataunit,STATUS='keep')

if (verbose.eq.'y') write (*,nml=Stereogram)

nump = npix + np + 1
bad = 0
nquats = 0

allocate(dict)
dict%Num_of_init = 3 
dict%Num_of_iterations = 30
dict%pgnum = pgnum
call DI_Init(dict,'VMF')

! generate a random distribution with 10^9 points and store the resulting data files
! in each of the standard projection modes
if (FZmode.eq.'Random') then
  one = 1
  seed = 432432
  tpi = 2.0*sngl(cPi)
  nquats = 1000000000
  sh = float(npix) / (3.0*sngl(cPi)/4.0)**0.3333333   ! scaled radius of homochoric ball
  sc = float(npix) / sngl(0.5D0 * LPs%ap)    ! scaled semi-edge length of cubochoric cube

! use the random quaternions to populate the various cubes
  allocate(cubeSP(-nump:nump,-nump:nump,-nump:nump))
  cubeSP = 0
  allocate(cubeHO(-nump:nump,-nump:nump,-nump:nump))
  cubeHO = 0
  allocate(cubeCU(-nump:nump,-nump:nump,-nump:nump))
  cubeCU = 0
  j = (npix+np)/2+1
  allocate(cubeEU(-nump:nump,-j:j,-nump:nump)) 
  cubeEU = 0

  do i=1,nquats
! we'll use the Marsaglia approach to generate a random quaternion
! and write them all to a file, so that this file can then be read 
! by this program again to generate the various arrays for MoR normalization
    loop1: do
      call R8VEC_uniform_01(2,seed,x1y1)
      x1y1 = 2.D0*x1y1-1.D0
      s1 = sum(x1y1*x1y1)
      if (s1.lt.1.D0) EXIT loop1
    end do loop1
  
    loop2: do
      call R8VEC_uniform_01(2,seed,x2y2)
      x2y2 = 2.D0*x2y2-1.D0
      s2 = sum(x2y2*x2y2)
      if (s2.lt.1.D0) EXIT loop2
    end do loop2
    s = dsqrt((1.D0-s1)/s2)
    qu = (/ x1y1(1), x1y1(2), x2y2(1)*s, x2y2(2)*s /)
    s = 1.D0/dsqrt(sum(qu*qu))
    qu = qu*s
    if (qu(1).lt.0.0) qu = -qu

! Euler projection
    eu = qu2eu(qu)
    if(eu(1).lt.0.0) eu(1) = eu(1) + tpi
    if(eu(2).lt.0.0) eu(2) = eu(2) + tpi
    if(eu(3).lt.0.0) eu(3) = eu(3) + tpi
    x = npix * (eu(1)-cPi)/cPi
    y = (npix/2) * (eu(2)-cPi/2.0)/(cPi/2.0)
    z = npix * (eu(3)-cPi)/cPi
    if ((abs(x).le.npix).and.(abs(y).le.npix/2).and.(abs(z).le.npix)) then
        cubeEU(nint(x),nint(y),nint(z)) = cubeEU(nint(x),nint(y),nint(z)) + 1
    end if 

! homochoric projection
    ho = qu2ho(qu)*sh
    x = ho(1)
    y = ho(2)
    z = ho(3)
    cubeHO(nint(x),nint(y),nint(z)) = cubeHO(nint(x),nint(y),nint(z)) + 1

! cubochoric projection
    cu = qu2cu(qu)*sc
    x = cu(1)
    y = cu(2)
    z = cu(3)
    cubeCU(nint(x),nint(y),nint(z)) = cubeCU(nint(x),nint(y),nint(z)) + 1

! stereographic projection
    s = float(npix)/(1.0+qu(1))
    x = qu(2)*s
    y = qu(3)*s
    z = qu(4)*s
    cubeSP(nint(x),nint(y),nint(z)) = cubeSP(nint(x),nint(y),nint(z)) + 1

    if (mod(i,1000000).eq.0) write (*,*) 'Number of quaternions generated : ',i
  end do

! finally, save all the cubes in separate files
  open(UNIT=dataunit,file='cuberandom_128_SP.raw',status='unknown',form='unformatted')
  write (dataunit) shape(cubeSP)
  write (dataunit) one, nquats
  write (dataunit) cubeSP
  close(UNIT=dataunit,status='keep')

  open(UNIT=dataunit,file='cuberandom_128_HO.raw',status='unknown',form='unformatted')
  write (dataunit) shape(cubeHO)
  write (dataunit) one, nquats
  write (dataunit) cubeHO
  close(UNIT=dataunit,status='keep')

  open(UNIT=dataunit,file='cuberandom_128_CU.raw',status='unknown',form='unformatted')
  write (dataunit) shape(cubeCU)
  write (dataunit) one, nquats
  write (dataunit) cubeCU
  close(UNIT=dataunit,status='keep')

  open(UNIT=dataunit,file='cuberandom_128_EU.raw',status='unknown',form='unformatted')
  write (dataunit) shape(cubeEU)
  write (dataunit) one, nquats
  write (dataunit) cubeEU
  close(UNIT=dataunit,status='keep')
end if


if ((FZmode.ne.'DrawFZ').and.(FZmode.ne.'Random')) then 
  open(unit=dataunit,file=eulerdatafile,form='unformatted',status='old')
  read(dataunit) nums
  allocate(eulers(3,nums))
  read(dataunit) eulers
  close(unit=dataunit,status='keep')
! first, we need to figure out how many points there are inside the Rodrigues FZ
  FZtype = FZtarray(dict%pgnum)
  FZorder = FZoarray(dict%pgnum)
  if (verbose.eq.'y') write (*,*) 'number of symmetry operators : ',dict%Nqsym, FZtype, FZorder
end if

if (FZmode.ne.'Random') nquats = nums*dict%Nqsym

! We'll create a binary raw data file with a 3D array of 256^3 points that 
! represents the 3D stereographic projection as a normalized histogram

if (FZmode.eq.'DrawFZ') then 
  numt = 0
  delta = dble(npix)**2
  do i = -npix,npix
   do j = -npix,npix
    do k = -npix,npix
       s = dble(i)**2+dble(j)**2+dble(k)**2
       if (s.lt.delta) then 
         qq = (/ delta-s, 2.D0*dble(i*npix), 2.D0*dble(j*npix), 2.D0*dble(k*npix) /) 
         qq = qq/(delta + s)
         if (IsinsideFZ(qu2ro(qq),FZtype,FZorder)) then 
          numt = numt + 1
          cube(i,j,k) = 1.0
         end if
       end if
     end do
   end do
  end do
  
  do i=0,720 
    c = float(npix)*cos(float(i)*cPi/360.0)
    s = float(npix)*sin(float(i)*cPi/360.0)
    cube(0,nint(c),nint(s)) = 1.0
    cube(nint(c),0,nint(s)) = 1.0
    cube(nint(c),nint(s),0) = 1.0
  end do
  cube(npix,-3:3,-3:3) = 1.0
  cube(-2:2,npix,-2:2) = 1.0
  cube(-1:1,-1:1,npix) = 1.0
end if


if (FZmode.eq.'EulerS') then ! Euler space representation
  i = (npix+np)/2+1
  allocate(cube(-nump:nump,-i:i,-nump:nump)) 
  cube = 0.0
  if (verbose.eq.'y') write (*,*) 'starting creation of euler space'
   numt = 0
   bad = 0
   tpi = 2.0*sngl(cPi)
   do i=1,nums
    qu = eu2qu(eulers(1:3,i))
    if (qu(1).lt.0.0) qu = -qu
    if ((abs(eulers(1,i))+abs(eulers(2,i))+abs(eulers(3,i))).ne.0.0) then 
     do j=1,dict%Nqsym
      qq = quat_mult(qu,dict%Pm(1:4,j))
      eu = qu2eu(qq)
      if(eu(1).lt.0.0) eu(1) = eu(1) + tpi
      if(eu(2).lt.0.0) eu(2) = eu(2) + tpi
      if(eu(3).lt.0.0) eu(3) = eu(3) + tpi
      x = npix * (eu(1)-cPi)/cPi
      y = (npix/2) * (eu(2)-cPi/2.0)/(cPi/2.0)
      z = npix * (eu(3)-cPi)/cPi
      if ((abs(x).le.npix).and.(abs(y).le.npix/2).and.(abs(z).le.npix)) then
        cube(nint(x),nint(y),nint(z)) = cube(nint(x),nint(y),nint(z)) + 1.0
        numt = numt+1
      else
        bad = bad+1
      end if
     end do    
    end if
   end do
! and put some identifiers in
   i = (npix+np)/2+1
   cube(-nump:(-nump+2),-i:(-i+2),-nump:(-nump+2)) = -10
   cube(-nump:(-nump+1),-i:(-i+1),(nump-1):nump) = -10
   cube(-nump:(-nump+1),(i-1):i,-nump:(-nump+1)) = -10
   cube(-nump:(-nump+1),(i-1):i,(nump-1):nump) = -10
   cube((nump-1):nump,-i:(-i+1),-nump:(-nump+1)) = -10
   cube((nump-1):nump,-i:(-i+1),(nump-1):nump) = -10
   cube((nump-1):nump,(i-1):i,-nump:(-nump+1)) = -10
   cube((nump-1):nump,(i-1):i,(nump-1):nump) = -10
end if

if (FZmode.eq.'HomocS') then ! homochoric ball
   s = float(npix) / (3.0*sngl(cPi)/4.0)**0.3333333   ! scaled radius of homochoric ball
   allocate(cube(-nump:nump,-nump:nump,-nump:nump))
   cube = 0.0
   numt = 0
   do i=1,nums
    qu = eu2qu(eulers(1:3,i))
    if (qu(1).lt.0.0) qu = -qu
    if ((abs(eulers(1,i))+abs(eulers(2,i))+abs(eulers(3,i))).ne.0.0) then 
     do j=1,dict%Nqsym
      qq = quat_mult(qu,dict%Pm(1:4,j))
      if (qq(1).lt.0.0) qq = -qq
      ho = qu2ho(qq)*s
      numt = numt + 1
      x = ho(1)
      y = ho(2)
      z = ho(3)
      cube(nint(x),nint(y),nint(z)) = cube(nint(x),nint(y),nint(z)) + 1.0
     end do
    end if
   end do
end if

if (FZmode.eq.'CubocS') then ! cubochoric cube
   s = float(npix) / sngl(0.5D0 * LPs%ap)    ! scaled semi-edge length of cubochoric cube
   allocate(cube(-nump:nump,-nump:nump,-nump:nump))
   cube = 0.0
   numt = 0
   bad = 0
   do i=1,nums
    qu = eu2qu(eulers(1:3,i))
    if (qu(1).lt.0.0) qu = -qu
    if ((abs(eulers(1,i))+abs(eulers(2,i))+abs(eulers(3,i))).ne.0.0) then 
     do j=1,dict%Nqsym
      qq = quat_mult(qu,dict%Pm(1:4,j))
      if (qq(1).lt.0.0) qq = -qq
      cu = qu2cu(qq)*s
      x = cu(1)
      y = cu(2)
      z = cu(3)
      if ((abs(x).lt.npix).and.(abs(y).le.npix).and.(abs(z).le.npix)) then
        cube(nint(x),nint(y),nint(z)) = cube(nint(x),nint(y),nint(z)) + 1.0
        numt = numt + 1
      else
        bad = bad+1
      end if
     end do
    end if
   end do
end if

if (FZmode.eq.'FZonly') then ! standard stereographic projection, Rodrigues FZ only
   allocate(cube(-nump:nump,-nump:nump,-nump:nump))
   cube = 0.0
   numt = 0
   do i=1,nums
    qu = eu2qu(eulers(1:3,i))
    if (qu(1).lt.0.0) qu = -qu
    if ((abs(eulers(1,i))+abs(eulers(2,i))+abs(eulers(3,i))).ne.0.0) then 
     do j=1,dict%Nqsym
      qq = quat_mult(qu,dict%Pm(1:4,j))
      if (qq(1).lt.0.0) qq = -qq
      if (IsinsideFZ(qu2ro(qq),FZtype,FZorder)) then 
        numt = numt + 1
        s = float(npix)/(1.0+qq(1))
        x = qq(2)*s
        y = qq(3)*s
        z = qq(4)*s
        cube(nint(x),nint(y),nint(z)) = cube(nint(x),nint(y),nint(z)) + 1.0
      end if
     end do
    end if
   end do
end if

if (FZmode.eq.'FullSP') then ! standard stereographic projection
   allocate(cube(-nump:nump,-nump:nump,-nump:nump))
   cube = 0.0
   numt = 0
   do i=1,nums
    qu = eu2qu(eulers(1:3,i))
    if (qu(1).lt.0.0) qu = -qu
    if ((abs(eulers(1,i))+abs(eulers(2,i))+abs(eulers(3,i))).ne.0.0) then 
     do j=1,dict%Nqsym
      qq = quat_mult(qu,dict%Pm(1:4,j))
      if (qq(1).lt.0.0) qq = -qq
        numt = numt + 1
        s = float(npix)/(1.0+qq(1))
        x = qq(2)*s
        y = qq(3)*s
        z = qq(4)*s
        cube(nint(x),nint(y),nint(z)) = cube(nint(x),nint(y),nint(z)) + 1.0
     end do
    end if
   end do
end if

if (FZmode.ne.'Random') then
 if (verbose.eq.'y') then
  write (*,*) 'number of starting Euler triplets : ',nums
  write (*,*) 'total number of possible points   : ',nums*dict%Nqsym
  write (*,*) 'actual number of good points      : ',numt
  if (bad.ne.0) write (*,*) 'number of bad points              : ',bad
 end if

 open(UNIT=dataunit,file=cubefile,status='unknown',form='unformatted')
 write (dataunit) shape(cube)
 write (dataunit) dict%Nqsym, nquats
 write (dataunit) cube
 close(UNIT=dataunit,status='keep')
end if
end program
