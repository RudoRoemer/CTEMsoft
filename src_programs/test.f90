!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CTEMsoft:  SRCBED.f90                                                        !
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
! -Fortran-90 Source Code-
! ###################################################################
!  Introduction to Conventional Transmission Electron Microscopy
!
!  by Marc De Graef
!
!  Cambridge University Press
!  SBN 0521629950 (Paperback)
!  SBN 0521620066 (Hardback)
! 
!  FILE: "SRCBED.f90"
!                                    created: 6/7/01  {9:29:46 AM} 
!                                last update: 6/8/2001 {7:42:24 AM} 
!  Author: Marc De Graef
!  E-mail: degraef@cmu.edu
!    mail: Carnegie Mellon
!          5000 Forbes Avenue, Department MSE
!     www: http://neon.mems.cmu.edu/degraef
!  
!  Description: Systematic Row Convergent Beam Electron Diffraction 
!               using scattering matrix approach
! 
!  History
! 
!  modified by  rev reason
!  -------- --- --- -----------
!  6/7/01   MDG 1.0 original
! ###################################################################
program SRCBED

use local
use io
use symmetryvars
use symmetry
use crystal
use files
use diffraction
use postscript

 progname = 'SRCBED.f90'
 progdesc = 'Systematic row convergent beam pattern using scattering matrix'
 call CTEMsoft
 
! first get the crystal data and microscope voltage
 SG % SYM_reduce=.TRUE.
 call CrystalData
 call GetVoltage
 call GetCameraLength(camlen)
! generate all atom positions
 call CalcPositions('v')

! open PostScript file
 write (*,"(1x,' PostScript file name : ',$)")
 read (*,"(A20)") PS % psname
 call PS_openfile
 pspage = 0
 imanum = 0
! generate a set of systematic row CBED patterns
 call SRCBEDPage
! close Postscript file
 call PS_closefile
end program

! ###################################################################
!
!  subroutine SRCBEDPage
!
!  Author: Marc De Graef
!
!  Description: draw systematic row convergent beam electron diffraction patterns
!
!  History
!
!  modified by  rev reason
!  -------- --- --- -----------
!   6/7/01  MDG 1.0 original
! ###################################################################
subroutine SRCBEDpage

use local
use constants
use postscript
use crystal
use crystalvars
use diffraction
use symmetry
use math
use io
use dynamical

real                :: laL,kt,z0,alp,thc,thb,omega_c,omega_min,omega_max,hkl(3),ind(3), &
                       dom,glen,qr,qi,bg,xgpz
integer             :: g(3),ira,dpcnt,ppi
character(1)        :: ans,c
complex             :: czero
logical             :: overlap,np,first
real,allocatable    :: row(:,:)
real,parameter      :: xoff(0:5)=(/0.0,3.3125,0.0,3.3125,0.0,3.3125/),yoff(0:5)=(/6.0,6.0,3.0,3.0,0.0,0.0/)
complex,allocatable :: SMz(:,:,:),disk(:,:,:),p(:),Sphiz(:,:)
complex,allocatable :: q(:,:),qin(:,:),qout(:,:),r(:,:)

! normal aborption factor
 call CalcUcg((/0,0,0/))
 xgpz= aimag(rlp%qg)
! camera length
 laL = sngl(cell%mLambda) * camlen
 mess = 'wavelength [nm] = '; oi_real(1) = sngl(cell%mLambda); call WriteReal(1,"(F10.6)")
 mess = ' L         [mm] = '; oi_real(1) = camlen; call WriteReal(1,"(f10.2)")
 mess = 'camera length lambda*L [mm nm] = '; oi_real(1) = laL; call WriteReal(1,"(f10.5)")
! the number of beams is the same for all patterns 
 mess = 'Maximum multiple of g to consider : '; call GetInt(1)
 ira=io_int(1)
! total number of beams
 nn = 2*ira+1
 izero = (nn-1)/2+1
 ans = 'y'
 dpcnt = 0
 imanum = 0
! main loop over patterns
 do while ((ans.eq.'y').or.(ans.eq.'Y'))
  dpcnt=dpcnt+1
! ask for parameters for this pattern
  mess = 'Diffraction vector           :'; call Message("(A)")
  call GetIndex(g,'r')
  mess = 'Enter k_t in units of |g|    : '; call GetReal(1)
  kt   = io_real(1)
  mess = 'Foil thickness        [nm]   : '; call GetReal(1)
  z0   = io_real(1)
! convert these parameters to various other parameters
! first get the Bragg angle for g
  thb = CalcDiffAngle(g(1),g(2),g(3))*0.5
  mess = 'Bragg angle [mrad] : '; oi_real(1) =thb*1000.0; call WriteReal(1,"(F10.5)")
  mess = 'Beam divergence angle [mrad] : '; call GetReal(1)
  thc  = io_real(1)/1000.0  ! convert to radians
! convert k_t to the alp and omega angles
  glen = CalcLength(float(g),'r')
  alp = -2.0*kt*thb
  omega_c = cPi*0.5+alp
  omega_min = omega_c - thc
  omega_max = omega_c + thc
! determine the number of pixels for this particular diffraction disk
  gc = thc/cell%mLambda     ! radius of disk in nm^-1
! scale bar (sc is the conversion factor from nm-1 to inches)
  sc = laL/25.4
  gci = gc*sc                 ! disk radius in inches
! define the number of pixels per inch for the computation
  ppi = 300
! for 300 dpi output we need this many pixels per disk diameter (odd number)
  npix = int(ppi * 2.0 * gci)
  if (mod(npix,2).eq.0) npix=npix+1
  mess = 'Number of pixels across disk diameter : '; oi_int(1)=npix; call WriteInt(1,"(I4)")
! step size in omega angle
  dom = (omega_max - omega_min)/float(npix-1)
! allocate dynamical variables 
  allocate(SMz(npix,nn,nn))
  allocate(Sphiz(npix,nn))
  allocate(Az(nn,nn))
  allocate(p(nn))
  allocate(DHWMz(nn,nn))
! compute the complex DHW matrix
  ind = float(g)
  do i=1,nn
   do j=1,nn
    hkl=(-ira+i-1)*ind-(-ira+j-1)*ind
    if (i.ne.j) then
     call CalcUcg(int(hkl))
     DHWMz(i,j) = cPi*cmplx(-aimag(rlp%qg),real(rlp%qg),dbl)
    else
     DHWMz(i,j) = cmplx(0.0,0.0)
    endif
   end do
  end do
  mess = 'Darwin-Howie-Whelan matrix initialized'; call Message("(A)")
! loop over all incident beam directions
  numi = 0
  allocate(q(nn,nn))
  allocate(qin(nn,nn))
  allocate(qout(nn,nn))
  allocate(r(nn,nn))
  do j=1,npix
! set omega angle
   omega = omega_min+float(j-1)*dom
   do i=1,nn
    n = -ira+i-1
! exer = excitation error
    exer = -n*glen*cos(omega)-(1.0-sqrt(1.0-(n*cell%mLambda*glen*sin(omega))**2))/cell%mLambda
    DHWMz(i,i)=2.0*cPi*cmplx(0.0,exer)
   end do
! compute the first slice scattering matrix from Taylor expansion
   sl = z0/100.0
   nsl = 100
! make sure slice thickness does not exceed 1 nm
   do while (sl.gt.0.4)
     sl = sl/2.0
     nsl = nsl*2
   end do
! use the series expansion for the scattering matrix (eq. 5.27)
   do i=1,nn
    Az(i,1:nn) = cmplx(0.0,0.0)
    Az(i,i) = cmplx(1.0,0.0)
   end do
   q = DHWMz*cmplx(sl,0.0)
   qout = q
   thr  = 1.0e-8
   zmax = 1.0
   i = 2
   do while (zmax.gt.thr)
    Az = Az + qout
    r = q/cmplx(i,0.0)
    qin = matmul(qout,r)
    zmax = maxval(cabs(qout-qin))
    qout = qin
    i = i+1
   end do
   Az = Az + qout
   numi = numi + i
   do i=1,nn
    SMz(j,i,1:nn) = Az(i,1:nn)
   end do
  end do
  mess = 'Average number of terms in series expansion : '; oi_real(1)=float(numi)/float(npix)
  call WriteReal(1,"(F8.3)")
  deallocate(q)
  deallocate(qin)
  deallocate(qout)
  deallocate(r)
! this completes the first slice scattering matrix; next multiply
! it with itself to the desired thickness
  do i=1,npix
   do j=1,nn
! the izero column of SMz is the actual initial wavefunction at z=sl
    Sphiz(i,j)=SMz(i,j,izero)
   end do
  end do
! loop to maximum thickness
  do l=2,nsl
   do ll=1,npix
    p(1:nn)=cmplx(0.0,0.0)
    do i=1,nn
     p(i)=p(i)+sum(SMz(ll,i,1:nn)*Sphiz(ll,1:nn))
    end do
    Sphiz(ll,1:nn)=p(1:nn)
   end do
  end do
! these still need to be multiplied by the normal absorption exponential
! amplitude, not the square !!!
  att = exp(-cPi*z0*xgpz)
write (*,*) 'attenuation = ',att
  Sphiz = Sphiz*att 
! deallocate all variables that are no longer needed
  deallocate(DHWMz)
  deallocate(Az)
  deallocate(SMz)
  deallocate(p)
! fill the disk variable with amplitudes
  allocate(disk(nn,npix,npix))
  do j=1,npix
   do i=1,nn
    disk(i,j,1:npix) = Sphiz(j,i)
   end do
  end do
! and get rid of phiz
  deallocate(Sphiz)
! next we apply a circular mask (using fourfold symmetry)
  czero = cmplx(0.0,0.0)
  np2 = npix/2
  nps = np2**2
  do i=1,npix/2
   is = (i-np2)**2
   do j=1,npix/2
    iq = is+(j-np2)**2
    if (iq.gt.nps) then 
     do k=1,nn
      disk(k,i,j) = czero
      disk(k,npix+1-i,j) = czero
      disk(k,i,npix+1-j) = czero
      disk(k,npix+1-i,npix+1-j) = czero
     end do
    end if
   end do
  end do
! open(unit=15,file='disk.out',status='unknown',form='unformatted')
! write (15) nn,npix
! write (15) disk
! close(unit=15,status='keep')
! allocate the image variable row (900 by whatever the vertical size is)
  npx = int(3.0 * 300)
  npy = npix
  allocate(row(npx,npy))
! next check for overlapping disks;  if they do overlap,
! then ask for coherent or incoherent illumination, otherwise
! simply take the modulus squared of the disk variable and
! place all disks in the row array
  overlap = .FALSE.
  if (thc.ge.thb) then   
   mess = 'The diffraction disks appear to overlap.'; call Message("(A)")
   mess = 'Coherent (c) or incoherent (i) illumination ? '; call Message("(A)")
   read (*,"(A1)") c
   if ((c.eq.'c').or.(c.eq.'C')) overlap = .TRUE.
  end if
  if (overlap) then
   write (*,*) 'not yet implemented'
  else
! the disks may overlap, but we are adding intensities, not amplitudes
! so the overlap does not matter.
! first the zero order beam
   row(1:900,1:npix) = cabs(czero)**2
   ipos = 450+np2+1
   do i=1,npix
    row(ipos-i,1:npix) = cabs(disk(izero,i,1:npix))**2
   end do
! negative disks
   do j=1,nn/2
     ipos = 450 - np2 - 1 + int(glen*(-ira+j-1)*ppi*sc)
     if (ipos.ge.1) then   ! the disk lies completely inside the row array
      do i=1,npix
       row(ipos+i,1:npix) = row(ipos+i,1:npix) + cabs(disk(j,i,1:npix))**2
      end do
     else    ! the disk lies partially or completely outside the visible area
      if (ipos+npix.gt.1) then   ! partially inside the visible area
       istart = 1-ipos
       do i=istart,npix-1
        row(ipos+i,1:npix) = row(ipos+i,1:npix) + cabs(disk(j,i+1,1:npix))**2
       end do
      end if   ! else is completely outside so we can skip it
     end if
   end do
! positive disks
   do j=nn/2+2,nn
    ipos = 450 - np2 - 1 + int(glen*(-ira+j-1)*ppi*sc)
    if (ipos+npix.lt.npx) then   ! the disk lies completely inside the row array
     do i=1,npix
      row(ipos+i,1:npix) = row(ipos+i,1:npix) + cabs(disk(j,i,1:npix))**2
     end do
    else    ! the disk lies partially or completely outside the visible area
     if (ipos.lt.npx) then   ! partially inside the visible area
      istop = npx-ipos+1
      do i=1,istop
       row(ipos+i-1,1:npix) = row(ipos+i-1,1:npix) + cabs(disk(j,i,1:npix))**2
      end do
     end if   ! else is completely outside so we can skip it
    end if
   end do
  end if
! deallocate the disk variable
  where(row.lt.0.0001) row = 0.0001
  row = alog10(1.0/row)
  rowmax = maxval(row)
  row = row/rowmax 
  row = 1.0-row
  deallocate(disk)
! finally, dump the row image to the Postscript file
! and create output in 2 columns, 3 rows 
  imo = mod(dpcnt-1,6)
  if (imo.eq.0) then 
   np=.TRUE.
  else
   np=.FALSE.
  endif
  mess = 'Creating systematic row CBED pattern ';
  oi_int(1:3) = g(1:3)
  call WriteInt(3,"('(',3i3,')  ')")
  if (.not.allocated(imaint)) then 
   allocate(imaint(npx,npy))
  end if
  write (*,*) xoff(imo),yoff(imo),g,np,npx,npy,laL,z0
  call DumpSRCBED(xoff(imo),yoff(imo),row,g,np,npx,npy,laL,z0)
  deallocate(row)
! another one ?
  mess = 'Another pattern ? (1/0) '; call GetInt(1)
  if (io_int(1).eq.1) then  
   ans='y' 
  else 
   ans='n'
  end if
 end do  ! main loop over patterns
end subroutine

! ###################################################################
!
!  subroutine DumpSRCBED
!
!  Author: Marc De Graef
!
!  Description: draw a single systematic row CBED pattern
!
!  History
!
!  modified by  rev reason
!  -------- --- --- -----------
!   6/ 7/01 MDG 1.0 f90
! ###################################################################
subroutine DumpSRCBED(xo,yo,row,g,np,ns,nt,laL,thick)

use local
use io
use postscript
use crystal
use crystalvars
use doublediff
use error

integer            :: i,ui,vi,wi,g(3)
real               :: xo,yo,laL,gmax,thick,row(ns,nt)
logical            :: np,first
real,parameter     :: le=3.25,he=2.9375


! do page preamble stuff if this is a new page
! [This assumes that the PostScript file has already been opened]
 if (np) then
  call PS_newpage(.FALSE.,'Systematic Row CBED Patterns')
  call PS_text(5.25,-0.05,'scale bar in reciprocal nm')
  gmax = laL
  call PS_textvar(5.25,PS % psfigheight+0.02,'Camera Constant [nm mm]',gmax)
  call PS_setfont(PSfonts(2),0.15)
  call PS_text(-0.25,PS % psfigheight+0.02,'Structure File : '//cell % fname)
 end if
! draw frame and related stuff
 call PS_setlinewidth(0.012)
 call PS_balloon(xo,yo,le,he,0.0312)
 call PS_setlinewidth(0.001)
! systematic row vector
 call PS_setfont(PSfonts(2),0.12)
 call PS_text(xo+0.05,yo+he-0.15,'Systematic row vector : ')
 ui=g(1)
 vi=g(2)
 wi=g(3)
 call PrintIndices('r',ui,vi,wi,xo+1.3,yo+he-0.15)
! foil thickness
 call PS_setfont(PSfonts(2),0.12)
 pp=p
 call PS_textvar(xo+0.05,yo+he-0.30,'Thickness ',thick)
! convergence angle

! scale bar (sc is the conversion factor from nm-1 to inches)
 sc = laL/25.4
 call PS_setlinewidth(0.020)
 call PS_line(xo+0.05,yo+0.06,xo+0.05+5.0*sc,yo+0.06)
 call PS_setfont(PSfonts(2),0.15)
 call PS_text(xo+0.05+2.5*sc,yo+0.10,'5 ')
! next dump the image
 scl = 3.0
 do i=1,ns
  do ii=1,nt
   imaint(i,ii)=int(255.0*row(i,ii))
  end do
 end do
! and dump the image
 call PS_DumpImage(xo+(le-3.0)/2.0,yo+1.5-float(nt)/2.0/300.0,ns,nt,scl)
end subroutine
