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
! CTEMsoft2013: CTEMSRCBED.f90
!--------------------------------------------------------------------------
!
! PROGRAM: CTEMSRCBED 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Systematic Row Convergent Beam Electron Diffraction 
! > using scattering matrix approach
!
!> @todo implement coherent vs. incoherent overlap of reflections
! 
!> @date  6/7/01   MDG 1.0 original
!> @date 4/18/13 MDG 2.0 rewrite
!--------------------------------------------------------------------------
program CTEMSRCBED

use local
use io
use symmetryvars
use symmetry
use crystal
use files
use diffraction
use postscript

IMPLICIT NONE

real(kind=sgl)		:: io_real(1)

 progname = 'CTEMSRCBED.f90'
 progdesc = 'Systematic row convergent beam pattern using scattering matrix'
 call CTEMsoft
 
! first get the crystal data and microscope voltage
 SG % SYM_reduce=.TRUE.
 call CrystalData
 call GetVoltage
 call ReadValue(' Enter the diffraction camera length L [mm,R]: ', io_real, 1)
 camlen = io_real(1)

! generate all atom positions
 call CalcPositions('v')

! open PostScript file
! write (*,"(1x,' PostScript file name : ',$)")
! read (*,"(A20)") PS % psname
 call PS_openfile
 PS % pspage = 0
 imanum = 0

! generate a set of systematic row CBED patterns
 call SRCBEDPage

! close Postscript file
 call PS_closefile

end program CTEMSRCBED

!--------------------------------------------------------------------------
!
! PROGRAM: SRCBEDPage 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief draw systematic row convergent beam electron diffraction patterns
! 
!> @todo split this humongously long subroutine into smaller reusable chunks
!
!> @date  6/7/01   MDG 1.0 original
!> @date 4/18/13 MDG 2.0 rewrite
!--------------------------------------------------------------------------
subroutine SRCBEDPage

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

IMPLICIT NONE

real(kind=sgl)                			:: laL,kt,z0,alp,thc,thb,omega_c,omega_min,omega_max,hkl(3),ind(3), &
                       					  dom,glen,xgpz, io_real(1),sc,omega,exer,sl,nsl,thr,zmax,att,gc,gci
integer(kind=irg)             		:: g(3),ira,dpcnt,ppi,io_int(1),nn,izero,npix,i,j,numi,n,l,ll,np2,nps,&
							  is,iq,npx,npy,ipos,istart,istop,rowmax,imo,k
character(1)        				:: ans,c
complex(kind=sgl)   			:: czero
logical             					:: overlap,np
real(kind=sgl),allocatable    		:: row(:,:)
real(kind=sgl),parameter      		:: xoff(0:5)=(/0.0,3.3125,0.0,3.3125,0.0,3.3125/),yoff(0:5)=(/6.0,6.0,3.0,3.0,0.0,0.0/)
complex(kind=sgl),allocatable 	:: SMz(:,:,:),disk(:,:,:),p(:),Sphiz(:,:)
complex(kind=sgl),allocatable 	:: q(:,:),qin(:,:),qout(:,:),r(:,:)

! normal aborption factor
 call CalcUcg((/0,0,0/))
 xgpz= aimag(rlp%qg)

! camera length
 laL = sngl(mLambda) * camlen
 io_real(1) = sngl(mLambda)
 call WriteValue(' wavelength [nm] = ', io_real, 1, "(F10.6)")
 io_real(1) = camlen
 call WriteValue('  L         [mm] = ', io_real, 1, "(f10.2)")
 io_real(1) = laL
 call WriteValue(' camera length lambda*L [mm nm] = ', io_real, 1, "(f10.5)")

! the number of beams is the same for all patterns 
 call ReadValue(' Maximum multiple of g to consider : ', io_int, 1)
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
  mess = ' Diffraction vector           :'; call Message("(A)")
  call GetIndex(g,'r')
  call ReadValue(' Enter k_t in units of |g|    : ', io_real, 1)
  kt   = io_real(1)
  call ReadValue(' Foil thickness        [nm]   : ', io_real, 1)
  z0   = io_real(1)

! convert these parameters to various other parameters
! first get the Bragg angle for g
  thb = CalcDiffAngle(g(1),g(2),g(3))*0.5
  io_real(1) = thb*1000.0
  call WriteValue(' Bragg angle [mrad] : ', io_real, 1, "(F10.5)")
  call ReadValue(' Beam divergence angle [mrad] : ', io_real, 1)
  thc  = io_real(1)/1000.0  ! convert to radians

! convert k_t to the alp and omega angles
  glen = CalcLength(float(g),'r')
  alp = -2.0*kt*thb
  omega_c = cPi*0.5+alp
  omega_min = omega_c - thc
  omega_max = omega_c + thc

! determine the number of pixels for this particular diffraction disk
  gc = thc/mLambda     ! radius of disk in nm^-1

! scale bar (sc is the conversion factor from nm-1 to inches)
  sc = laL/25.4
  gci = gc*sc                 ! disk radius in inches
write (*,*) 'disk radius = ',gci
! define the number of pixels per inch for the computation
  ppi = 300

! for 300 dpi output we need this many pixels per disk diameter (odd number)
  npix = int(ppi * 2.0 * gci)
  if (mod(npix,2).eq.0) npix=npix+1
  io_int(1) = npix
  call WriteValue(' Number of pixels across disk diameter : ', io_int, 1, "(I4)")

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
    exer = -n*glen*cos(omega)-(1.0-sqrt(1.0-(n*mLambda*glen*sin(omega))**2))/mLambda
    DHWMz(i,i)=2.0*cPi*cmplx(0.0,exer)
   end do

! compute the first slice scattering matrix from Taylor expansion
   sl = z0/100.0
   nsl = 100

! make sure slice thickness does not exceed 0.4 nm (somewhat arbitrary)
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

  io_real(1) = float(numi)/float(npix)
  call WriteValue(' Average number of terms in series expansion : ', io_real, 1, "(F8.3)")

  deallocate(q)
  deallocate(qin)
  deallocate(qout)
  deallocate(r)

! this completes the first slice scattering matrix; next multiply
! it with itself to the desired thickness
! the izero column of SMz is the actual initial wavefunction at z=sl
  Sphiz(1:npix,1:nn)=SMz(1:npix,1:nn,izero)

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

! allocate the image variable row (900x900 for historical reasons)
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
  imanum = 0
  imo = mod(dpcnt-1,6)
  if (imo.eq.0) then 
   np=.TRUE.
  else
   np=.FALSE.
  endif
  io_int(1:3) = g(1:3)
  call WriteValue(' Creating systematic row CBED pattern ', io_int, 3, "('(',3i3,')  ')")

  if (.not.allocated(imaint)) then 
   allocate(imaint(npx,npy))
  end if
  
  call DumpSRCBED(xoff(imo),yoff(imo),row,g,np,npx,npy,laL,z0)
  
  deallocate(row)

! another one ?
  call ReadValue(' Another pattern ? (1/0) ', io_int, 1)
  if (io_int(1).eq.1) then  
   ans='y' 
  else 
   ans='n'
  end if

 end do  ! main loop over patterns

end subroutine SRCBEDPage

!--------------------------------------------------------------------------
!
! SUBROUTINE: DumpSRCBED 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief draw a single systematic row CBED pattern
!
!> @date  6/7/01   MDG 1.0 original
!> @date 4/18/13 MDG 2.0 rewrite
!--------------------------------------------------------------------------
subroutine DumpSRCBED(xo,yo,row,g,np,ns,nt,laL,thick)

use local
use io
use postscript
use crystal
use crystalvars
use doublediff
use error

IMPLICIT NONE

real(kind=sgl),INTENT(IN)	:: xo
real(kind=sgl),INTENT(IN)	:: yo
real(kind=sgl),INTENT(IN)	:: row(ns,nt)
integer(kind=irg),INTENT(IN)	:: g(3)
logical,INTENT(IN)			:: np
integer(kind=irg),INTENT(IN)	:: ns
integer(kind=irg),INTENT(IN)	:: nt
real(kind=sgl),INTENT(IN)	:: laL
real(kind=sgl),INTENT(IN)	:: thick


integer(kind=irg)            		:: i,ui,vi,wi,ind(3),ii
real(kind=sgl)               		:: gmax,pp,p,sc,scl
real(kind=sgl),parameter   	:: le=3.25,he=2.9375


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
imanum = 0
 scl = 3.0
 do i=1,ns
  do ii=1,nt
   imaint(i,ii)=int(255.0*row(i,ii))
  end do
 end do

! and dump the image
 call PS_DumpImageDistort(xo+(le-3.0)/2.0,yo+1.5-float(nt)/2.0/300.0,ns,nt,scl,scl)

end subroutine DumpSRCBED
